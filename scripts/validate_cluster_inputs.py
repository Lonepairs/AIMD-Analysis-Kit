#!/usr/bin/env python3
"""Read-only validation of the eight cluster AIMD/Bader input cases.

The report directory must be outside ``--data-root``.  The script never edits
the source tree and writes only portable relative paths plus SHA-256 hashes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aimd_analysis.inventory import (  # noqa: E402
    DEFAULT_BADER_STEPS,
    DEFAULT_CASES,
    InventoryError,
    inspect_xdatcar,
    parse_acf,
)
from aimd_analysis.potcar import PotcarError, validate_potcar_species  # noqa: E402
from aimd_analysis.structure import parse_poscar  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def _read_variant_manifest(case_dir: Path, override: Path | None = None) -> tuple[dict[str, object] | None, str | None]:
    candidates = [override] if override is not None else [case_dir / name for name in ("analysis_input_manifest.json", "input_manifest.json", "manifest.json")]
    for path in candidates:
        if path is None:
            continue
        if not path.is_file():
            continue
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            return None, f"cannot read {name}: {exc}"
        if not isinstance(value, dict):
            return None, f"{name} must contain a JSON object"
        return value, None
    return None, None


def _check_manifest(
    manifest: dict[str, object] | None,
    case_dir: Path,
    root: Path,
    xdatcar: Path,
    poscar: Path,
) -> list[str]:
    if manifest is None:
        return ["updated neutral trajectory is configured by case matrix; no checksum manifest supplied"]
    errors: list[str] = []
    variant = manifest.get("trajectory_variant")
    if variant != "updated_neutral":
        errors.append(f"CuxPd1 neutral manifest trajectory_variant={variant!r}, expected 'updated_neutral'")
    for key, path in (("xdatcar_sha256", xdatcar), ("poscar_sha256", poscar)):
        expected = manifest.get(key)
        if expected is not None and expected != _sha256(path):
            errors.append(f"{key} does not match manifest")
    matching = manifest.get("matching_inputs", {})
    if matching is not None and not isinstance(matching, dict):
        errors.append("matching_inputs must be an object")
    if isinstance(matching, dict):
        for label, raw_paths in matching.items():
            values = raw_paths if isinstance(raw_paths, list) else [raw_paths]
            for raw_path in values:
                path = Path(str(raw_path))
                if not path.is_absolute():
                    path = case_dir / path
                try:
                    path.resolve().relative_to(root.resolve())
                except ValueError:
                    errors.append(f"matching {label} input is outside data-root")
                    continue
                if not path.is_file():
                    errors.append(f"matching {label} input missing: {_relative(path, root)}")
    return errors


def _validate_case(root: Path, system: str, ensemble: str, start: int, end: int, bader_steps: tuple[int, ...], cuxpd1_manifest: Path | None = None) -> dict[str, object]:
    case_dir = root / system / ensemble
    xdatcar = case_dir / "XDATCAR"
    poscar = case_dir / "POSCAR"
    row: dict[str, object] = {
        "case": f"{system}/{ensemble}",
        "system": system,
        "ensemble": ensemble,
        "trajectory_variant": "updated_neutral" if system == "CuxPd1" and ensemble == "aimd1" else "k_containing_2k" if ensemble == "aimd1_e" else "neutral",
        "xdatcar": _relative(xdatcar, root),
        "poscar": _relative(poscar, root),
        "potcar": _relative(case_dir / "POTCAR", root),
        "xdatcar_sha256": "",
        "poscar_sha256": "",
        "potcar_sha256": "",
        "species": "",
        "counts": "",
        "frame_count": "",
        "requested_range_complete": False,
        "bader_valid": 0,
        "bader_expected": len(bader_steps),
        "status": "INVALID",
        "messages": [],
    }
    messages: list[str] = []
    structure = None
    if not poscar.is_file():
        messages.append("missing POSCAR")
    else:
        try:
            structure = parse_poscar(poscar)
            row.update({
                "poscar_sha256": _sha256(poscar),
                "species": json.dumps(list(structure.species), separators=(",", ":")),
                "counts": json.dumps(list(structure.counts), separators=(",", ":")),
            })
        except (OSError, ValueError) as exc:
            messages.append(f"POSCAR: {exc}")
    if not xdatcar.is_file():
        messages.append("missing XDATCAR")
    else:
        try:
            trajectory = inspect_xdatcar(xdatcar, start, end)
            row.update({
                "xdatcar_sha256": _sha256(xdatcar),
                "frame_count": trajectory["frame_count"],
                "requested_range_complete": trajectory["requested_range_complete"],
            })
            if structure is not None and (tuple(trajectory["species"]), tuple(trajectory["counts"])) != (structure.species, structure.counts):
                messages.append("XDATCAR header species/count order differs from POSCAR")
            if not trajectory["requested_range_complete"]:
                messages.append(f"requested frame range {start}-{end} is incomplete")
        except (OSError, InventoryError) as exc:
            messages.append(f"XDATCAR: {exc}")
    potcar = case_dir / "POTCAR"
    if structure is None:
        messages.append("POTCAR cannot be checked without a valid POSCAR")
    elif not potcar.is_file():
        messages.append("missing POTCAR")
    else:
        try:
            validated = validate_potcar_species(potcar, structure.species)
            row["potcar_sha256"] = _sha256(potcar)
            row["zval"] = json.dumps(validated["zval"], sort_keys=True, separators=(",", ":"))
            if ensemble == "aimd1_e" and validated["zval"].get("K") != 7.0:
                messages.append(f"K ZVAL is {validated['zval'].get('K')!r}, expected 7.000")
            for element in ("Cu", "Pd", "O", "H"):
                if element in structure.atom_species and element not in validated["zval"]:
                    messages.append(f"POTCAR has no ZVAL for {element}")
        except (OSError, PotcarError, ValueError) as exc:
            messages.append(f"POTCAR: {exc}")

    bader_dir = root / system / ("bader_e" if ensemble == "aimd1_e" else "bader")
    valid_bader = 0
    for step in bader_steps:
        snapshot = bader_dir / f"bader_{system}_{step}"
        snapshot_poscar = snapshot / "POSCAR"
        acf = snapshot / "ACF.dat"
        if not snapshot_poscar.is_file() or not acf.is_file():
            continue
        try:
            snapshot_structure = parse_poscar(snapshot_poscar)
            acf_info = parse_acf(acf)
            if structure is None or (snapshot_structure.species, snapshot_structure.counts) != (structure.species, structure.counts):
                continue
            if acf_info["atom_rows"] != snapshot_structure.atom_count:
                continue
            valid_bader += 1
        except (OSError, ValueError, InventoryError):
            continue
    row["bader_valid"] = valid_bader
    if valid_bader != len(bader_steps):
        messages.append(f"Bader snapshots complete {valid_bader}/{len(bader_steps)}")

    if system == "CuxPd1" and ensemble == "aimd1" and xdatcar.is_file() and poscar.is_file():
        manifest, manifest_error = _read_variant_manifest(case_dir, cuxpd1_manifest)
        if manifest_error:
            messages.append(manifest_error)
        else:
            messages.extend(_check_manifest(manifest, case_dir, root, xdatcar, poscar))
    row["messages"] = messages
    row["status"] = "OK" if not messages else "INVALID"
    return row


def _check_compositions(rows: list[dict[str, object]]) -> None:
    by_system = {str(row["system"]): row for row in rows if row["ensemble"] == "aimd1"}
    by_e = {str(row["system"]): row for row in rows if row["ensemble"] == "aimd1_e"}
    for system in sorted(set(by_system) & set(by_e)):
        neutral = json.loads(str(by_system[system]["counts"])) if by_system[system]["counts"] else {}
        electro = json.loads(str(by_e[system]["counts"])) if by_e[system]["counts"] else {}
        n_species = json.loads(str(by_system[system]["species"])) if by_system[system]["species"] else []
        e_species = json.loads(str(by_e[system]["species"])) if by_e[system]["species"] else []
        n = dict(zip(n_species, neutral))
        e = dict(zip(e_species, electro))
        if e.get("K", 0) != 2:
            by_e[system]["messages"].append(f"aimd1_e must contain exactly 2 K atoms, found {e.get('K', 0)}")
        if e.get("O", 0) != n.get("O", 0) - 4 or e.get("H", 0) != n.get("H", 0) - 8:
            by_e[system]["messages"].append("aimd1_e must contain four fewer H2O molecules (O -4, H -8)")
        for element in set(n) | set(e):
            if element not in {"O", "H", "K"} and n.get(element, 0) != e.get(element, 0):
                by_e[system]["messages"].append(f"non-water composition changed for {element}: {n.get(element, 0)} -> {e.get(element, 0)}")
        for row in (by_system[system], by_e[system]):
            row["status"] = "OK" if not row["messages"] else "INVALID"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--start", type=int, default=5001)
    parser.add_argument("--end", type=int, default=9000)
    parser.add_argument("--bader-steps", type=int, nargs="+", default=list(DEFAULT_BADER_STEPS))
    parser.add_argument("--case", action="append", dest="cases", default=None)
    parser.add_argument("--cuxpd1-neutral-manifest", type=Path, default=None, help="Override the CuxPd1/aimd1 input manifest used for updated-neutral verification")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    root = args.data_root.expanduser().resolve()
    report_dir = args.report_dir.expanduser().resolve()
    if report_dir == root or root in report_dir.parents:
        parser.error("report-dir must be outside data-root; production validation cannot write into source data")
    if not root.is_dir():
        parser.error(f"data-root is not a directory: {root}")
    if report_dir.exists() and any(report_dir.iterdir()) and not args.overwrite:
        parser.error(f"report directory is non-empty: {report_dir}; use --overwrite or a new staging path")
    steps = tuple(args.bader_steps)
    if not steps or tuple(sorted(set(steps))) != steps or any(step <= 0 for step in steps):
        parser.error("bader steps must be unique, positive, and increasing")
    rows: list[dict[str, object]] = []
    for case in (args.cases or DEFAULT_CASES):
        parts = case.replace("\\", "/").split("/")
        if len(parts) != 2:
            parser.error(f"invalid case: {case}")
        rows.append(_validate_case(root, parts[0], parts[1], args.start, args.end, steps, args.cuxpd1_neutral_manifest))
    _check_compositions(rows)
    report_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": 1,
        "data_root_label": "data_root",
        "requested_frame_range": {"start": args.start, "end": args.end},
        "bader_steps": list(steps),
        "case_count": len(rows),
        "ok_count": sum(row["status"] == "OK" for row in rows),
        "invalid_count": sum(row["status"] != "OK" for row in rows),
        "cases": rows,
        "notes": [
            "No PZC or capacitance parameters are requested or validated.",
            "CuxPd1/aimd1 requires a checksum manifest to prove which updated neutral trajectory was selected; without one the report records a warning.",
        ],
    }
    (report_dir / "cluster_input_validation.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    fields = ["case", "system", "ensemble", "trajectory_variant", "xdatcar", "poscar", "potcar", "xdatcar_sha256", "poscar_sha256", "potcar_sha256", "species", "counts", "zval", "frame_count", "requested_range_complete", "bader_valid", "bader_expected", "status", "messages"]
    with (report_dir / "cluster_input_validation.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: " | ".join(row[field]) if field == "messages" else row.get(field, "") for field in fields})
    print(f"Validated {report['ok_count']}/{report['case_count']} cases; report: {report_dir}")
    return 0 if report["invalid_count"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
