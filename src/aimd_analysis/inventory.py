"""Read-only inventory and validation for AIMD case inputs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, TextIO

from .io import prepare_output_dir


DEFAULT_CASES = (
    "Cu/aimd1",
    "Cu/aimd1_e",
    "CuxPd1/aimd1",
    "CuxPd1/aimd1_e",
    "Pd/aimd1",
    "Pd/aimd1_e",
    "PdxCu1/aimd1",
    "PdxCu1/aimd1_e",
)
DEFAULT_BADER_STEPS = tuple(range(5500, 9001, 500))
CONFIGURATION_RE = re.compile(r"(?:direct|cartesian)?\s*configuration\s*=\s*(\d+)", re.I)


class InventoryError(ValueError):
    """Raised when a source file cannot be safely inventoried."""


@dataclass(frozen=True)
class StructureInfo:
    species: tuple[str, ...]
    counts: tuple[int, ...]
    atom_count: int
    coordinate_mode: str

    @property
    def composition(self) -> dict[str, int]:
        return dict(zip(self.species, self.counts))


def _read_header(lines: list[str], path: Path) -> tuple[tuple[str, ...], tuple[int, ...], int]:
    if len(lines) < 7:
        raise InventoryError(f"{path.name}: fewer than 7 header lines")
    try:
        scale = float(lines[1].split()[0])
        lattice = [[float(value) for value in lines[i].split()[:3]] for i in range(2, 5)]
    except (ValueError, IndexError) as exc:
        raise InventoryError(f"{path.name}: malformed scale or lattice") from exc
    if not math.isfinite(scale) or scale == 0 or any(len(row) != 3 for row in lattice):
        raise InventoryError(f"{path.name}: invalid scale or lattice")
    if any(not math.isfinite(value) for row in lattice for value in row):
        raise InventoryError(f"{path.name}: non-finite lattice value")
    determinant = (
        lattice[0][0] * (lattice[1][1] * lattice[2][2] - lattice[1][2] * lattice[2][1])
        - lattice[0][1] * (lattice[1][0] * lattice[2][2] - lattice[1][2] * lattice[2][0])
        + lattice[0][2] * (lattice[1][0] * lattice[2][1] - lattice[1][1] * lattice[2][0])
    )
    if abs(determinant) < 1e-12:
        raise InventoryError(f"{path.name}: singular lattice")

    species_fields = lines[5].split()
    count_fields = lines[6].split()
    try:
        counts = tuple(int(value) for value in count_fields)
    except ValueError as exc:
        raise InventoryError(f"{path.name}: VASP4/count-only headers are unsupported") from exc
    if not species_fields or len(species_fields) != len(counts):
        raise InventoryError(f"{path.name}: species and count columns do not match")
    if any(count <= 0 for count in counts):
        raise InventoryError(f"{path.name}: atom counts must be positive")
    if any(not re.fullmatch(r"[A-Z][a-z]?", species) for species in species_fields):
        raise InventoryError(f"{path.name}: invalid element symbol in species line")
    return tuple(species_fields), counts, sum(counts)


def _coordinate_mode_and_start(lines: list[str], path: Path) -> tuple[str, int]:
    cursor = 7
    if lines[cursor].strip().lower().startswith("s"):
        cursor += 1
    if cursor >= len(lines):
        raise InventoryError(f"{path.name}: missing coordinate mode")
    mode_text = lines[cursor].strip().lower()
    if mode_text.startswith("d"):
        mode = "Direct"
    elif mode_text.startswith(("c", "k")):
        mode = "Cartesian"
    else:
        raise InventoryError(f"{path.name}: unknown coordinate mode {lines[cursor]!r}")
    return mode, cursor + 1


def parse_poscar(path: Path) -> StructureInfo:
    """Read and validate a POSCAR header and all atom coordinates."""
    try:
        lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise InventoryError(f"cannot read {path.name}: {exc}") from exc
    species, counts, atom_count = _read_header(lines, path)
    mode, coordinate_start = _coordinate_mode_and_start(lines, path)
    _validate_coordinate_rows(iter(lines[coordinate_start:]), atom_count, path.name)
    return StructureInfo(species, counts, atom_count, mode)


def _validate_coordinate_rows(lines: Iterable[str], atom_count: int, label: str) -> None:
    found = 0
    for line in lines:
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < 3:
            raise InventoryError(f"{label}: malformed coordinate row {found + 1}")
        try:
            xyz = [float(value) for value in fields[:3]]
        except ValueError as exc:
            raise InventoryError(f"{label}: non-numeric coordinate row {found + 1}") from exc
        if not all(math.isfinite(value) for value in xyz):
            raise InventoryError(f"{label}: non-finite coordinate row {found + 1}")
        found += 1
        if found == atom_count:
            return
    raise InventoryError(f"{label}: expected {atom_count} coordinate rows, found {found}")


def _parse_xdatcar_header(stream: TextIO, path: Path) -> tuple[tuple[str, ...], tuple[int, ...], int]:
    header = [stream.readline() for _ in range(7)]
    if any(line == "" for line in header):
        raise InventoryError(f"{path.name}: incomplete XDATCAR header")
    species, counts, atom_count = _read_header([line.rstrip("\r\n") for line in header], path)
    return species, counts, atom_count


def inspect_xdatcar(path: Path, start: int, end: int) -> dict[str, object]:
    """Stream-validate frame markers, coordinates, monotonicity, and requested range."""
    if start < 1 or end < start:
        raise InventoryError("frame range must satisfy 1 <= start <= end")
    try:
        with path.open("r", encoding="utf-8", errors="strict") as stream:
            species, counts, atom_count = _parse_xdatcar_header(stream, path)
            frame_numbers: list[int] = []
            previous: int | None = None
            line_number = 7
            while True:
                line = stream.readline()
                if line == "":
                    break
                line_number += 1
                if not line.strip():
                    continue
                match = CONFIGURATION_RE.search(line)
                if match is None:
                    raise InventoryError(
                        f"{path.name}:{line_number}: expected configuration marker, got {line.strip()!r}"
                    )
                frame_number = int(match.group(1))
                if previous is not None and frame_number <= previous:
                    raise InventoryError(
                        f"{path.name}:{line_number}: configuration numbers are not strictly increasing"
                    )
                previous = frame_number
                frame_numbers.append(frame_number)
                _validate_coordinate_rows(stream, atom_count, f"{path.name} configuration {frame_number}")
                line_number += atom_count
    except (OSError, UnicodeError) as exc:
        raise InventoryError(f"cannot read {path.name}: {exc}") from exc

    if not frame_numbers:
        raise InventoryError(f"{path.name}: no configurations found")
    requested = [frame for frame in frame_numbers if start <= frame <= end]
    expected_count = end - start + 1
    missing = _missing_ranges(start, end, requested)
    return {
        "species": list(species),
        "counts": list(counts),
        "atom_count": atom_count,
        "frame_count": len(frame_numbers),
        "first_configuration": frame_numbers[0],
        "last_configuration": frame_numbers[-1],
        "requested_start": start,
        "requested_end": end,
        "requested_frame_count": len(requested),
        "expected_requested_frame_count": expected_count,
        "missing_requested_ranges": missing,
        "requested_range_complete": len(requested) == expected_count and not missing,
    }


def _missing_ranges(start: int, end: int, present: list[int]) -> list[list[int]]:
    present_set = set(present)
    ranges: list[list[int]] = []
    missing_start: int | None = None
    for frame in range(start, end + 1):
        if frame not in present_set and missing_start is None:
            missing_start = frame
        elif frame in present_set and missing_start is not None:
            ranges.append([missing_start, frame - 1])
            missing_start = None
    if missing_start is not None:
        ranges.append([missing_start, end])
    return ranges


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_acf(path: Path) -> dict[str, object]:
    """Count indexed ACF rows and validate their atom indices and charge fields."""
    rows: list[int] = []
    try:
        with path.open("r", encoding="utf-8", errors="strict") as stream:
            for line in stream:
                fields = line.split()
                if not fields:
                    continue
                if not fields[0].isdigit():
                    continue
                if len(fields) < 5:
                    raise InventoryError(f"{path.name}: incomplete ACF atom row")
                try:
                    charge = float(fields[4])
                except ValueError as exc:
                    raise InventoryError(f"{path.name}: invalid charge for atom {fields[0]}") from exc
                if not math.isfinite(charge):
                    raise InventoryError(f"{path.name}: non-finite charge for atom {fields[0]}")
                rows.append(int(fields[0]))
    except (OSError, UnicodeError) as exc:
        raise InventoryError(f"cannot read {path.name}: {exc}") from exc
    expected_indices = list(range(1, len(rows) + 1))
    if rows != expected_indices:
        raise InventoryError(f"{path.name}: atom indices are not contiguous from 1")
    if not rows:
        raise InventoryError(f"{path.name}: no indexed atom rows found")
    return {"atom_rows": len(rows), "indices_contiguous": True}


def _relative(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def _case_parts(case: str) -> tuple[str, str]:
    parts = case.replace("\\", "/").split("/")
    valid_component = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]*")
    if len(parts) != 2 or any(not valid_component.fullmatch(part) for part in parts):
        raise InventoryError(f"invalid case {case!r}; expected SYSTEM/ENSEMBLE")
    return parts[0], parts[1]


def _inspect_bader_case(
    root: Path,
    system: str,
    ensemble: str,
    reference: StructureInfo,
    steps: tuple[int, ...],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    is_k_case = ensemble == "aimd1_e"
    bader_dir = "bader_e" if is_k_case else "bader"
    snapshot_rows: list[dict[str, object]] = []
    valid_steps: list[int] = []
    messages: list[str] = []
    for step in steps:
        snapshot_dir = root / system / bader_dir / f"bader_{system}_{step}"
        poscar_path = snapshot_dir / "POSCAR"
        acf_path = snapshot_dir / "ACF.dat"
        row: dict[str, object] = {
            "system": system,
            "ensemble": ensemble,
            "step": step,
            "snapshot_path": _relative(snapshot_dir, root),
            "poscar_path": _relative(poscar_path, root),
            "acf_path": _relative(acf_path, root),
            "status": "MISSING",
            "atom_count": None,
            "acf_atom_rows": None,
            "poscar_sha256": None,
            "acf_sha256": None,
            "message": "",
        }
        if not poscar_path.is_file() or not acf_path.is_file():
            absent = [name for name, path in (("POSCAR", poscar_path), ("ACF.dat", acf_path)) if not path.is_file()]
            row["message"] = f"missing {', '.join(absent)}"
            snapshot_rows.append(row)
            messages.append(f"{step}: {row['message']}")
            continue
        try:
            snapshot = parse_poscar(poscar_path)
            acf = parse_acf(acf_path)
            if (snapshot.species, snapshot.counts) != (reference.species, reference.counts):
                raise InventoryError("Bader POSCAR species/count order differs from trajectory reference")
            if acf["atom_rows"] != snapshot.atom_count:
                raise InventoryError(
                    f"ACF has {acf['atom_rows']} atoms but POSCAR has {snapshot.atom_count}"
                )
            row.update(
                {
                    "status": "OK",
                    "atom_count": snapshot.atom_count,
                    "acf_atom_rows": acf["atom_rows"],
                    "poscar_sha256": sha256_file(poscar_path),
                    "acf_sha256": sha256_file(acf_path),
                    "message": "",
                }
            )
            valid_steps.append(step)
        except InventoryError as exc:
            row["status"] = "INVALID"
            row["message"] = str(exc)
            messages.append(f"{step}: {exc}")
        snapshot_rows.append(row)

    coverage = {
        "expected_steps": list(steps),
        "valid_steps": valid_steps,
        "missing_or_invalid_steps": [step for step in steps if step not in valid_steps],
        "valid_count": len(valid_steps),
        "expected_count": len(steps),
        "complete": len(valid_steps) == len(steps),
        "messages": messages,
    }
    return snapshot_rows, coverage


def inventory_cases(
    data_root: Path,
    cases: Iterable[str] = DEFAULT_CASES,
    start: int = 5001,
    end: int = 9000,
    bader_steps: Iterable[int] = DEFAULT_BADER_STEPS,
) -> dict[str, object]:
    """Inspect trajectory/reference/Bader inputs and return a portable manifest."""
    root = data_root.expanduser().resolve(strict=True)
    if not root.is_dir():
        raise InventoryError(f"data root is not a directory: {root}")
    if start < 1 or end < start:
        raise InventoryError("frame range must satisfy 1 <= start <= end")
    selected_cases = list(cases)
    steps = tuple(int(step) for step in bader_steps)
    if not selected_cases:
        raise InventoryError("at least one case is required")
    if len(set(selected_cases)) != len(selected_cases):
        raise InventoryError("duplicate case IDs are not allowed")
    if any(step <= 0 for step in steps) or tuple(sorted(set(steps))) != steps:
        raise InventoryError("Bader steps must be unique, positive, and increasing")

    case_rows: list[dict[str, object]] = []
    bader_rows: list[dict[str, object]] = []
    for case in selected_cases:
        system, ensemble = _case_parts(case)
        case_dir = root / system / ensemble
        xdatcar = case_dir / "XDATCAR"
        reference_path = case_dir / "POSCAR"
        errors: list[str] = []
        row: dict[str, object] = {
            "case_id": f"{system}__{ensemble}",
            "system": system,
            "ensemble": ensemble,
            "trajectory_variant": (
                "updated_neutral" if system == "CuxPd1" and ensemble == "aimd1"
                else "k_containing_2k" if ensemble == "aimd1_e"
                else "neutral"
            ),
            "xdatcar_path": _relative(xdatcar, root),
            "reference_path": _relative(reference_path, root),
            "xdatcar_sha256": None,
            "reference_sha256": None,
            "species": None,
            "counts": None,
            "composition": None,
            "atom_count": None,
            "frame_count": None,
            "first_configuration": None,
            "last_configuration": None,
            "requested_start": start,
            "requested_end": end,
            "requested_frame_count": None,
            "requested_range_complete": False,
            "bader_valid_snapshots": None,
            "bader_expected_snapshots": len(steps),
            "status": "INVALID",
            "messages": [],
        }
        reference: StructureInfo | None = None
        if not reference_path.is_file():
            errors.append(f"missing reference POSCAR: {row['reference_path']}")
        else:
            try:
                reference = parse_poscar(reference_path)
                row.update(
                    {
                        "reference_sha256": sha256_file(reference_path),
                        "species": list(reference.species),
                        "counts": list(reference.counts),
                        "composition": reference.composition,
                        "atom_count": reference.atom_count,
                    }
                )
            except InventoryError as exc:
                errors.append(f"reference POSCAR: {exc}")
        if not xdatcar.is_file():
            errors.append(f"missing XDATCAR: {row['xdatcar_path']}")
        else:
            try:
                trajectory = inspect_xdatcar(xdatcar, start, end)
                row.update(
                    {
                        "xdatcar_sha256": sha256_file(xdatcar),
                        "frame_count": trajectory["frame_count"],
                        "first_configuration": trajectory["first_configuration"],
                        "last_configuration": trajectory["last_configuration"],
                        "requested_frame_count": trajectory["requested_frame_count"],
                        "requested_range_complete": trajectory["requested_range_complete"],
                    }
                )
                if reference is not None and (
                    trajectory["species"], trajectory["counts"]
                ) != (list(reference.species), list(reference.counts)):
                    errors.append("XDATCAR header species/count order differs from reference POSCAR")
                if not trajectory["requested_range_complete"]:
                    errors.append(
                        f"requested configurations {start}-{end} incomplete; "
                        f"missing ranges: {trajectory['missing_requested_ranges']}"
                    )
            except InventoryError as exc:
                errors.append(f"XDATCAR: {exc}")
        if reference is not None:
            snapshots, coverage = _inspect_bader_case(root, system, ensemble, reference, steps)
            bader_rows.extend(snapshots)
            row["bader_valid_snapshots"] = coverage["valid_count"]
            row["bader_expected_snapshots"] = coverage["expected_count"]
            if not coverage["complete"]:
                errors.append(f"Bader snapshot coverage {coverage['valid_count']}/{coverage['expected_count']}")
        row["messages"] = errors
        row["status"] = "OK" if not errors else "INCOMPLETE"
        case_rows.append(row)

    return {
        "schema_version": "1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "data_root_label": "data_root",
        "requested_frame_range": {"start": start, "end": end},
        "bader_steps": list(steps),
        "case_count": len(case_rows),
        "ok_count": sum(row["status"] == "OK" for row in case_rows),
        "incomplete_count": sum(row["status"] != "OK" for row in case_rows),
        "cases": case_rows,
        "bader_snapshots": bader_rows,
    }


def write_inventory(report: dict[str, object], output_dir: Path, *, overwrite: bool = False) -> list[Path]:
    """Write JSON and two flat CSV reports atomically."""
    output_dir = prepare_output_dir(output_dir, overwrite=overwrite)
    cases = report["cases"]
    bader_snapshots = report["bader_snapshots"]
    assert isinstance(cases, list) and isinstance(bader_snapshots, list)
    outputs = [
        output_dir / "input_inventory.json",
        output_dir / "case_inventory.csv",
        output_dir / "bader_snapshot_coverage.csv",
    ]
    _atomic_write_text(outputs[0], json.dumps(report, indent=2, sort_keys=True) + "\n")
    _atomic_write_csv(
        outputs[1],
        [
            "case_id", "system", "ensemble", "trajectory_variant", "xdatcar_path",
            "reference_path", "xdatcar_sha256", "reference_sha256", "species", "counts",
            "composition", "atom_count", "frame_count", "first_configuration",
            "last_configuration", "requested_start", "requested_end", "requested_frame_count",
            "requested_range_complete", "bader_valid_snapshots", "bader_expected_snapshots",
            "status", "messages",
        ],
        [
            {
                **row,
                "species": json.dumps(row["species"], separators=(",", ":")),
                "counts": json.dumps(row["counts"], separators=(",", ":")),
                "composition": json.dumps(row["composition"], sort_keys=True, separators=(",", ":")),
                "messages": " | ".join(row["messages"]),
            }
            for row in cases
        ],
    )
    _atomic_write_csv(
        outputs[2],
        [
            "system", "ensemble", "step", "snapshot_path", "poscar_path", "acf_path",
            "status", "atom_count", "acf_atom_rows", "poscar_sha256", "acf_sha256", "message",
        ],
        bader_snapshots,
    )
    return outputs


def _atomic_write_text(path: Path, content: str) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(content, encoding="utf-8", newline="\n")
    temporary.replace(path)


def _atomic_write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)
