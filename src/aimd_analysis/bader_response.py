"""Collect neutral/2K Bader surface summaries into response tables."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
from pathlib import Path
from typing import Iterable

from . import __version__
from .bader import SYSTEMS
from .density import _git_revision, _portable_input_identifier, _sha256
from .io import prepare_output_dir


def _optional_float(value: object) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _write_csv(path: Path, header: list[str], rows: Iterable[Iterable[object]]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    temporary.replace(path)


def _parse_case_output(value: str) -> tuple[str, str, Path]:
    if "=" not in value:
        raise ValueError("case output must be SYSTEM/ENSEMBLE=PATH")
    case, raw_path = value.split("=", 1)
    parts = case.replace("\\", "/").split("/")
    if len(parts) != 2 or parts[1] not in {"aimd1", "aimd1_e"} or parts[0] not in SYSTEMS:
        raise ValueError(f"invalid case output identifier: {case!r}")
    return parts[0], parts[1], Path(raw_path).expanduser()


def _read_summary(path: Path) -> dict[str, str]:
    source = path / "bader_surface_metal_charge_summary.csv"
    if not source.is_file():
        raise FileNotFoundError(f"Bader surface summary not found: {source}")
    with source.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError(f"expected one summary row in {source}")
    return rows[0]


def _read_snapshots(path: Path) -> list[dict[str, str]]:
    source = path / "bader_surface_metal_charge.csv"
    if not source.is_file():
        raise FileNotFoundError(f"Bader surface snapshot table not found: {source}")
    with source.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return rows


def _sem(values: list[float]) -> float | None:
    if len(values) < 2:
        return None
    mean = sum(values) / len(values)
    return math.sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - 1)) / math.sqrt(len(values))


def collect_bader_response(args: argparse.Namespace) -> dict:
    case_paths = [_parse_case_output(value) for value in args.case_output]
    by_system: dict[str, dict[str, tuple[Path, dict[str, str], list[dict[str, str]]]]] = {}
    snapshot_rows: list[tuple[object, ...]] = []
    for system, ensemble, path in case_paths:
        if ensemble in by_system.get(system, {}):
            raise ValueError(f"duplicate Bader output for {system}/{ensemble}")
        summary = _read_summary(path)
        snapshots = _read_snapshots(path)
        if not snapshots:
            raise ValueError(f"no Bader snapshots found for {system}/{ensemble}")
        if summary.get("system") != system or summary.get("ensemble") != ensemble:
            raise ValueError(f"Bader summary identity does not match {system}/{ensemble}")
        if int(summary["n_snapshots"]) != len(snapshots):
            raise ValueError(f"Bader summary snapshot count does not match snapshot rows for {system}/{ensemble}")
        if any(row.get("system") != system or row.get("ensemble") != ensemble for row in snapshots):
            raise ValueError(f"Bader snapshot identity does not match {system}/{ensemble}")
        by_system.setdefault(system, {})[ensemble] = (path, summary, snapshots)
        for row in snapshots:
            snapshot_rows.append(
                (
                    system, ensemble, int(row["step"]), float(row["net_surface_metal_charge_q"]),
                    float(row["sigma_uC_cm-2"]), _optional_float(row.get("delta_phi_V")), _optional_float(row.get("phi_V_vs_SHE")),
                    int(row["n_surface_metal_atoms"]),
                )
            )
    missing_pairs = [system for system, cases in by_system.items() if "aimd1" not in cases or "aimd1_e" not in cases]
    if missing_pairs:
        raise ValueError("both aimd1 and aimd1_e outputs are required for: " + ", ".join(sorted(missing_pairs)))
    surface_summary_rows = []
    response_rows = []
    for system in sorted(by_system):
        _, neutral, neutral_snapshots = by_system[system]["aimd1"]
        _, electro, electro_snapshots = by_system[system]["aimd1_e"]
        neutral_by_step = {int(row["step"]): row for row in neutral_snapshots}
        electro_by_step = {int(row["step"]): row for row in electro_snapshots}
        if set(neutral_by_step) != set(electro_by_step):
            raise ValueError(f"neutral and aimd1_e snapshot grids differ for {system}")
        if len(neutral_by_step) != len(neutral_snapshots) or len(electro_by_step) != len(electro_snapshots):
            raise ValueError(f"duplicate Bader snapshot steps for {system}")
        if int(neutral["n_surface_metal_atoms"]) != int(electro["n_surface_metal_atoms"]):
            raise ValueError(f"neutral and aimd1_e surface group sizes differ for {system}")
        neutral_q = float(neutral["mean_net_surface_metal_charge_q"])
        electro_q = float(electro["mean_net_surface_metal_charge_q"])
        neutral_sigma = float(neutral["mean_sigma_uC_cm-2"])
        electro_sigma = float(electro["mean_sigma_uC_cm-2"])
        neutral_phi = _optional_float(neutral.get("mean_phi_V_vs_SHE"))
        electro_phi = _optional_float(electro.get("mean_phi_V_vs_SHE"))
        neutral_delta = _optional_float(neutral.get("mean_delta_phi_V"))
        electro_delta = _optional_float(electro.get("mean_delta_phi_V"))
        neutral_phi_values = [value for value in (_optional_float(row.get("phi_V_vs_SHE")) for row in neutral_snapshots) if value is not None]
        electro_phi_values = [value for value in (_optional_float(row.get("phi_V_vs_SHE")) for row in electro_snapshots) if value is not None]
        neutral_phi_sem = _sem(neutral_phi_values) if len(neutral_phi_values) == len(neutral_snapshots) else None
        electro_phi_sem = _sem(electro_phi_values) if len(electro_phi_values) == len(electro_snapshots) else None
        response = electro_phi - neutral_phi if neutral_phi is not None and electro_phi is not None else None
        response_sem = math.hypot(neutral_phi_sem, electro_phi_sem) if neutral_phi_sem is not None and electro_phi_sem is not None else None
        response_rows.append(
            (
                system, neutral_q, electro_q, neutral_sigma, electro_sigma, neutral_delta, electro_delta,
                neutral_phi, electro_phi, response, "" if response_sem is None else response_sem,
            )
        )
        for ensemble, summary in (("aimd1", neutral), ("aimd1_e", electro)):
            snapshot_sem_q = summary["snapshot_sem_net_surface_metal_charge_q"]
            surface_summary_rows.append(
                (
                    system, ensemble, int(summary["n_snapshots"]), int(summary["n_surface_metal_atoms"]),
                    float(summary["mean_net_surface_metal_charge_q"]), _optional_float(summary.get("snapshot_std_net_surface_metal_charge_q")),
                    "" if snapshot_sem_q == "" else float(snapshot_sem_q), float(summary["mean_sigma_uC_cm-2"]),
                    _optional_float(summary.get("mean_delta_phi_V")), _optional_float(summary.get("mean_phi_V_vs_SHE")),
                )
            )
    output_dir = prepare_output_dir(Path(args.output_dir).expanduser(), overwrite=getattr(args, "overwrite", False))
    _write_csv(
        output_dir / "bader_snapshot_potential_response.csv",
        ["system", "ensemble", "step", "net_surface_metal_charge_q", "sigma_uC_cm-2", "delta_phi_V", "phi_V_vs_SHE", "n_surface_metal_atoms"],
        sorted(snapshot_rows, key=lambda row: (row[0], row[1], row[2])),
    )
    _write_csv(
        output_dir / "bader_surface_charge_summary.csv",
        ["system", "ensemble", "n_snapshots", "n_surface_metal_atoms", "mean_net_surface_metal_charge_q", "snapshot_std_q", "snapshot_sem_q", "mean_sigma_uC_cm-2", "mean_delta_phi_V", "mean_phi_V_vs_SHE"],
        surface_summary_rows,
    )
    _write_csv(
        output_dir / "bader_potential_response_summary.csv",
        ["system", "Q_neutral_e", "Q_electro_e", "sigma_neutral_uC_cm-2", "sigma_electro_uC_cm-2", "delta_phi_neutral_V", "delta_phi_electro_V", "phi_neutral_V", "phi_electro_V", "potential_response_V", "potential_response_independent_sem_V"],
        response_rows,
    )
    metadata = {
        "schema_version": 1,
        "analysis": "Bader neutral/2K surface-charge response collection",
        "package_version": __version__,
        "git_revision": _git_revision(),
        "runtime": {"python": platform.python_version(), "operating_system": platform.system()},
        "equations": {"sigma": "Q_surface / A * 1602.176634", "potential_response": "optional; requires supplied phi_V_vs_SHE fields"},
        "uncertainty": "Independent-sample SEM propagation applies only when both input outputs contain finite potential fields; otherwise potential response fields remain blank.",
        "systems": sorted(by_system),
        "potential_fields": "optional; PZC and capacitance are not supplied by this kit",
        "input_case_outputs": [
            {"system": system, "ensemble": ensemble, "path": _portable_input_identifier(path, args.data_root)[0], "absolute_path_redacted": _portable_input_identifier(path, args.data_root)[1], "summary_sha256": _sha256(path / "bader_surface_metal_charge_summary.csv")}
            for system, ensemble, path in case_paths
        ],
        "output_files": ["bader_snapshot_potential_response.csv", "bader_surface_charge_summary.csv", "bader_potential_response_summary.csv", "analysis_metadata.json"],
    }
    metadata_path = output_dir / "analysis_metadata.json"
    temporary_metadata = metadata_path.with_name(metadata_path.name + ".tmp")
    temporary_metadata.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    temporary_metadata.replace(metadata_path)
    return metadata


def add_bader_response_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--case-output", action="append", required=True, metavar="SYSTEM/ENSEMBLE=PATH", help="Per-case Bader output directory; repeat for neutral and aimd1_e")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true", help="Allow reuse of a non-empty output directory")
