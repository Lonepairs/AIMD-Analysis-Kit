#!/usr/bin/env python3
"""Prepare portable Bader snapshot job directories and a manifest."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aimd_analysis.structure import parse_poscar
from aimd_analysis.potcar import PotcarError, validate_potcar_species


SYSTEMS = ("Cu", "CuxPd1", "Pd", "PdxCu1")
ENSEMBLES = ("aimd1", "aimd1_e")
REQUIRED = ("POSCAR", "CHGCAR", "AECCAR0", "AECCAR2", "POTCAR")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="Source data root")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--system", action="append", choices=SYSTEMS)
    parser.add_argument("--ensemble", action="append", choices=ENSEMBLES)
    parser.add_argument("--start-step", type=int, default=5500)
    parser.add_argument("--end-step", type=int, default=9000)
    parser.add_argument("--step-stride", type=int, default=500)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.start_step < 1 or args.end_step < args.start_step or args.step_stride < 1:
        parser.error("invalid Bader step grid")
    root = args.root.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()
    systems = args.system or list(SYSTEMS)
    ensembles = args.ensemble or list(ENSEMBLES)
    steps = list(range(args.start_step, args.end_step + 1, args.step_stride))
    rows = []
    for system in systems:
        for ensemble in ensembles:
            bader_dir = "bader_e" if ensemble == "aimd1_e" else "bader"
            case_reference = root / system / ensemble / "POSCAR"
            neutral_reference = root / system / "aimd1" / "POSCAR"
            for reference_source, reference_target in (
                (case_reference, output / system / ensemble / "POSCAR"),
                (neutral_reference, output / system / "aimd1" / "POSCAR"),
            ):
                if reference_source.is_file() and not reference_target.exists():
                    reference_target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(reference_source, reference_target)
            case_potcar = root / system / ensemble / "POTCAR"
            staged_potcar = output / system / ensemble / "POTCAR"
            if case_potcar.is_file() and (args.force or not staged_potcar.exists()):
                staged_potcar.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(case_potcar, staged_potcar)
            elif not case_potcar.is_file():
                rows.append({
                    "system": system, "ensemble": ensemble, "step": "",
                    "source_path": case_potcar.relative_to(root).as_posix(),
                    "job_path": staged_potcar.relative_to(output).as_posix(),
                    "expected_atom_count": "", "status": "MISSING_SOURCE",
                    "message": "missing POTCAR",
                })
                continue
            try:
                validate_potcar_species(case_potcar, parse_poscar(case_reference).species)
            except (OSError, PotcarError, ValueError) as exc:
                rows.append({
                    "system": system, "ensemble": ensemble, "step": "",
                    "source_path": case_potcar.relative_to(root).as_posix(),
                    "job_path": staged_potcar.relative_to(output).as_posix(),
                    "expected_atom_count": "", "status": "INVALID",
                    "message": f"POTCAR validation: {exc}",
                })
                continue
            for step in steps:
                source = root / system / bader_dir / f"bader_{system}_{step}"
                destination = output / system / bader_dir / f"bader_{system}_{step}"
                existing_acf = destination / "ACF.dat"
                status = "MISSING_SOURCE"
                message = ""
                present = [name for name in REQUIRED if (source / name).is_file()]
                missing = [name for name in REQUIRED if name not in present]
                if missing:
                    message = "missing " + ", ".join(missing)
                elif existing_acf.is_file() and not args.force:
                    status = "EXISTS"
                    message = "completed ACF.dat preserved"
                else:
                    destination.mkdir(parents=True, exist_ok=True)
                    for name in REQUIRED:
                        target = destination / name
                        if args.force or not target.exists():
                            shutil.copy2(source / name, target)
                    status = "PREPARED"
                expected_atom_count = ""
                source_poscar = source / "POSCAR"
                if source_poscar.is_file():
                    try:
                        expected_atom_count = parse_poscar(source_poscar).atom_count
                    except ValueError:
                        expected_atom_count = ""
                rows.append(
                    {
                        "system": system,
                        "ensemble": ensemble,
                        "step": step,
                        "source_path": source.relative_to(root).as_posix(),
                        "job_path": destination.relative_to(output).as_posix(),
                        "expected_atom_count": expected_atom_count,
                        "status": status,
                        "message": message,
                    }
                )
    output.mkdir(parents=True, exist_ok=True)
    with (output / "bader_manifest.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Prepared {sum(row['status'] == 'PREPARED' for row in rows)} jobs; manifest: {output / 'bader_manifest.csv'}")
    return 0 if all(row["status"] in {"PREPARED", "EXISTS"} for row in rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
