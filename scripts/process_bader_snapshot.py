#!/usr/bin/env python3
"""Run an approved Bader command in one prepared snapshot directory."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

_PATCH_PROBE = None


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aimd_analysis.bader import (
    ELECTRON_COUNT_TOLERANCE,
    parse_acf,
    parse_nelect,
    validate_density_grids,
)
from aimd_analysis.structure import parse_poscar


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job-dir", type=Path, required=True)
    parser.add_argument("--command", nargs="+", required=True, help="Executable and arguments; no shell expansion is used")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    job_dir = args.job_dir.expanduser().resolve()
    acf = job_dir / "ACF.dat"

    def validate_completed_output() -> tuple[bool, str | None]:
        try:
            parsed_acf = parse_acf(acf)
            structure = parse_poscar(job_dir / "POSCAR")
        except (OSError, ValueError) as exc:
            return False, f"invalid Bader output ({exc})"
        if len(parsed_acf) != structure.atom_count:
            return False, "ACF/POSCAR atom count mismatch"
        density_paths = [
            job_dir / name
            for name in ("CHGCAR", "CHGCAR_sum", "AECCAR0", "AECCAR2")
            if (job_dir / name).is_file()
        ]
        if len(density_paths) >= 2:
            try:
                validate_density_grids(density_paths)
            except ValueError as exc:
                return False, str(exc)
        outcar = job_dir / "OUTCAR"
        if outcar.is_file():
            try:
                expected = parse_nelect(outcar)
            except ValueError as exc:
                return False, str(exc)
            total = sum(float(row["electron_population"]) for row in parsed_acf)
            residual = total - expected
            if abs(residual) > ELECTRON_COUNT_TOLERANCE:
                return False, (
                    f"ACF/NELECT residual {residual:.6g} e exceeds "
                    f"{ELECTRON_COUNT_TOLERANCE:.6g} e"
                )
        return True, None

    if acf.is_file() and not args.force:
        valid, error = validate_completed_output()
        if valid:
            print(f"SKIP {job_dir}: validated ACF.dat already exists")
            return 0
        print(f"FAILED {job_dir}: {error}; rerun with --force")
        return 2
    if not job_dir.is_dir():
        parser.error(f"job directory not found: {job_dir}")
    log = job_dir / "bader_command.log"
    with log.open("w", encoding="utf-8") as stream:
        stream.write("$ " + " ".join(args.command) + "\n")
        result = subprocess.run(args.command, cwd=job_dir, stdout=stream, stderr=subprocess.STDOUT, check=False)
        stream.write(f"\nexit_code={result.returncode}\n")
    if result.returncode:
        print(f"FAILED {job_dir}; see {log}")
        return result.returncode
    if not acf.is_file():
        print(f"FAILED {job_dir}: command returned zero but did not create ACF.dat")
        return 2
    valid, error = validate_completed_output()
    if not valid:
        print(f"FAILED {job_dir}: {error}")
        return 2
    print(f"OK {job_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
