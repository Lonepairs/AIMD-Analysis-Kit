#!/usr/bin/env python3
"""Run an approved Bader command in one prepared snapshot directory."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aimd_analysis.bader import parse_acf
from aimd_analysis.structure import parse_poscar


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job-dir", type=Path, required=True)
    parser.add_argument("--command", nargs="+", required=True, help="Executable and arguments; no shell expansion is used")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    job_dir = args.job_dir.expanduser().resolve()
    acf = job_dir / "ACF.dat"
    if acf.is_file() and not args.force:
        try:
            parsed_acf = parse_acf(acf)
            structure = parse_poscar(job_dir / "POSCAR")
            if len(parsed_acf) == structure.atom_count:
                print(f"SKIP {job_dir}: validated ACF.dat already exists")
                return 0
            print(f"FAILED {job_dir}: existing ACF.dat does not match POSCAR; rerun with --force")
            return 2
        except (OSError, ValueError) as exc:
            print(f"FAILED {job_dir}: existing ACF.dat is invalid ({exc}); rerun with --force")
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
    try:
        parsed_acf = parse_acf(acf)
        structure = parse_poscar(job_dir / "POSCAR")
    except (OSError, ValueError) as exc:
        print(f"FAILED {job_dir}: invalid Bader output: {exc}")
        return 2
    if len(parsed_acf) != structure.atom_count:
        print(f"FAILED {job_dir}: ACF/POSCAR atom count mismatch")
        return 2
    print(f"OK {job_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
