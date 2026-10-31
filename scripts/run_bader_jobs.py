#!/usr/bin/env python3
"""Select and optionally execute prepared Bader jobs from a manifest."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def _selection(raw: str, count: int) -> list[int]:
    value = raw.strip().casefold()
    if value == "all":
        return list(range(count))
    if value == "none":
        return []
    result = []
    for token in value.replace(",", " ").split():
        index = int(token) - 1
        if index < 0 or index >= count:
            raise ValueError(f"selection {index + 1} is outside 1-{count}")
        if index not in result:
            result.append(index)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--select", default=None, help="number, list, all, or none")
    parser.add_argument("--command", nargs="+", help="Optional command template; {job_dir} is replaced")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--yes", action="store_true", help="Skip interactive confirmation")
    parser.add_argument("--stop-on-error", action="store_true")
    args = parser.parse_args()
    manifest = args.manifest.expanduser()
    with manifest.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        parser.error("manifest has no jobs")
    raw = args.select if args.select is not None else input("Select jobs (number/list/all/none): ")
    try:
        selected = _selection(raw, len(rows))
    except ValueError as exc:
        parser.error(str(exc))
    if not selected:
        print("No Bader jobs selected.")
        return 0
    print("Selected jobs:")
    for index in selected:
        print(f"{index + 1}: {rows[index]['system']}/{rows[index]['ensemble']} step {rows[index]['step']} [{rows[index]['status']}]")
    if not args.yes and not args.dry_run:
        if input("Submit selected jobs? [y/N] ").strip().casefold() not in {"y", "yes"}:
            print("Cancelled.")
            return 0
    failures = 0
    for index in selected:
        row = rows[index]
        job_dir = manifest.parent / row["job_path"]
        if row["status"] == "EXISTS":
            print(f"SKIP {index + 1}: completed output preserved")
            continue
        if row["status"] != "PREPARED":
            print(f"FAILED {index + 1}: manifest status is {row['status']}")
            failures += 1
            if args.stop_on_error:
                break
            continue
        if args.dry_run or not args.command:
            print(f"DRY-RUN {index + 1}: {job_dir}")
            continue
        command = [part.replace("{job_dir}", str(job_dir)) for part in args.command]
        result = subprocess.run(command, cwd=job_dir, check=False)
        if result.returncode:
            failures += 1
            if args.stop_on_error:
                break
    return 2 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
