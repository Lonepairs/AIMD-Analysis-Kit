#!/usr/bin/env python3
"""Run configured analysis jobs with an auditable status CSV."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aimd_analysis.cli import main as cli_main
from aimd_analysis.workflow import interactive_selection, run_workflow


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--select", default=None, help="analysis choices 1-5, all, or none")
    parser.add_argument("--status", type=Path, default=None, help="Status CSV; defaults beside the config")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--continue-on-error", action="store_true")
    args = parser.parse_args()
    selection = args.select if args.select is not None else interactive_selection()
    if selection.casefold() == "none":
        print("No analyses selected.")
        return 0
    status = args.status or args.config.with_name("batch_status.csv")
    rows = run_workflow(
        args.config,
        cli_main,
        selection=selection,
        status_path=status,
        dry_run=args.dry_run,
        stop_on_error=not args.continue_on_error,
    )
    for row in rows:
        print(f"{row['job']}: {row['status']}" + (f" ({row['message']})" if row["message"] else ""))
    return 0 if all(row["status"] in {"OK", "PLANNED"} for row in rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
