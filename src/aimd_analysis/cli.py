"""Command-line entry point for AIMD Analysis Kit."""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .inventory import (
    DEFAULT_BADER_STEPS,
    DEFAULT_CASES,
    InventoryError,
    inventory_cases,
    write_inventory,
)
from .density import add_density_arguments, analyze_density
from .polar_angle import add_polar_angle_arguments, analyze_polar_angle
from .hydrogen_bonds import add_hydrogen_bond_arguments, analyze_hydrogen_bonds
from .rdf import add_rdf_arguments, analyze_rdf
from .bader import add_bader_arguments, analyze_bader
from .workflow import run_workflow
from .bader_response import add_bader_response_arguments, collect_bader_response

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="aimd-analysis",
        description="Reproducible analysis tools for VASP AIMD trajectories.",
        epilog="Use named commands for scripts and scripts/run_analysis.py for the numbered compatibility interface.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    density = subparsers.add_parser("density", help="Analyze slice-based water density.")
    add_density_arguments(density)
    polar_angle = subparsers.add_parser("polar-angle", help="Analyze water orientation and rho*cos(phi).")
    add_polar_angle_arguments(polar_angle)
    hydrogen_bonds = subparsers.add_parser("hydrogen-bonds", help="Analyze hydrogen bonds and donor-water counts.")
    add_hydrogen_bond_arguments(hydrogen_bonds)
    rdf = subparsers.add_parser("rdf", help="Analyze radial distribution and coordination numbers.")
    add_rdf_arguments(rdf)
    bader = subparsers.add_parser("bader", help="Analyze Bader charge and surface charge density.")
    add_bader_arguments(bader)
    bader_response = subparsers.add_parser("bader-response", help="Collect neutral/2K Bader surface-charge tables.")
    add_bader_response_arguments(bader_response)
    batch = subparsers.add_parser("batch", help="Run analyses from a YAML or JSON workflow configuration.")
    batch.add_argument("--config", type=Path, required=True)
    batch.add_argument("--select", default="all", help="analysis choices 1-5, all, or none")
    batch.add_argument("--status", type=Path, default=None)
    batch.add_argument("--dry-run", action="store_true")
    batch.add_argument("--continue-on-error", action="store_true")

    inventory = subparsers.add_parser(
        "inventory", help="Validate case inputs and write a portable inventory manifest."
    )
    inventory.add_argument("--data-root", type=Path, required=True)
    inventory.add_argument("--output-dir", type=Path, required=True)
    inventory.add_argument("--case", action="append", dest="cases", metavar="SYSTEM/ENSEMBLE")
    inventory.add_argument("--start", type=int, default=5001)
    inventory.add_argument("--end", type=int, default=9000)
    inventory.add_argument("--bader-steps", type=int, nargs="+", default=list(DEFAULT_BADER_STEPS))
    inventory.add_argument("--overwrite", action="store_true", help="Allow reuse of a non-empty output directory")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "inventory":
        try:
            report = inventory_cases(
                data_root=args.data_root,
                cases=args.cases or DEFAULT_CASES,
                start=args.start,
                end=args.end,
                bader_steps=args.bader_steps,
            )
            outputs = write_inventory(report, args.output_dir, overwrite=args.overwrite)
        except (InventoryError, OSError, ValueError) as exc:
            parser.error(str(exc))
        print(
            f"Inventory: {report['ok_count']}/{report['case_count']} cases complete; "
            f"reports written: {', '.join(path.name for path in outputs)}"
        )
        return 0 if report["incomplete_count"] == 0 else 2
    if args.command == "density":
        try:
            metadata = analyze_density(args)
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        print(
            json.dumps(
                {
                    "output_dir": str(Path(args.output_dir)),
                    "frames": metadata["frames_analyzed"],
                    "profile_range_A": metadata["profile_range_A"],
                    "integrated_oxygen_equivalent_count": metadata["integrated_oxygen_equivalent_count"],
                    "warnings": metadata["warnings"],
                },
                indent=2,
            )
        )
        return 0
    if args.command == "polar-angle":
        try:
            metadata = analyze_polar_angle(args)
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        print(
            json.dumps(
                {
                    "output_dir": str(Path(args.output_dir)),
                    "frames": metadata["frames_analyzed"],
                    "valid_water_observations": metadata["valid_water_observations"],
                    "mean_alpha_deg": metadata["mean_alpha_deg"],
                    "layer_analysis_enabled": metadata["layer_analysis_enabled"],
                },
                indent=2,
            )
        )
        return 0
    if args.command == "hydrogen-bonds":
        try:
            metadata = analyze_hydrogen_bonds(args)
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        summary_name = "ndonor_waters_summary.csv" if args.traversal == "donor" else "nacceptor_waters_summary.csv"
        print(
            json.dumps(
                {
                    "output_dir": str(Path(args.output_dir)),
                    "frames": metadata["frames_analyzed"],
                    "regions": len(metadata["regions"]),
                    "summary_csv": str(Path(args.output_dir) / summary_name),
                },
                indent=2,
            )
        )
        return 0
    if args.command == "rdf":
        try:
            metadata = analyze_rdf(args)
        except (ImportError, OSError, ValueError) as exc:
            parser.error(str(exc))
        print(
            json.dumps(
                {
                    "output_dir": str(Path(args.output_dir)),
                    "pair": metadata["pair"],
                    "frames": metadata["frames_analyzed"],
                    "block_count": metadata["block_count"],
                    "origin_output": next((name for name in metadata["output_files"] if name.endswith("_origin.csv")), None),
                },
                indent=2,
            )
        )
        return 0
    if args.command == "bader":
        try:
            metadata = analyze_bader(args)
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        print(
            json.dumps(
                {
                    "output_dir": str(Path(args.output_dir)),
                    "system": metadata["system"],
                    "ensemble": metadata["ensemble"],
                    "snapshots": len(metadata["step_grid"]),
                    "top_layer_atom_count": metadata["top_layer_atom_count"],
                },
                indent=2,
            )
        )
        return 0
    if args.command == "bader-response":
        try:
            metadata = collect_bader_response(args)
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        print(json.dumps({"output_dir": str(Path(args.output_dir)), "systems": len(metadata["systems"])}, indent=2))
        return 0
    if args.command == "batch":
        if args.select.casefold() == "none":
            print("No analyses selected.")
            return 0
        status = args.status or args.config.with_name("batch_status.csv")
        rows = run_workflow(
            args.config,
            main,
            selection=args.select,
            status_path=status,
            dry_run=args.dry_run,
            stop_on_error=not args.continue_on_error,
        )
        for row in rows:
            print(f"{row['job']}: {row['status']}" + (f" ({row['message']})" if row["message"] else ""))
        return 0 if all(row["status"] in {"OK", "PLANNED"} for row in rows) else 2
    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
