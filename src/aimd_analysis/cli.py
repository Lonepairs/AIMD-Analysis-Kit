"""Command-line entry point for AIMD Analysis Kit."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from . import __version__

_PLANNED_COMMANDS = {
    "density": "Analyze slice-based water density.",
    "polar-angle": "Analyze water orientation and rho*cos(phi).",
    "hydrogen-bonds": "Analyze hydrogen bonds and donor-water counts.",
    "rdf": "Analyze radial distribution and coordination numbers.",
    "bader": "Analyze Bader charge and charge-derived potential response.",
    "batch": "Run analyses from a YAML or JSON workflow configuration.",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="aimd-analysis",
        description="Reproducible analysis tools for VASP AIMD trajectories.",
        epilog="Analysis subcommands are being added incrementally.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    for name, help_text in _PLANNED_COMMANDS.items():
        subparsers.add_parser(name, help=help_text)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command is not None:
        parser.error(
            f"the '{args.command}' command is not implemented yet; "
            "the package scaffold is in place for the analysis migration"
        )
    parser.print_help()
    return 0
