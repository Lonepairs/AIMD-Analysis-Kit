#!/usr/bin/env python3
"""Numbered compatibility interface for the five analysis commands."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aimd_analysis.cli import main


CHOICES = {"1": "density", "2": "polar-angle", "3": "hydrogen-bonds", "4": "rdf", "5": "bader"}


def run() -> int:
    if len(sys.argv) == 2 and sys.argv[1] in {"-h", "--help"}:
        print("usage: python scripts/run_analysis.py 1|2|3|4|5 [analysis options]")
        print("1  water density")
        print("2  polar angle and rho*cos(phi)")
        print("3  hydrogen bonds and donor counts")
        print("4  RDF and coordination number")
        print("5  Bader charge and potential response")
        return 0
    if len(sys.argv) < 2 or sys.argv[1] not in CHOICES:
        print("usage: python scripts/run_analysis.py 1|2|3|4|5 [analysis options]", file=sys.stderr)
        return 2
    return main([CHOICES[sys.argv[1]], *sys.argv[2:]])


if __name__ == "__main__":
    raise SystemExit(run())
