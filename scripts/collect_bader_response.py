#!/usr/bin/env python3
"""Collect paired neutral/2K Bader response tables."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aimd_analysis.cli import main


if __name__ == "__main__":
    raise SystemExit(main(["bader-response", *sys.argv[1:]]))
