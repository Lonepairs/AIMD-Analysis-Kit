#!/usr/bin/env python3
"""Compatibility wrapper for ``aimd-analysis bader``."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aimd_analysis.cli import main


if __name__ == "__main__":
    raise SystemExit(main(["bader", *sys.argv[1:]]))
