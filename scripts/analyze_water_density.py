#!/usr/bin/env python3
"""Compatibility wrapper for ``aimd-analysis density``."""

from __future__ import annotations

import sys
from pathlib import Path

# Allow this compatibility wrapper to run from a source checkout before the
# package has been installed; installed users should use ``aimd-analysis``.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aimd_analysis.cli import main


if __name__ == "__main__":
    raise SystemExit(main(["density", *sys.argv[1:]]))
