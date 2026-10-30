# AIMD Analysis Kit

Toolkit for reproducible analysis of ab initio molecular dynamics trajectories.

## Status

Package scaffold v0.1.0. The scientific analysis commands are being migrated
from the research scripts and are not available yet.

## Installation

Python 3.10 or newer is required. Install the core trajectory-analysis
dependencies with:

```bash
python -m pip install .
```

Install plotting support with `python -m pip install '.[plot]'`. Development
tools are available with `python -m pip install '.[dev]'`.

The `aimd-analysis --help` command lists the planned analysis interfaces. The
analysis commands will become available as their scientific modules are
migrated and validated.

## Skills

The repository includes a portable `document-ingestion` Skill for routing local
PDF, image, DOCX, PPTX, HTML, text, and table files to an appropriate parser.
It uses local MinerU when available, preserves a machine-readable result
manifest, and validates generated artifacts with a standard-library Python
script. The Skill is under `skills/document-ingestion/` and follows the open
Agent Skills directory layout.
