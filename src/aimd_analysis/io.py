"""Shared filesystem guards for generated analysis outputs."""

from __future__ import annotations

from pathlib import Path


def prepare_output_dir(path: Path, *, overwrite: bool = False) -> Path:
    """Create an output directory, refusing reuse of non-empty directories by default."""
    output_dir = Path(path).expanduser()
    if output_dir.exists() and not output_dir.is_dir():
        raise FileExistsError(f"output path is not a directory: {output_dir}")
    if output_dir.is_dir() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(f"output directory is not empty: {output_dir}; choose a new path or pass --overwrite")
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir
