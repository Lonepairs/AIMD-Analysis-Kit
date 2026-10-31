"""Streaming XDATCAR parsing."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from .structure import Structure, StructureParseError, parse_xdatcar_header


CONFIGURATION_RE = re.compile(r"(?:direct|cartesian)?\s*configuration\s*=\s*(\d+)", re.I)


@dataclass(frozen=True)
class TrajectoryFrame:
    """One XDATCAR configuration stored as fractional coordinates."""

    configuration: int
    positions_frac: np.ndarray

    def __post_init__(self) -> None:
        positions = np.asarray(self.positions_frac, dtype=float)
        if positions.ndim != 2 or positions.shape[1] != 3 or not np.all(np.isfinite(positions)):
            raise StructureParseError("trajectory frame has invalid coordinates")
        positions.setflags(write=False)
        object.__setattr__(self, "configuration", int(self.configuration))
        object.__setattr__(self, "positions_frac", positions)


def _configuration_marker(line: str, path: Path) -> tuple[int, str]:
    match = CONFIGURATION_RE.search(line)
    if match is None:
        raise StructureParseError(f"{path}: expected configuration marker, got {line.strip()!r}")
    mode = "Cartesian" if line.lstrip().lower().startswith("cartesian") else "Direct"
    return int(match.group(1)), mode


def _coordinate(line: str, path: Path, configuration: int, atom_index: int) -> np.ndarray:
    fields = line.split()
    if len(fields) < 3:
        raise StructureParseError(
            f"{path}: configuration {configuration} atom {atom_index}: expected three coordinates"
        )
    try:
        values = np.asarray([float(value) for value in fields[:3]], dtype=float)
    except ValueError as exc:
        raise StructureParseError(
            f"{path}: configuration {configuration} atom {atom_index}: invalid coordinate"
        ) from exc
    if not np.all(np.isfinite(values)):
        raise StructureParseError(
            f"{path}: configuration {configuration} atom {atom_index}: non-finite coordinate"
        )
    return values


def iter_xdatcar(
    path: str | Path,
    start: int | None = None,
    end: int | None = None,
    stride: int = 1,
) -> Iterable[TrajectoryFrame]:
    """Yield frames in the requested configuration range without retaining the trajectory."""
    if stride < 1:
        raise StructureParseError("trajectory stride must be positive")
    if start is not None and start < 1:
        raise StructureParseError("trajectory start must be positive")
    if end is not None and (end < 1 or (start is not None and end < start)):
        raise StructureParseError("trajectory end must be >= start")
    source = Path(path)
    try:
        with source.open("r", encoding="utf-8", errors="strict") as stream:
            structure = parse_xdatcar_header(stream, str(source))
            previous: int | None = None
            in_range_count = 0
            while True:
                marker = stream.readline()
                if marker == "":
                    break
                if not marker.strip():
                    continue
                configuration, mode = _configuration_marker(marker, source)
                if previous is not None and configuration <= previous:
                    raise StructureParseError(f"{source}: configurations are not strictly increasing")
                previous = configuration
                if end is not None and configuration > end:
                    break
                selected_by_range = start is None or configuration >= start
                selected_by_stride = selected_by_range and (in_range_count % stride == 0)
                if selected_by_range:
                    in_range_count += 1
                rows = []
                for atom_index in range(structure.atom_count):
                    line = stream.readline()
                    if line == "":
                        raise StructureParseError(
                            f"{source}: configuration {configuration} has only {atom_index} coordinate rows"
                        )
                    coordinates = _coordinate(line, source, configuration, atom_index + 1)
                    if selected_by_stride:
                        rows.append(coordinates)
                if not selected_by_stride:
                    continue
                coordinates = np.asarray(rows, dtype=float)
                if mode == "Cartesian":
                    coordinates = (coordinates * structure.scale) @ np.linalg.inv(structure.lattice)
                yield TrajectoryFrame(configuration, coordinates)
    except (OSError, UnicodeError) as exc:
        raise StructureParseError(f"cannot read {source}: {exc}") from exc


def read_xdatcar_header(path: str | Path) -> Structure:
    """Read an XDATCAR header without loading any frames."""
    source = Path(path)
    try:
        with source.open("r", encoding="utf-8", errors="strict") as stream:
            return parse_xdatcar_header(stream, str(source))
    except (OSError, UnicodeError) as exc:
        raise StructureParseError(f"cannot read {source}: {exc}") from exc
