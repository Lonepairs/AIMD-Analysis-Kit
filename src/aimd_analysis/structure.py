"""VASP structure parsing and shared lattice geometry."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math
from pathlib import Path
import re
from typing import Sequence

import numpy as np


class StructureParseError(ValueError):
    """Raised when a POSCAR or XDATCAR structure header is malformed."""


def _as_lattice(rows: Sequence[Sequence[float]], scale: float) -> np.ndarray:
    raw = np.asarray(rows, dtype=float)
    if raw.shape != (3, 3) or not np.all(np.isfinite(raw)):
        raise StructureParseError("lattice must contain three finite 3-vectors")
    raw_volume = float(np.linalg.det(raw))
    if abs(raw_volume) < 1e-12:
        raise StructureParseError("lattice vectors are singular")
    if not np.isfinite(scale) or scale == 0:
        raise StructureParseError("POSCAR scale must be finite and nonzero")
    effective_scale = scale if scale > 0 else (abs(scale) / abs(raw_volume)) ** (1.0 / 3.0)
    lattice = raw * effective_scale
    if np.linalg.det(lattice) <= 0:
        raise StructureParseError("lattice determinant must be positive")
    return lattice


def _validate_species_counts(species: Sequence[str], counts: Sequence[int]) -> tuple[tuple[str, ...], tuple[int, ...]]:
    if not species or len(species) != len(counts):
        raise StructureParseError("species and atom-count fields must have equal nonzero length")
    if any(not isinstance(value, (int, np.integer)) or int(value) <= 0 for value in counts):
        raise StructureParseError("atom counts must be positive integers")
    return tuple(str(value) for value in species), tuple(int(value) for value in counts)


def _parse_numeric_row(line: str, label: str) -> np.ndarray:
    fields = line.split()
    if len(fields) < 3:
        raise StructureParseError(f"{label}: expected at least three coordinate values")
    try:
        values = np.asarray([float(value) for value in fields[:3]], dtype=float)
    except ValueError as exc:
        raise StructureParseError(f"{label}: non-numeric coordinates") from exc
    if not np.all(np.isfinite(values)):
        raise StructureParseError(f"{label}: non-finite coordinates")
    return values


def _parse_species_counts(lines: Sequence[str], path_label: str) -> tuple[tuple[str, ...], tuple[int, ...], int]:
    if len(lines) < 7:
        raise StructureParseError(f"{path_label}: incomplete VASP header")
    try:
        counts = [int(value) for value in lines[6].split()]
    except ValueError as exc:
        raise StructureParseError(
            f"{path_label}: VASP 4 count-only headers are unsupported; species are required"
        ) from exc
    species, counts_tuple = _validate_species_counts(lines[5].split(), counts)
    return species, counts_tuple, sum(counts_tuple)


def _parse_header(lines: Sequence[str], path_label: str) -> tuple[str, float, np.ndarray, tuple[str, ...], tuple[int, ...], int]:
    if len(lines) < 7:
        raise StructureParseError(f"{path_label}: incomplete VASP header")
    try:
        scale = float(lines[1].split()[0])
        raw_lattice = [[float(value) for value in lines[index].split()[:3]] for index in range(2, 5)]
    except (IndexError, ValueError) as exc:
        raise StructureParseError(f"{path_label}: malformed scale or lattice") from exc
    lattice = _as_lattice(raw_lattice, scale)
    raw_volume = abs(float(np.linalg.det(np.asarray(raw_lattice, dtype=float))))
    effective_scale = scale if scale > 0 else (abs(scale) / raw_volume) ** (1.0 / 3.0)
    species, counts, atom_count = _parse_species_counts(lines, path_label)
    return lines[0].rstrip("\r\n"), effective_scale, lattice, species, counts, atom_count


def _coordinate_start(lines: Sequence[str], atom_count: int, path_label: str) -> tuple[str, int]:
    cursor = 7
    if cursor < len(lines) and lines[cursor].strip().lower().startswith("s"):
        cursor += 1
    if cursor >= len(lines):
        raise StructureParseError(f"{path_label}: missing coordinate mode")
    mode_text = lines[cursor].strip().lower()
    if mode_text.startswith("d"):
        return "Direct", cursor + 1
    if mode_text.startswith(("c", "k")):
        return "Cartesian", cursor + 1
    raise StructureParseError(f"{path_label}: unknown coordinate mode {lines[cursor]!r}")


@dataclass(frozen=True)
class Structure:
    """A VASP structure with fractional coordinates and row-vector lattice."""

    comment: str
    scale: float
    lattice: np.ndarray
    species: tuple[str, ...]
    counts: tuple[int, ...]
    positions_frac: np.ndarray
    coordinate_mode: str = "Direct"

    def __post_init__(self) -> None:
        lattice = np.asarray(self.lattice, dtype=float)
        positions = np.asarray(self.positions_frac, dtype=float)
        if lattice.shape != (3, 3) or positions.ndim != 2 or positions.shape[1] != 3:
            raise StructureParseError("invalid lattice or coordinate array shape")
        if positions.shape[0] != self.atom_count:
            raise StructureParseError("coordinate count does not match species counts")
        if not np.all(np.isfinite(lattice)) or not np.all(np.isfinite(positions)):
            raise StructureParseError("structure contains non-finite values")
        lattice.setflags(write=False)
        positions.setflags(write=False)
        object.__setattr__(self, "lattice", lattice)
        object.__setattr__(self, "positions_frac", positions)

    @property
    def atom_count(self) -> int:
        return sum(self.counts)

    @property
    def volume(self) -> float:
        return float(np.linalg.det(self.lattice))

    @property
    def atom_species(self) -> tuple[str, ...]:
        return tuple(element for element, count in zip(self.species, self.counts) for _ in range(count))

    @property
    def positions_cart(self) -> np.ndarray:
        return self.positions_frac @ self.lattice

    def indices(self, element: str) -> np.ndarray:
        return np.asarray([index for index, value in enumerate(self.atom_species) if value == element], dtype=int)

    def with_positions(self, positions_frac: np.ndarray) -> "Structure":
        return Structure(
            self.comment,
            self.scale,
            self.lattice,
            self.species,
            self.counts,
            np.asarray(positions_frac, dtype=float),
            "Direct",
        )


def parse_poscar(path: str | Path) -> Structure:
    """Parse a VASP 5/6 POSCAR or CONTCAR into a :class:`Structure`."""
    source = Path(path)
    try:
        lines = source.read_text(encoding="utf-8", errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise StructureParseError(f"cannot read {source}: {exc}") from exc
    comment, scale, lattice, species, counts, atom_count = _parse_header(lines, str(source))
    mode, cursor = _coordinate_start(lines, atom_count, str(source))
    coordinates = []
    for index in range(atom_count):
        if cursor + index >= len(lines):
            raise StructureParseError(f"{source}: expected {atom_count} coordinate rows")
        coordinates.append(_parse_numeric_row(lines[cursor + index], f"{source} atom {index + 1}"))
    coordinates_array = np.asarray(coordinates, dtype=float)
    if mode == "Cartesian":
        coordinates_array = coordinates_array * scale
        coordinates_array = coordinates_array @ np.linalg.inv(lattice)
    return Structure(comment, scale, lattice, species, counts, coordinates_array, mode)


def parse_xdatcar_header(stream, path_label: str = "XDATCAR") -> Structure:
    """Read the seven-line XDATCAR header and return an empty-coordinate structure."""
    lines = [stream.readline() for _ in range(7)]
    if any(line == "" for line in lines):
        raise StructureParseError(f"{path_label}: incomplete XDATCAR header")
    comment, scale, lattice, species, counts, atom_count = _parse_header(
        [line.rstrip("\r\n") for line in lines], path_label
    )
    return Structure(comment, scale, lattice, species, counts, np.zeros((atom_count, 3)), "Direct")


def minimum_image_fractional(delta_frac: np.ndarray) -> np.ndarray:
    """Wrap fractional displacement components to ``[-0.5, 0.5)``."""
    delta = np.asarray(delta_frac, dtype=float)
    return delta - np.floor(delta + 0.5)


def fractional_to_cartesian(lattice: np.ndarray, vector_frac: np.ndarray) -> np.ndarray:
    """Convert one or more row-vector fractional coordinates to Cartesian."""
    return np.asarray(vector_frac, dtype=float) @ np.asarray(lattice, dtype=float)


def minimum_image_vector(lattice: np.ndarray, first_frac: np.ndarray, second_frac: np.ndarray) -> np.ndarray:
    """Return the minimum-image Cartesian vector from first to second."""
    vectors = np.asarray(lattice, dtype=float)
    delta = np.asarray(second_frac, dtype=float) - np.asarray(first_frac, dtype=float)
    if vectors.shape != (3, 3) or delta.shape != (3,):
        raise ValueError("minimum_image_vector expects a 3x3 lattice and two 3-vectors")
    center = np.rint(delta).astype(int)
    best = (delta - center) @ vectors
    best_norm_sq = float(np.dot(best, best))
    sigma_min = float(np.linalg.svd(vectors, compute_uv=False)[-1])
    if sigma_min <= 0:
        raise StructureParseError("lattice is singular")
    radius = math.sqrt(best_norm_sq) / sigma_min
    bounds = [range(math.floor(value - radius) - 1, math.ceil(value + radius) + 2) for value in delta]
    for image in product(*bounds):
        candidate = (delta - np.asarray(image, dtype=float)) @ vectors
        norm_sq = float(np.dot(candidate, candidate))
        if norm_sq < best_norm_sq:
            best = candidate
            best_norm_sq = norm_sq
    return best


def minimum_image_distance(lattice: np.ndarray, first_frac: np.ndarray, second_frac: np.ndarray) -> float:
    """Return the minimum-image distance between two fractional coordinates."""
    return float(np.linalg.norm(minimum_image_vector(lattice, first_frac, second_frac)))


def surface_normal(lattice: np.ndarray) -> np.ndarray:
    """Return the unit normal to lattice vectors a and b oriented with c."""
    vectors = np.asarray(lattice, dtype=float)
    normal = np.cross(vectors[0], vectors[1])
    length = np.linalg.norm(normal)
    if length == 0:
        raise StructureParseError("surface normal is undefined for collinear in-plane vectors")
    normal = normal / length
    if np.dot(normal, vectors[2]) < 0:
        normal = -normal
    return normal
