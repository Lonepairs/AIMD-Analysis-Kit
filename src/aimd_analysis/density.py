"""Slice-based water density analysis for VASP XDATCAR trajectories."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

from . import __version__
from .io import prepare_output_dir
from .structure import Structure, parse_poscar, surface_normal
from .trajectory import TrajectoryFrame, iter_xdatcar, read_xdatcar_header


AVOGADRO = 6.02214076e23
OXYGEN_MOLAR_MASS_G_MOL = 15.999
WATER_MOLAR_MASS_G_MOL = 18.015
NON_SURFACE_SPECIES = frozenset({"o", "h", "he", "k"})
CONSERVATION_TOLERANCE = 1.0e-8


def cell_geometry(lattice: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Return outward surface normal, normal cell period, and lateral area."""
    vectors = np.asarray(lattice, dtype=float)
    normal = surface_normal(vectors)
    period = float(np.dot(vectors[2], normal))
    area = float(np.linalg.norm(np.cross(vectors[0], vectors[1])))
    if period <= 0 or area <= 0:
        raise ValueError("cell period and lateral area must be positive")
    return normal, period, area


def select_surface_indices(
    labels: tuple[str, ...] | list[str], surface_species: str | None = None
) -> tuple[np.ndarray, list[str]]:
    """Select substrate atoms; by default exclude O, H, He, and K."""
    if surface_species is None:
        selected = {label.casefold() for label in labels if label.casefold() not in NON_SURFACE_SPECIES}
    else:
        selected = {item.strip().casefold() for item in surface_species.split(",") if item.strip()}
        if not selected:
            raise ValueError("--surface-species must contain at least one species")

    available = {label.casefold() for label in labels}
    unknown = sorted(selected - available)
    if unknown:
        raise ValueError("surface species not present in XDATCAR: " + ", ".join(unknown))
    indices = np.asarray(
        [index for index, label in enumerate(labels) if label.casefold() in selected], dtype=int
    )
    if indices.size == 0:
        raise ValueError("no atoms match the selected surface species")
    names = [label for label in dict.fromkeys(labels) if label.casefold() in selected]
    return indices, names


def surface_relative_heights(
    frame: TrajectoryFrame,
    structure: Structure,
    normal: np.ndarray,
    period: float,
    surface_indices: np.ndarray,
    oxygen_indices: np.ndarray,
    top_layer_tolerance: float,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Return instantaneous surface height, wrapped O heights, and raw heights."""
    projected = frame.positions_frac @ structure.lattice @ normal
    surface_heights = projected[surface_indices]
    top = surface_heights >= float(np.max(surface_heights)) - top_layer_tolerance
    if not np.any(top):
        raise ValueError("no substrate atoms satisfy the top-layer surface criterion")
    surface_z = float(np.mean(surface_heights[top]))
    raw_heights = projected[oxygen_indices] - surface_z
    return surface_z, np.mod(raw_heights, period), raw_heights


def slice_intersections(oxygen_heights: np.ndarray, edges: np.ndarray, radius: float) -> np.ndarray:
    """Return an oxygen-by-slice mask for nonzero ball/slice overlap."""
    centers = np.asarray(oxygen_heights, dtype=float)[:, None]
    edges = np.asarray(edges, dtype=float)
    tolerance = max(1.0e-12, np.finfo(float).eps * max(1.0, float(np.max(np.abs(edges)))) * 16)
    return (centers + radius > edges[:-1][None, :] + tolerance) & (
        centers - radius < edges[1:][None, :] - tolerance
    )


def normalized_slice_contributions(
    oxygen_heights: np.ndarray, edges: np.ndarray, radius: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return raw hits, normalized oxygen equivalents per slice, and denominators."""
    intersections = slice_intersections(oxygen_heights, edges, radius)
    denominators = intersections.sum(axis=1).astype(float)
    contributions = np.divide(
        intersections,
        denominators[:, None],
        out=np.zeros_like(intersections, dtype=float),
        where=denominators[:, None] > 0,
    )
    return intersections.sum(axis=0), contributions.sum(axis=0), denominators


def rebin_slice_values(
    values: np.ndarray, source_edges: np.ndarray, target_width: float
) -> tuple[np.ndarray, np.ndarray]:
    """Rebin slice totals for plotting by uniform overlap, preserving their sum."""
    values = np.asarray(values, dtype=float)
    source_edges = np.asarray(source_edges, dtype=float)
    if values.ndim != 1 or source_edges.ndim != 1 or len(source_edges) != len(values) + 1:
        raise ValueError("values and source_edges have incompatible dimensions")
    if target_width <= 0 or np.any(np.diff(source_edges) <= 0):
        raise ValueError("source edges must increase and target width must be positive")

    count = int(math.ceil(source_edges[-1] / target_width - 1.0e-12))
    target_edges = np.minimum(np.arange(count + 1, dtype=float) * target_width, source_edges[-1])
    target_edges[-1] = source_edges[-1]
    rebinned = np.zeros(count, dtype=float)
    for index, value in enumerate(values):
        lower, upper = source_edges[index : index + 2]
        first = max(0, int(math.floor(lower / target_width)))
        last = min(count - 1, int(math.ceil(upper / target_width - 1.0e-12) - 1))
        if last < first:
            continue
        target_indices = np.arange(first, last + 1)
        overlaps = np.maximum(
            0.0,
            np.minimum(upper, target_edges[target_indices + 1])
            - np.maximum(lower, target_edges[target_indices]),
        )
        rebinned[target_indices] += value * overlaps / (upper - lower)
    return target_edges, rebinned


def find_local_minima(
    x: np.ndarray, y: np.ndarray, *, atol: float = 1.0e-12
) -> list[dict[str, float | int]]:
    """Return interior local minima; represent a flat minimum at its midpoint."""
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    if x_values.ndim != 1 or y_values.ndim != 1 or len(x_values) != len(y_values):
        raise ValueError("x and y must be one-dimensional arrays of equal length")
    if len(x_values) < 3:
        return []
    if not np.all(np.isfinite(x_values)) or not np.all(np.isfinite(y_values)):
        raise ValueError("x and y must contain only finite values")
    if np.any(np.diff(x_values) <= 0) or atol < 0:
        raise ValueError("x must strictly increase and atol must be nonnegative")

    minima: list[dict[str, float | int]] = []
    index = 0
    while index < len(y_values):
        plateau_start = index
        value = y_values[index]
        index += 1
        while index < len(y_values) and abs(y_values[index] - value) <= atol:
            index += 1
        plateau_end = index - 1
        if plateau_start == 0 or plateau_end == len(y_values) - 1:
            continue
        left, right = y_values[plateau_start - 1], y_values[plateau_end + 1]
        if value <= left + atol and value <= right + atol and (value < left - atol or value < right - atol):
            representative = (plateau_start + plateau_end) // 2
            minima.append(
                {
                    "index": representative,
                    "x": float(x_values[representative]),
                    "y": float(value),
                    "plateau_start_index": plateau_start,
                    "plateau_end_index": plateau_end,
                    "left_y": float(left),
                    "right_y": float(right),
                }
            )
    return minima


def _write_csv(path: Path, header: list[str], rows: Iterable[Iterable[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _portable_input_identifier(path: Path, data_root: Path | None) -> tuple[str, bool]:
    supplied = path.expanduser()
    if not supplied.is_absolute():
        return supplied.as_posix(), False
    if data_root is not None:
        try:
            return supplied.resolve().relative_to(data_root.expanduser().resolve()).as_posix(), False
        except ValueError:
            pass
    return supplied.name, True


def _git_revision() -> str | None:
    package_root = Path(__file__).resolve().parents[2]
    try:
        result = subprocess.run(
            ["git", "-C", str(package_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=2,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout.strip() or None


def _validate_args(args: argparse.Namespace) -> None:
    if args.start < 1 or args.end < args.start:
        raise ValueError("configuration range must satisfy 1 <= start <= end")
    if args.stride < 1 or args.block_frames < 1:
        raise ValueError("stride and block-frames must be positive")
    if args.oxygen_radius <= 0 or args.slice_width <= 0 or args.top_layer_tolerance <= 0:
        raise ValueError("oxygen-radius, slice-width, and top-layer-tolerance must be positive")
    if args.origin is not None and args.origin <= 0:
        raise ValueError("origin width must be positive")
    if args.zmax is not None and args.zmax <= 0:
        raise ValueError("zmax must be positive")


def analyze_density(args: argparse.Namespace) -> dict:
    """Analyze XDATCAR water density and write detailed plus plotting CSVs."""
    _validate_args(args)
    trajectory = Path(args.xdatcar).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    if not trajectory.is_file():
        raise FileNotFoundError(f"XDATCAR not found: {trajectory}")

    structure = read_xdatcar_header(trajectory)
    reference_path = getattr(args, "reference", None)
    reference = None
    if reference_path is not None:
        reference = parse_poscar(Path(reference_path).expanduser())
        if reference.species != structure.species or reference.counts != structure.counts:
            raise ValueError("reference POSCAR species/counts do not match XDATCAR")
    labels = structure.atom_species
    oxygen_indices = np.asarray(
        [index for index, label in enumerate(labels) if label.casefold() == "o"], dtype=int
    )
    if oxygen_indices.size == 0:
        raise ValueError("XDATCAR does not contain oxygen atoms labelled O")
    surface_indices, selected_surface_species = select_surface_indices(labels, args.surface_species)
    normal, period, area = cell_geometry(structure.lattice)

    frame_steps: list[int] = []
    frame_surfaces: list[float] = []
    frame_heights: list[np.ndarray] = []
    frame_raw_heights: list[np.ndarray] = []
    max_ball_edge = 0.0
    expected_steps = list(range(args.start, args.end + 1, args.stride))
    for frame in iter_xdatcar(trajectory, args.start, args.end, args.stride):
        if frame.positions_frac.shape[0] != structure.atom_count:
            raise ValueError(f"configuration {frame.configuration} atom count changed")
        surface_z, heights, raw_heights = surface_relative_heights(
            frame,
            structure,
            normal,
            period,
            surface_indices,
            oxygen_indices,
            args.top_layer_tolerance,
        )
        frame_steps.append(frame.configuration)
        frame_surfaces.append(surface_z)
        frame_heights.append(heights)
        frame_raw_heights.append(raw_heights)
        max_ball_edge = max(max_ball_edge, float(np.max(heights)) + args.oxygen_radius)

    if not frame_steps:
        raise ValueError(f"no XDATCAR configurations found in {args.start}-{args.end}")
    if frame_steps != expected_steps:
        missing = sorted(set(expected_steps) - set(frame_steps))
        extra = sorted(set(frame_steps) - set(expected_steps))
        detail = f"missing configurations: {missing[:10]}" if missing else f"unexpected configurations: {extra[:10]}"
        raise ValueError(f"requested frame grid is incomplete ({detail})")

    requested_zmax = args.zmax if args.zmax is not None else max_ball_edge
    zmax = math.ceil(requested_zmax / args.slice_width - 1.0e-12) * args.slice_width
    if zmax <= 0:
        raise ValueError("profile upper edge must be positive")
    slice_count = int(round(zmax / args.slice_width))
    edges = np.arange(slice_count + 1, dtype=float) * args.slice_width
    edges[-1] = zmax
    centers = 0.5 * (edges[:-1] + edges[1:])

    frame_hits: list[np.ndarray] = []
    frame_equivalents: list[np.ndarray] = []
    all_denominators: list[np.ndarray] = []
    diagnostics: list[tuple[object, ...]] = []
    lower_clipped = upper_clipped = zero_intersections = 0
    oxygen_instances = len(frame_steps) * oxygen_indices.size
    for step, surface_z, heights, raw_heights in zip(
        frame_steps, frame_surfaces, frame_heights, frame_raw_heights
    ):
        hits, equivalents, denominators = normalized_slice_contributions(
            heights, edges, args.oxygen_radius
        )
        frame_hits.append(hits)
        frame_equivalents.append(equivalents)
        all_denominators.append(denominators)
        lower_mask = heights - args.oxygen_radius < edges[0] - 1.0e-12
        upper_mask = heights + args.oxygen_radius > edges[-1] + 1.0e-12
        lower_clipped += int(np.count_nonzero(lower_mask))
        upper_clipped += int(np.count_nonzero(upper_mask))
        zero_intersections += int(np.count_nonzero(denominators == 0))
        wrapped_count = int(np.count_nonzero((raw_heights < 0.0) | (raw_heights >= period)))
        diagnostics.append(
            (
                step,
                surface_z,
                float(np.min(heights)),
                float(np.max(heights)),
                wrapped_count,
                wrapped_count / oxygen_indices.size,
                int(np.count_nonzero(lower_mask)),
                int(np.count_nonzero(upper_mask)),
            )
        )

    hits_array = np.asarray(frame_hits, dtype=float)
    equivalents_array = np.asarray(frame_equivalents, dtype=float)
    denominators_array = np.asarray(all_denominators, dtype=float)
    frame_count = len(frame_steps)
    oxygen_equivalent = equivalents_array.mean(axis=0)
    oxygen_density_cm3 = oxygen_equivalent * OXYGEN_MOLAR_MASS_G_MOL / AVOGADRO / (area * args.slice_width) * 1.0e24
    oxygen_density_a3 = oxygen_equivalent * OXYGEN_MOLAR_MASS_G_MOL / AVOGADRO / (area * args.slice_width)
    water_density_cm3 = oxygen_equivalent * WATER_MOLAR_MASS_G_MOL / AVOGADRO / (area * args.slice_width) * 1.0e24
    water_density_a3 = oxygen_equivalent * WATER_MOLAR_MASS_G_MOL / AVOGADRO / (area * args.slice_width)
    integrated_equivalent = float(np.sum(oxygen_equivalent))
    expected_oxygen_count = int(oxygen_indices.size)
    conservation_error = integrated_equivalent - expected_oxygen_count
    minima = find_local_minima(centers, water_density_cm3)
    warnings: list[str] = []
    if abs(conservation_error) > CONSERVATION_TOLERANCE:
        warnings.append("integrated oxygen-equivalent count differs from the trajectory oxygen count")
    if lower_clipped:
        warnings.append("some oxygen-ball support is clipped at the lower profile boundary")
    if upper_clipped:
        warnings.append("some oxygen-ball support is clipped at the upper profile boundary")
    if zero_intersections:
        warnings.append("some oxygen balls have no intersection with the requested profile")

    output_dir = prepare_output_dir(output_dir, overwrite=getattr(args, "overwrite", False))
    output_files = [
        "water_density_profile.csv",
        "average_density_delta_z.csv",
        "density_local_minima.csv",
        "surface_diagnostics.csv",
        "water_density_blocks.csv",
    ]
    _write_csv(
        output_dir / "water_density_profile.csv",
        [
            "slice_lower_A", "slice_upper_A", "slice_center_A", "mean_intersected_oxygen_balls",
            "mean_oxygen_equivalent_count", "oxygen_mass_density_g_cm-3", "oxygen_mass_density_g_A-3",
            "water_mass_density_g_cm-3", "water_mass_density_g_A-3", "oxygen_radius_A", "slice_width_A",
        ],
        (
            (
                edges[index], edges[index + 1], centers[index], hits_array[:, index].mean(),
                oxygen_equivalent[index], oxygen_density_cm3[index], oxygen_density_a3[index],
                water_density_cm3[index], water_density_a3[index], args.oxygen_radius, args.slice_width,
            )
            for index in range(len(centers))
        ),
    )
    _write_csv(
        output_dir / "average_density_delta_z.csv",
        [
            "delta_z_A", "average_oxygen_mass_density_g_cm-3", "average_oxygen_mass_density_g_A-3",
            "average_water_mass_density_g_cm-3", "average_water_mass_density_g_A-3",
        ],
        zip(centers, oxygen_density_cm3, oxygen_density_a3, water_density_cm3, water_density_a3),
    )
    _write_csv(
        output_dir / "density_local_minima.csv",
        [
            "minimum_index", "delta_z_A", "water_mass_density_g_cm-3", "plateau_start_index",
            "plateau_end_index", "left_density_g_cm-3", "right_density_g_cm-3",
        ],
        (
            (
                item["index"], item["x"], item["y"], item["plateau_start_index"],
                item["plateau_end_index"], item["left_y"], item["right_y"],
            )
            for item in minima
        ),
    )
    _write_csv(
        output_dir / "surface_diagnostics.csv",
        [
            "configuration", "surface_z_A", "minimum_oxygen_delta_z_A", "maximum_oxygen_delta_z_A",
            "wrapped_oxygen_count", "wrapped_oxygen_fraction", "lower_clipped_oxygen_count",
            "upper_clipped_oxygen_count",
        ],
        diagnostics,
    )

    block_rows = []
    for block_index, block_start in enumerate(range(0, frame_count, args.block_frames), start=1):
        block_stop = min(block_start + args.block_frames, frame_count)
        block_values = equivalents_array[block_start:block_stop].mean(axis=0)
        block_oxygen = block_values * OXYGEN_MOLAR_MASS_G_MOL / AVOGADRO / (area * args.slice_width) * 1.0e24
        block_water = block_values * WATER_MOLAR_MASS_G_MOL / AVOGADRO / (area * args.slice_width) * 1.0e24
        block_rows.extend(
            (
                block_index, block_stop - block_start, edges[index], edges[index + 1], centers[index],
                block_values[index], block_oxygen[index], block_water[index],
            )
            for index in range(len(centers))
        )
    _write_csv(
        output_dir / "water_density_blocks.csv",
        [
            "block", "frames", "slice_lower_A", "slice_upper_A", "slice_center_A",
            "mean_oxygen_equivalent_count", "oxygen_mass_density_g_cm-3", "water_mass_density_g_cm-3",
        ],
        block_rows,
    )

    origin_name = None
    if args.origin is not None:
        target_edges, rebinned = rebin_slice_values(oxygen_equivalent, edges, args.origin)
        target_widths = np.diff(target_edges)
        origin_centers = 0.5 * (target_edges[:-1] + target_edges[1:])
        origin_oxygen = rebinned * OXYGEN_MOLAR_MASS_G_MOL / AVOGADRO / (area * target_widths) * 1.0e24
        origin_water = rebinned * WATER_MOLAR_MASS_G_MOL / AVOGADRO / (area * target_widths) * 1.0e24
        origin_name = f"water_density_{format(args.origin, 'g')}A_origin.csv"
        _write_csv(
            output_dir / origin_name,
            ["delta_z_A", "oxygen_mass_density_g_cm-3", "water_mass_density_g_cm-3"],
            zip(origin_centers, origin_oxygen, origin_water),
        )
        output_files.append(origin_name)

    input_identifier, redacted = _portable_input_identifier(trajectory, args.data_root)
    metadata = {
        "schema_version": 1,
        "analysis": "slice-based oxygen-ball water density",
        "package_version": __version__,
        "git_revision": _git_revision(),
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "operating_system": platform.system(),
        },
        "input": {
            "xdatcar": input_identifier,
            "absolute_path_redacted": redacted,
            "size_bytes": trajectory.stat().st_size,
            "sha256": _sha256(trajectory),
        },
        "reference_title": structure.comment,
        "reference_poscar": (
            {
                "path": _portable_input_identifier(Path(reference_path).expanduser(), args.data_root)[0],
                "absolute_path_redacted": _portable_input_identifier(
                    Path(reference_path).expanduser(), args.data_root
                )[1],
                "size_bytes": Path(reference_path).expanduser().stat().st_size,
                "sha256": _sha256(Path(reference_path).expanduser()),
            }
            if reference_path is not None
            else None
        ),
        "species": list(structure.species),
        "counts": list(structure.counts),
        "atom_count": structure.atom_count,
        "oxygen_atom_count": expected_oxygen_count,
        "frames_analyzed": frame_count,
        "first_configuration": frame_steps[0],
        "last_configuration": frame_steps[-1],
        "configuration_stride": args.stride,
        "surface_species": selected_surface_species,
        "surface_definition": "instantaneous mean projected height of selected substrate atoms with z_i >= z_max(surface) - top_layer_tolerance",
        "top_layer_tolerance_A": args.top_layer_tolerance,
        "profile_direction": "+ surface normal from instantaneous surface; oxygen heights mapped to [0, cell_period)",
        "surface_normal": normal.tolist(),
        "cell_period_normal_A": period,
        "lateral_area_A2": area,
        "oxygen_radius_A": args.oxygen_radius,
        "slice_width_A": args.slice_width,
        "slice_rule": "Each upward slice with nonzero overlap receives 1/n_i oxygen equivalents from an oxygen ball; edge tangencies are excluded.",
        "oxygen_molar_mass_g_mol": OXYGEN_MOLAR_MASS_G_MOL,
        "water_molar_mass_g_mol": WATER_MOLAR_MASS_G_MOL,
        "water_density_model": "Each oxygen equivalent represents one H2O molecule; hydrogen atoms are not independently spatially resolved.",
        "profile_range_A": [float(edges[0]), float(edges[-1])],
        "requested_zmax_A": args.zmax,
        "integrated_oxygen_equivalent_count": integrated_equivalent,
        "expected_oxygen_atom_count": expected_oxygen_count,
        "integrated_oxygen_equivalent_error": conservation_error,
        "oxygen_equivalent_conservation_tolerance": CONSERVATION_TOLERANCE,
        "oxygen_slice_intersection_count": {
            "minimum": float(np.min(denominators_array)),
            "maximum": float(np.max(denominators_array)),
            "mean": float(np.mean(denominators_array)),
            "zero_intersection_count": zero_intersections,
        },
        "profile_edge_diagnostics": {
            "oxygen_frame_instances": oxygen_instances,
            "lower_clipped_count": lower_clipped,
            "lower_clipped_fraction": lower_clipped / oxygen_instances,
            "upper_clipped_count": upper_clipped,
            "upper_clipped_fraction": upper_clipped / oxygen_instances,
            "wrapped_height_count": int(sum(row[4] for row in diagnostics)),
            "wrapped_height_fraction": float(sum(row[4] for row in diagnostics) / oxygen_instances),
        },
        "block_frames": args.block_frames,
        "local_minima_curve": "average_water_mass_density_g_cm-3",
        "local_minima_definition": "Interior sampled minima; plateaus represented once at midpoint; endpoints excluded.",
        "local_minima_count": len(minima),
        "origin_bin_width_A": args.origin,
        "warnings": warnings,
        "output_files": output_files + ["analysis_metadata.json"],
    }
    (output_dir / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8"
    )
    return metadata


def add_density_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the public density command options to an argparse parser."""
    parser.add_argument("--xdatcar", type=Path, required=True, help="Input VASP XDATCAR")
    parser.add_argument("--reference", type=Path, default=None, help="Optional POSCAR used to verify species/order")
    parser.add_argument("--output-dir", "--output_dir", type=Path, required=True, help="Output directory")
    parser.add_argument("--start", type=int, required=True, help="First configuration number")
    parser.add_argument("--end", type=int, required=True, help="Last configuration number")
    parser.add_argument("--stride", "--frame-stride", dest="stride", type=int, default=1, help="Configuration stride (default: 1)")
    parser.add_argument("--oxygen-radius", type=float, default=0.2, metavar="A")
    parser.add_argument("--slice-width", type=float, default=0.1, metavar="A")
    parser.add_argument("--top-layer-tolerance", "--top_layer_tolerance", type=float, default=0.45, metavar="A")
    parser.add_argument(
        "--surface-species", default=None, metavar="SPECIES[,SPECIES... ]",
        help="Override substrate atom selection; default excludes O, H, He, and K",
    )
    parser.add_argument("--origin", type=float, default=None, metavar="A", help="Also write a rebinned Origin CSV")
    parser.add_argument("--zmax", type=float, default=None, metavar="A", help="Optional upper profile edge above the surface")
    parser.add_argument("--block-frames", type=int, default=500, help="Frames per uncertainty block")
    parser.add_argument("--data-root", type=Path, default=None, help="Optional root used to store a relative input identifier")
    parser.add_argument("--overwrite", action="store_true", help="Allow reuse of a non-empty output directory")
