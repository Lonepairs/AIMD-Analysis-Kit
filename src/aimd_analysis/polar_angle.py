"""Water orientation and density-weighted polar-angle analysis."""

from __future__ import annotations

import argparse
import csv
import json
from itertools import product
import math
import platform
from pathlib import Path
from typing import Iterable

import numpy as np

from . import __version__
from .io import prepare_output_dir
from .density import (
    _git_revision,
    _portable_input_identifier,
    _sha256,
    cell_geometry,
    select_surface_indices,
    surface_relative_heights,
)
from .structure import parse_poscar
from .trajectory import TrajectoryFrame, iter_xdatcar, read_xdatcar_header


AVOGADRO = 6.02214076e23
WATER_MOLAR_MASS_G_MOL = 18.015


def _write_csv(path: Path, header: list[str], rows: Iterable[Iterable[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def minimum_image_vectors(lattice: np.ndarray, deltas: np.ndarray) -> np.ndarray:
    """Apply the minimum-image calculation to many vectors without per-pair SVDs."""
    values = np.asarray(deltas, dtype=float)
    if values.shape[-1] != 3:
        raise ValueError("vectors must end in a length-three coordinate")
    flat = values.reshape(-1, 3)
    rounded = np.rint(flat)
    base = flat - rounded
    gram = lattice @ lattice.T
    off_diagonal = gram - np.diag(np.diag(gram))
    if np.max(np.abs(off_diagonal)) <= 1.0e-10 * max(1.0, float(np.max(np.diag(gram)))):
        wrapped = base @ lattice
        return wrapped.reshape(values.shape)

    best = None
    best_sq = np.full(len(flat), np.inf)
    for image in product((-1, 0, 1), repeat=3):
        candidate = (base + np.asarray(image, dtype=float)) @ lattice
        squared = np.einsum("ij,ij->i", candidate, candidate)
        improve = squared < best_sq
        if np.any(improve):
            if best is None:
                best = np.zeros_like(candidate)
            best[improve] = candidate[improve]
            best_sq[improve] = squared[improve]
    wrapped = best
    return wrapped.reshape(values.shape)


def assign_hydrogen_parents(
    positions_frac: np.ndarray,
    lattice: np.ndarray,
    oxygen_indices: np.ndarray,
    hydrogen_indices: np.ndarray,
    max_distance: float,
) -> dict[int, tuple[int, int]]:
    """Assign two unique H atoms per O by global greedy minimum distance."""
    if len(hydrogen_indices) != 2 * len(oxygen_indices):
        raise ValueError(
            f"expected two H atoms per O; found {len(oxygen_indices)} O and {len(hydrogen_indices)} H"
        )
    deltas = positions_frac[hydrogen_indices][None, :, :] - positions_frac[oxygen_indices][:, None, :]
    vectors = minimum_image_vectors(lattice, deltas)
    distances = np.linalg.norm(vectors, axis=-1)
    pairs = sorted(
        (float(distances[oxygen_position, hydrogen_position]), int(oxygen), int(hydrogen))
        for oxygen_position, oxygen in enumerate(oxygen_indices)
        for hydrogen_position, hydrogen in enumerate(hydrogen_indices)
    )
    assigned: dict[int, list[int]] = {int(oxygen): [] for oxygen in oxygen_indices}
    used_hydrogens: set[int] = set()
    for distance, oxygen, hydrogen in pairs:
        if distance > max_distance:
            break
        if len(assigned[oxygen]) < 2 and hydrogen not in used_hydrogens:
            assigned[oxygen].append(hydrogen)
            used_hydrogens.add(hydrogen)
            if len(used_hydrogens) == len(hydrogen_indices):
                break
    if len(used_hydrogens) != len(hydrogen_indices) or any(len(values) != 2 for values in assigned.values()):
        raise ValueError(
            f"could not assign two unique H atoms per O within {max_distance:.3f} A"
        )
    return {oxygen: tuple(sorted(values)) for oxygen, values in assigned.items()}


def _profile_edges(period: float, width: float) -> np.ndarray:
    count = max(1, int(math.ceil(period / width - 1.0e-12)))
    edges = np.arange(count + 1, dtype=float) * width
    edges[-1] = period
    return edges


def _isotropic_bin_reference(edges: np.ndarray) -> np.ndarray:
    """Return the normalized isotropic probability in each angle bin."""
    lower = np.radians(edges[:-1])
    upper = np.radians(edges[1:])
    return (np.cos(lower) - np.cos(upper)) / 2.0


def _layer_index(delta_z: float, boundaries: tuple[float, float] | None) -> int | None:
    if boundaries is None:
        return None
    first, second = boundaries
    if delta_z < first:
        return 1
    if delta_z < second:
        return 2
    return 3


def _validate_args(args: argparse.Namespace) -> None:
    if args.start < 1 or args.end < args.start:
        raise ValueError("configuration range must satisfy 1 <= start <= end")
    if args.stride < 1 or args.block_frames < 1:
        raise ValueError("stride and block-frames must be positive")
    if args.parent_oh_max <= 0 or args.oh_min < 0 or args.oh_max <= args.oh_min:
        raise ValueError("require parent-oh-max > 0 and 0 <= oh-min < oh-max")
    if args.alpha_bins < 1 or args.profile_slice_width <= 0:
        raise ValueError("alpha-bins must be positive and profile-slice-width must be positive")
    if args.origin is not None and (args.origin <= 0 or args.origin > 180):
        raise ValueError("origin angle width must be in (0, 180]")
    if args.boundaries is not None and (
        len(args.boundaries) != 2 or args.boundaries[0] <= 0 or args.boundaries[1] <= args.boundaries[0]
    ):
        raise ValueError("boundaries must satisfy 0 < first < second")


def analyze_polar_angle(args: argparse.Namespace) -> dict:
    """Analyze water polar angle, orientation density, and optional layers."""
    _validate_args(args)
    trajectory = Path(args.xdatcar).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    if not trajectory.is_file():
        raise FileNotFoundError(f"XDATCAR not found: {trajectory}")
    structure = read_xdatcar_header(trajectory)
    reference_path = getattr(args, "reference", None)
    if reference_path is not None:
        reference = parse_poscar(Path(reference_path).expanduser())
        if reference.species != structure.species or reference.counts != structure.counts:
            raise ValueError("reference POSCAR species/counts do not match XDATCAR")

    labels = structure.atom_species
    oxygen_indices = np.asarray([i for i, label in enumerate(labels) if label.casefold() == "o"], dtype=int)
    hydrogen_indices = np.asarray([i for i, label in enumerate(labels) if label.casefold() == "h"], dtype=int)
    if oxygen_indices.size == 0 or hydrogen_indices.size == 0:
        raise ValueError("XDATCAR must contain atoms labelled O and H")
    surface_indices, surface_species = select_surface_indices(labels, args.surface_species)
    normal, period, area = cell_geometry(structure.lattice)
    profile_edges = _profile_edges(period, args.profile_slice_width)
    profile_widths = np.diff(profile_edges)
    profile_centers = 0.5 * (profile_edges[:-1] + profile_edges[1:])
    expected_steps = list(range(args.start, args.end + 1, args.stride))
    parent_mapping: dict[int, tuple[int, int]] | None = None
    steps: list[int] = []
    surfaces: list[float] = []
    all_angles: list[float] = []
    layer_angles: dict[int, list[float]] | None = {1: [], 2: [], 3: []} if args.boundaries else None
    frame_rows: list[tuple[object, ...]] = []
    frame_angles: list[list[float]] = []
    point_rows: list[tuple[object, ...]] = []
    profile_oxygen = []
    profile_valid = []
    profile_cosine = []
    block_data: list[dict[str, object]] = []

    for frame in iter_xdatcar(trajectory, args.start, args.end, args.stride):
        if frame.positions_frac.shape[0] != structure.atom_count:
            raise ValueError(f"configuration {frame.configuration} atom count changed")
        if parent_mapping is None:
            parent_mapping = assign_hydrogen_parents(
                frame.positions_frac,
                structure.lattice,
                oxygen_indices,
                hydrogen_indices,
                args.parent_oh_max,
            )
        surface_z, oxygen_heights, _ = surface_relative_heights(
            frame,
            structure,
            normal,
            period,
            surface_indices,
            oxygen_indices,
            args.top_layer_tolerance,
        )
        steps.append(frame.configuration)
        surfaces.append(surface_z)
        vector_rows = []
        for oxygen in oxygen_indices:
            hydrogens = parent_mapping[int(oxygen)]
            hydrogen_cart = frame.positions_frac[list(hydrogens)]
            vector_rows.append(
                minimum_image_vectors(
                    structure.lattice,
                    hydrogen_cart - frame.positions_frac[int(oxygen)],
                )
            )
        oh_vectors = np.asarray(vector_rows, dtype=float)
        oh_lengths = np.linalg.norm(oh_vectors, axis=2)
        valid_bonds = (oh_lengths >= args.oh_min) & (oh_lengths < args.oh_max)
        unit_vectors = np.divide(
            oh_vectors,
            oh_lengths[:, :, None],
            out=np.zeros_like(oh_vectors),
            where=oh_lengths[:, :, None] > 0,
        )
        bisectors = unit_vectors.sum(axis=1)
        bisector_norms = np.linalg.norm(bisectors, axis=1)
        valid = np.all(valid_bonds, axis=1) & (bisector_norms > 1.0e-12)
        cosine = np.full(oxygen_indices.size, np.nan)
        angles = np.full(oxygen_indices.size, np.nan)
        normalized = bisectors[valid] / bisector_norms[valid, None]
        cosine[valid] = np.clip(normalized @ normal, -1.0, 1.0)
        angles[valid] = np.degrees(np.arccos(cosine[valid]))
        all_angles.extend(float(value) for value in angles[valid])
        frame_angles.append([float(value) for value in angles[valid]])
        if layer_angles is not None:
            for oxygen_position in np.flatnonzero(valid):
                layer = _layer_index(float(oxygen_heights[oxygen_position]), args.boundaries)
                layer_angles[layer].append(float(angles[oxygen_position]))
        oxygen_counts, _ = np.histogram(oxygen_heights, bins=profile_edges)
        valid_counts, _ = np.histogram(oxygen_heights[valid], bins=profile_edges)
        cosine_sums, _ = np.histogram(oxygen_heights[valid], bins=profile_edges, weights=cosine[valid])
        profile_oxygen.append(oxygen_counts.astype(float))
        profile_valid.append(valid_counts.astype(float))
        profile_cosine.append(cosine_sums.astype(float))
        frame_mean = float(np.mean(angles[valid])) if np.any(valid) else ""
        frame_std = float(np.std(angles[valid])) if np.any(valid) else ""
        frame_median = float(np.median(angles[valid])) if np.any(valid) else ""
        frame_rows.append((frame.configuration, frame_mean, frame_std, frame_median, int(np.count_nonzero(valid)), int(oxygen_indices.size)))
        for position in np.flatnonzero(valid):
            hydrogens = parent_mapping[int(oxygen_indices[position])]
            point_rows.append(
                (
                    frame.configuration,
                    int(position),
                    int(oxygen_indices[position]),
                    int(hydrogens[0]),
                    int(hydrogens[1]),
                    float(oxygen_heights[position]),
                    float(angles[position]),
                    float(cosine[position]),
                    float(oh_lengths[position, 0]),
                    float(oh_lengths[position, 1]),
                    _layer_index(float(oxygen_heights[position]), args.boundaries) or "",
                )
            )

        block_index = (len(steps) - 1) // args.block_frames
        if block_index == len(block_data):
            block_data.append(
                {
                    "frames": 0,
                    "oxygen": np.zeros(len(profile_widths)),
                    "valid": np.zeros(len(profile_widths)),
                    "cosine": np.zeros(len(profile_widths)),
                }
            )
        block_data[block_index]["frames"] = int(block_data[block_index]["frames"]) + 1
        block_data[block_index]["oxygen"] += oxygen_counts
        block_data[block_index]["valid"] += valid_counts
        block_data[block_index]["cosine"] += cosine_sums

    if not steps:
        raise ValueError(f"no XDATCAR configurations found in {args.start}-{args.end}")
    if steps != expected_steps:
        missing = sorted(set(expected_steps) - set(steps))
        raise ValueError(f"requested frame grid is incomplete; missing configurations: {missing[:10]}")
    if parent_mapping is None:
        raise ValueError("no parent mapping was established")

    profile_oxygen_array = np.asarray(profile_oxygen, dtype=float)
    profile_valid_array = np.asarray(profile_valid, dtype=float)
    profile_cosine_array = np.asarray(profile_cosine, dtype=float)
    frame_count = len(steps)
    oxygen_mean = profile_oxygen_array.mean(axis=0)
    valid_mean = profile_valid_array.mean(axis=0)
    cosine_mean = profile_cosine_array.mean(axis=0)
    volume = area * profile_widths * frame_count
    oxygen_density = oxygen_mean / (area * profile_widths)
    valid_density = valid_mean / (area * profile_widths)
    mean_cosine = np.divide(cosine_mean, valid_mean, out=np.full_like(cosine_mean, np.nan), where=valid_mean > 0)
    rho_cosine = cosine_mean / (area * profile_widths)
    water_mass_density = oxygen_density * WATER_MOLAR_MASS_G_MOL / AVOGADRO * 1.0e24
    valid_mass_density = valid_density * WATER_MOLAR_MASS_G_MOL / AVOGADRO * 1.0e24
    rho_mass_cosine = rho_cosine * WATER_MOLAR_MASS_G_MOL / AVOGADRO * 1.0e24
    valid_fraction = np.divide(valid_mean, oxygen_mean, out=np.full_like(valid_mean, np.nan), where=oxygen_mean > 0)

    output_dir = prepare_output_dir(output_dir, overwrite=getattr(args, "overwrite", False))
    output_files = [
        "polar_angle_per_frame.csv",
        "polar_angle_points.csv",
        "polar_angle_distribution.csv",
        "polar_angle_blocks.csv",
        "polar_orientation_density_profile.csv",
        "polar_orientation_density_blocks.csv",
    ]
    _write_csv(
        output_dir / "polar_angle_per_frame.csv",
        ["configuration", "mean_alpha_deg", "std_alpha_deg", "median_alpha_deg", "valid_waters", "total_waters"],
        frame_rows,
    )
    _write_csv(
        output_dir / "polar_angle_points.csv",
        [
            "configuration", "water_index", "oxygen_index_zero_based", "hydrogen_1_index_zero_based",
            "hydrogen_2_index_zero_based", "delta_z_A", "alpha_deg", "cos_phi", "oh_1_A", "oh_2_A", "layer",
        ],
        point_rows,
    )
    angle_edges = np.linspace(0.0, 180.0, args.alpha_bins + 1)
    histogram, _ = np.histogram(all_angles, bins=angle_edges)
    total_angles = int(histogram.sum())
    isotropic = _isotropic_bin_reference(angle_edges)
    _write_csv(
        output_dir / "polar_angle_distribution.csv",
        ["alpha_lower_deg", "alpha_upper_deg", "alpha_center_deg", "count", "fraction", "isotropic_reference"],
        (
            (
                angle_edges[i], angle_edges[i + 1], 0.5 * (angle_edges[i] + angle_edges[i + 1]),
                int(histogram[i]), float(histogram[i] / total_angles) if total_angles else 0.0, isotropic[i],
            )
            for i in range(args.alpha_bins)
        ),
    )
    _write_csv(
        output_dir / "polar_angle_blocks.csv",
        ["block", "frames", "alpha_lower_deg", "alpha_upper_deg", "alpha_center_deg", "count", "fraction", "isotropic_reference"],
        (
            row
            for block_index, block_start in enumerate(range(0, frame_count, args.block_frames), start=1)
            for row in _angle_block_rows(
                block_index,
                frame_angles[block_start : block_start + args.block_frames],
                angle_edges,
            )
        ),
    )
    _write_csv(
        output_dir / "polar_orientation_density_profile.csv",
        [
            "delta_z_lower_A", "delta_z_upper_A", "delta_z_A", "mean_water_oxygen_count_per_frame",
            "water_number_density_A-3", "mean_valid_orientation_count_per_frame", "valid_orientation_number_density_A-3",
            "valid_orientation_fraction", "water_mass_density_g_cm-3", "valid_orientation_mass_density_g_cm-3",
            "mean_cos_phi", "rho_H2O_cos_phi_A-3", "rho_H2O_cos_phi_g_cm-3",
        ],
        (
            (
                profile_edges[i], profile_edges[i + 1], profile_centers[i], oxygen_mean[i], oxygen_density[i],
                valid_mean[i], valid_density[i], valid_fraction[i] if np.isfinite(valid_fraction[i]) else "",
                water_mass_density[i], valid_mass_density[i], mean_cosine[i] if np.isfinite(mean_cosine[i]) else "",
                rho_cosine[i], rho_mass_cosine[i],
            )
            for i in range(len(profile_centers))
        ),
    )
    profile_block_rows = []
    for block_index, block in enumerate(block_data, start=1):
        frames_in_block = int(block["frames"])
        oxygen_block = np.asarray(block["oxygen"], dtype=float)
        valid_block = np.asarray(block["valid"], dtype=float)
        cosine_block = np.asarray(block["cosine"], dtype=float)
        block_oxygen_density = oxygen_block / (frames_in_block * area * profile_widths)
        block_valid_density = valid_block / (frames_in_block * area * profile_widths)
        block_mean_cosine = np.divide(cosine_block, valid_block, out=np.full_like(cosine_block, np.nan), where=valid_block > 0)
        block_rho_cosine = cosine_block / (frames_in_block * area * profile_widths)
        profile_block_rows.extend(
            (
                block_index, frames_in_block, profile_edges[i], profile_edges[i + 1], profile_centers[i],
                block_oxygen_density[i], block_valid_density[i],
                block_oxygen_density[i] * WATER_MOLAR_MASS_G_MOL / AVOGADRO * 1.0e24,
                block_valid_density[i] * WATER_MOLAR_MASS_G_MOL / AVOGADRO * 1.0e24,
                block_mean_cosine[i] if np.isfinite(block_mean_cosine[i]) else "", block_rho_cosine[i],
                block_rho_cosine[i] * WATER_MOLAR_MASS_G_MOL / AVOGADRO * 1.0e24,
            )
            for i in range(len(profile_centers))
        )
    _write_csv(
        output_dir / "polar_orientation_density_blocks.csv",
        [
            "block", "frames", "delta_z_lower_A", "delta_z_upper_A", "delta_z_A", "water_number_density_A-3",
            "valid_orientation_number_density_A-3", "water_mass_density_g_cm-3", "valid_orientation_mass_density_g_cm-3",
            "mean_cos_phi", "rho_H2O_cos_phi_A-3", "rho_H2O_cos_phi_g_cm-3",
        ],
        profile_block_rows,
    )

    origin_name = None
    if args.origin is not None:
        origin_edges = np.arange(0.0, 180.0, float(args.origin))
        if origin_edges.size == 0 or not math.isclose(origin_edges[-1], 180.0, abs_tol=1.0e-12):
            origin_edges = np.append(origin_edges, 180.0)
        origin_centers = 0.5 * (origin_edges[:-1] + origin_edges[1:])
        origin_name = f"polar_angle_distribution_{format(args.origin, 'g')}deg"
        if layer_angles is None:
            origin_histogram, _ = np.histogram(all_angles, bins=origin_edges)
            _write_csv(
                output_dir / f"{origin_name}.csv",
                ["alpha_center_deg", "fraction", "isotropic_reference"],
                (
                    (origin_centers[i], origin_histogram[i] / total_angles if total_angles else 0.0, _isotropic_bin_reference(origin_edges)[i])
                    for i in range(len(origin_centers))
                ),
            )
        else:
            layer_histograms = {layer: np.histogram(values, bins=origin_edges)[0] for layer, values in layer_angles.items()}
            layer_totals = {layer: int(values.sum()) for layer, values in layer_histograms.items()}
            _write_csv(
                output_dir / f"{origin_name}_by_layer.csv",
                ["alpha_center_deg", "Layer_I", "Layer_II", "Layer_III", "isotropic_reference"],
                (
                    (
                        origin_centers[i],
                        layer_histograms[1][i] / layer_totals[1] if layer_totals[1] else 0.0,
                        layer_histograms[2][i] / layer_totals[2] if layer_totals[2] else 0.0,
                        layer_histograms[3][i] / layer_totals[3] if layer_totals[3] else 0.0,
                        _isotropic_bin_reference(origin_edges)[i],
                    )
                    for i in range(len(origin_centers))
                ),
            )
        output_files.append(f"{origin_name}{'_by_layer' if layer_angles is not None else ''}.csv")

    if layer_angles is not None:
        layer_summary_name = "polar_angle_layer_summary.csv"
        _write_csv(
            output_dir / layer_summary_name,
            ["layer", "lower_bound_A", "upper_bound_A", "valid_water_observations", "mean_alpha_deg", "sem_alpha_deg"],
            (
                (
                    f"Layer_{['I', 'II', 'III'][layer - 1]}",
                    0.0 if layer == 1 else args.boundaries[layer - 2],
                    args.boundaries[layer - 1] if layer < 3 else "",
                    len(values),
                    float(np.mean(values)) if values else "",
                    float(np.std(values, ddof=1) / math.sqrt(len(values))) if len(values) > 1 else "",
                )
                for layer, values in layer_angles.items()
            ),
        )
        output_files.append(layer_summary_name)

    input_identifier, redacted = _portable_input_identifier(trajectory, args.data_root)
    metadata = {
        "schema_version": 1,
        "analysis": "water polar angle and orientation-weighted density",
        "package_version": __version__,
        "git_revision": _git_revision(),
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "operating_system": platform.system()},
        "input": {
            "xdatcar": input_identifier,
            "absolute_path_redacted": redacted,
            "size_bytes": trajectory.stat().st_size,
            "sha256": _sha256(trajectory),
        },
        "reference_poscar": (
            {
                "path": _portable_input_identifier(Path(reference_path).expanduser(), args.data_root)[0],
                "absolute_path_redacted": _portable_input_identifier(Path(reference_path).expanduser(), args.data_root)[1],
                "size_bytes": Path(reference_path).expanduser().stat().st_size,
                "sha256": _sha256(Path(reference_path).expanduser()),
            }
            if reference_path is not None
            else None
        ),
        "species": list(structure.species),
        "counts": list(structure.counts),
        "frames_analyzed": frame_count,
        "first_configuration": steps[0],
        "last_configuration": steps[-1],
        "configuration_stride": args.stride,
        "oxygen_count": int(oxygen_indices.size),
        "hydrogen_count": int(hydrogen_indices.size),
        "surface_species": surface_species,
        "surface_definition": "instantaneous mean projected height of selected substrate atoms with z_i >= z_max(surface) - top_layer_tolerance",
        "top_layer_tolerance_A": args.top_layer_tolerance,
        "parent_assignment": "Global greedy minimum-distance O-H matching in the first selected frame; fixed for all selected frames.",
        "parent_assignment_max_OH_A": args.parent_oh_max,
        "parent_hydrogen_indices_zero_based": [
            {"oxygen": int(oxygen), "hydrogens": [int(hydrogen) for hydrogen in hydrogens]}
            for oxygen, hydrogens in parent_mapping.items()
        ],
        "indicator": {
            "name": "phi",
            "definition": "angle between the sum of separately normalized parent O-H vectors and outward surface normal",
            "range_deg": [0.0, 180.0],
            "periodicity": "full-cell minimum-image convention",
        },
        "valid_OH_range_A": [args.oh_min, args.oh_max],
        "valid_water_rule": "Both parent O-H lengths are in [oh_min, oh_max) and the bisector norm is nonzero.",
        "valid_water_observations": total_angles,
        "mean_alpha_deg": float(np.mean(all_angles)) if all_angles else None,
        "orientation_density_profile": {
            "coordinate": "O-atom height relative to instantaneous surface, wrapped into [0, cell_period)",
            "surface_normal": normal.tolist(),
            "cell_period_A": period,
            "lateral_area_A2": area,
            "profile_range_A": [float(profile_edges[0]), float(profile_edges[-1])],
            "bin_width_A": args.profile_slice_width,
            "rho_H2O_cos_phi_definition": "water mass density weighted by the mean cosine of valid parent-OH bisectors",
            "water_molar_mass_g_mol": WATER_MOLAR_MASS_G_MOL,
        },
        "block_frames": args.block_frames,
        "origin_bin_width_deg": args.origin,
        "layer_analysis_enabled": args.boundaries is not None,
        "layer_boundaries_A": list(args.boundaries) if args.boundaries else None,
        "layer_definitions": (
            {"Layer_I": "[0, zmin)", "Layer_II": "[zmin, zmax)", "Layer_III": "[zmax, infinity)"}
            if args.boundaries
            else None
        ),
        "isotropic_reference": "sin(phi), normalized over 0-180 degrees; bin reference is [cos(lower)-cos(upper)]/2",
        "output_files": output_files + ["analysis_metadata.json"],
    }
    (output_dir / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    return metadata


def _angle_block_rows(
    block_index: int,
    block_angle_lists: list[list[float]],
    edges: np.ndarray,
) -> list[tuple[object, ...]]:
    """Return a normalized histogram for one actual frame block."""
    values = [angle for angles in block_angle_lists for angle in angles]
    histogram, _ = np.histogram(values, bins=edges)
    total = int(histogram.sum())
    reference = _isotropic_bin_reference(edges)
    return [
        (
            block_index,
            len(block_angle_lists),
            edges[index],
            edges[index + 1],
            0.5 * (edges[index] + edges[index + 1]),
            int(histogram[index]),
            histogram[index] / total if total else 0.0,
            reference[index],
        )
        for index in range(len(edges) - 1)
    ]


def add_polar_angle_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--xdatcar", type=Path, required=True, help="Input VASP XDATCAR")
    parser.add_argument("--reference", type=Path, default=None, help="Optional POSCAR used to verify species/order")
    parser.add_argument("--output-dir", "--output_dir", type=Path, required=True, help="Output directory")
    parser.add_argument("--start", type=int, required=True, help="First configuration number")
    parser.add_argument("--end", type=int, required=True, help="Last configuration number")
    parser.add_argument("--stride", "--frame-stride", dest="stride", type=int, default=1, help="Configuration stride")
    parser.add_argument("--parent-oh-max", type=float, default=1.35, metavar="A")
    parser.add_argument("--oh-min", type=float, default=0.50, metavar="A")
    parser.add_argument("--oh-max", type=float, default=1.25, metavar="A")
    parser.add_argument("--alpha-bins", type=int, default=36)
    parser.add_argument("--block-frames", type=int, default=500)
    parser.add_argument("--profile-slice-width", type=float, default=0.1, metavar="A")
    parser.add_argument("--origin", type=float, default=None, metavar="DEG", help="Origin angle-bin width")
    parser.add_argument("--boundaries", type=float, nargs=2, default=None, metavar=("ZMIN", "ZMAX"))
    parser.add_argument("--top-layer-tolerance", type=float, default=0.45, metavar="A")
    parser.add_argument("--surface-species", default=None, metavar="SPECIES[,SPECIES...]")
    parser.add_argument("--data-root", type=Path, default=None, help="Optional root used for portable metadata paths")
    parser.add_argument("--overwrite", action="store_true", help="Allow reuse of a non-empty output directory")
