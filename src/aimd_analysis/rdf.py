"""Periodic radial distribution functions and cumulative coordination numbers."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
from pathlib import Path
from typing import Iterable

import numpy as np

from . import __version__
from .density import _git_revision, _portable_input_identifier, _sha256, cell_geometry
from .io import prepare_output_dir
from .structure import parse_poscar, surface_normal
from .trajectory import iter_xdatcar, read_xdatcar_header


DEFAULT_CN_CUTOFF_A = 2.5


def _write_csv(path: Path, header: list[str], rows: Iterable[Iterable[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def cell_dimensions(lattice: np.ndarray) -> np.ndarray:
    vectors = np.asarray(lattice, dtype=float)
    lengths = np.linalg.norm(vectors, axis=1)
    if np.any(lengths <= 0):
        raise ValueError("cell vectors must have nonzero length")

    def angle(first: np.ndarray, second: np.ndarray) -> float:
        cosine = np.dot(first, second) / (np.linalg.norm(first) * np.linalg.norm(second))
        return math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))

    return np.asarray(
        [*lengths, angle(vectors[1], vectors[2]), angle(vectors[0], vectors[2]), angle(vectors[0], vectors[1])],
        dtype=np.float32,
    )


def shortest_lattice_translation(lattice: np.ndarray) -> float:
    """Find the shortest nonzero lattice translation in a 3D periodic cell."""
    vectors = np.asarray(lattice, dtype=float)
    best = math.inf
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                if i == j == k == 0:
                    continue
                best = min(best, float(np.linalg.norm(np.asarray([i, j, k]) @ vectors)))
    return best


def _extrema(x: np.ndarray, y: np.ndarray, kind: str) -> list[dict[str, float | int]]:
    candidates: list[dict[str, float | int]] = []
    for index in range(1, len(y) - 1):
        if kind == "maximum" and y[index] >= y[index - 1] and y[index] > y[index + 1]:
            candidates.append({"index": index, "r_A": float(x[index]), "value": float(y[index])})
        if kind == "minimum" and y[index] <= y[index - 1] and y[index] < y[index + 1]:
            candidates.append({"index": index, "r_A": float(x[index]), "value": float(y[index])})
    return candidates


def _rebin_histogram_counts(
    counts: np.ndarray, source_edges: np.ndarray, target_width: float
) -> tuple[np.ndarray, np.ndarray]:
    """Rebin RDF pair counts by radial overlap for a plotting export."""
    count = int(math.ceil((source_edges[-1] - source_edges[0]) / target_width - 1.0e-12))
    target_edges = source_edges[0] + np.arange(count + 1, dtype=float) * target_width
    target_edges[-1] = source_edges[-1]
    rebinned = np.zeros(count, dtype=float)
    for index, value in enumerate(counts):
        lower, upper = source_edges[index : index + 2]
        first = max(0, int(math.floor((lower - source_edges[0]) / target_width)))
        last = min(count - 1, int(math.ceil((upper - source_edges[0]) / target_width - 1.0e-12) - 1))
        if last < first:
            continue
        target_indices = np.arange(first, last + 1)
        overlap = np.maximum(
            0.0,
            np.minimum(upper, target_edges[target_indices + 1])
            - np.maximum(lower, target_edges[target_indices]),
        )
        rebinned[target_indices] += value * overlap / (upper - lower)
    return target_edges, rebinned


def _complete_block_sem(values: list[np.ndarray], frame_counts: list[int], block_frames: int) -> np.ndarray:
    """Estimate SEM from equal-length blocks, excluding a short final block."""
    complete = [value for value, count in zip(values, frame_counts) if count == block_frames]
    if len(complete) < 2:
        width = len(values[0]) if values else 0
        return np.full(width, np.nan)
    return np.std(np.asarray(complete), axis=0, ddof=1) / math.sqrt(len(complete))


def _validate_args(args: argparse.Namespace) -> None:
    if args.start < 1 or args.end < args.start or args.frame_stride < 1:
        raise ValueError("configuration range/frame stride is invalid")
    if args.nbins < 1 or args.block_frames < 1:
        raise ValueError("nbins and block-frames must be positive")
    if args.rmin < 0 or args.rmax <= args.rmin:
        raise ValueError("require 0 <= rmin < rmax")
    if args.x_scope not in {"all", "top-layer"}:
        raise ValueError("x-scope must be all or top-layer")
    if args.top_layer_tolerance <= 0:
        raise ValueError("top-layer-tolerance must be positive")
    if args.origin is not None and args.origin != -1.0 and args.origin <= 0:
        raise ValueError("numeric origin width must be positive")
    if any(cutoff <= args.rmin or cutoff > args.rmax for cutoff in args.cn_cutoff):
        raise ValueError("each CN cutoff must satisfy rmin < cutoff <= rmax")


def analyze_rdf(args: argparse.Namespace) -> dict:
    try:
        import MDAnalysis as mda
        from MDAnalysis.lib.distances import distance_array
    except ImportError as exc:
        raise ImportError("MDAnalysis is required for RDF analysis; install the package dependencies") from exc
    _validate_args(args)
    trajectory = Path(args.xdatcar).expanduser()
    reference_path = Path(args.reference).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    if not trajectory.is_file():
        raise FileNotFoundError(f"XDATCAR not found: {trajectory}")
    if not reference_path.is_file():
        raise FileNotFoundError(f"reference structure not found: {reference_path}")
    structure = parse_poscar(reference_path)
    labels = structure.atom_species
    available = list(dict.fromkeys(structure.species))
    by_lower = {element.casefold(): element for element in available}
    for label in (args.element_x, args.element_y):
        if label.casefold() not in by_lower:
            raise ValueError(f"element {label!r} not found in reference; available: {', '.join(available)}")
    element_x = by_lower[args.element_x.casefold()]
    element_y = by_lower[args.element_y.casefold()]
    y_indices = np.asarray([i for i, label in enumerate(labels) if label.casefold() == element_y.casefold()], dtype=int)
    x_all = np.asarray([i for i, label in enumerate(labels) if label.casefold() == element_x.casefold()], dtype=int)
    x_indices = x_all.copy()
    x_selection = "all atoms of selected element"
    if args.x_scope == "top-layer":
        normal = surface_normal(structure.lattice)
        surface_elements = {label.casefold() for label in args.surface_species.split(",") if label.strip()} if args.surface_species else {label.casefold() for label in labels if label.casefold() not in {"o", "h", "he", "k"}}
        unknown = surface_elements - {label.casefold() for label in labels}
        if unknown:
            raise ValueError("surface species not present in reference: " + ", ".join(sorted(unknown)))
        substrate_indices = np.asarray([i for i, label in enumerate(labels) if label.casefold() in surface_elements], dtype=int)
        if not substrate_indices.size:
            raise ValueError("top-layer selection has no surface atoms")
        projected = structure.positions_cart @ normal
        maximum = float(np.max(projected[substrate_indices]))
        top_mask = projected[x_all] >= maximum - args.top_layer_tolerance
        x_indices = x_all[top_mask]
        x_selection = "top-layer subset of X selected from reference POSCAR using projected height tolerance"
    if not x_indices.size or not y_indices.size:
        raise ValueError("selected reference and neighbor groups must both contain at least one atom")
    lattice = structure.lattice
    dimensions = cell_dimensions(lattice)
    inscribed_radius = shortest_lattice_translation(lattice) / 2.0
    if args.rmax > inscribed_radius + 1.0e-8:
        raise ValueError(f"rmax={args.rmax:g} A exceeds half the shortest lattice translation ({inscribed_radius:.6g} A)")
    trajectory_structure = read_xdatcar_header(trajectory)
    if structure.species != trajectory_structure.species or structure.counts != trajectory_structure.counts:
        raise ValueError("reference POSCAR species/counts do not match XDATCAR")
    if not np.allclose(structure.lattice, trajectory_structure.lattice, rtol=1.0e-6, atol=1.0e-6):
        raise ValueError("reference POSCAR lattice does not match XDATCAR")
    if len(labels) != len(trajectory_structure.atom_species):
        raise ValueError("reference POSCAR atom count does not match XDATCAR")
    x_overlap = len(set(x_indices.tolist()) & set(y_indices.tolist()))
    effective_pairs_per_center = len(y_indices) - (x_overlap / len(x_indices))
    if effective_pairs_per_center <= 0:
        raise ValueError("selected groups contain no unique X-Y pairs after self exclusion")

    edges = np.linspace(args.rmin, args.rmax, args.nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    shell_volumes = 4.0 * np.pi / 3.0 * (edges[1:] ** 3 - edges[:-1] ** 3)
    box = dimensions
    total_counts = np.zeros(args.nbins, dtype=float)
    block_count_arrays: list[np.ndarray] = []
    block_steps: list[list[int]] = []
    steps_in_block: list[int] = []
    counts_in_block = np.zeros(args.nbins, dtype=float)
    steps: list[int] = []
    for frame in iter_xdatcar(trajectory, args.start, args.end, args.frame_stride):
        if frame.positions_frac.shape[0] != len(labels):
            raise ValueError(f"configuration {frame.configuration} atom count does not match reference")
        positions = frame.positions_frac @ lattice
        distances = distance_array(positions[x_indices], positions[y_indices], box=box)
        if x_overlap:
            y_position_by_atom = {int(atom): position for position, atom in enumerate(y_indices)}
            for x_position, atom in enumerate(x_indices):
                y_position = y_position_by_atom.get(int(atom))
                if y_position is not None:
                    distances[x_position, y_position] = np.inf
        counts, _ = np.histogram(distances[np.isfinite(distances)], bins=edges)
        total_counts += counts
        counts_in_block += counts
        steps.append(frame.configuration)
        steps_in_block.append(frame.configuration)
        if len(steps_in_block) == args.block_frames:
            block_count_arrays.append(counts_in_block.copy())
            block_steps.append(steps_in_block.copy())
            steps_in_block = []
            counts_in_block.fill(0.0)
    if steps_in_block:
        block_steps.append(steps_in_block)
        block_count_arrays.append(counts_in_block.copy())
    if not steps:
        raise ValueError(f"no XDATCAR configurations found in {args.start}-{args.end}")
    expected_steps = list(range(args.start, args.end + 1, args.frame_stride))
    if steps != expected_steps:
        missing = sorted(set(expected_steps) - set(steps))
        raise ValueError(f"requested frame grid is incomplete; missing configurations: {missing[:10]}")

    frames = len(steps)
    volume = structure.volume
    expected_histogram = frames * len(x_indices) * effective_pairs_per_center * shell_volumes / volume
    rdf = np.divide(total_counts, expected_histogram, out=np.zeros_like(total_counts), where=expected_histogram > 0)
    coordination = np.cumsum(total_counts) / (frames * len(x_indices))
    block_rdf: list[np.ndarray] = []
    block_cn: list[np.ndarray] = []
    block_rows = []
    for block_number, (step_group, count_array) in enumerate(zip(block_steps, block_count_arrays), start=1):
        block_frames = len(step_group)
        expected_block = block_frames * len(x_indices) * effective_pairs_per_center * shell_volumes / volume
        values = np.divide(count_array, expected_block, out=np.zeros_like(count_array), where=expected_block > 0)
        cn = np.cumsum(count_array) / (block_frames * len(x_indices))
        block_rdf.append(values)
        block_cn.append(cn)
        block_rows.extend(
            (block_number, block_frames, step_group[0], step_group[-1], centers[i], values[i], cn[i])
            for i in range(args.nbins)
        )
    block_frame_counts = [len(step_group) for step_group in block_steps]
    rdf_se = _complete_block_sem(block_rdf, block_frame_counts, args.block_frames)
    cn_se = _complete_block_sem(block_cn, block_frame_counts, args.block_frames)

    output_dir = prepare_output_dir(output_dir, overwrite=getattr(args, "overwrite", False))
    pair_name = f"{element_x}-{element_y}"
    detailed_name = f"rdf_{pair_name}.csv"
    blocks_name = f"rdf_{pair_name}_blocks.csv"
    metadata_name = "rdf_metadata.json"
    analysis_metadata_name = "analysis_metadata.json"
    _write_csv(
        output_dir / detailed_name,
        ["r_A", "RDF", "RDF_block_sem", "pair_count_all_frames", "mean_pairs_per_frame", "CN", "CN_block_sem"],
        (
            (centers[i], rdf[i], rdf_se[i] if np.isfinite(rdf_se[i]) else "", total_counts[i], total_counts[i] / frames, coordination[i], cn_se[i] if np.isfinite(cn_se[i]) else "")
            for i in range(args.nbins)
        ),
    )
    _write_csv(
        output_dir / blocks_name,
        ["block", "frames", "first_configuration", "last_configuration", "r_A", "RDF", "CN"],
        block_rows,
    )
    origin_name = None
    if args.origin is not None:
        origin_name = f"rdf_{pair_name}_origin.csv"
        if isinstance(args.origin, (float, int)) and not isinstance(args.origin, bool) and args.origin > 0:
            origin_edges, origin_counts = _rebin_histogram_counts(total_counts, edges, float(args.origin))
            origin_centers = 0.5 * (origin_edges[:-1] + origin_edges[1:])
            origin_shells = 4.0 * np.pi / 3.0 * (origin_edges[1:] ** 3 - origin_edges[:-1] ** 3)
            origin_expected = frames * len(x_indices) * effective_pairs_per_center * origin_shells / volume
            origin_rdf = np.divide(origin_counts, origin_expected, out=np.zeros_like(origin_counts), where=origin_expected > 0)
            origin_cn = np.cumsum(origin_counts) / (frames * len(x_indices))
        else:
            origin_centers, origin_rdf, origin_cn = centers, rdf, coordination
        _write_csv(output_dir / origin_name, ["r_A", "RDF", "CN"], zip(origin_centers, origin_rdf, origin_cn))

    peak_candidates = _extrema(centers, rdf, "maximum")
    minimum_candidates = _extrema(centers, rdf, "minimum")
    named_cn = {}
    for cutoff in args.cn_cutoff:
        index = int(np.searchsorted(centers, cutoff, side="right") - 1)
        index = max(0, min(index, len(centers) - 1))
        named_cn[format(cutoff, "g")] = {
            "cutoff_A": cutoff,
            "CN_per_reference_atom": float(coordination[index]),
            "CN_block_sem": float(cn_se[index]) if np.isfinite(cn_se[index]) else None,
            "mean_neighbor_count_per_reference_atom": float(coordination[index]),
            "raw_pair_count_all_frames": float(np.sum(total_counts[centers <= cutoff])),
        }
    input_identifier, redacted = _portable_input_identifier(trajectory, args.data_root)
    ref_identifier, ref_redacted = _portable_input_identifier(reference_path, args.data_root)
    metadata = {
        "schema_version": 1,
        "analysis": "periodic radial distribution and cumulative coordination number",
        "package_version": __version__,
        "git_revision": _git_revision(),
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "mdanalysis": mda.__version__, "operating_system": platform.system()},
        "inputs": {
            "xdatcar": {"path": input_identifier, "absolute_path_redacted": redacted, "size_bytes": trajectory.stat().st_size, "sha256": _sha256(trajectory)},
            "reference": {"path": ref_identifier, "absolute_path_redacted": ref_redacted, "size_bytes": reference_path.stat().st_size, "sha256": _sha256(reference_path)},
        },
        "reference_title": structure.comment,
        "pair": {"X": element_x, "Y": element_y},
        "atom_counts": {"X": len(x_indices), "Y": len(y_indices), "all_X_element": len(x_all), "overlap_excluded_self_pairs": x_overlap},
        "atom_indices_zero_based": {"X": x_indices.tolist(), "Y": y_indices.tolist()},
        "reference_group_scope": "top_layer_metal" if args.x_scope == "top-layer" else "all_surface_and_subsurface_metal" if element_x.casefold() in {"cu", "pd"} else "all_atoms_of_element",
        "reference_group_selection": x_selection,
        "surface_species": args.surface_species,
        "top_layer_tolerance_A": args.top_layer_tolerance if args.x_scope == "top-layer" else None,
        "frames_analyzed": frames,
        "first_configuration": steps[0],
        "last_configuration": steps[-1],
        "frame_stride": args.frame_stride,
        "cell_dimensions_A_degrees": dimensions.tolist(),
        "cell_volume_A3": volume,
        "minimum_image_rmax_limit_A": inscribed_radius,
        "method": "MDAnalysis.lib.distances.distance_array with orthogonal/triclinic minimum-image box",
        "radial_bins": {"nbins": args.nbins, "range_A": [args.rmin, args.rmax], "bin_width_A": (args.rmax - args.rmin) / args.nbins},
        "normalization": "RDF = pair histogram / (frames * N_X * effective_N_Y_per_X * spherical_shell_volume / cell_volume); periodic full-cell normalization, no slab accessible-volume correction.",
        "coordination_definition": "Cumulative minimum-image pair count from rmin through each bin divided by analyzed frames and N_X; CN is neighbors per reference atom.",
        "coordination_cutoffs": named_cn,
        "first_rdf_peak_candidates": peak_candidates,
        "first_rdf_minimum_candidates": minimum_candidates,
        "block_frames": args.block_frames,
        "block_count": len(block_steps),
        "uncertainty_block_count": sum(count == args.block_frames for count in block_frame_counts),
        "partial_block_frames_excluded_from_uncertainty": sum(count for count in block_frame_counts if count != args.block_frames),
        "block_uncertainty_method": "SEM across equal-length complete blocks; a shorter trailing block remains in the mean and block CSV but is excluded from SEM.",
        "origin_export": {
            "enabled": args.origin is not None,
            "requested_bin_width_A": args.origin if isinstance(args.origin, (float, int)) and not isinstance(args.origin, bool) and args.origin > 0 else None,
            "rebin_method": "uniform overlap distribution of detailed-bin pair counts" if isinstance(args.origin, (float, int)) and not isinstance(args.origin, bool) and args.origin > 0 else "detailed RDF bins",
        },
        "scope_note": "For slab cells containing vacuum this full-cell RDF need not approach one at large r.",
        "output_files": [detailed_name, blocks_name, *([origin_name] if origin_name else []), metadata_name, analysis_metadata_name],
    }
    serialized_metadata = json.dumps(metadata, indent=2, sort_keys=True)
    (output_dir / metadata_name).write_text(serialized_metadata, encoding="utf-8")
    (output_dir / analysis_metadata_name).write_text(serialized_metadata, encoding="utf-8")
    return metadata


def add_rdf_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--xdatcar", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", "--output_dir", type=Path, required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--element-x", required=True)
    parser.add_argument("--element-y", required=True)
    parser.add_argument("--nbins", type=int, default=210)
    parser.add_argument("--rmin", type=float, default=0.0, metavar="A")
    parser.add_argument("--rmax", type=float, default=4.2, metavar="A")
    parser.add_argument("--frame-stride", type=int, default=1)
    parser.add_argument("--block-frames", type=int, default=500)
    parser.add_argument("--origin", nargs="?", const=-1.0, default=None, type=float, metavar="WIDTH_A", help="Write compact r_A,RDF,CN CSV; optionally rebin to this radial width")
    parser.add_argument("--overwrite", action="store_true", help="Allow reuse of a non-empty output directory")
    parser.add_argument("--cn-cutoff", type=float, action="append", default=[DEFAULT_CN_CUTOFF_A], metavar="A", help=f"Named CN cutoff; may be repeated (default: {DEFAULT_CN_CUTOFF_A:g} A)")
    parser.add_argument("--x-scope", choices=("all", "top-layer"), default="all", help="Reference element group scope")
    parser.add_argument("--surface-species", default=None, metavar="SPECIES[,SPECIES...]")
    parser.add_argument("--top-layer-tolerance", type=float, default=0.45, metavar="A")
    parser.add_argument("--data-root", type=Path, default=None)
