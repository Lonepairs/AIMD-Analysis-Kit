"""Hydrogen-bond triplets and unique donor-water statistics."""

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
from .io import prepare_output_dir
from .density import _git_revision, _portable_input_identifier, _sha256, cell_geometry, select_surface_indices, surface_relative_heights
from .polar_angle import assign_hydrogen_parents, minimum_image_vectors
from .structure import parse_poscar
from .trajectory import iter_xdatcar, read_xdatcar_header


def _write_csv(path: Path, header: list[str], rows: Iterable[Iterable[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def find_hydrogen_bonds(
    positions_frac: np.ndarray,
    lattice: np.ndarray,
    oxygen_indices: np.ndarray,
    parent_mapping: dict[int, tuple[int, int]],
    selected_donors: set[int] | None = None,
    oo_cutoff: float = 3.5,
    angle_max_deg: float = 30.0,
) -> list[tuple[int, int, int, float, float]]:
    """Return one best acceptor triplet per qualifying donor H."""
    selected = set(int(value) for value in oxygen_indices) if selected_donors is None else set(selected_donors)
    cos_cutoff = math.cos(math.radians(angle_max_deg))
    triplets: list[tuple[int, int, int, float, float]] = []
    for donor in sorted(selected):
        donor_position = positions_frac[donor]
        acceptor_vectors = minimum_image_vectors(
            lattice, positions_frac[oxygen_indices] - donor_position
        )
        acceptor_distances = np.linalg.norm(acceptor_vectors, axis=1)
        candidate_mask = (acceptor_distances < oo_cutoff) & (oxygen_indices != donor)
        if not np.any(candidate_mask):
            continue
        for hydrogen in parent_mapping[donor]:
            donor_h_vector = minimum_image_vectors(
                lattice, positions_frac[hydrogen] - donor_position
            )
            oh_distance = float(np.linalg.norm(donor_h_vector))
            if oh_distance <= 0:
                continue
            denominator = acceptor_distances * oh_distance
            usable = candidate_mask & (denominator > 0)
            cosine = np.zeros_like(acceptor_distances)
            cosine[usable] = np.einsum("ij,j->i", acceptor_vectors[usable], donor_h_vector) / denominator[usable]
            qualified = usable & (cosine > cos_cutoff)
            candidates: list[tuple[float, float, int, int]] = []
            for position in np.flatnonzero(qualified):
                angle = math.degrees(math.acos(float(np.clip(cosine[position], -1.0, 1.0))))
                acceptor = int(oxygen_indices[position])
                candidates.append((angle, float(acceptor_distances[position]), acceptor, int(hydrogen)))
            if candidates:
                angle, distance, acceptor, hydrogen_index = min(candidates, key=lambda item: (item[0], item[1], item[2]))
                triplets.append((donor, hydrogen_index, acceptor, distance, angle))
    return triplets


def _regions(boundaries: tuple[float, float] | None, zmin: float | None, zmax: float | None, period: float) -> list[tuple[str, float, float | None]]:
    if boundaries is not None:
        first, second = boundaries
        if not (0 < first < second):
            raise ValueError("boundaries must satisfy 0 < first < second")
        return [("Layer I", 0.0, first), ("Layer II", first, second), ("Layer III", second, None)]
    if zmin is None:
        raise ValueError("provide --boundaries or --zmin")
    if zmin < 0 or zmin > period or (zmax is not None and zmax <= zmin) or (zmax is not None and zmax > period):
        raise ValueError("z bounds must be ordered and lie within the cell period")
    return [("Selected region", zmin, zmax)]


def _sem(values: np.ndarray) -> float:
    return float(np.std(values, ddof=1) / math.sqrt(len(values))) if len(values) > 1 else 0.0


def _validate_args(args: argparse.Namespace) -> None:
    if args.start < 1 or args.end < args.start or args.stride < 1:
        raise ValueError("configuration range/stride is invalid")
    if args.oo_cutoff <= 0 or args.parent_oh_max <= 0 or args.angle_max_deg <= 0 or args.angle_max_deg > 180:
        raise ValueError("distance cutoffs must be positive and angle must be in (0, 180]")
    if args.block_frames < 1:
        raise ValueError("block-frames must be positive")
    if args.boundaries is None and args.zmin is None:
        raise ValueError("provide --boundaries or --zmin")
    if args.boundaries is not None and (args.zmin is not None or args.zmax is not None):
        raise ValueError("use --boundaries or --zmin/--zmax, not both")


def analyze_hydrogen_bonds(args: argparse.Namespace) -> dict:
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
    oxygen_indices = np.asarray([index for index, label in enumerate(labels) if label.casefold() == "o"], dtype=int)
    hydrogen_indices = np.asarray([index for index, label in enumerate(labels) if label.casefold() == "h"], dtype=int)
    if oxygen_indices.size == 0 or hydrogen_indices.size == 0:
        raise ValueError("XDATCAR must contain atoms labelled O and H")
    surface_indices, surface_species = select_surface_indices(labels, args.surface_species)
    normal, period, area = cell_geometry(structure.lattice)
    regions = _regions(args.boundaries, args.zmin, args.zmax, period)
    parent_mapping: dict[int, tuple[int, int]] | None = None
    frame_rows: list[tuple[object, ...]] = []
    triplet_rows: list[tuple[object, ...]] = []
    frame_steps: list[int] = []
    region_frame_rows: dict[str, list[tuple[object, ...]]] = {name: [] for name, _, _ in regions}
    per_frame_profile: list[dict[str, np.ndarray | int]] = []

    for frame in iter_xdatcar(trajectory, args.start, args.end, args.stride):
        if parent_mapping is None:
            parent_mapping = assign_hydrogen_parents(
                frame.positions_frac, structure.lattice, oxygen_indices, hydrogen_indices, args.parent_oh_max
            )
        _, oxygen_heights, _ = surface_relative_heights(
            frame, structure, normal, period, surface_indices, oxygen_indices, args.top_layer_tolerance
        )
        height_by_oxygen = {int(oxygen): float(height) for oxygen, height in zip(oxygen_indices, oxygen_heights)}
        all_triplets = find_hydrogen_bonds(
            frame.positions_frac, structure.lattice, oxygen_indices, parent_mapping,
            selected_donors=set(int(value) for value in oxygen_indices),
            oo_cutoff=args.oo_cutoff, angle_max_deg=args.angle_max_deg,
        )
        triplets_by_donor: dict[int, list[tuple[int, int, int, float, float]]] = {}
        for triplet in all_triplets:
            triplets_by_donor.setdefault(triplet[0], []).append(triplet)
            triplet_rows.append((frame.configuration, *triplet, height_by_oxygen[triplet[0]]))
        for region_name, lower, upper in regions:
            selected = {
                oxygen for oxygen, height in height_by_oxygen.items()
                if height >= lower and (upper is None or height < upper)
            }
            donor_counts = {oxygen: len(triplets_by_donor.get(oxygen, [])) for oxygen in selected}
            donor_waters = sum(count > 0 for count in donor_counts.values())
            donor_bonds = sum(donor_counts.values())
            water_count = len(selected)
            row = (
                frame.configuration, region_name, lower, "" if upper is None else upper,
                water_count, donor_waters, donor_bonds,
                donor_waters / water_count if water_count else "",
                donor_bonds / water_count if water_count else "",
                sum(count == 0 for count in donor_counts.values()),
                sum(count == 1 for count in donor_counts.values()),
                sum(count == 2 for count in donor_counts.values()),
            )
            frame_rows.append(row)
            region_frame_rows[region_name].append(row)
        frame_steps.append(frame.configuration)

    if not frame_steps or parent_mapping is None:
        raise ValueError(f"no XDATCAR configurations found in {args.start}-{args.end}")
    expected_steps = list(range(args.start, args.end + 1, args.stride))
    if frame_steps != expected_steps:
        missing = sorted(set(expected_steps) - set(frame_steps))
        raise ValueError(f"requested frame grid is incomplete; missing configurations: {missing[:10]}")
    output_dir = prepare_output_dir(output_dir, overwrite=getattr(args, "overwrite", False))
    _write_csv(
        output_dir / "hydrogen_bond_triplets.csv",
        ["configuration", "donor_oxygen_index_zero_based", "donor_hydrogen_index_zero_based", "acceptor_oxygen_index_zero_based", "oo_distance_A", "donor_angle_deg", "donor_delta_z_A"],
        triplet_rows,
    )
    _write_csv(
        output_dir / "hydrogen_bonds_per_frame.csv",
        ["configuration", "total_hydrogen_bonds"],
        ((step, sum(1 for row in triplet_rows if row[0] == step)) for step in frame_steps),
    )
    _write_csv(
        output_dir / "ndonor_waters_per_frame.csv",
        ["configuration", "region", "zmin_A", "zmax_A", "water_count", "N_donor_waters", "N_donor_hbonds", "donor_water_fraction", "donor_hbonds_per_water", "N_water_0_donor_hbonds", "N_water_1_donor_hbond", "N_water_2_donor_hbonds"],
        frame_rows,
    )
    summary_rows: list[tuple[object, ...]] = []
    block_rows: list[tuple[object, ...]] = []
    for region_name, lower, upper in regions:
        rows = region_frame_rows[region_name]
        donor_values = np.asarray([row[5] for row in rows], dtype=float)
        water_values = np.asarray([row[4] for row in rows], dtype=float)
        donor_bond_values = np.asarray([row[6] for row in rows], dtype=float)
        total_water = float(water_values.sum())
        summary_rows.append(
            (
                region_name, lower, "" if upper is None else upper, len(rows),
                donor_values.mean(), donor_values.std(ddof=1) if len(rows) > 1 else 0.0, _sem(donor_values),
                water_values.mean(), np.mean([row[7] for row in rows if row[7] != ""]) if total_water else "",
                donor_bond_values.mean(), donor_bond_values.sum() / total_water if total_water else "",
                sum(row[9] for row in rows) / total_water if total_water else "",
                sum(row[10] for row in rows) / total_water if total_water else "",
                sum(row[11] for row in rows) / total_water if total_water else "",
            )
        )
        for block_index, start in enumerate(range(0, len(rows), args.block_frames), start=1):
            block = rows[start : start + args.block_frames]
            block_water = sum(row[4] for row in block)
            values = np.asarray([row[5] for row in block], dtype=float)
            bonds = np.asarray([row[6] for row in block], dtype=float)
            block_rows.append(
                (
                    block_index, region_name, len(block), values.mean(), values.std(ddof=1) if len(values) > 1 else 0.0,
                    bonds.mean(), bonds.sum() / block_water if block_water else "",
                )
            )
    _write_csv(
        output_dir / "ndonor_waters_summary.csv",
        ["region", "zmin_A", "zmax_A", "frames", "mean_N_donor_waters", "std_N_donor_waters", "sem_N_donor_waters", "mean_water_count", "mean_donor_water_fraction", "mean_N_donor_hbonds", "mean_donor_hbonds_per_water", "fraction_0_donor_hbonds", "fraction_1_donor_hbond", "fraction_2_donor_hbonds"],
        summary_rows,
    )
    _write_csv(
        output_dir / "ndonor_waters_blocks.csv",
        ["block", "region", "frames", "mean_N_donor_waters", "std_N_donor_waters", "mean_N_donor_hbonds", "mean_donor_hbonds_per_water"],
        block_rows,
    )
    origin_rows = []
    for region_name, _, _ in regions:
        rows = region_frame_rows[region_name]
        total_water = sum(int(row[4]) for row in rows)
        total_hbonds = sum(int(row[6]) for row in rows)
        origin_rows.append(
            (
                region_name,
                total_hbonds / total_water if total_water else 0.0,
                sum(int(row[5]) for row in rows) / len(rows),
                total_hbonds / len(rows),
            )
        )
    _write_csv(
        output_dir / "ndonor_waters_origin.csv",
        ["region", "Hbond_per_water", "mean_N_donor_waters", "mean_N_donor_hbonds"],
        origin_rows,
    )

    input_identifier, redacted = _portable_input_identifier(trajectory, args.data_root)
    metadata = {
        "schema_version": 1,
        "analysis": "hydrogen bonds and unique donor-water statistics",
        "package_version": __version__,
        "git_revision": _git_revision(),
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "operating_system": platform.system()},
        "input": {"xdatcar": input_identifier, "absolute_path_redacted": redacted, "size_bytes": trajectory.stat().st_size, "sha256": _sha256(trajectory)},
        "reference_poscar": (
            {"path": _portable_input_identifier(Path(reference_path).expanduser(), args.data_root)[0], "absolute_path_redacted": _portable_input_identifier(Path(reference_path).expanduser(), args.data_root)[1], "size_bytes": Path(reference_path).expanduser().stat().st_size, "sha256": _sha256(Path(reference_path).expanduser())}
            if reference_path is not None else None
        ),
        "species": list(structure.species), "counts": list(structure.counts), "frames_analyzed": len(frame_steps),
        "first_configuration": frame_steps[0], "last_configuration": frame_steps[-1], "configuration_stride": args.stride,
        "surface_species": surface_species, "surface_normal": normal.tolist(), "cell_period_normal_A": period, "lateral_area_A2": area,
        "regions": [{"name": name, "zmin_A": lower, "zmax_A": upper} for name, lower, upper in regions],
        "region_rule": "Oxygen height is measured from the instantaneous mean top surface; intervals are [zmin, zmax), and an omitted upper bound is unbounded.",
        "criterion": {
            "oo_distance": f"minimum-image O_d-O_a distance < {args.oo_cutoff} A",
            "angle": f"angle O_a-O_d-H < {args.angle_max_deg} degrees",
            "parent_oh_max_A": args.parent_oh_max,
            "hydrogen_parent_assignment": "Global greedy minimum-distance O-H matching in the first selected frame; fixed for all selected frames.",
            "acceptor_selection": "Each donor H selects at most one acceptor, preferring smallest angle, then distance, then atom index.",
            "one_water_count_rule": "A water oxygen is counted once if at least one parent H forms a qualifying hydrogen bond.",
            "donor_hbond_count_rule": "Each water contributes 0, 1, or 2 donor H bonds according to its two parent H atoms.",
        },
        "output_files": ["hydrogen_bond_triplets.csv", "hydrogen_bonds_per_frame.csv", "ndonor_waters_per_frame.csv", "ndonor_waters_summary.csv", "ndonor_waters_blocks.csv", "ndonor_waters_origin.csv", "analysis_metadata.json"],
    }
    (output_dir / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    return metadata


def add_hydrogen_bond_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--xdatcar", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=None)
    parser.add_argument("--output-dir", "--output_dir", type=Path, required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--stride", "--frame-stride", dest="stride", type=int, default=1)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--boundaries", type=float, nargs=2, metavar=("ZMIN", "ZMAX"))
    group.add_argument("--zmin", "-zmin", type=float, metavar="A")
    parser.add_argument("--zmax", "-zmax", type=float, default=None, metavar="A")
    parser.add_argument("--surface-species", default=None, metavar="SPECIES[,SPECIES...]")
    parser.add_argument("--top-layer-tolerance", type=float, default=0.45, metavar="A")
    parser.add_argument("--oo-cutoff", type=float, default=3.5, metavar="A")
    parser.add_argument("--angle-max-deg", type=float, default=30.0)
    parser.add_argument("--parent-oh-max", type=float, default=1.35, metavar="A")
    parser.add_argument("--block-frames", type=int, default=500)
    parser.add_argument("--data-root", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true", help="Allow reuse of a non-empty output directory")
