#!/usr/bin/env python3
"""Export a surface-relative H-bond-per-water profile for Origin.

Example
-------
python scripts/analyze_hbond_z_profile.py \
    --xdatcar data/Cu/aimd1/XDATCAR --reference data/Cu/aimd1/POSCAR \
    --output-dir _smoke/raw-data/Cu-aimd1/hbond-z-profile-5001-5003 \
    --start 5001 --end 5003 --bin-width 0.3 --surface-species Cu \
    --traversal acceptor
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aimd_analysis.density import cell_geometry, select_surface_indices, surface_relative_heights
from aimd_analysis.hydrogen_bonds import find_hydrogen_bonds, find_hydrogen_bonds_acceptor_centered
from aimd_analysis.io import prepare_output_dir
from aimd_analysis.polar_angle import assign_hydrogen_parents
from aimd_analysis.structure import parse_poscar
from aimd_analysis.trajectory import iter_xdatcar, read_xdatcar_header


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xdatcar", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=None, help="Optional POSCAR used to verify species/counts")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--bin-width", type=float, default=0.3, metavar="A")
    parser.add_argument("--surface-species", default=None, metavar="SPECIES[,SPECIES...]")
    parser.add_argument("--top-layer-tolerance", type=float, default=0.45, metavar="A")
    parser.add_argument("--oo-cutoff", type=float, default=3.5, metavar="A")
    parser.add_argument("--angle-max-deg", type=float, default=30.0)
    parser.add_argument("--parent-oh-max", type=float, default=1.35, metavar="A")
    parser.add_argument(
        "--traversal", choices=("donor", "acceptor"), default="acceptor",
        help="Count donor waters or unique acceptor-centered waters (default: acceptor)",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.start < 1 or args.end < args.start or args.stride < 1:
        parser.error("configuration range and stride must be positive and ordered")
    if args.bin_width <= 0 or args.top_layer_tolerance <= 0 or args.oo_cutoff <= 0 or args.parent_oh_max <= 0:
        parser.error("bin width and distance tolerances must be positive")
    if not 0 < args.angle_max_deg <= 180:
        parser.error("angle maximum must lie in (0, 180]")
    return args


def analyze(args: argparse.Namespace) -> dict[str, object]:
    xdatcar = args.xdatcar.expanduser().resolve()
    if not xdatcar.is_file():
        raise FileNotFoundError("XDATCAR not found: %s" % xdatcar)
    structure = read_xdatcar_header(xdatcar)
    if args.reference is not None:
        reference = parse_poscar(args.reference.expanduser())
        if reference.species != structure.species or reference.counts != structure.counts:
            raise ValueError("reference POSCAR species/counts do not match XDATCAR")

    labels = structure.atom_species
    oxygen_indices = np.asarray([i for i, label in enumerate(labels) if label.casefold() == "o"], dtype=int)
    hydrogen_indices = np.asarray([i for i, label in enumerate(labels) if label.casefold() == "h"], dtype=int)
    if oxygen_indices.size == 0 or hydrogen_indices.size == 0:
        raise ValueError("XDATCAR must contain atoms labelled O and H")
    surface_indices, surface_species = select_surface_indices(labels, args.surface_species)
    normal, period, _ = cell_geometry(structure.lattice)
    bin_count = int(math.ceil(period / args.bin_width - 1.0e-12))
    water_observations = np.zeros(bin_count, dtype=np.int64)
    hbond_observations = np.zeros(bin_count, dtype=np.int64)
    parent_mapping = None
    frame_steps: list[int] = []
    total_triplets = 0

    for frame in iter_xdatcar(xdatcar, args.start, args.end, args.stride):
        if args.traversal == "donor" and parent_mapping is None:
            parent_mapping = assign_hydrogen_parents(
                frame.positions_frac, structure.lattice, oxygen_indices, hydrogen_indices, args.parent_oh_max
            )

        _, heights, _ = surface_relative_heights(
            frame, structure, normal, period, surface_indices, oxygen_indices, args.top_layer_tolerance
        )
        bin_by_oxygen = np.minimum((heights / args.bin_width).astype(int), bin_count - 1)
        water_observations += np.bincount(bin_by_oxygen, minlength=bin_count)

        if args.traversal == "acceptor":
            triplets = find_hydrogen_bonds_acceptor_centered(
                frame.positions_frac,
                structure.lattice,
                oxygen_indices,
                hydrogen_indices,
                oo_cutoff=args.oo_cutoff,
                parent_oh_max=args.parent_oh_max,
                angle_max_deg=args.angle_max_deg,
            )
            center_column = 2
        else:
            assert parent_mapping is not None
            triplets = find_hydrogen_bonds(
                frame.positions_frac,
                structure.lattice,
                oxygen_indices,
                parent_mapping,
                oo_cutoff=args.oo_cutoff,
                angle_max_deg=args.angle_max_deg,
            )
            center_column = 0

        oxygen_bin_by_index = {int(oxygen): int(bin_by_oxygen[pos]) for pos, oxygen in enumerate(oxygen_indices)}
        hbond_bins = [oxygen_bin_by_index[triplet[center_column]] for triplet in triplets]
        hbond_observations += np.bincount(hbond_bins, minlength=bin_count)
        total_triplets += len(triplets)
        frame_steps.append(frame.configuration)

    expected_steps = list(range(args.start, args.end + 1, args.stride))
    if frame_steps != expected_steps:
        missing = sorted(set(expected_steps) - set(frame_steps))
        raise ValueError("requested configurations were not all found; missing: %s" % missing[:10])
    if not frame_steps:
        raise ValueError("no XDATCAR configurations were selected")
    if int(water_observations.sum()) != len(oxygen_indices) * len(frame_steps):
        raise RuntimeError("oxygen-to-bin assignment did not conserve water observations")
    if int(hbond_observations.sum()) != total_triplets:
        raise RuntimeError("z-binned H-bond count does not match the triplet count")

    edges = np.minimum(np.arange(bin_count + 1, dtype=float) * args.bin_width, period)
    edges[-1] = period
    centers = (edges[:-1] + edges[1:]) / 2.0
    hbond_per_water = np.divide(
        hbond_observations,
        water_observations,
        out=np.full(bin_count, np.nan),
        where=water_observations > 0,
    )
    cumulative_hbonds_per_frame = np.cumsum(hbond_observations) / len(frame_steps)

    output_dir = prepare_output_dir(args.output_dir.expanduser(), overwrite=args.overwrite)
    profile_path = output_dir / "hbond_z_profile.csv"
    with profile_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["delta_z", "Hbond_per_water", "total_Hbond", "water_observations", "hbond_count"])
        for center, ratio, total, waters, hbonds in zip(
            centers, hbond_per_water, cumulative_hbonds_per_frame, water_observations, hbond_observations
        ):
            writer.writerow((
                f"{center:.10g}",
                "" if not np.isfinite(ratio) else f"{ratio:.10g}",
                f"{total:.10g}",
                int(waters),
                int(hbonds),
            ))

    metadata = {
        "analysis": "surface-relative hydrogen bonds per water profile",
        "trajectory": str(xdatcar),
        "reference_poscar": str(args.reference.expanduser().resolve()) if args.reference else None,
        "species": list(structure.species),
        "counts": list(structure.counts),
        "surface_species": surface_species,
        "surface_definition": "instantaneous mean projected height of selected top-layer substrate atoms",
        "surface_normal": normal.tolist(),
        "cell_period_normal_A": period,
        "frames_analyzed": len(frame_steps),
        "first_configuration": frame_steps[0],
        "last_configuration": frame_steps[-1],
        "configuration_stride": args.stride,
        "oxygen_count": int(oxygen_indices.size),
        "water_observations": int(water_observations.sum()),
        "total_hydrogen_bonds": total_triplets,
        "global_hydrogen_bonds_per_water": total_triplets / int(water_observations.sum()),
        "traversal": args.traversal,
        "counted_oxygen": "acceptor O_a" if args.traversal == "acceptor" else "donor O_d",
        "criterion": {
            "oo_cutoff_A": args.oo_cutoff,
            "angle_max_deg_strict": args.angle_max_deg,
            "parent_oh_max_A": args.parent_oh_max,
            "acceptor_rule": (
                "For each O_a, inspect O_b within the O-O cutoff, choose the two nearest H atoms to each O_b, "
                "and retain at most one qualifying O_b-H witness per O_a."
                if args.traversal == "acceptor" else
                "Each donor H selects its best qualifying acceptor; each donor water contributes 0-2 bonds."
            ),
        },
        "bin_width_A": args.bin_width,
        "bin_count": bin_count,
        "last_bin_width_A": round(float(edges[-1] - edges[-2]), 10),
        "empty_bins_are_blank": True,
        "profile_definition": "sum of binned H-bond counts divided by sum of binned oxygen observations over selected frames",
        "output_csv": profile_path.name,
    }
    (output_dir / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def main() -> int:
    args = _parse_args()
    try:
        metadata = analyze(args)
    except (OSError, ValueError, IndexError, KeyError) as exc:
        print("Error: %s" % exc, file=sys.stderr)
        return 2
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
