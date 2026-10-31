from __future__ import annotations

import argparse
import csv
import json
import unittest
from pathlib import Path

from aimd_analysis.polar_angle import analyze_polar_angle, assign_hydrogen_parents
from aimd_analysis.trajectory import read_xdatcar_header, iter_xdatcar


class PolarAngleTests(unittest.TestCase):
    def _args(self, output: Path, boundaries=None) -> argparse.Namespace:
        return argparse.Namespace(
            xdatcar=Path("examples/minimal_case/XDATCAR"),
            reference=Path("examples/minimal_case/POSCAR"),
            output_dir=output,
            start=1,
            end=3,
            stride=1,
            parent_oh_max=1.35,
            oh_min=0.50,
            oh_max=1.25,
            alpha_bins=36,
            block_frames=2,
            profile_slice_width=0.2,
            origin=5.0,
            boundaries=boundaries,
            top_layer_tolerance=0.45,
            surface_species=None,
            data_root=Path("examples"),
            overwrite=True,
        )

    def test_parent_assignment_uses_two_unique_hydrogens(self) -> None:
        structure = read_xdatcar_header("examples/minimal_case/XDATCAR")
        frame = next(iter_xdatcar("examples/minimal_case/XDATCAR", 1, 1))
        oxygen = [index for index, label in enumerate(structure.atom_species) if label == "O"]
        hydrogen = [index for index, label in enumerate(structure.atom_species) if label == "H"]
        mapping = assign_hydrogen_parents(frame.positions_frac, structure.lattice, oxygen, hydrogen, 1.35)
        self.assertEqual(sorted(sum((list(values) for values in mapping.values()), [])), hydrogen)
        self.assertTrue(all(len(values) == 2 for values in mapping.values()))

    def test_unlayered_run_has_no_layer_summary(self) -> None:
        output = Path("_smoke") / "polar-test"
        metadata = analyze_polar_angle(self._args(output))
        self.assertEqual(metadata["valid_water_observations"], 6)
        self.assertFalse(metadata["layer_analysis_enabled"])
        self.assertFalse((output / "polar_angle_layer_summary.csv").exists())
        with (output / "polar_angle_distribution.csv").open(newline="", encoding="utf-8") as stream:
            header = next(csv.reader(stream))
        self.assertIn("isotropic_reference", header)

    def test_boundaries_emit_layer_summary(self) -> None:
        output = Path("_smoke") / "polar-layer-test"
        analyze_polar_angle(self._args(output, boundaries=(1.0, 3.0)))
        text = (output / "polar_angle_layer_summary.csv").read_text(encoding="utf-8")
        self.assertIn("Layer_III", text)
        metadata = json.loads((output / "analysis_metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["layer_boundaries_A"], [1.0, 3.0])


if __name__ == "__main__":
    unittest.main()
