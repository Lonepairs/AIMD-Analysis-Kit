from __future__ import annotations

import argparse
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from aimd_analysis.density import (
    analyze_density,
    find_local_minima,
    normalized_slice_contributions,
    rebin_slice_values,
    select_surface_indices,
)


def _xdatcar() -> str:
    header = (
        "minimal density fixture\n1.0\n"
        "10 0 0\n0 10 0\n0 0 20\n"
        "Cu O H K\n2 2 4 1\n"
    )
    atoms = (
        "0.1 0.1 0.10\n0.2 0.2 0.10\n"
        "0.3 0.3 0.30\n0.4 0.4 0.40\n"
        "0.30 0.30 0.30\n0.31 0.30 0.30\n"
        "0.39 0.40 0.40\n0.40 0.40 0.40\n0.7 0.7 0.20\n"
    )
    return header + "".join(
        f"Direct configuration= {step}\n{atoms}" for step in (1, 2, 3)
    )


class DensityTests(unittest.TestCase):
    def test_each_oxygen_contributes_one_equivalent(self) -> None:
        edges = np.arange(0.0, 2.1, 0.1)
        hits, equivalents, denominators = normalized_slice_contributions(
            np.asarray([0.55, 1.25]), edges, 0.2
        )
        self.assertEqual(hits.shape, (20,))
        self.assertTrue(np.all(denominators >= 1))
        self.assertAlmostEqual(float(equivalents.sum()), 2.0)

    def test_rebinning_conserves_profile(self) -> None:
        source_edges = np.arange(0.0, 1.01, 0.1)
        values = np.arange(1.0, 11.0)
        _, rebinned = rebin_slice_values(values, source_edges, 0.3)
        self.assertAlmostEqual(float(rebinned.sum()), float(values.sum()))

    def test_surface_selection_excludes_k(self) -> None:
        indices, names = select_surface_indices(("Cu", "Pd", "O", "H", "K"))
        self.assertEqual(indices.tolist(), [0, 1])
        self.assertEqual(names, ["Cu", "Pd"])

    def test_minima_excludes_endpoints_and_constant_tails(self) -> None:
        minima = find_local_minima(np.arange(6.0), np.asarray([0.0, 2.0, 1.0, 2.0, 1.0, 1.0]))
        self.assertEqual([item["index"] for item in minima], [2])

    def test_density_writes_portable_outputs_and_conserves_oxygen(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            xdatcar = root / "case" / "XDATCAR"
            output = root / "results"
            xdatcar.parent.mkdir()
            xdatcar.write_text(_xdatcar(), encoding="utf-8")
            args = argparse.Namespace(
                xdatcar=xdatcar,
                output_dir=output,
                start=1,
                end=3,
                stride=1,
                oxygen_radius=0.2,
                slice_width=0.1,
                top_layer_tolerance=0.45,
                surface_species=None,
                origin=0.3,
                zmax=None,
                block_frames=2,
                data_root=root,
            )
            metadata = analyze_density(args)
            self.assertAlmostEqual(metadata["integrated_oxygen_equivalent_count"], 2.0)
            self.assertAlmostEqual(metadata["integrated_oxygen_equivalent_error"], 0.0)
            self.assertTrue((output / "water_density_0.3A_origin.csv").is_file())
            self.assertTrue((output / "analysis_metadata.json").is_file())
            encoded = json.dumps(json.loads((output / "analysis_metadata.json").read_text()))
            self.assertNotIn(str(root), encoded)
            header = (output / "average_density_delta_z.csv").read_text().splitlines()[0]
            self.assertEqual(header, "delta_z_A,average_oxygen_mass_density_g_cm-3,average_oxygen_mass_density_g_A-3,average_water_mass_density_g_cm-3,average_water_mass_density_g_A-3")


if __name__ == "__main__":
    unittest.main()
