from __future__ import annotations

import argparse
import csv
import json
import unittest
from pathlib import Path

import numpy as np

from aimd_analysis.rdf import _complete_block_sem, _rebin_histogram_counts, analyze_rdf, shortest_lattice_translation


class RDFTests(unittest.TestCase):
    def test_shortest_lattice_translation(self) -> None:
        lattice = np.asarray([[10.0, 0.0, 0.0], [4.0, 9.0, 0.0], [0.0, 0.0, 12.0]])
        self.assertAlmostEqual(shortest_lattice_translation(lattice), 9.8488578018, places=8)

    def test_origin_rebin_conserves_pair_counts_with_partial_final_bin(self) -> None:
        edges = np.asarray([0.0, 0.25, 0.5, 0.75, 1.0])
        counts = np.asarray([2.0, 3.0, 5.0, 7.0])
        target_edges, rebinned = _rebin_histogram_counts(counts, edges, 0.3)

        self.assertAlmostEqual(float(np.sum(rebinned)), float(np.sum(counts)))
        self.assertAlmostEqual(float(target_edges[-1]), 1.0)
        self.assertAlmostEqual(float(target_edges[-1] - target_edges[-2]), 0.1)

    def test_uncertainty_uses_only_complete_equal_length_blocks(self) -> None:
        samples = [np.asarray([1.0]), np.asarray([3.0]), np.asarray([100.0])]
        sem = _complete_block_sem(samples, [500, 500, 12], 500)

        self.assertAlmostEqual(float(sem[0]), 1.0)

    def test_rdf_cn_monotonic_and_origin_csv(self) -> None:
        output = Path("_smoke") / "rdf-test"
        args = argparse.Namespace(
            xdatcar=Path("examples/minimal_case/XDATCAR"),
            reference=Path("examples/minimal_case/POSCAR"),
            output_dir=output,
            start=1,
            end=3,
            element_x="O",
            element_y="O",
            nbins=30,
            rmin=0.0,
            rmax=4.5,
            frame_stride=1,
            block_frames=2,
            origin=True,
            cn_cutoff=[3.0],
            x_scope="all",
            surface_species=None,
            top_layer_tolerance=0.45,
            data_root=Path("examples"),
            overwrite=True,
        )
        metadata = analyze_rdf(args)
        self.assertEqual(metadata["atom_counts"]["X"], 2)
        self.assertAlmostEqual(metadata["coordination_cutoffs"]["3"]["CN_per_reference_atom"], 1.0)
        with (output / "rdf_O-O_origin.csv").open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            self.assertEqual(reader.fieldnames, ["r_A", "RDF", "CN"])
            cn = [float(row["CN"]) for row in reader]
        self.assertTrue(all(right >= left for left, right in zip(cn, cn[1:])))
        encoded = json.dumps(json.loads((output / "rdf_metadata.json").read_text(encoding="utf-8")))
        self.assertNotIn(str(Path("examples").resolve()), encoded)


if __name__ == "__main__":
    unittest.main()
