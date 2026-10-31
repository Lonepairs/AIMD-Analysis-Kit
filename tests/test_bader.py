from __future__ import annotations

import argparse
import csv
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from aimd_analysis.bader import analyze_bader, parse_acf


def _poscar() -> str:
    return (
        "bader fixture\n1.0\n10 0 0\n0 10 0\n0 0 20\n"
        "Cu O H\n2 1 2\nDirect\n"
        "0.10 0.10 0.90\n0.20 0.20 0.90\n0.30 0.30 0.30\n"
        "0.31 0.30 0.30\n0.29 0.30 0.30\n"
    )


def _acf(populations: list[float]) -> str:
    rows = [" # X Y Z CHARGE MIN DIST ATOMIC VOL", " -----------------------------"]
    rows.extend(f"{index} 0 0 0 {value:.6f} 0.5 1.0" for index, value in enumerate(populations, start=1))
    return "\n".join(rows) + "\nVACUUM CHARGE: 0\n"


class BaderTests(unittest.TestCase):
    def test_parse_acf_rejects_noncontiguous_indices(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "ACF.dat"
            path.write_text(_acf([10.0, 5.0, 1.0, 1.0]).replace("2 0", "3 0", 1), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "contiguous"):
                parse_acf(path)

    def test_bader_charge_and_potential_conversion(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            case = root / "Cu" / "aimd1"
            snapshot = root / "Cu" / "bader" / "bader_Cu_5500"
            case.mkdir(parents=True)
            snapshot.mkdir(parents=True)
            (case / "POSCAR").write_text(_poscar(), encoding="utf-8")
            (snapshot / "POSCAR").write_text(_poscar(), encoding="utf-8")
            (snapshot / "ACF.dat").write_text(_acf([10.5, 10.5, 5.5, 0.5, 0.5]), encoding="utf-8")
            output = root / "results"
            args = argparse.Namespace(
                root=root,
                system="Cu",
                ensemble="aimd1",
                output_dir=output,
                start_step=5500,
                end_step=5500,
                step_stride=500,
                reference=None,
                surface_reference=None,
                surface_metadata=None,
                surface_species=None,
                top_layer_tolerance=0.45,
                bader_dir=None,
                valence=[],
            )
            metadata = analyze_bader(args)
            self.assertEqual(metadata["top_layer_atom_indices_zero_based"], [0, 1])
            self.assertAlmostEqual(metadata["step_grid"][0], 5500)
            with (output / "bader_surface_metal_charge.csv").open(newline="", encoding="utf-8") as stream:
                row = next(csv.DictReader(stream))
            self.assertAlmostEqual(float(row["net_surface_metal_charge_q"]), 1.0)
            self.assertAlmostEqual(float(row["sigma_uC_cm-2"]), 17.698, places=2)
            with (output / "bader_surface_metal_charge_summary.csv").open(newline="", encoding="utf-8") as stream:
                summary = next(csv.DictReader(stream))
            self.assertEqual(summary["snapshot_std_net_surface_metal_charge_q"], "")
            self.assertEqual(summary["snapshot_sem_net_surface_metal_charge_q"], "")
            encoded = json.dumps(json.loads((output / "analysis_metadata.json").read_text(encoding="utf-8")))
            self.assertNotIn(str(root), encoded)

    def test_k_case_requires_explicit_valence(self) -> None:
        with self.assertRaisesRegex(ValueError, "K-containing"):
            analyze_bader(
                argparse.Namespace(
                    root=Path("."), system="Cu", ensemble="aimd1_e", output_dir=Path("out"),
                    start_step=5500, end_step=5500, step_stride=500, reference=None,
                    surface_reference=None, surface_metadata=None, surface_species=None,
                    top_layer_tolerance=0.45, bader_dir=None, valence=[],
                )
            )


if __name__ == "__main__":
    unittest.main()
