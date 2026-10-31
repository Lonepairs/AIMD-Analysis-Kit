from __future__ import annotations

import argparse
import csv
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from aimd_analysis.bader_response import collect_bader_response


def _write_case(path: Path, ensemble: str, charge: float, phis: list[float]) -> None:
    path.mkdir(parents=True)
    mean_phi = sum(phis) / len(phis)
    (path / "bader_surface_metal_charge_summary.csv").write_text(
        "system,ensemble,n_snapshots,n_surface_metal_atoms,mean_net_surface_metal_charge_q,snapshot_std_net_surface_metal_charge_q,snapshot_sem_net_surface_metal_charge_q,mean_sigma_uC_cm-2,mean_delta_phi_V,mean_phi_V_vs_SHE\n"
        f"Cu,{ensemble},{len(phis)},2,{charge},0,0,{charge * 10},0,{mean_phi}\n",
        encoding="utf-8",
    )
    rows = "".join(
        f"Cu,{ensemble},{5500 + index * 500},2,{charge},{charge / 2},{charge * 10},0,{phi},{{}}\n"
        for index, phi in enumerate(phis)
    )
    (path / "bader_surface_metal_charge.csv").write_text(
        "system,ensemble,step,n_surface_metal_atoms,net_surface_metal_charge_q,mean_surface_metal_charge_q,sigma_uC_cm-2,delta_phi_V,phi_V_vs_SHE, surface_charge_by_element_json\n"
        + rows,
        encoding="utf-8",
    )


class BaderResponseTests(unittest.TestCase):
    def test_collects_paired_response(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            _write_case(root / "neutral", "aimd1", -0.2, [-0.11, 0.09])
            _write_case(root / "electro", "aimd1_e", -1.2, [-0.60, -0.40])
            output = root / "response"
            metadata = collect_bader_response(
                argparse.Namespace(
                    case_output=[f"Cu/aimd1={root / 'neutral'}", f"Cu/aimd1_e={root / 'electro'}"],
                    output_dir=output,
                    data_root=root,
                )
            )
            with (output / "bader_potential_response_summary.csv").open(newline="", encoding="utf-8") as stream:
                row = next(csv.DictReader(stream))
            self.assertAlmostEqual(float(row["potential_response_V"]), -0.49)
            self.assertAlmostEqual(float(row["potential_response_independent_sem_V"]), 2**0.5 * 0.1)
            self.assertEqual(metadata["output_files"][-1], "analysis_metadata.json")
            self.assertNotIn(str(root), json.dumps(metadata))

    def test_one_snapshot_reports_undefined_response_sem(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            _write_case(root / "neutral", "aimd1", -0.2, [-0.01])
            _write_case(root / "electro", "aimd1_e", -1.2, [-0.60])
            output = root / "response"
            collect_bader_response(argparse.Namespace(
                case_output=[f"Cu/aimd1={root / 'neutral'}", f"Cu/aimd1_e={root / 'electro'}"],
                output_dir=output,
                data_root=root,
            ))
            with (output / "bader_potential_response_summary.csv").open(newline="", encoding="utf-8") as stream:
                row = next(csv.DictReader(stream))
            self.assertEqual(row["potential_response_independent_sem_V"], "")


if __name__ == "__main__":
    unittest.main()
