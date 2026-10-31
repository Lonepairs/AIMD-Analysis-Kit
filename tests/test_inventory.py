from __future__ import annotations

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from aimd_analysis.inventory import (
    InventoryError,
    inspect_xdatcar,
    inventory_cases,
    parse_acf,
    parse_poscar,
    write_inventory,
)


def _poscar(composition: str = "Cu O H He\n1 1 2 1") -> str:
    return (
        "fixture\n1.0\n"
        "10 0 0\n0 10 0\n0 0 20\n"
        f"{composition}\n"
        "Direct\n"
        "0 0 0\n0.1 0.1 0.1\n0.2 0.1 0.1\n0.3 0.1 0.1\n0.4 0.1 0.1\n"
    )


def _xdatcar(configurations: tuple[int, ...] = (1, 2, 3)) -> str:
    header = "fixture\n1.0\n10 0 0\n0 10 0\n0 0 20\nCu O H He\n1 1 2 1\n"
    atoms = "0 0 0\n0.1 0.1 0.1\n0.2 0.1 0.1\n0.3 0.1 0.1\n0.4 0.1 0.1\n"
    return header + "".join(f"Direct configuration= {frame}\n{atoms}" for frame in configurations)


def _acf(atom_count: int = 5) -> str:
    rows = [" # X Y Z CHARGE MIN DIST ATOMIC VOL", " -----------------------------"]
    rows.extend(f"{i} 0 0 0 1.0 0.5 1.0" for i in range(1, atom_count + 1))
    rows.append("VACUUM CHARGE: 0.0")
    return "\n".join(rows) + "\n"


class InventoryTests(unittest.TestCase):
    def test_poscar_and_xdatcar_parse_species_and_frames(self) -> None:
        with self.subTest("parsers"):
            with TemporaryDirectory() as directory:
                tmp_path = Path(directory)
                poscar = tmp_path / "POSCAR"
                xdatcar = tmp_path / "XDATCAR"
                poscar.write_text(_poscar(), encoding="utf-8")
                xdatcar.write_text(_xdatcar(), encoding="utf-8")

                structure = parse_poscar(poscar)
                trajectory = inspect_xdatcar(xdatcar, 1, 3)

                self.assertEqual(structure.composition, {"Cu": 1, "O": 1, "H": 2, "He": 1})
                self.assertEqual(structure.atom_count, 5)
                self.assertEqual(trajectory["frame_count"], 3)
                self.assertTrue(trajectory["requested_range_complete"])

    def test_xdatcar_reports_missing_requested_frames(self) -> None:
        with TemporaryDirectory() as directory:
            xdatcar = Path(directory) / "XDATCAR"
            xdatcar.write_text(_xdatcar((1, 3)), encoding="utf-8")
            result = inspect_xdatcar(xdatcar, 1, 3)

        self.assertFalse(result["requested_range_complete"])
        self.assertEqual(result["missing_requested_ranges"], [[2, 2]])

    def test_xdatcar_rejects_nonmonotonic_frame_numbers(self) -> None:
        with TemporaryDirectory() as directory:
            xdatcar = Path(directory) / "XDATCAR"
            xdatcar.write_text(_xdatcar((1, 3, 2)), encoding="utf-8")
            with self.assertRaisesRegex(InventoryError, "strictly increasing"):
                inspect_xdatcar(xdatcar, 1, 3)

    def test_acf_requires_contiguous_indices(self) -> None:
        with TemporaryDirectory() as directory:
            acf = Path(directory) / "ACF.dat"
            acf.write_text(_acf().replace("3 0", "4 0", 1), encoding="utf-8")
            with self.assertRaisesRegex(InventoryError, "not contiguous"):
                parse_acf(acf)

    def test_inventory_records_only_relative_paths(self) -> None:
        with TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            case_dir = tmp_path / "Cu" / "aimd1"
            snapshot = tmp_path / "Cu" / "bader" / "bader_Cu_5500"
            case_dir.mkdir(parents=True)
            snapshot.mkdir(parents=True)
            (case_dir / "POSCAR").write_text(_poscar(), encoding="utf-8")
            (case_dir / "XDATCAR").write_text(_xdatcar(), encoding="utf-8")
            (snapshot / "POSCAR").write_text(_poscar(), encoding="utf-8")
            (snapshot / "ACF.dat").write_text(_acf(), encoding="utf-8")

            report = inventory_cases(tmp_path, ["Cu/aimd1"], 1, 3, [5500])
            output_paths = write_inventory(report, tmp_path / "reports")
            manifest = json.loads(output_paths[0].read_text(encoding="utf-8"))
            encoded = json.dumps(manifest)

            self.assertEqual(manifest["cases"][0]["status"], "OK")
            self.assertEqual(manifest["cases"][0]["xdatcar_path"], "Cu/aimd1/XDATCAR")
            self.assertNotIn(str(tmp_path), encoded)

    def test_inventory_flags_bader_atom_count_mismatch(self) -> None:
        with TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            case_dir = tmp_path / "Cu" / "aimd1"
            snapshot = tmp_path / "Cu" / "bader" / "bader_Cu_5500"
            case_dir.mkdir(parents=True)
            snapshot.mkdir(parents=True)
            (case_dir / "POSCAR").write_text(_poscar(), encoding="utf-8")
            (case_dir / "XDATCAR").write_text(_xdatcar(), encoding="utf-8")
            (snapshot / "POSCAR").write_text(_poscar(), encoding="utf-8")
            (snapshot / "ACF.dat").write_text(_acf(4), encoding="utf-8")

            report = inventory_cases(tmp_path, ["Cu/aimd1"], 1, 3, [5500])

            self.assertEqual(report["cases"][0]["status"], "INCOMPLETE")
            self.assertEqual(report["bader_snapshots"][0]["status"], "INVALID")
            self.assertIn(
                "ACF has 4 atoms but POSCAR has 5",
                report["bader_snapshots"][0]["message"],
            )


if __name__ == "__main__":
    unittest.main()
