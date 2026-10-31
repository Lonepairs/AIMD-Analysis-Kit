from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from aimd_analysis.potcar import PotcarError, parse_potcar, validate_potcar_species


POTCAR = """TITEL  = PAW_PBE Cu 06Sep2000
ZVAL   = 11.000
TITEL  = PAW_PBE Pd 06Sep2000
ZVAL   = 10.000
TITEL  = PAW_PBE O 08Apr2002
ZVAL   = 6.000
TITEL  = PAW_PBE H 15Jun2001
ZVAL   = 1.000
"""


class PotcarTests(unittest.TestCase):
    def test_reads_ordered_zval(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "POTCAR"
            path.write_text(POTCAR, encoding="ascii")
            self.assertEqual(parse_potcar(path), {"Cu": 11.0, "Pd": 10.0, "O": 6.0, "H": 1.0})
            result = validate_potcar_species(path, ("Cu", "Pd", "O", "H"))
            self.assertEqual(result["zval"]["Pd"], 10.0)

    def test_rejects_species_order_mismatch(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "POTCAR"
            path.write_text(POTCAR, encoding="ascii")
            with self.assertRaisesRegex(PotcarError, "does not match"):
                validate_potcar_species(path, ("Pd", "Cu", "O", "H"))


if __name__ == "__main__":
    unittest.main()
