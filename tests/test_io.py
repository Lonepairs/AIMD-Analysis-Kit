from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from aimd_analysis.io import prepare_output_dir


class OutputDirectoryTests(unittest.TestCase):
    def test_non_empty_output_requires_explicit_overwrite(self) -> None:
        with TemporaryDirectory() as directory:
            output_dir = Path(directory) / "results"
            prepare_output_dir(output_dir)
            marker = output_dir / "existing.csv"
            marker.write_text("preserve", encoding="utf-8")

            with self.assertRaisesRegex(FileExistsError, "--overwrite"):
                prepare_output_dir(output_dir)

            prepare_output_dir(output_dir, overwrite=True)
            self.assertEqual(marker.read_text(encoding="utf-8"), "preserve")


if __name__ == "__main__":
    unittest.main()
