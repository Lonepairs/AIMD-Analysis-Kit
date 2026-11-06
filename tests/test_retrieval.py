from __future__ import annotations

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from aimd_analysis.retrieval import build_index, format_results, load_index, search_index, tokenize


class RetrievalTests(unittest.TestCase):
    def test_tokenize_keeps_identifiers_and_chinese_terms_searchable(self) -> None:
        tokens = tokenize("rho_H2O*cos_phi 氢键")
        self.assertIn("rho_h2o", tokens)
        self.assertIn("cos_phi", tokens)
        self.assertIn("氢", tokens)
        self.assertIn("键", tokens)

    def test_build_excludes_generated_directories_and_searches_chunks(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "docs").mkdir()
            (root / "_staging").mkdir()
            (root / "docs" / "methods.md").write_text(
                "Hydrogen bond donor waters\nLayer I uses density boundaries.\n", encoding="utf-8"
            )
            (root / "_staging" / "secret.md").write_text("hydrogen bond should be excluded", encoding="utf-8")
            index_path = root / ".retrieval" / "index.json"
            metadata = build_index(root, index_path, chunk_lines=2, overlap_lines=0)
            self.assertEqual(metadata["document_count"], 1)
            index = load_index(index_path)
            results = search_index(index, "hydrogen bond Layer I", top_k=3)
            self.assertEqual(results[0].path, "docs/methods.md")
            self.assertIn("Hydrogen bond", results[0].snippet)
            self.assertNotIn("secret.md", format_results(results))

    def test_json_output_is_machine_readable(self) -> None:
        index = {
            "schema_version": 1,
            "chunks": [
                {
                    "path": "README.md",
                    "start_line": 1,
                    "end_line": 1,
                    "text": "water density",
                    "term_frequency": {"water": 1, "density": 1},
                    "length": 2,
                }
            ],
            "document_frequency": {"water": 1, "density": 1},
        }
        output = format_results(search_index(index, "water", top_k=1), "json")
        self.assertEqual(json.loads(output)[0]["path"], "README.md")


if __name__ == "__main__":
    unittest.main()
