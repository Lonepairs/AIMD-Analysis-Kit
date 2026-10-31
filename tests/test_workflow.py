from __future__ import annotations

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from aimd_analysis.workflow import _expand_case_matrix, _parse_analysis_selection, job_argv, load_workflow, run_workflow
from aimd_analysis.cli import build_parser


class WorkflowTests(unittest.TestCase):
    def test_selection_supports_all_none_and_lists(self) -> None:
        self.assertEqual(_parse_analysis_selection("1, 3"), {"1", "3"})
        self.assertIsNone(_parse_analysis_selection("all"))
        self.assertEqual(_parse_analysis_selection("none"), set())
        with self.assertRaises(ValueError):
            _parse_analysis_selection("6")

    def test_job_argv_uses_named_command_and_flags(self) -> None:
        choice, argv = job_argv({"choice": "1", "args": {"xdatcar": "in/XDATCAR", "origin": 0.3, "dry": False}})
        self.assertEqual(choice, "1")
        self.assertEqual(argv[:2], ["density", "--xdatcar"])
        self.assertIn("--origin", argv)
        _, bader_argv = job_argv({"choice": "5", "args": {"root": "data", "system": "Cu", "ensemble": "aimd1_e", "output_dir": "out", "valence": ["K=7", "O=6"]}})
        self.assertEqual(bader_argv.count("--valence"), 2)

    def test_none_does_not_write_status(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "workflow.json"
            config.write_text(json.dumps({"jobs": [{"choice": "1", "args": {"xdatcar": "missing"}}]}), encoding="utf-8")
            status = root / "status.csv"
            rows = run_workflow(config, lambda argv: 0, selection="none", status_path=status)
            self.assertEqual(rows, [])
            self.assertFalse(status.exists())

    def test_dry_run_writes_planned_status(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "workflow.json"
            input_file = root / "XDATCAR"
            input_file.write_text("fixture", encoding="utf-8")
            config.write_text(json.dumps({"jobs": [{"name": "density", "choice": "1", "args": {"xdatcar": str(input_file), "output_dir": "out", "start": 1, "end": 1}}]}), encoding="utf-8")
            status = root / "status.csv"
            rows = run_workflow(config, lambda argv: 0, selection="all", status_path=status, dry_run=True)
            self.assertEqual(rows[0]["status"], "PLANNED")
            self.assertTrue(status.is_file())
            self.assertIn("density", status.read_text(encoding="utf-8"))

    def test_relative_paths_resolve_from_config_directory(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            input_file = root / "inputs" / "XDATCAR"
            input_file.parent.mkdir()
            input_file.write_text("fixture", encoding="utf-8")
            config = root / "workflows" / "workflow.json"
            config.parent.mkdir()
            config.write_text(
                json.dumps({"jobs": [{"choice": "1", "args": {
                    "xdatcar": "../inputs/XDATCAR", "output_dir": "../results/density",
                    "start": 1, "end": 1,
                }}]}),
                encoding="utf-8",
            )
            rows = run_workflow(config, lambda argv: 0, dry_run=True)

            self.assertEqual(rows[0]["status"], "PLANNED")
            self.assertIn(str(input_file.resolve()), rows[0]["message"])
            self.assertIn(str((root / "results" / "density").resolve()), rows[0]["message"])

    def test_case_matrix_expands_cases_and_rdf_pairs(self) -> None:
        jobs = _expand_case_matrix({
            "cases": [
                {"system": "Cu", "ensemble": "aimd1", "rdf_pairs": [{"x": "O", "y": "O"}, {"x": "Cu", "y": "O", "x_scope": "top-layer"}]},
                {"system": "Pd", "ensemble": "aimd1_e", "rdf_pairs": [{"x": "O", "y": "O"}]},
            ],
            "templates": [
                {"name": "{system}-{ensemble}", "choice": "1", "args": {"output_dir": "results/{system}/{ensemble}"}},
                {"name": "{system}-{pair}", "choice": "4", "expand_pairs": "rdf_pairs", "args": {"element_x": "{element_x}", "element_y": "{element_y}", "x_scope": "{x_scope}"}},
            ],
        })

        self.assertEqual([job["name"] for job in jobs], ["Cu-aimd1", "Cu-O-O", "Cu-Cu-O", "Pd-aimd1_e", "Pd-O-O"])
        self.assertEqual(jobs[2]["args"]["x_scope"], "top-layer")

    def test_eight_case_config_expands_to_parseable_jobs(self) -> None:
        config_path = Path("workflows/example_all_cases.json")
        config = load_workflow(config_path)
        jobs = _expand_case_matrix(config["case_matrix"])
        parser = build_parser()
        parsed_jobs = [parser.parse_args(job_argv(job)[1]) for job in jobs]

        self.assertEqual(len(config["case_matrix"]["cases"]), 8)
        self.assertEqual({namespace.command for namespace in parsed_jobs}, {
            "density", "polar-angle", "hydrogen-bonds", "rdf", "bader",
        })
        self.assertEqual(sum(namespace.command == "rdf" for namespace in parsed_jobs), 40)
        self.assertEqual(sum(namespace.command == "bader" for namespace in parsed_jobs), 8)

    def test_analysis_selection_runs_all_matching_configured_jobs(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            input_file = root / "XDATCAR"
            input_file.write_text("fixture", encoding="utf-8")
            config = root / "workflow.json"
            jobs = [
                {"name": "density-a", "choice": "1", "args": {"xdatcar": str(input_file), "output_dir": str(root / "density-a"), "start": 1, "end": 1}},
                {"name": "polar-a", "choice": "2", "args": {"xdatcar": str(input_file), "output_dir": str(root / "polar-a"), "start": 1, "end": 1}},
                {"name": "density-b", "choice": "1", "args": {"xdatcar": str(input_file), "output_dir": str(root / "density-b"), "start": 1, "end": 1}},
            ]
            config.write_text(json.dumps({"jobs": jobs}), encoding="utf-8")
            received = []
            rows = run_workflow(config, lambda argv: received.append(argv) or 0, selection="1", dry_run=True)
            self.assertEqual([row["job"] for row in rows], ["density-a", "density-b"])


if __name__ == "__main__":
    unittest.main()
