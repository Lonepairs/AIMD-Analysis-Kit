"""Portable numbered and config-driven workflow dispatch."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Iterable


CHOICES = {
    "1": "density",
    "2": "polar-angle",
    "3": "hydrogen-bonds",
    "4": "rdf",
    "5": "bader",
}
CHOICE_LABELS = {
    "1": "water density",
    "2": "polar angle and rho*cos(phi)",
    "3": "hydrogen bonds and donor counts",
    "4": "RDF and coordination number",
    "5": "Bader charge and potential response",
}
REPEATABLE_OPTIONS = {"valence", "cn_cutoff", "case_output", "case"}
PATH_OPTIONS = {"xdatcar", "reference", "output_dir", "root", "data_root", "surface_reference", "surface_metadata", "potcar"}


def load_workflow(path: Path) -> dict[str, Any]:
    """Load a JSON or optional YAML workflow file."""
    source = path.expanduser()
    if not source.is_file():
        raise FileNotFoundError(f"workflow config not found: {source}")
    text = source.read_text(encoding="utf-8")
    if source.suffix.casefold() in {".yaml", ".yml"}:
        try:
            import yaml
        except ImportError as exc:
            raise ValueError("YAML workflow requires the optional cluster dependency PyYAML") from exc
        config = yaml.safe_load(text)
    else:
        config = json.loads(text)
    if not isinstance(config, dict):
        raise ValueError("workflow root must be an object")
    jobs = config.get("jobs")
    matrix = config.get("case_matrix")
    if jobs is not None and not isinstance(jobs, list):
        raise ValueError("workflow jobs must be a list")
    if matrix is not None and not isinstance(matrix, dict):
        raise ValueError("workflow case_matrix must be an object")
    if jobs is None and matrix is None:
        raise ValueError("workflow must contain a jobs list or case_matrix")
    return config


def _format_value(value: Any, fields: dict[str, Any]) -> Any:
    if isinstance(value, str):
        if value.startswith("{") and value.endswith("}"):
            field_name = value[1:-1]
            if field_name in fields and not isinstance(fields[field_name], str):
                return fields[field_name]
        return value.format_map(fields)
    if isinstance(value, list):
        return [_format_value(item, fields) for item in value]
    if isinstance(value, dict):
        return {key: _format_value(item, fields) for key, item in value.items()}
    return value


def _expand_case_matrix(matrix: dict[str, Any]) -> list[dict[str, Any]]:
    cases = matrix.get("cases")
    templates = matrix.get("templates")
    if not isinstance(cases, list) or not cases or not isinstance(templates, list) or not templates:
        raise ValueError("case_matrix requires non-empty cases and templates lists")
    jobs: list[dict[str, Any]] = []
    seen_cases: set[tuple[str, str]] = set()
    for case_number, case in enumerate(cases, start=1):
        if not isinstance(case, dict) or not case.get("system") or not case.get("ensemble"):
            raise ValueError(f"case_matrix case {case_number} requires system and ensemble")
        identity = (str(case["system"]), str(case["ensemble"]))
        if identity in seen_cases:
            raise ValueError(f"duplicate case_matrix case: {identity[0]}/{identity[1]}")
        seen_cases.add(identity)
        for template_number, template in enumerate(templates, start=1):
            if not isinstance(template, dict):
                raise ValueError(f"case_matrix template {template_number} must be an object")
            pair_field = template.get("expand_pairs")
            pair_values = case.get(pair_field, [None]) if pair_field else [None]
            if not isinstance(pair_values, list) or not pair_values:
                raise ValueError(f"case {identity[0]}/{identity[1]} has no values for {pair_field}")
            for pair in pair_values:
                fields = {**case, "pair": ""}
                if pair is not None:
                    if not isinstance(pair, dict) or not pair.get("x") or not pair.get("y"):
                        raise ValueError("expanded RDF pair entries require x and y")
                    fields.update({"element_x": pair["x"], "element_y": pair["y"], "x_scope": pair.get("x_scope", "all")})
                    fields["pair"] = f"{pair['x']}-{pair['y']}"
                name = str(_format_value(template.get("name", f"case-{case_number}-job-{template_number}"), fields))
                job = {
                    key: _format_value(value, fields)
                    for key, value in template.items()
                    if key not in {"name", "expand_pairs"}
                }
                job["name"] = name
                jobs.append(job)
    return jobs


def _resolve_path_arguments(argv: list[str], base_dir: Path) -> list[str]:
    """Resolve relative path-valued CLI options from the workflow file location."""
    resolved = list(argv)
    index = 1
    while index < len(resolved):
        token = resolved[index]
        if not token.startswith("--"):
            index += 1
            continue
        option, separator, attached = token.partition("=")
        key = option[2:].replace("-", "_")
        if key not in PATH_OPTIONS:
            index += 1
            continue
        if separator:
            raw_path = Path(attached).expanduser()
            if not raw_path.is_absolute():
                resolved[index] = f"{option}={(base_dir / raw_path).resolve()}"
            index += 1
            continue
        if index + 1 < len(resolved):
            raw_path = Path(resolved[index + 1]).expanduser()
            if not raw_path.is_absolute():
                resolved[index + 1] = str((base_dir / raw_path).resolve())
            index += 2
            continue
        index += 1
    return resolved


def _choice_for_job(job: dict[str, Any]) -> str:
    raw = str(job.get("choice", job.get("analysis", ""))).casefold()
    if raw in CHOICES:
        return raw
    for choice, name in CHOICES.items():
        if raw == name:
            return choice
    raise ValueError(f"invalid workflow analysis choice: {raw!r}")


def _option_name(key: str) -> str:
    return "--" + key.replace("_", "-")


def _value_tokens(value: Any) -> list[str]:
    if isinstance(value, bool):
        return [] if not value else ["__FLAG__"]
    if isinstance(value, (list, tuple)):
        tokens: list[str] = []
        for item in value:
            tokens.extend(_value_tokens(item))
        return tokens
    if value is None:
        return []
    return [str(value)]


def job_argv(job: dict[str, Any]) -> tuple[str, list[str]]:
    """Convert a config job into the named CLI command and argument vector."""
    choice = _choice_for_job(job)
    command = CHOICES[choice]
    values = dict(job.get("args", {}))
    values.update({key: value for key, value in job.items() if key not in {"choice", "analysis", "args", "name"}})
    argv = [command]
    for key, value in values.items():
        if value is None or value is False:
            continue
        flag = _option_name(key)
        if key in REPEATABLE_OPTIONS and isinstance(value, (list, tuple)):
            for item in value:
                argv.extend((flag, str(item)))
            continue
        tokens = _value_tokens(value)
        if tokens == ["__FLAG__"]:
            argv.append(flag)
        elif tokens:
            argv.append(flag)
            argv.extend(tokens)
    return choice, argv


def _write_status(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["job", "choice", "analysis", "status", "message"]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def run_workflow(
    config_path: Path,
    main_callable,
    *,
    selection: str = "all",
    status_path: Path | None = None,
    dry_run: bool = False,
    stop_on_error: bool = True,
) -> list[dict[str, Any]]:
    """Run selected config jobs through the package CLI callable."""
    config = load_workflow(config_path)
    jobs = config.get("jobs")
    if jobs is None:
        jobs = _expand_case_matrix(config["case_matrix"])
    if selection.casefold() == "none":
        return []
    selected_choices = _parse_analysis_selection(selection)
    selected_indices = [
        index for index, job in enumerate(jobs)
        if selected_choices is None or _choice_for_job(job) in selected_choices
    ]
    if selected_choices is not None and not selected_indices:
        raise ValueError("selected analyses have no matching jobs in the workflow")
    planned: list[tuple[str, list[str]]] = []
    from .cli import build_parser

    command_parser = build_parser()
    for index in selected_indices:
        job = jobs[index]
        if not isinstance(job, dict):
            raise ValueError(f"workflow job {index + 1} must be an object")
        choice, argv = job_argv(job)
        argv = _resolve_path_arguments(argv, config_path.expanduser().resolve().parent)
        parsed = command_parser.parse_args(argv)
        _preflight_inputs(choice, parsed, index)
        planned.append((choice, argv))
    rows: list[dict[str, Any]] = []
    for index, (choice, argv) in zip(selected_indices, planned):
        job = jobs[index]
        row = {"job": job.get("name", f"job-{index + 1}"), "choice": choice, "analysis": CHOICE_LABELS[choice], "status": "PLANNED", "message": ""}
        rows.append(row)
        if dry_run:
            row["message"] = " ".join(argv)
            continue
        try:
            exit_code = main_callable(argv)
            if exit_code:
                raise RuntimeError(f"command returned exit code {exit_code}")
            row["status"] = "OK"
        except Exception as exc:  # noqa: BLE001 - status table must retain the job failure
            row["status"] = "FAILED"
            row["message"] = str(exc)
            if stop_on_error:
                break
    if status_path is not None:
        _write_status(Path(status_path), rows)
    return rows


def _parse_analysis_selection(selection: str) -> set[str] | None:
    value = selection.strip().casefold()
    if value == "all":
        return None
    if value == "none":
        return set()
    selected: set[str] = set()
    for token in value.replace(",", " ").split():
        if token not in CHOICES:
            raise ValueError("analysis selection must be 1-5, all, or none")
        selected.add(token)
    return selected


def _preflight_inputs(choice: str, args: argparse.Namespace, index: int) -> None:
    if choice in {"1", "2", "3", "4"}:
        trajectory = Path(args.xdatcar).expanduser()
        if not trajectory.is_file():
            raise FileNotFoundError(f"workflow job {index + 1}: XDATCAR not found: {trajectory}")
    if choice in {"1", "2", "3", "4"} and getattr(args, "reference", None) is not None:
        reference = Path(args.reference).expanduser()
        if not reference.is_file():
            raise FileNotFoundError(f"workflow job {index + 1}: reference structure not found: {reference}")
    if choice == "4":
        reference = Path(args.reference).expanduser()
        if not reference.is_file():
            raise FileNotFoundError(f"workflow job {index + 1}: reference structure not found: {reference}")
    if choice == "5" and not Path(args.root).expanduser().is_dir():
        raise FileNotFoundError(f"workflow job {index + 1}: Bader root not found: {args.root}")


def interactive_selection() -> str:
    print("1  water density")
    print("2  polar angle and rho*cos(phi)")
    print("3  hydrogen bonds and donor counts")
    print("4  RDF and coordination number")
    print("5  Bader charge and potential response")
    print("all run all five analyses")
    print("none exit without running")
    return input("Select analyses: ").strip()
