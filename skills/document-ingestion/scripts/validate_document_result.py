#!/usr/bin/env python3
"""Validate a document-ingestion manifest using only the Python standard library."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


SUPPORTED_FORMATS = {
    "pdf",
    "png",
    "jpg",
    "jpeg",
    "webp",
    "bmp",
    "tif",
    "tiff",
    "docx",
    "pptx",
    "html",
    "htm",
    "md",
    "txt",
    "rst",
    "xlsx",
    "xls",
    "csv",
    "tsv",
    "doc",
    "ppt",
}

MODE_ARTIFACTS = {
    "markdown": "markdown",
    "middle_json": "structured",
    "structured": "structured",
    "json": "structured",
    "zip": "archive",
}


def add_check(checks: list[dict[str, Any]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": passed, "detail": detail})


def resolve_path(value: str, manifest_path: Path) -> Path:
    candidate = Path(value).expanduser()
    return candidate if candidate.is_absolute() else (manifest_path.parent / candidate)


def validate(manifest_path: Path) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    warnings: list[str] = []
    errors: list[str] = []

    try:
        raw = json.loads(manifest_path.read_text(encoding="utf-8"))
        add_check(checks, "manifest_json", True, "Manifest is valid JSON.")
    except FileNotFoundError:
        return {
            "status": "error",
            "manifest": str(manifest_path),
            "passed": False,
            "errors": [f"Manifest does not exist: {manifest_path}"],
            "checks": [],
            "warnings": [],
        }
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        return {
            "status": "error",
            "manifest": str(manifest_path),
            "passed": False,
            "errors": [f"Manifest cannot be read as UTF-8 JSON: {exc}"],
            "checks": [{"name": "manifest_json", "passed": False, "detail": str(exc)}],
            "warnings": [],
        }

    if not isinstance(raw, dict):
        errors.append("Manifest root must be a JSON object.")
        add_check(checks, "manifest_object", False, errors[-1])
        return _result(manifest_path, checks, errors, warnings)

    for field in ("schema_version", "source", "backend", "artifacts"):
        present = field in raw and raw[field] is not None
        add_check(checks, f"required_{field}", present, "Present." if present else f"Missing {field}.")
        if not present:
            errors.append(f"Missing required field: {field}")

    source = raw.get("source") if isinstance(raw.get("source"), dict) else {}
    backend = raw.get("backend") if isinstance(raw.get("backend"), dict) else {}
    artifacts = raw.get("artifacts") if isinstance(raw.get("artifacts"), dict) else {}

    backend_name = backend.get("name")
    backend_name_ok = isinstance(backend_name, str) and bool(backend_name.strip())
    add_check(checks, "backend_name", backend_name_ok, str(backend_name or "missing"))
    if not backend_name_ok:
        errors.append("backend.name must be a non-empty string.")

    source_path_value = source.get("path")
    if isinstance(source_path_value, str) and source_path_value.strip():
        source_path = resolve_path(source_path_value, manifest_path)
        exists = source_path.is_file()
        add_check(checks, "source_exists", exists, str(source_path))
        if not exists:
            errors.append(f"Source file does not exist: {source_path}")
        suffix = source_path.suffix.lower().lstrip(".")
        declared_format = str(source.get("format", suffix)).lower().lstrip(".")
        format_ok = declared_format in SUPPORTED_FORMATS
        add_check(checks, "source_format", format_ok, declared_format)
        if not format_ok:
            errors.append(f"Unsupported source format: {declared_format}")
        if declared_format and suffix and declared_format != suffix:
            warnings.append(f"Declared format {declared_format!r} differs from suffix {suffix!r}.")
    else:
        errors.append("source.path must be a non-empty string.")
        add_check(checks, "source_path", False, errors[-1])

    mode = str(backend.get("mode", "")).lower()
    mode_key = MODE_ARTIFACTS.get(mode)
    add_check(checks, "backend_mode", mode_key is not None, mode or "missing")
    if mode and mode_key is None:
        errors.append(f"Unsupported backend mode: {mode}")
    elif not mode:
        errors.append("backend.mode must be provided.")

    artifact_values = [value for value in artifacts.values() if isinstance(value, str) and value.strip()]
    if not artifact_values:
        errors.append("artifacts must contain at least one path.")
        add_check(checks, "artifact_paths", False, errors[-1])
    else:
        add_check(checks, "artifact_paths", True, f"{len(artifact_values)} path(s) declared.")
        for name, value in artifacts.items():
            if not isinstance(value, str) or not value.strip():
                errors.append(f"Artifact {name!r} must be a non-empty path string.")
                add_check(checks, f"artifact_{name}", False, "Missing or invalid path.")
                continue
            artifact_path = resolve_path(value, manifest_path)
            exists = artifact_path.exists()
            add_check(checks, f"artifact_{name}", exists, str(artifact_path))
            if not exists:
                errors.append(f"Declared artifact does not exist: {artifact_path}")

    if mode_key:
        required_artifact = artifacts.get(mode_key)
        if not isinstance(required_artifact, str) or not required_artifact.strip():
            errors.append(f"Mode {mode!r} requires artifacts.{mode_key}.")
            add_check(checks, "mode_artifact_declared", False, errors[-1])
        else:
            artifact_path = resolve_path(required_artifact, manifest_path)
            exists = artifact_path.exists()
            add_check(checks, "mode_artifact_exists", exists, str(artifact_path))
            if not exists:
                errors.append(f"Declared {mode_key} artifact does not exist: {artifact_path}")
            elif not artifact_path.is_file():
                errors.append(f"Declared {mode_key} artifact must be a file: {artifact_path}")
                add_check(checks, "mode_artifact_is_file", False, str(artifact_path))
            else:
                add_check(checks, "mode_artifact_is_file", True, str(artifact_path))
                if artifact_path.stat().st_size == 0:
                    errors.append(f"Declared {mode_key} artifact is empty: {artifact_path}")

            if mode_key == "structured" and exists and artifact_path.is_file():
                try:
                    json.loads(artifact_path.read_text(encoding="utf-8"))
                    add_check(checks, "structured_json", True, "Structured artifact is valid JSON.")
                except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                    errors.append(f"Structured artifact is not valid UTF-8 JSON: {exc}")
                    add_check(checks, "structured_json", False, str(exc))

    counts = raw.get("counts", {})
    if counts is not None and not isinstance(counts, dict):
        errors.append("counts must be a JSON object when present.")
        add_check(checks, "counts_object", False, errors[-1])
    elif isinstance(counts, dict):
        for name, value in counts.items():
            valid = isinstance(value, int) and not isinstance(value, bool) and value >= 0
            add_check(checks, f"count_{name}", valid, str(value))
            if not valid:
                errors.append(f"Count {name!r} must be a non-negative integer.")

    return _result(manifest_path, checks, errors, warnings)


def _result(
    manifest_path: Path,
    checks: list[dict[str, Any]],
    errors: list[str],
    warnings: list[str],
) -> dict[str, Any]:
    return {
        "status": "ok" if not errors else "error",
        "manifest": str(manifest_path),
        "passed": not errors,
        "errors": errors,
        "checks": checks,
        "warnings": warnings,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate a document-ingestion JSON manifest.")
    parser.add_argument("manifest", type=Path, help="Path to the JSON manifest")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print the JSON report")
    args = parser.parse_args()

    report = validate(args.manifest.resolve())
    print(json.dumps(report, ensure_ascii=True, indent=2 if args.pretty else None, sort_keys=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
