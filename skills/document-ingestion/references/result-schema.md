# Result Manifest Schema

Each ingestion run should write one JSON manifest. Paths may be absolute or relative to the manifest directory; relative paths are preferred for portable results.

```json
{
  "schema_version": "1",
  "source": {
    "path": "input/report.pdf",
    "format": "pdf",
    "sha256": "optional"
  },
  "backend": {
    "name": "mineru",
    "mode": "markdown",
    "tier": "standard",
    "version": "optional"
  },
  "parameters": {
    "pages": "optional",
    "language": "optional"
  },
  "artifacts": {
    "markdown": "report.md",
    "structured": "report.json",
    "archive": "optional.zip"
  },
  "counts": {
    "pages": 12,
    "slides": 0,
    "text_blocks": 0,
    "tables": 0,
    "images": 0
  },
  "validation": {
    "passed": true,
    "checks": [],
    "warnings": []
  }
}
```

Required fields are `schema_version`, `source.path`, `source.format`, `backend.name`, `backend.mode`, and `artifacts`. A validator may accept a manifest before the producer has filled `validation`; it must then calculate and report its own result.

Artifact values must point to existing files or directories. A mode of `markdown` requires a Markdown artifact; `middle_json` requires a structured JSON artifact; `zip` requires an archive artifact. Counts are optional, but when present they must be non-negative integers.

The validation report should identify:

- missing source or artifact paths;
- unsupported or inconsistent format/mode combinations;
- empty Markdown or structured JSON files;
- invalid JSON manifests or structured outputs;
- declared counts that are negative or not integers;
- warnings that do not prevent downstream reading.
