# Contributing

Thank you for contributing to AIMD Analysis Kit. Keep changes reproducible,
small, and tied to a documented analysis definition.

## Development environment

Python 3.10 or newer is required. The cluster-supported baseline is Python
3.10.20.

```bash
python -m pip install -e ".[dev,plot,cluster]"
python -m unittest discover -s tests -v
python -m compileall -q src scripts tests
python scripts/run_batch.py \
  --config workflows/smoke_test.json \
  --select all \
  --status _smoke/ci_status.csv
```

The smoke workflow uses synthetic or repository-provided lightweight fixtures.
Do not add production `XDATCAR`, `CHGCAR`, `AECCAR*`, `WAVECAR`, `POTCAR`, or
full Bader output to a pull request.

## Scientific changes

When changing a cutoff, denominator, traversal, layer rule, normalization, or
charge/potential conversion, update all of the following in the same pull
request:

1. implementation and command-line help;
2. `docs/scientific_definitions.md` and the relevant README section;
3. deterministic tests covering the changed behavior;
4. `CHANGELOG.md`;
5. output metadata and representative compact examples, when applicable.

Do not silently overwrite an existing production result. Use a new output
directory or method/version label and record the repository commit, command,
input identifiers, and dependency versions.

## Pull requests

Use a focused branch and open a pull request against `main`. Include the
motivation, changed scientific definition, input/output schema changes, and
the exact validation commands used. GitHub Actions must pass before merge.

Avoid unrelated formatting or generated-file changes. Review public diffs for
local paths, credentials, private data, and accidental large files.

## Release process

1. Move the completed `Unreleased` entries in `CHANGELOG.md` to a dated
   version section.
2. Update `project.version` in `pyproject.toml` and `version` in `CITATION.cff`.
3. Run the full test and smoke commands above.
4. Create an annotated tag such as `v0.2.0` and publish GitHub release notes.
5. State the exact commit, environment, and analysis definitions used for any
   published scientific result.
