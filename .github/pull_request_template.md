## Change

Describe the user-visible or scientific-method change.

## Scientific definition impact

- [ ] No scientific definition changed.
- [ ] If a definition changed, `README.md` and `docs/scientific_definitions.md` were updated.
- [ ] Cutoffs, traversal, layer boundaries, normalization, and output schema are stated.

## Validation

- [ ] `python -m unittest discover -s tests -v`
- [ ] `python -m compileall -q src scripts tests`
- [ ] `python scripts/run_batch.py --config workflows/smoke_test.json --select all --status _smoke/ci_status.csv`
- [ ] Public diff checked for credentials, private data, local paths, and large files.

## Data and compatibility

- [ ] No production trajectories or large VASP files were added.
- [ ] Output metadata records the method and input provenance.
- [ ] Python 3.10 compatibility was considered.
