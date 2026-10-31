# AIMD Analysis Kit

Toolkit for reproducible analysis of ab initio molecular dynamics trajectories.

## Status

Package v0.1.0. The five named analyses, numbered compatibility interface,
config-driven batch runner, and paired Bader response collector are
implemented and covered by deterministic fixture tests.

## Installation

Python 3.10 or newer is required. Install the core trajectory-analysis
dependencies with:

```bash
python -m pip install .
```

Install plotting support with `python -m pip install '.[plot]'`. Development
tools are available with `python -m pip install '.[dev]'`.

The `aimd-analysis --help` command lists the available interfaces.
Use a relative or user-supplied path for each input and output so the same
command can run on a workstation or cluster.

## Quick Start

Each analysis has its own command and output directory:

```bash
aimd-analysis density --xdatcar DATA/Cu/aimd1/XDATCAR --output-dir results/Cu/aimd1/density --start 5001 --end 9000
aimd-analysis polar-angle --xdatcar DATA/Cu/aimd1/XDATCAR --output-dir results/Cu/aimd1/polar --start 5001 --end 9000
aimd-analysis hydrogen-bonds --xdatcar DATA/Cu/aimd1/XDATCAR --output-dir results/Cu/aimd1/hbond --start 5001 --end 9000 --boundaries 4.35 10.95
aimd-analysis rdf --xdatcar DATA/Cu/aimd1/XDATCAR --reference DATA/Cu/aimd1/POSCAR --output-dir results/Cu/aimd1/rdf/Cu-O --start 5001 --end 9000 --element-x Cu --element-y O --origin
aimd-analysis bader --root DATA --system Cu --ensemble aimd1 --output-dir results/Cu/aimd1/bader --potcar DATA/Cu/aimd1/POTCAR
```

Bader reads ordered `ZVAL` values from the case POTCAR. An explicit
`--valence ELEMENT=Z` is only an override and is checked against POTCAR;
`aimd1_e` validation requires `K=7.000`.

## Water density

The density command uses the validated discrete oxygen-ball/slice method. By
default it uses a 0.2 A oxygen radius, 0.1 A slices, and a 0.45 A tolerance
for selecting the instantaneous top substrate layer. Oxygen heights are
relative to that surface and wrapped into the normal cell period; K, O, H,
and He are excluded from surface selection by default.

```bash
python -m aimd_analysis.cli density \
  --xdatcar data/Cu/aimd1/XDATCAR --output-dir results/Cu/aimd1/density \
  --start 5001 --end 9000 --origin 0.3 --data-root data
```

The command writes detailed and compact CSVs, 500-frame block profiles,
surface diagnostics, local-minimum candidates, and `analysis_metadata.json`.
The metadata records the input checksum, frame grid, physical constants,
surface rule, clipping diagnostics, and oxygen-equivalent conservation check.
`--origin WIDTH` adds `water_density_<WIDTH>A_origin.csv` for plotting; it
does not change the analysis grid. Local-minimum candidates are diagnostics,
not automatic layer assignments. A nonzero first slice is not by itself
evidence that water touches the metal because the oxygen ball has finite
radius and the profile is periodic.

## Polar angle and orientation profile

The polar-angle command assigns two unique parent H atoms to each O in the
first selected frame, reuses that mapping, and accepts an orientation only
when both O-H distances are in `[0.50, 1.25) A`. The angle is measured between
the normalized O-H bisector and the outward surface normal. The distribution
CSV includes the normalized isotropic `sin(phi)` reference. The unlayered
orientation profile always includes water density, valid-orientation density,
mean `cos(phi)`, and `rho_H2O*cos(phi)`.

```bash
python -m aimd_analysis.cli polar-angle \
  --xdatcar data/Cu/aimd1/XDATCAR --reference data/Cu/aimd1/POSCAR \
  --output-dir results/Cu/aimd1/polar --start 5001 --end 9000 \
  --origin 5 --data-root data
```

Layer files are written only when density-derived boundaries are supplied
explicitly. The convention is Layer I `[0, zmin)`, Layer II `[zmin, zmax)`,
and Layer III `[zmax, infinity)`:

```bash
python -m aimd_analysis.cli polar-angle \
  --xdatcar data/Cu/aimd1/XDATCAR --reference data/Cu/aimd1/POSCAR \
  --output-dir results/Cu/aimd1/polar --start 5001 --end 9000 \
  --boundaries 4.35 10.95 --origin 5 --data-root data
```

## Hydrogen bonds and donor waters

Hydrogen-bond analysis requires explicit layer boundaries or one selected
region. It uses minimum-image `O_d-O_a < 3.5 A` and a strict donor angle
`O_a-O_d-H < 30 degrees`; each parent H selects at most one acceptor, by
smallest angle, then shortest distance, then atom index. Results distinguish
unique donor-water counts from donor-H-bond counts (0, 1, or 2 per water).

```bash
python -m aimd_analysis.cli hydrogen-bonds \
  --xdatcar data/Cu/aimd1/XDATCAR --reference data/Cu/aimd1/POSCAR \
  --output-dir results/Cu/aimd1/hydrogen-bonds --start 5001 --end 9000 \
  --boundaries 4.35 10.95 --data-root data
```

The detailed triplets, per-frame region counts, block summaries, aggregate
summary, and metadata are written separately. `N_donor_waters` counts unique
donor-water oxygens with at least one qualifying H; `N_donor_hbonds` counts
qualifying donor-H/acceptor triplets, so one water contributes 0, 1, or 2.
No legacy half-split is used.

## RDF and coordination number

RDF analysis requires a matching POSCAR/CONTCAR and explicit element labels.
It streams frames through MDAnalysis' periodic distance kernel, reports a
full-cell normalized RDF and cumulative CN per reference atom, and warns about
the half-shortest-lattice-translation limit. For slab cells with vacuum the
full-cell RDF need not approach one at large distance. The default X group
contains all atoms of the selected element; use `--x-scope top-layer` for a
top-layer subset selected from the reference POSCAR.

```bash
python -m aimd_analysis.cli rdf \
  --xdatcar data/Cu/aimd1/XDATCAR --reference data/Cu/aimd1/POSCAR \
  --output-dir results/Cu/aimd1/rdf/Cu-O --start 5001 --end 9000 \
  --element-x Cu --element-y O --origin --cn-cutoff 2.5 --data-root data
```

`--cn-cutoff` may be repeated. Detailed CSVs include block SEM for RDF and CN;
metadata records atom groups, normalization, named-cutoff CN values, and local
peak/minimum candidates. The Origin CSV has exactly `r_A,RDF,CN` columns.

## Numbered and batch workflows

The stable named commands are also available through the requested numbered
compatibility interface:

```bash
python scripts/run_analysis.py 1 --xdatcar data/Cu/aimd1/XDATCAR \
  --reference data/Cu/aimd1/POSCAR --output-dir results/Cu/aimd1/density \
  --start 5001 --end 9000
```

The choices are `1` density, `2` polar angle, `3` hydrogen bonds, `4` RDF/CN,
and `5` Bader. For repeatable automation, use a JSON or YAML config. JSON
needs no extra dependency; YAML requires the optional `cluster` extra.

Batch selection accepts analysis choices `1` through `5`, `all`, or `none`;
selecting a number runs that analysis for every matching case and RDF pair in
the workflow. `none` exits without creating a status file or partial outputs.
Relative input and output paths in a workflow are resolved from the workflow
file's directory, so the same config works from any current directory. Each
run writes a status CSV; `--continue-on-error` records later jobs after a
failure instead of stopping at the first failed job.
Analysis commands refuse to write into non-empty output directories unless
`--overwrite` is supplied. For production reruns, use a new results directory
or set `overwrite: true` in the selected workflow template after reviewing
the existing outputs.

The synthetic smoke config is separate from the eight-case production matrix:

```bash
python scripts/run_batch.py --config workflows/smoke_test.json --select all
```

For the production matrix, place inputs under `data/<system>/<ensemble>/`
using the layout in `docs/cluster_workflow.md`, then run a dry check before
starting analyses:

```bash
python scripts/run_batch.py --config workflows/example_all_cases.json --select all --dry-run
python scripts/run_batch.py --config workflows/example_all_cases.json --select all --continue-on-error
```

The matrix contains all eight system/ensemble cases, the recorded layer
boundaries, O-O/O-H RDFs, applicable top-layer metal-water pairs, and the
default eight-step Bader grid. Metal RDF jobs explicitly use
`x_scope=top-layer`. The `CuxPd1/aimd1` case is labeled as the updated neutral
trajectory; use its matching POSCAR, RDF, hydrogen-bond, and Bader inputs in
the staged data tree. Run `scripts/validate_cluster_inputs.py` before
production analysis.
Use `--select 1`, for example, to run density for all eight cases. The Bader
potential-response collector runs after both ensembles have completed.

## Input inventory

Before migrating analyses, validate the available case inputs and produce a
portable inventory report. The command checks POSCAR/XDATCAR species and atom
counts, streams configuration markers and coordinates, hashes input files,
and checks Bader snapshot coverage and ACF/POSCAR atom counts. Paths in the
report are relative to `--data-root`; the local root path is not stored.

```bash
aimd-analysis inventory --data-root /path/to/data_root \
  --output-dir ./inventory --start 5001 --end 9000
```

Use `--case SYSTEM/ENSEMBLE` one or more times to select cases. The default is
the four systems across `aimd1` and `aimd1_e`. The default Bader grid is 5500
through 9000 by 500; override it with `--bader-steps` when checking a smaller
fixture. The command returns exit code 2 if any case or required Bader snapshot
is incomplete, while still writing the diagnostic report.

For paired neutral/2K response tables, first run the per-case `bader` command
for both ensembles, then collect only those generated summaries:

```bash
python scripts/collect_bader_response.py \
  --case-output Cu/aimd1=results/Cu/aimd1/bader \
  --case-output Cu/aimd1_e=results/Cu/aimd1_e/bader \
  --output-dir results/bader-response
```

The collector writes snapshot response rows and surface-charge summaries
without recalculating Bader charges. PZC and capacitance are not provided, so
potential response columns are blank unless both input summaries explicitly
contain finite potential fields. Its optional uncertainty is propagated from
independent neutral and `aimd1_e` snapshot SEMs; matching step labels do not
make samples from distinct trajectories physically paired.

Cluster preparation, approved-command execution, packaging, checksum
verification, and `scp`/resumable `rsync` examples are documented in
[`docs/cluster_workflow.md`](docs/cluster_workflow.md).

## Skills

The repository includes a portable `document-ingestion` Skill for routing local
PDF, image, DOCX, PPTX, HTML, text, and table files to an appropriate parser.
It uses local MinerU when available, preserves a machine-readable result
manifest, and validates generated artifacts with a standard-library Python
script. The Skill is under `skills/document-ingestion/` and follows the open
Agent Skills directory layout.
