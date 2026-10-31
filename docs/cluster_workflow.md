# Cluster Bader Workflow

The scheduler is not required for local analysis. On a cluster, prepare a
portable job tree and manifest from a user-supplied data root:

```bash
python scripts/prepare_bader_jobs.py \
  --root /path/to/data_root --output-dir bader_jobs
```

The default static grid is 5500 through 9000 by 500. Existing `ACF.dat` files
are preserved unless `--force` is supplied. The manifest records case, step,
source path relative to the root, job path relative to the output directory,
and preparation status. Scheduler submission remains site-specific; submit
the rows selected by your cluster wrapper and record its command output.
The output tree mirrors the analyzer's expected `system/bader[_e]/bader_SYSTEM_STEP`
layout and includes the matching ensemble POSCAR plus the neutral surface
reference POSCAR and the matching case POTCAR. After jobs complete, use the staging directory itself as
`--root` for Bader aggregation.

Before submitting jobs, validate all eight source cases into a separate report
directory. This is read-only with respect to the source tree:

```bash
python scripts/validate_cluster_inputs.py \
  --data-root /path/to/data_root \
  --report-dir ./validation/inputs \
  --start 5001 --end 9000 \
  --bader-steps 5500 6000 6500 7000 7500 8000 8500 9000
```

The validator checks every `XDATCAR/POSCAR` pair, ordered POTCAR species and
`ZVAL`, `K=7.000` for each `aimd1_e` case, the two-K/four-water composition
change, all Bader snapshot `POSCAR`/`ACF.dat` pairs, and the configured
`CuxPd1/aimd1` updated-neutral variant. To prove the exact updated neutral
trajectory and matching RDF/hydrogen-bond/Bader inputs, place an
`analysis_input_manifest.json` in that case directory with the trajectory
hashes and `trajectory_variant: "updated_neutral"`; without this manifest the
case is reported as incomplete and production validation returns a nonzero
status.

Minimal manifest example (paths are relative to `CuxPd1/aimd1/`):

```json
{
  "trajectory_variant": "updated_neutral",
  "xdatcar_sha256": "...",
  "poscar_sha256": "...",
  "matching_inputs": {
    "rdf": ["POSCAR", "XDATCAR"],
    "hydrogen_bonds": ["POSCAR", "XDATCAR"],
    "bader": ["POTCAR"]
  }
}
```

Select jobs explicitly before submission. `all` and `none` are supported, as
are space/comma-separated table numbers. `--dry-run` prints the selection and
does not execute anything:

```bash
python scripts/run_bader_jobs.py --manifest bader_jobs/bader_manifest.csv \
  --select all --dry-run
```

An approved site wrapper can be supplied as a list of arguments. The
`{job_dir}` placeholder is replaced with the selected prepared directory and
no shell is invoked. Its Bader reference-density preparation must match the
site's VASP/Bader toolchain:

```bash
python scripts/run_bader_jobs.py --manifest bader_jobs/bader_manifest.csv \
  --select 1,2 --command ./site_bader_wrapper.sh {job_dir} --yes
```

After Bader completes, aggregate a case with:

```bash
python -m aimd_analysis.cli bader \
  --root bader_jobs --system Cu --ensemble aimd1 \
  --output-dir results/Cu/aimd1/bader \
  --start-step 5500 --end-step 9000 --step-stride 500
```

For `aimd1_e`, the analyzer reads and checks the case POTCAR automatically:

```bash
python -m aimd_analysis.cli bader \
  --root /path/to/data_root --system Cu --ensemble aimd1_e \
  --output-dir results/Cu/aimd1_e/bader \
  --potcar /path/to/data_root/Cu/aimd1_e/POTCAR
```

`--valence ELEMENT=Z` remains an explicit override for exceptional staging
cases and is rejected when it disagrees with POTCAR. Bader output reports atom
charge and `sigma_uC_cm-2`. PZC, capacitance, and potential response are not
assumed; the potential columns remain blank unless an external workflow has
supplied them.

Package compact results with:

```bash
python scripts/package_results.py --root results --output aimd-analysis-results.tgz
sha256sum -c aimd-analysis-results.tgz.sha256
tar -xzf aimd-analysis-results.tgz
```

Use `--include-snapshots` only when selected `ACF.dat` and matching POSCAR
files are needed. Large VASP files remain excluded. To download from a site
specific cluster host, use either:

```bash
scp user@cluster-host:/path/to/aimd-analysis-results.tgz* .
rsync -avP user@cluster-host:/path/to/aimd-analysis-results.tgz* .
```

Replace `cluster-host` with the hostname, IP, or SSH config alias provided by
your cluster. Verify the checksum before extraction.
