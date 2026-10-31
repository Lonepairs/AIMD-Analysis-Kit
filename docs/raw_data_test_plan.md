# Raw-data integration test

The repository includes a local integration workflow at
`workflows/raw_data_test.json`. It exercises the five named analyses against
one real neutral trajectory (`Cu/aimd1`) and one real Bader snapshot.

The expected local input layout is:

```text
data/
  Cu/
    aimd1/
      POSCAR
      XDATCAR
    bader/
      bader_Cu_5500/
        POSCAR
        ACF.dat
```

The raw files are intentionally ignored by Git. Keep them in the local
checkout or provision them through a separate data-transfer step; do not add
production trajectories or Bader snapshots to the source repository.

Run the integration test from the repository root:

```bash
python scripts/run_batch.py \
  --config workflows/raw_data_test.json \
  --select all \
  --status _smoke/raw-data/status.csv
```

The workflow reads trajectory frames 5001-5003 for density, polar angle,
hydrogen bonds, and O-O RDF. Bader reads the single 5500 snapshot. Its RDF
`rmax` is 4.2 A because the Cu cell's half-shortest lattice translation is
about 4.4271 A; a larger value must be rejected by the periodic-distance
guard.

The test is successful only when all five jobs report `OK` and the outputs
contain the expected invariants:

- density preserves 50 oxygen equivalents across the three selected frames;
- polar-angle reports 150 valid water observations;
- hydrogen-bond output contains the three configured regions;
- RDF writes the O-O origin table and metadata for three frames;
- Bader processes one snapshot and identifies the Cu top layer.

This test complements, rather than replaces, the deterministic unit suite.
The unit suite checks geometry, parsing, selection, conservation, error
handling, and workflow expansion. After the raw-data integration test passes,
run `scripts/validate_cluster_inputs.py` on the selected production cases.
That production check additionally requires complete frame ranges, POTCAR
ZVAL records, and the configured Bader snapshot grid.
