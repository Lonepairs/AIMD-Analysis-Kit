# Troubleshooting

`XDATCAR not found` means the command received a path that does not exist in
the current environment. Supply `--xdatcar` explicitly; the kit does not use
research-workspace defaults.

If a reference POSCAR reports a species/count or lattice mismatch, confirm
that it belongs to the same trajectory and that atom order is unchanged. RDF,
polar-angle, and hydrogen-bond analyses use the reference order for atom
labels and parent assignments.

For Bader analysis, every requested snapshot needs both `POSCAR` and
`ACF.dat`. ACF rows must be contiguous from atom 1 and must match the snapshot
POSCAR order. The production grid is 5500, 6000, ..., 9000. Missing snapshots
are reported before aggregation. The analyzer reads `ZVAL` from the matching
POTCAR and checks `K=7.000` for `aimd1_e`; explicit valence overrides are
accepted only when they agree with POTCAR. PZC and capacitance are not
assumed, so potential fields are blank.

Python 3.10 or newer is required by the package. The original cluster Python
3.6.8 environment cannot install this repository unchanged; run the kit in a
Python 3.10+ environment or use the legacy scripts there. MDAnalysis is needed
only for RDF/CN and is installed by the core package dependency set.

The RDF implementation uses full-cell periodic normalization. A slab cell with
vacuum therefore need not approach `g(r)=1` at large distance. CN is cumulative
neighbors per selected reference atom, not a bond count. The default named CN
cutoff is 2.5 A and can be replaced or supplemented with `--cn-cutoff`.

If an output contains an unexpected nonzero first density slice, inspect
`analysis_metadata.json` for profile clipping, periodic wrapping, the selected
surface species, and the oxygen-ball radius. The first slice is not direct
evidence that water touches the metal.

For stale or old `CuxPd1` results, rerun the analysis from the current
trajectory and matching reference/Bader inputs. The kit records input SHA256
checksums and does not import legacy plotted images or private default paths.
