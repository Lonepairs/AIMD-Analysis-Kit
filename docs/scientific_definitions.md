# Scientific Definitions

## Slice-based water density

For each trajectory frame, the surface normal is the unit vector along
`a x b`, oriented so that its projection on `c` is positive. The instantaneous
surface is the mean projected height of selected substrate atoms satisfying

```text
z_i >= max(z_surface_atoms) - 0.45 A
```

The default substrate selection excludes `O`, `H`, `He`, and `K`. Oxygen
heights are measured from this instantaneous surface and wrapped into
`[0, cell_period)`.

An oxygen is represented by a sphere of radius `0.2 A`. For a profile with
slice width `0.1 A`, every slice with nonzero geometric overlap receives
`1/n_i` oxygen equivalents, where `n_i` is the number of intersected slices
for that oxygen. Tangencies at a slice edge are excluded. Thus each oxygen
contributes one total oxygen equivalent over its represented slices. This is a
discrete intersection model, rather than a kernel-density estimate or an
analytic sphere-volume projection.

For lateral cell area `A`, slice width `Delta z`, and mean oxygen-equivalent
count `N_O(z)`, the reported densities are

```text
rho_O(z)   = N_O(z) * 15.999 / (N_A * A * Delta z)
rho_H2O(z) = N_O(z) * 18.015 / (N_A * A * Delta z)
```

The factor `1e24` is applied when writing `g cm-3` from Angstrom-based cell
dimensions. `--origin WIDTH` only rebins the validated profile for plotting;
it does not change the analysis grid.

The output metadata reports the integrated oxygen-equivalent count, expected
oxygen count, clipping at both profile edges, and periodic wrapping. Local
minima are interior sampled candidates for human boundary selection; they do
not automatically define physical water layers. A nonzero first slice can
come from the finite oxygen-ball radius or periodic mapping and is not, by
itself, evidence that water touches the metal.

## Polar angle and orientation-weighted density

Assign two unique parent H atoms to each O in the first selected frame by
global greedy minimum-distance matching with a maximum O-H distance of
`1.35 A`. Keep this parent map fixed during the selected trajectory. For each
frame, compute minimum-image O-to-H vectors, normalize each vector separately,
and sum them to form the water bisector. The polar angle `phi` is the angle
between this bisector and the outward surface normal, in `[0, 180]` degrees.
Only waters with both O-H lengths in `[0.50, 1.25) A` and a nonzero bisector
are included in orientation statistics.

The isotropic angle-bin reference is the normalized integral of `sin(phi)`:

```text
P_iso([phi_1, phi_2]) = [cos(phi_1) - cos(phi_2)] / 2
```

The unlayered profile includes all oxygen-centered water counts, the valid
orientation fraction, mean `cos(phi)`, and water mass density weighted by
`cos(phi)`. This weighted density is not a direct potential, dipole moment, or
polarizability. Layer assignment is disabled unless the user supplies density
derived `zmin` and `zmax`; the intervals are `[0,zmin)`, `[zmin,zmax)`, and
`[zmax,infinity)`.

## Hydrogen-bond donor statistics

A directed donor triplet `(O_d, H, O_a)` qualifies when the minimum-image
`O_d-O_a` distance is strictly less than `3.5 A` and the angle
`O_a-O_d-H`, measured at `O_d`, is strictly less than 30 degrees. Parent H
assignment uses the same fixed first-frame global matching with maximum
O-H distance `1.35 A`. Each parent H selects at most one acceptor, preferring
the smallest angle, then shortest O-O distance, then lowest acceptor index.

`N_donor_waters` is the number of unique donor-water oxygen atoms for which at
least one of the two parent H atoms has a qualifying acceptor. A water is
counted once even when both H atoms donate. `N_donor_hbonds` is the number of
qualifying directed `(O_d, H, O_a)` triplets, so one water contributes 0, 1,
or 2. The final cutoffs are configurable CLI parameters and default to
`O_d-O_a < 3.5 A` and `O_a-O_d-H < 30 degrees`. Layer intervals follow the
explicit polar-angle convention and are never inferred by a half split.

## RDF and coordination number

Select reference atoms X and neighbor atoms Y by element from the matching
POSCAR/CONTCAR. The optional `top-layer` X scope selects X atoms whose projected
height is within the declared tolerance of the top substrate atom in that
reference structure. Metal RDF jobs in the production matrix use this scope
explicitly; they do not use all Cu/Pd atoms. Periodic pair distances use
MDAnalysis' orthogonal or triclinic minimum-image distance kernel.

For shell `k`, the full-cell RDF is normalized as

```text
g_XY(k) = n_pairs(k) * V / (N_frames * N_X * N_Y,eff * shell_volume(k))
```

where `N_Y,eff` excludes self pairs when X and Y overlap. No accessible-slab
volume correction is applied, so a slab cell with vacuum need not have
`g(r) -> 1` at large distance. The cumulative coordination number is

```text
CN(r) = cumulative_pair_count(rmin..r) / (N_frames * N_X)
```

It therefore reports mean neighbors per reference atom, not a bond count.
Named cutoff values, block uncertainty, first RDF peak candidates, and
first-minimum candidates are recorded as diagnostics. The default named CN
cutoff is `2.5 A`; repeat `--cn-cutoff` to use a different value. The radial
limit must not exceed half the shortest periodic lattice translation.

## Bader charge and surface charge density

For atom `i`, convert the ACF Bader electron population to charge using

```text
q_i = Z_val,i - N_Bader,i
```

Positive charge means electron depletion under this convention. The surface
charge `Q_surface` is the sum over top-layer substrate metal atoms selected
from the neutral reference structure. The same zero-based indices are reused
for the matching K-containing structure; the top layer is not reselected from
the K-containing geometry.

`Z_val,i` is read from the matching POTCAR and its species order must match the
POSCAR. Explicit `--valence ELEMENT=Z` values are accepted only as overrides
and are checked against POTCAR when POTCAR is available. This is especially
important for `aimd1_e`, where validation requires `K ZVAL = 7.000`.

The reported surface-density conversion is

```text
sigma    = Q_surface / A * 1602.176634  [uC/cm2]
```

`A` is the surface cell area derived from the neutral reference POSCAR. PZC,
capacitance, `delta_phi`, and `phi_vs_SHE` are not supplied by this kit, so
potential columns are left blank and metadata records why. The response
collector always combines charge and `sigma`; it calculates a potential
response only when finite potential fields were explicitly supplied in both
input summaries. Snapshot SEM values are descriptive because static AIMD
snapshots are time correlated. The neutral and `aimd1_e` trajectories are
independent; any optional potential-response SEM is propagated as
`sqrt(SEM_neutral^2 + SEM_aimd1_e^2)`.
