# Production donor-centered H-bond z profiles

These compact CSVs are the November 2026 production results for the four
metal-water systems, with neutral (`aimd1`) and 2K (`aimd1_e`) ensembles.
They are intended for Origin or other plotting tools; each file contains
`delta_z`, `Hbond_per_water`, and `total_Hbond`.

The analysis used:

- configurations 5001-9000, stride 1;
- surface-relative z bins of 0.3 A;
- donor-centered traversal;
- minimum-image `O-O < 3.5 A`;
- strict `angle O_a-O_d-H < 30 degrees`;
- parent `O-H <= 1.35 A`.

`Hbond_per_water` is the number of directed donor-H-bond triplets divided by
donor-water oxygen observations in the same z bin. A donor water may contribute
zero, one, or two bonds, so this quantity can be greater than one. The raw
XDATCAR/POSCAR files and the full per-frame tables are intentionally not part
of this repository.

The summary below records the global ratio from the same run. The `0K` label
corresponds to `aimd1`; `2K` corresponds to `aimd1_e`.
