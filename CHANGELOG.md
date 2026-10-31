# Changelog

## 0.1.0

- Added VASP POSCAR/XDATCAR parsing, inventory validation, and shared geometry.
- Added slice-based water density, polar angle with `rho_H2O*cos(phi)`,
  hydrogen-bond donor statistics, RDF/CN, and Bader snapshot aggregation.
- Added numbered and config-driven batch interfaces, neutral/2K Bader response
  collection, and portable cluster preparation/packaging helpers.
- Added deterministic synthetic fixtures and scientific method tests.
- Added an eight-case production workflow matrix with portable paths and kept
  the synthetic smoke workflow separate.
- Corrected Bader response uncertainty to independent-ensemble propagation;
  undefined single-snapshot SEMs are left blank. RDF SEM now excludes a short
  trailing block instead of weighting it like a full block.
- Analysis and inventory commands now refuse to reuse non-empty output
  directories unless `--overwrite` is explicitly passed.
- Bader valences now come from ordered POTCAR `ZVAL` records, with checked
  overrides and explicit K=7.000 validation for `aimd1_e`.
- Bader output now reports charge and surface charge density without assumed
  PZC/capacitance; potential fields remain optional. Added read-only
  `scripts/validate_cluster_inputs.py` for all eight cases and composition
  checks.
- Metal RDF matrix jobs use `top-layer` scope and record a default adjustable
  CN cutoff of 2.5 A.
