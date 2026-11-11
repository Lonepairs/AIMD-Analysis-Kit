"""Bader snapshot aggregation and surface charge-density analysis."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import re
import warnings as warnings_module
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np

from . import __version__
from .density import _git_revision, _portable_input_identifier, _sha256
from .io import prepare_output_dir
from .potcar import PotcarError, parse_potcar, validate_potcar_species
from .structure import Structure, parse_poscar, surface_normal


SYSTEMS = ("Cu", "CuxPd1", "Pd", "PdxCu1")
# Kept only as a compatibility fallback for old neutral fixture calls that do
# not contain a POTCAR. Production validation and every K-containing case
# require POTCAR values or explicit, checked overrides.
LEGACY_FIXTURE_VALENCE = {"Cu": 11.0, "Pd": 10.0, "O": 6.0, "H": 1.0, "He": 2.0, "K": 7.0}
LEGACY_FIXTURE_SURFACE_AREA_A2 = {"Cu": 90.52, "CuxPd1": 90.52, "Pd": 104.87, "PdxCu1": 104.87}
ACF_ROW = re.compile(
    r"^\s*(\d+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$"
)
NELECT_ROW = re.compile(r"\bNELECT\s*=\s*([-+0-9.Ee]+)")
NONMETAL = frozenset({"O", "H", "He", "K"})
ELECTRON_COUNT_TOLERANCE = 1.0e-3


def _write_csv(path: Path, header: list[str], rows: Iterable[Iterable[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def parse_acf(path: Path) -> list[dict[str, float | int]]:
    rows = []
    with path.open("r", encoding="utf-8", errors="strict") as stream:
        for line in stream:
            match = ACF_ROW.match(line)
            if not match:
                continue
            atom, x, y, z, population, min_dist, volume = match.groups()
            values = [float(x), float(y), float(z), float(population), float(min_dist), float(volume)]
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"non-finite ACF value for atom {atom} in {path.name}")
            rows.append(
                {
                    "atom_index": int(atom), "x": values[0], "y": values[1], "z": values[2],
                    "electron_population": values[3], "min_dist_A": values[4], "atomic_volume_A3": values[5],
                }
            )
    if not rows:
        raise ValueError(f"no atom rows found in {path}")
    if [row["atom_index"] for row in rows] != list(range(1, len(rows) + 1)):
        raise ValueError(f"ACF atom indices are not contiguous from 1 in {path}")
    return rows


def parse_chgcar_grid(path: Path) -> tuple[int, int, int]:
    """Read the first VASP charge-density grid dimensions."""

    source = Path(path)
    try:
        lines = source.read_text(encoding="utf-8", errors="strict").splitlines()
        structure = parse_poscar(source)
    except (OSError, UnicodeError, ValueError) as exc:
        raise ValueError(f"cannot read CHGCAR header from {source}: {exc}") from exc
    cursor = 7
    if cursor < len(lines) and lines[cursor].strip().lower().startswith("s"):
        cursor += 1
    cursor += 1  # coordinate mode
    grid_index = cursor + structure.atom_count
    while grid_index < len(lines) and not lines[grid_index].strip():
        grid_index += 1
    if grid_index >= len(lines):
        raise ValueError(f"no charge-density grid found in {source}")
    fields = lines[grid_index].split()
    if len(fields) != 3:
        raise ValueError(f"invalid charge-density grid in {source}: {lines[grid_index]!r}")
    try:
        grid = tuple(int(value) for value in fields)
    except ValueError as exc:
        raise ValueError(f"invalid charge-density grid in {source}: {lines[grid_index]!r}") from exc
    if any(value <= 0 for value in grid):
        raise ValueError(f"charge-density grid must be positive in {source}")
    return grid


def validate_density_grids(paths: Iterable[Path]) -> tuple[int, int, int]:
    """Require all supplied charge-density files to use the same grid."""

    sources = [Path(path) for path in paths]
    if len(sources) < 2:
        raise ValueError("at least two charge-density files are required for a grid check")
    grids = [(source, parse_chgcar_grid(source)) for source in sources]
    reference_path, reference_grid = grids[0]
    mismatches = [
        f"{source.name}={grid} (expected {reference_path.name}={reference_grid})"
        for source, grid in grids[1:]
        if grid != reference_grid
    ]
    if mismatches:
        raise ValueError("charge-density grid mismatch: " + "; ".join(mismatches))
    return reference_grid


def parse_nelect(path: Path) -> float:
    """Read the first ``NELECT =`` value from an OUTCAR."""

    source = Path(path)
    try:
        with source.open("r", encoding="utf-8", errors="strict") as stream:
            for line in stream:
                match = NELECT_ROW.search(line)
                if match:
                    value = float(match.group(1))
                    if not math.isfinite(value) or value < 0:
                        break
                    return value
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"cannot read NELECT from {source}: {exc}") from exc
    raise ValueError(f"no finite NELECT value found in {source}")


def select_top_layer(structure: Structure, tolerance: float = 0.45, surface_species: str | None = None) -> list[int]:
    labels = structure.atom_species
    if surface_species:
        selected_species = {value.strip() for value in surface_species.split(",") if value.strip()}
        unknown = selected_species - set(labels)
        if unknown:
            raise ValueError("surface species absent from reference POSCAR: " + ", ".join(sorted(unknown)))
        selected_species -= NONMETAL
    else:
        selected_species = set(labels) - NONMETAL
    if not selected_species:
        raise ValueError("no substrate metal species selected for surface detection")
    normal = surface_normal(structure.lattice)
    projected = structure.positions_cart @ normal
    substrate = np.asarray([i for i, label in enumerate(labels) if label in selected_species], dtype=int)
    maximum = float(np.max(projected[substrate]))
    top = substrate[projected[substrate] >= maximum - tolerance]
    if not len(top):
        raise ValueError("no substrate atoms satisfy the top-layer criterion")
    return [int(value) for value in top]


def _load_top_layer_metadata(path: Path, structure: Structure) -> list[int] | None:
    if not path.is_file():
        return None
    metadata = json.loads(path.read_text(encoding="utf-8"))
    if metadata.get("atom_count", structure.atom_count) != structure.atom_count:
        raise ValueError("surface metadata atom count does not match reference structure")
    raw = metadata.get("top_layer_atom_indices_zero_based")
    if raw is None:
        return None
    indices = sorted(set(int(value) for value in raw))
    if not indices or indices[0] < 0 or indices[-1] >= structure.atom_count:
        raise ValueError("surface metadata contains invalid atom indices")
    if any(structure.atom_species[index] in NONMETAL for index in indices):
        raise ValueError("surface metadata includes a non-substrate atom")
    return indices


def _charge_to_sigma(charge_e: float, area_A2: float) -> float:
    """Convert elementary-charge/Angstrom-squared to microcoulomb/cm2."""

    return charge_e / area_A2 * 1602.176634


def _parse_valences(items: list[str]) -> dict[str, float]:
    values: dict[str, float] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"valence must be ELEMENT=Z, got {item!r}")
        element, raw = item.split("=", 1)
        try:
            value = float(raw)
        except ValueError as exc:
            raise ValueError(f"invalid valence {item!r}") from exc
        if not element or not math.isfinite(value) or value <= 0:
            raise ValueError(f"valence must be finite and positive in {item!r}")
        values[element] = value
    return values


def _resolve_potcar_path(args: argparse.Namespace, case_dir: Path, snapshots: list[tuple[int, Path]] | None = None) -> Path | None:
    explicit = getattr(args, "potcar", None)
    if explicit is not None:
        return Path(explicit).expanduser()
    candidates = [case_dir / "POTCAR"]
    if snapshots:
        candidates.append(snapshots[0][1] / "POTCAR")
    candidates.append(case_dir.parent / "aimd1" / "POTCAR")
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def _resolve_valences(
    args: argparse.Namespace,
    reference: Structure,
    case_dir: Path,
    snapshots: list[tuple[int, Path]] | None = None,
) -> tuple[dict[str, float], Path | None, str]:
    """Resolve valences from POTCAR, checking explicit overrides if present."""

    overrides = _parse_valences(getattr(args, "valence", []))
    potcar_path = _resolve_potcar_path(args, case_dir, snapshots)
    if potcar_path is not None:
        try:
            parsed = validate_potcar_species(potcar_path, reference.species)
        except PotcarError as exc:
            raise ValueError(str(exc)) from exc
        valence = {str(key): float(value) for key, value in parsed["zval"].items()}
        mismatches = [
            f"{element}: override {overrides[element]:.6f} != POTCAR ZVAL {valence[element]:.6f}"
            for element in overrides
            if element in valence and not math.isclose(overrides[element], valence[element], rel_tol=0.0, abs_tol=1e-6)
        ]
        if mismatches:
            raise ValueError("POTCAR valence mismatch: " + "; ".join(mismatches))
        missing = sorted(set(reference.atom_species) - set(valence))
        if missing:
            raise ValueError("POTCAR is missing valence reference(s): " + ", ".join(missing))
        return valence, potcar_path, "POTCAR ZVAL"
    if overrides:
        missing = sorted(set(reference.atom_species) - set(overrides))
        if missing:
            raise ValueError("missing valence reference(s) and no POTCAR: " + ", ".join(missing))
        if args.ensemble == "aimd1_e" and not math.isclose(overrides.get("K", -1.0), 7.0, rel_tol=0.0, abs_tol=1e-6):
            raise ValueError("aimd1_e K valence must be verified as POTCAR ZVAL K=7.000")
        return overrides, None, "explicit override (POTCAR not found)"
    if args.ensemble == "aimd1_e":
        raise ValueError("K-containing cases require a POTCAR with ZVAL K=7.000")
    # Compatibility for the old minimal neutral fixture. Production runs are
    # blocked by validate_cluster_inputs.py unless a matching POTCAR exists.
    missing = sorted(set(reference.atom_species) - set(LEGACY_FIXTURE_VALENCE))
    if missing:
        raise ValueError("missing POTCAR valence reference(s): " + ", ".join(missing))
    return {element: LEGACY_FIXTURE_VALENCE[element] for element in reference.atom_species}, None, "legacy fixture fallback (no POTCAR)"


def _std(values: list[float]) -> float | None:
    return float(np.std(values, ddof=1)) if len(values) > 1 else None


def _sem(values: list[float]) -> float | None:
    standard_deviation = _std(values)
    return standard_deviation / math.sqrt(len(values)) if standard_deviation is not None else None


def _csv_stat(value: float | None) -> float | str:
    return "" if value is None else value


def _validate_args(args: argparse.Namespace) -> None:
    if args.start_step < 1 or args.end_step < args.start_step or args.step_stride < 1:
        raise ValueError("Bader step grid is invalid")
    if args.top_layer_tolerance <= 0:
        raise ValueError("top-layer-tolerance must be positive")
    electron_tolerance = float(getattr(args, "electron_tolerance", ELECTRON_COUNT_TOLERANCE))
    if not math.isfinite(electron_tolerance) or electron_tolerance < 0:
        raise ValueError("electron-tolerance must be finite and nonnegative")
    if args.ensemble not in {"aimd1", "aimd1_e"}:
        raise ValueError("ensemble must be aimd1 or aimd1_e")
    if args.ensemble == "aimd1_e" and "K" not in _parse_valences(args.valence):
        root = Path(args.root).expanduser()
        case_dir = root / args.system / args.ensemble
        explicit = getattr(args, "potcar", None)
        candidates = [Path(explicit).expanduser()] if explicit is not None else [
            case_dir / "POTCAR",
            root / args.system / "bader_e" / f"bader_{args.system}_{args.start_step}" / "POTCAR",
        ]
        if not any(path.is_file() for path in candidates):
            raise ValueError("K-containing cases require a POTCAR with verified ZVAL K=7.000")


def analyze_bader(args: argparse.Namespace) -> dict:
    _validate_args(args)
    root = Path(args.root).expanduser()
    system = args.system
    ensemble = args.ensemble
    case_dir = root / system / ensemble
    reference_path = Path(args.reference).expanduser() if args.reference else case_dir / "POSCAR"
    surface_reference_path = Path(args.surface_reference).expanduser() if args.surface_reference else root / system / "aimd1" / "POSCAR"
    output_dir = Path(args.output_dir).expanduser()
    if not reference_path.is_file():
        raise FileNotFoundError(f"reference POSCAR not found for {system}/{ensemble}")
    if not surface_reference_path.is_file():
        raise FileNotFoundError("neutral surface-reference POSCAR is missing")
    reference = parse_poscar(reference_path)
    surface_reference = parse_poscar(surface_reference_path)

    top_layer = None
    top_source = "derived from neutral reference POSCAR"
    if args.surface_metadata:
        top_layer = _load_top_layer_metadata(Path(args.surface_metadata).expanduser(), surface_reference)
        if top_layer is not None:
            top_source = "user-supplied neutral surface metadata"
    if top_layer is None:
        top_layer = select_top_layer(surface_reference, args.top_layer_tolerance, args.surface_species)
    if any(index >= reference.atom_count for index in top_layer):
        raise ValueError("neutral top-layer indices exceed target structure atom count")
    for index in top_layer:
        if surface_reference.atom_species[index] != reference.atom_species[index]:
            raise ValueError(f"neutral/electrochemical substrate element order differs at atom {index + 1}")

    bader_subdir = args.bader_dir or ("bader_e" if ensemble == "aimd1_e" else "bader")
    steps = list(range(args.start_step, args.end_step + 1, args.step_stride))
    snapshot_root = root / system / bader_subdir
    snapshots = [(step, snapshot_root / f"bader_{system}_{step}") for step in steps]
    missing_dirs = [str(step) for step, path in snapshots if not (path / "ACF.dat").is_file() or not (path / "POSCAR").is_file()]
    if missing_dirs:
        raise FileNotFoundError("missing Bader snapshot files for steps: " + ", ".join(missing_dirs))
    valence, potcar_path, valence_source = _resolve_valences(args, reference, case_dir, snapshots)

    metal_elements = [element for element in dict.fromkeys(reference.atom_species) if element not in NONMETAL]
    element_counts = Counter(reference.atom_species)
    singleton_metals = {element for element in metal_elements if element_counts[element] == 1}
    groups = {}
    for index, element in enumerate(reference.atom_species):
        if element == "O":
            groups[index] = "water_O"
        elif element == "H":
            groups[index] = "water_H"
        elif element == "He":
            groups[index] = "He"
        elif element == "K":
            groups[index] = "counterion_K"
        elif element in singleton_metals:
            groups[index] = "dopant"
        elif index in set(top_layer):
            groups[index] = "top_host"
        else:
            groups[index] = "subsurface_metal"

    atom_rows = []
    group_snapshot_rows = []
    electron_rows = []
    surface_rows = []
    electron_checks = []
    density_grids = []
    analysis_warnings: list[str] = []
    common_density_grid: tuple[int, int, int] | None = None
    expected_valence_population = sum(valence[element] for element in reference.atom_species)
    electron_tolerance = float(getattr(args, "electron_tolerance", ELECTRON_COUNT_TOLERANCE))
    strict_electron_count = bool(getattr(args, "strict_electron_count", False))
    for step, snapshot_dir in snapshots:
        snapshot_structure = parse_poscar(snapshot_dir / "POSCAR")
        if (snapshot_structure.species, snapshot_structure.counts) != (reference.species, reference.counts):
            raise ValueError(f"snapshot POSCAR species/count order mismatch at step {step}")
        acf_rows = parse_acf(snapshot_dir / "ACF.dat")
        if len(acf_rows) != reference.atom_count:
            raise ValueError(f"ACF/POSCAR atom count mismatch at step {step}: {len(acf_rows)} vs {reference.atom_count}")
        density_paths = [
            snapshot_dir / name
            for name in ("CHGCAR", "CHGCAR_sum", "AECCAR0", "AECCAR2")
            if (snapshot_dir / name).is_file()
        ]
        if len(density_paths) >= 2:
            current_grid = validate_density_grids(density_paths)
            density_grids.append({"step": step, "dimensions": list(current_grid), "files": [path.name for path in density_paths]})
            if common_density_grid is None:
                common_density_grid = current_grid
            elif current_grid != common_density_grid:
                raise ValueError(
                    f"charge-density grid changes between snapshots at step {step}: "
                    f"{current_grid} vs {common_density_grid}"
                )
        elif density_paths:
            message = f"step {step}: cannot compare charge-density grids; only {density_paths[0].name} is present"
            analysis_warnings.append(message)
            warnings_module.warn(message, RuntimeWarning, stacklevel=2)
        if (snapshot_dir / "CHGCAR").is_file() and (snapshot_dir / "AECCAR0").is_file() and (snapshot_dir / "AECCAR2").is_file() and not (snapshot_dir / "CHGCAR_sum").is_file():
            message = f"step {step}: CHGCAR_sum is missing although CHGCAR, AECCAR0, and AECCAR2 are present"
            analysis_warnings.append(message)
            warnings_module.warn(message, RuntimeWarning, stacklevel=2)
        group_values: dict[str, list[float]] = {}
        group_populations: dict[str, list[float]] = {}
        for index, (element, acf_row) in enumerate(zip(reference.atom_species, acf_rows)):
            if acf_row["atom_index"] != index + 1:
                raise ValueError(f"ACF row order mismatch at step {step}, atom {index + 1}")
            population = float(acf_row["electron_population"])
            charge = valence[element] - population
            group = groups[index]
            group_values.setdefault(group, []).append(charge)
            group_populations.setdefault(group, []).append(population)
            atom_rows.append(
                (
                    system, ensemble, step, index + 1, element, group, population, valence[element], charge,
                    acf_row["x"], acf_row["y"], acf_row["z"], acf_row["min_dist_A"], acf_row["atomic_volume_A3"],
                )
            )
        for group, values in group_values.items():
            populations = group_populations[group]
            group_snapshot_rows.append(
                (system, ensemble, step, group, len(values), np.mean(values), _csv_stat(_std(values)), np.sum(values), np.mean(populations), np.sum(populations))
            )
        total_population = sum(float(row["electron_population"]) for row in acf_rows)
        outcar_path = snapshot_dir / "OUTCAR"
        electron_reference = expected_valence_population
        electron_reference_source = "sum of POTCAR ZVAL references"
        if outcar_path.is_file():
            electron_reference = parse_nelect(outcar_path)
            electron_reference_source = "OUTCAR NELECT"
        electron_residual = total_population - electron_reference
        electron_checks.append(
            {
                "step": step,
                "total_bader_electron_population": total_population,
                "reference_electron_population": electron_reference,
                "residual": electron_residual,
                "reference_source": electron_reference_source,
            }
        )
        if abs(electron_residual) > electron_tolerance:
            message = (
                f"step {step}: Bader electron count residual {electron_residual:.6g} e exceeds "
                f"tolerance {electron_tolerance:.6g} e ({electron_reference_source})"
            )
            analysis_warnings.append(message)
            if strict_electron_count:
                raise ValueError(message)
            warnings_module.warn(message, RuntimeWarning, stacklevel=2)
        electron_rows.append(
            (
                system, ensemble, step, len(acf_rows), total_population,
                expected_valence_population, total_population - expected_valence_population,
                electron_reference, electron_residual,
            )
        )
        surface_charge = sum(
            valence[reference.atom_species[index]] - float(acf_rows[index]["electron_population"])
            for index in top_layer
        )
        surface_area_A2 = float(np.linalg.norm(np.cross(surface_reference.lattice[0], surface_reference.lattice[1])))
        if valence_source.startswith("legacy fixture"):
            surface_area_A2 = LEGACY_FIXTURE_SURFACE_AREA_A2[system]
        if surface_area_A2 <= 0:
            raise ValueError("surface-reference POSCAR has invalid surface area")
        sigma = _charge_to_sigma(surface_charge, surface_area_A2)
        element_surface = {}
        for index in top_layer:
            element = reference.atom_species[index]
            element_surface[element] = element_surface.get(element, 0.0) + valence[element] - float(acf_rows[index]["electron_population"])
        surface_rows.append(
            (system, ensemble, step, len(top_layer), surface_charge, surface_charge / len(top_layer), sigma, "", "", json.dumps(element_surface, sort_keys=True))
        )

    output_dir = prepare_output_dir(output_dir, overwrite=getattr(args, "overwrite", False))
    _write_csv(
        output_dir / "bader_atom_charges.csv",
        ["system", "ensemble", "step", "atom_index", "element", "group", "bader_electron_population", "valence_electrons", "charge_q", "x_A", "y_A", "z_A", "min_dist_A", "atomic_volume_A3"],
        atom_rows,
    )
    _write_csv(
        output_dir / "bader_group_snapshot.csv",
        ["system", "ensemble", "step", "group", "n_atoms", "mean_charge_q", "std_charge_q_within_group", "total_charge_q", "mean_electron_population", "total_electron_population"],
        group_snapshot_rows,
    )
    _write_csv(
        output_dir / "bader_total_electrons.csv",
        [
            "system", "ensemble", "step", "n_atoms", "total_bader_electron_population",
            "expected_valence_population", "population_residual", "electron_count_reference",
            "electron_count_residual",
        ],
        electron_rows,
    )
    _write_csv(
        output_dir / "bader_surface_metal_charge.csv",
        ["system", "ensemble", "step", "n_surface_metal_atoms", "net_surface_metal_charge_q", "mean_surface_metal_charge_q", "sigma_uC_cm-2", "delta_phi_V", "phi_V_vs_SHE", "surface_charge_by_element_json"],
        surface_rows,
    )

    group_summary = []
    for group in sorted({row[3] for row in group_snapshot_rows}):
        rows = [row for row in group_snapshot_rows if row[3] == group]
        means = [float(row[5]) for row in rows]
        totals = [float(row[7]) for row in rows]
        group_summary.append(
            (system, ensemble, group, int(rows[0][4]), len(rows), np.mean(means), _csv_stat(_std(means)), _csv_stat(_sem(means)), np.mean(totals))
        )
    _write_csv(
        output_dir / "bader_group_summary.csv",
        ["system", "ensemble", "group", "n_atoms_per_snapshot", "n_snapshots", "mean_charge_q", "snapshot_std_q", "snapshot_sem_q", "mean_total_charge_q"],
        group_summary,
    )
    water_values = []
    for step in steps:
        water_values.append(sum(float(row[7]) for row in group_snapshot_rows if row[2] == step and row[3] in {"water_O", "water_H"}))
    water_mean = float(np.mean(water_values))
    oxygen_count = element_counts.get("O", 0)
    _write_csv(
        output_dir / "bader_water_charge_summary.csv",
        ["system", "ensemble", "n_snapshots", "mean_total_water_charge_q", "snapshot_std_total_water_charge_q", "snapshot_sem_total_water_charge_q", "mean_water_charge_per_molecule_q"],
        [(system, ensemble, len(water_values), water_mean, _csv_stat(_std(water_values)), _csv_stat(_sem(water_values)), water_mean / oxygen_count if oxygen_count else "")],
    )
    surface_values = [float(row[4]) for row in surface_rows]
    sigma_values = [float(row[6]) for row in surface_rows]
    surface_area_A2 = float(np.linalg.norm(np.cross(surface_reference.lattice[0], surface_reference.lattice[1])))
    if valence_source.startswith("legacy fixture"):
        surface_area_A2 = LEGACY_FIXTURE_SURFACE_AREA_A2[system]
    surface_summary = (
        system, ensemble, len(steps), len(top_layer), np.mean(surface_values), _csv_stat(_std(surface_values)),
        _csv_stat(_sem(surface_values)), np.mean(sigma_values), "", "",
    )
    _write_csv(
        output_dir / "bader_surface_metal_charge_summary.csv",
        ["system", "ensemble", "n_snapshots", "n_surface_metal_atoms", "mean_net_surface_metal_charge_q", "snapshot_std_net_surface_metal_charge_q", "snapshot_sem_net_surface_metal_charge_q", "mean_sigma_uC_cm-2", "mean_delta_phi_V", "mean_phi_V_vs_SHE"],
        [surface_summary],
    )

    root_path = root.resolve()
    reference_id, reference_redacted = _portable_input_identifier(reference_path, root_path)
    surface_reference_id, surface_reference_redacted = _portable_input_identifier(surface_reference_path, root_path)
    snapshot_ids = []
    for step, directory in snapshots:
        snapshot_ids.append({
            "step": step,
            "poscar": _portable_input_identifier(directory / "POSCAR", root_path)[0],
            "acf": _portable_input_identifier(directory / "ACF.dat", root_path)[0],
            "poscar_sha256": _sha256(directory / "POSCAR"),
            "acf_sha256": _sha256(directory / "ACF.dat"),
        })
    metadata = {
        "schema_version": 1,
        "analysis": "Bader charge and surface charge-density analysis",
        "package_version": __version__,
        "git_revision": _git_revision(),
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "operating_system": platform.system()},
        "system": system,
        "ensemble": ensemble,
        "step_grid": steps,
        "valence_references": valence,
        "valence_source": valence_source,
        "potcar": None if potcar_path is None else {
            "path": _portable_input_identifier(potcar_path, root_path)[0],
            "absolute_path_redacted": _portable_input_identifier(potcar_path, root_path)[1],
            "sha256": _sha256(potcar_path),
        },
        "charge_definition": "q_i = Z_val,i - N_Bader,i; positive q means electron depletion",
        "surface_charge_definition": "Sum of q over neutral-reference top-layer substrate atoms; indices are reused for aimd1_e.",
        "surface_reference": {"path": surface_reference_id, "absolute_path_redacted": surface_reference_redacted, "sha256": _sha256(surface_reference_path)},
        "trajectory_reference": {"path": reference_id, "absolute_path_redacted": reference_redacted, "sha256": _sha256(reference_path)},
        "top_layer_atom_indices_zero_based": top_layer,
        "top_layer_atom_count": len(top_layer),
        "top_layer_source": top_source,
        "surface_species": args.surface_species,
        "top_layer_tolerance_A": args.top_layer_tolerance,
        "snapshots": snapshot_ids,
        "density_grid_check": {
            "common_dimensions": None if common_density_grid is None else list(common_density_grid),
            "per_snapshot": density_grids,
            "files_checked": ["CHGCAR", "CHGCAR_sum", "AECCAR0", "AECCAR2"],
        },
        "electron_count_check": {
            "tolerance_e": electron_tolerance,
            "reference_rule": "OUTCAR NELECT when present, otherwise sum of POTCAR ZVAL references",
            "strict": strict_electron_count,
            "per_snapshot": electron_checks,
        },
        "warnings": analysis_warnings,
        "surface_charge_density": {
            "equation": "sigma_uC_cm-2 = Q_surface_e / area_A2 * 1602.176634",
            "conversion_uC_cm2_per_e_A2": 1602.176634,
            "surface_area_A2": surface_area_A2,
            "interpretation": "Surface charge density only; no PZC, capacitance, or potential is assumed.",
        },
        "potential_estimate": {
            "status": "not calculated",
            "reason": "PZC and capacitance were not supplied",
            "pzc_V": None,
            "capacitance_uF_cm2": None,
        },
        "expected_valence_population": expected_valence_population,
        "population_residuals_e": [float(row[6]) for row in electron_rows],
        "output_files": ["bader_atom_charges.csv", "bader_group_snapshot.csv", "bader_total_electrons.csv", "bader_surface_metal_charge.csv", "bader_group_summary.csv", "bader_water_charge_summary.csv", "bader_surface_metal_charge_summary.csv", "analysis_metadata.json", "bader_metadata.json"],
    }
    serialized_metadata = json.dumps(metadata, indent=2, sort_keys=True)
    (output_dir / "analysis_metadata.json").write_text(serialized_metadata, encoding="utf-8")
    (output_dir / "bader_metadata.json").write_text(serialized_metadata, encoding="utf-8")
    return metadata


def add_bader_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--root", type=Path, required=True, help="Portable data root with system/ensemble and Bader snapshot directories")
    parser.add_argument("--system", choices=SYSTEMS, required=True)
    parser.add_argument("--ensemble", choices=("aimd1", "aimd1_e"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start-step", type=int, default=5500)
    parser.add_argument("--end-step", type=int, default=9000)
    parser.add_argument("--step-stride", type=int, default=500)
    parser.add_argument("--reference", type=Path, default=None, help="Override POSCAR for the selected ensemble")
    parser.add_argument("--surface-reference", type=Path, default=None, help="Neutral POSCAR used for stable top-layer indices")
    parser.add_argument("--surface-metadata", type=Path, default=None, help="Optional metadata JSON with top_layer_atom_indices_zero_based")
    parser.add_argument("--surface-species", default=None, metavar="SPECIES[,SPECIES...]")
    parser.add_argument("--top-layer-tolerance", type=float, default=0.45, metavar="A")
    parser.add_argument("--bader-dir", default=None, help="Override snapshot subdirectory; default follows ensemble")
    parser.add_argument("--potcar", type=Path, default=None, help="POTCAR to validate; defaults to the selected case POTCAR")
    parser.add_argument("--valence", action="append", default=[], metavar="ELEMENT=Z", help="Explicit valence override; checked against POTCAR when present")
    parser.add_argument("--electron-tolerance", type=float, default=ELECTRON_COUNT_TOLERANCE, metavar="E", help="Allowed ACF/NELECT electron-count residual before warning")
    parser.add_argument("--strict-electron-count", action="store_true", help="Fail when the ACF electron-count residual exceeds --electron-tolerance")
    parser.add_argument("--overwrite", action="store_true", help="Allow reuse of a non-empty output directory")
