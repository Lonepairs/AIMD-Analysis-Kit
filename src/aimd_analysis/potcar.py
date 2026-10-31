"""Small, dependency-free POTCAR parser used for input validation.

Only the information needed by this project is read: the ordered potential
species and their ``ZVAL`` values.  A POTCAR may contain several concatenated
potential blocks, one per POSCAR species.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path


class PotcarError(ValueError):
    """Raised when a POTCAR cannot provide an unambiguous species/ZVAL list."""


@dataclass(frozen=True)
class PotcarBlock:
    """One concatenated POTCAR potential block."""

    element: str
    zval: float


# Symbols are used to recognize ``Cu_pv``/``K_sv`` title variants without
# mistaking functional names or dates for an element symbol.
ELEMENT_SYMBOLS = frozenset(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
)

_TITLE_RE = re.compile(r"^\s*TITEL\s*=\s*(?P<body>.+?)\s*$", re.IGNORECASE)
_ZVAL_RE = re.compile(r"^\s*ZVAL\s*=\s*(?P<value>[+\-0-9.eE]+)", re.IGNORECASE)
_SYMBOL_RE = re.compile(r"(?<![A-Za-z])(?P<symbol>[A-Z][a-z]?)(?:_[A-Za-z0-9]+)?(?![A-Za-z])")


def _element_from_title(body: str) -> str:
    for match in _SYMBOL_RE.finditer(body):
        symbol = match.group("symbol")
        if symbol in ELEMENT_SYMBOLS:
            return symbol
    raise PotcarError(f"cannot identify element in POTCAR TITEL line: {body!r}")


def parse_potcar_blocks(path: Path) -> list[PotcarBlock]:
    """Parse ordered potential blocks from *path*.

    ``TITEL`` and the following ``ZVAL`` are paired by block order.  A second
    ``TITEL`` before a ``ZVAL`` is rejected so a truncated or malformed file
    cannot silently supply the wrong valence.
    """

    titles: list[str] = []
    zvals: list[float] = []
    pending_title: str | None = None
    try:
        stream = path.open("r", encoding="utf-8", errors="strict")
    except OSError as exc:
        raise PotcarError(f"cannot read POTCAR {path}: {exc}") from exc
    try:
        with stream:
            for line in stream:
                title = _TITLE_RE.match(line)
                if title:
                    if pending_title is not None:
                        raise PotcarError(f"POTCAR block has no ZVAL before next TITEL in {path}")
                    pending_title = _element_from_title(title.group("body"))
                    continue
                zval = _ZVAL_RE.match(line)
                if zval and pending_title is not None:
                    value = float(zval.group("value"))
                    if not math.isfinite(value) or value <= 0:
                        raise PotcarError(f"invalid ZVAL {value!r} for {pending_title} in {path}")
                    titles.append(pending_title)
                    zvals.append(value)
                    pending_title = None
    except (OSError, UnicodeError, ValueError) as exc:
        if isinstance(exc, PotcarError):
            raise
        raise PotcarError(f"cannot parse POTCAR {path}: {exc}") from exc
    if pending_title is not None:
        raise PotcarError(f"POTCAR block for {pending_title} has no ZVAL in {path}")
    if not titles:
        raise PotcarError(f"no TITEL/ZVAL blocks found in POTCAR {path}")
    return [PotcarBlock(element, zval) for element, zval in zip(titles, zvals)]


def parse_potcar(path: Path) -> dict[str, float]:
    """Return ``element -> ZVAL`` and reject duplicate element blocks."""

    blocks = parse_potcar_blocks(path)
    values: dict[str, float] = {}
    for block in blocks:
        if block.element in values:
            raise PotcarError(f"duplicate POTCAR element {block.element!r} in {path}")
        values[block.element] = block.zval
    return values


def validate_potcar_species(path: Path, expected_species: tuple[str, ...] | list[str]) -> dict[str, object]:
    """Validate ordered POTCAR species against a POSCAR species sequence."""

    blocks = parse_potcar_blocks(path)
    expected = tuple(expected_species)
    actual = tuple(block.element for block in blocks)
    if actual != expected:
        raise PotcarError(
            f"POTCAR species order {list(actual)!r} does not match POSCAR species {list(expected)!r}"
        )
    return {
        "path": path,
        "species": list(actual),
        "zval": {block.element: block.zval for block in blocks},
    }
