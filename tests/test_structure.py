from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from aimd_analysis.structure import (
    StructureParseError,
    fractional_to_cartesian,
    minimum_image_distance,
    minimum_image_fractional,
    minimum_image_vector,
    parse_poscar,
    surface_normal,
)
from aimd_analysis.trajectory import iter_xdatcar


def poscar(scale: str = "1.0", mode: str = "Direct") -> str:
    return (
        "fixture\n"
        f"{scale}\n"
        "2 0 0\n0 3 0\n0 0 4\n"
        "Cu O H\n1 1 2\n"
        f"{mode}\n"
        "0 0 0\n0.5 0.5 0.5\n0.2 0.2 0.2\n0.3 0.2 0.2\n"
    )


def xdatcar(markers: tuple[int, ...] = (1, 2, 3)) -> str:
    header = "fixture\n1.0\n2 0 0\n0 3 0\n0 0 4\nCu O H\n1 1 2\n"
    rows = "0 0 0\n0.5 0.5 0.5\n0.2 0.2 0.2\n0.3 0.2 0.2\n"
    return header + "".join(f"Direct configuration= {marker}\n{rows}" for marker in markers)


class StructureTests(unittest.TestCase):
    def test_negative_scale_sets_target_volume(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "POSCAR"
            path.write_text(poscar("-24.0"), encoding="utf-8")
            structure = parse_poscar(path)
        self.assertAlmostEqual(structure.volume, 24.0)

    def test_negative_scale_is_applied_to_cartesian_positions(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "POSCAR"
            path.write_text(
                "fixture\n-192.0\n2 0 0\n0 3 0\n0 0 4\nCu O H\n1 1 2\n"
                "Cartesian\n0 0 0\n1 1.5 2\n0.4 0.6 0.8\n0.6 0.6 0.8\n",
                encoding="utf-8",
            )
            structure = parse_poscar(path)
        np.testing.assert_allclose(structure.positions_frac[1], [0.5, 0.5, 0.5])

    def test_cartesian_poscar_converts_to_fractional(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "POSCAR"
            path.write_text(poscar(mode="Cartesian"), encoding="utf-8")
            structure = parse_poscar(path)
        np.testing.assert_allclose(structure.positions_frac[1], [0.25, 1 / 6, 0.125])

    def test_xdatcar_streaming_range_and_stride(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "XDATCAR"
            path.write_text(xdatcar(), encoding="utf-8")
            frames = list(iter_xdatcar(path, start=1, end=3, stride=2))
        self.assertEqual([frame.configuration for frame in frames], [1, 3])
        self.assertEqual(frames[0].positions_frac.shape, (4, 3))

    def test_xdatcar_rejects_nonmonotonic_markers(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "XDATCAR"
            path.write_text(xdatcar((1, 3, 2)), encoding="utf-8")
            with self.assertRaisesRegex(StructureParseError, "strictly increasing"):
                list(iter_xdatcar(path))

    def test_minimum_image_uses_fractional_wrapping(self) -> None:
        lattice = np.diag([10.0, 10.0, 20.0])
        delta = minimum_image_fractional(np.array([0.9, -0.6, 0.49]))
        np.testing.assert_allclose(delta, [-0.1, 0.4, 0.49])
        np.testing.assert_allclose(
            minimum_image_vector(lattice, [0.95, 0.0, 0.0], [0.05, 0.0, 0.0]),
            [1.0, 0.0, 0.0],
        )
        self.assertAlmostEqual(minimum_image_distance(lattice, [0.95, 0, 0], [0.05, 0, 0]), 1.0)

    def test_minimum_image_handles_skewed_lattice(self) -> None:
        lattice = np.array([[1.0, 0.0, 0.0], [0.9, 0.1, 0.0], [0.0, 0.0, 1.0]])
        vector = minimum_image_vector(lattice, [0.0, 0.0, 0.0], [0.49, 0.49, 0.0])
        self.assertAlmostEqual(np.linalg.norm(vector), np.linalg.norm([0.031, -0.051, 0.0]), places=10)

    def test_geometry_helpers_and_surface_normal(self) -> None:
        lattice = np.diag([2.0, 3.0, 4.0])
        np.testing.assert_allclose(fractional_to_cartesian(lattice, [0.5, 0.5, 0.5]), [1.0, 1.5, 2.0])
        np.testing.assert_allclose(surface_normal(lattice), [0.0, 0.0, 1.0])


if __name__ == "__main__":
    unittest.main()
