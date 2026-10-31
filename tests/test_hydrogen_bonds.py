from __future__ import annotations

import unittest

import numpy as np

from aimd_analysis.hydrogen_bonds import find_hydrogen_bonds


class HydrogenBondTests(unittest.TestCase):
    def test_each_parent_h_selects_one_best_acceptor(self) -> None:
        lattice = np.diag([20.0, 20.0, 20.0])
        cart = np.asarray(
            [
                [5.0, 5.0, 5.0],  # donor O
                [7.0, 5.0, 5.0],  # closer x acceptor
                [8.0, 5.0, 5.0],  # farther collinear x acceptor
                [5.0, 7.0, 5.0],  # y acceptor
                [6.0, 5.0, 5.0],  # H along x
                [5.0, 6.0, 5.0],  # H along y
            ]
        )
        fractional = cart @ np.linalg.inv(lattice)
        triplets = find_hydrogen_bonds(
            fractional,
            lattice,
            np.asarray([0, 1, 2, 3]),
            {0: (4, 5)},
            selected_donors={0},
        )
        self.assertEqual(len(triplets), 2)
        self.assertEqual({triplet[1] for triplet in triplets}, {4, 5})
        self.assertEqual({triplet[2] for triplet in triplets}, {1, 3})

    def test_oo_cutoff_is_strict(self) -> None:
        lattice = np.diag([20.0, 20.0, 20.0])
        cart = np.asarray([[5.0, 5.0, 5.0], [7.5, 5.0, 5.0], [6.0, 5.0, 5.0], [5.0, 6.0, 5.0]])
        fractional = cart @ np.linalg.inv(lattice)
        triplets = find_hydrogen_bonds(
            fractional, lattice, np.asarray([0, 1]), {0: (2, 3)}, selected_donors={0}, oo_cutoff=2.5
        )
        self.assertEqual(triplets, [])

    def test_angle_cutoff_is_strict(self) -> None:
        lattice = np.diag([20.0, 20.0, 20.0])
        angle = np.radians(30.0)
        cart = np.asarray(
            [
                [5.0, 5.0, 5.0],
                [5.0 + 3.0 * np.cos(angle), 5.0 + 3.0 * np.sin(angle), 5.0],
                [6.0, 5.0, 5.0],
                [5.0, 6.0, 5.0],
            ]
        )
        fractional = cart @ np.linalg.inv(lattice)
        triplets = find_hydrogen_bonds(
            fractional,
            lattice,
            np.asarray([0, 1]),
            {0: (2, 3)},
            selected_donors={0},
            angle_max_deg=30.0,
        )
        self.assertEqual(triplets, [])


if __name__ == "__main__":
    unittest.main()
