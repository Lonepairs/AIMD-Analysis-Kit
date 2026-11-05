from __future__ import annotations

import unittest

import numpy as np

from aimd_analysis.hydrogen_bonds import find_hydrogen_bonds, find_hydrogen_bonds_acceptor_centered


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

    def test_acceptor_centered_counts_each_acceptor_once(self) -> None:
        lattice = np.diag([20.0, 20.0, 20.0])
        cart = np.asarray(
            [
                [5.0, 5.0, 5.0],  # O_a
                [7.0, 5.0, 5.0],  # closer O_b
                [8.0, 5.0, 5.0],  # second qualifying O_b
                [5.0, 5.8, 5.0],  # O_a H
                [5.0, 4.2, 5.0],  # O_a H
                [6.0, 5.0, 5.0],  # O_b 1 H
                [7.0, 6.0, 5.0],  # O_b 1 H
                [7.2, 5.0, 5.0],  # nearest H to both O_b atoms
                [8.0, 6.0, 5.0],  # O_b 2 H
            ]
        )
        fractional = cart @ np.linalg.inv(lattice)
        triplets = find_hydrogen_bonds_acceptor_centered(
            fractional, lattice, np.asarray([0, 1, 2]), np.asarray([3, 4, 5, 6, 7, 8])
        )
        acceptors = [triplet[2] for triplet in triplets]
        self.assertEqual(len(acceptors), len(set(acceptors)))
        target = [triplet for triplet in triplets if triplet[2] == 0]
        self.assertEqual(len(target), 1)
        self.assertEqual(target[0][0], 1)
        self.assertEqual(target[0][1], 5)

    def test_acceptor_centered_oo_cutoff_is_strict(self) -> None:
        lattice = np.diag([20.0, 20.0, 20.0])
        cart = np.asarray(
            [
                [5.0, 5.0, 5.0], [7.5, 5.0, 5.0],
                [5.0, 6.0, 5.0], [5.0, 4.0, 5.0],
                [6.5, 5.0, 5.0], [7.5, 6.0, 5.0],
            ]
        )
        fractional = cart @ np.linalg.inv(lattice)
        triplets = find_hydrogen_bonds_acceptor_centered(
            fractional, lattice, np.asarray([0, 1]), np.asarray([2, 3, 4, 5]), oo_cutoff=2.5
        )
        self.assertEqual(triplets, [])

    def test_acceptor_centered_angle_cutoff_is_strict(self) -> None:
        lattice = np.diag([20.0, 20.0, 20.0])
        angle = np.radians(30.0)
        cart = np.asarray(
            [
                [5.0 + 3.0 * np.cos(angle), 5.0 + 3.0 * np.sin(angle), 5.0],  # O_a
                [5.0, 5.0, 5.0],  # O_b
                [6.0, 5.0, 5.0],  # H along x, exactly 30 degrees
                [5.0, 6.0, 5.0],  # second nearest H does not qualify
                [8.5, 6.5, 5.0],  # O_a H
                [7.4, 7.3, 5.0],  # O_a H
            ]
        )
        fractional = cart @ np.linalg.inv(lattice)
        triplets = find_hydrogen_bonds_acceptor_centered(
            fractional, lattice, np.asarray([0, 1]), np.asarray([2, 3, 4, 5])
        )
        self.assertEqual(triplets, [])

    def test_acceptor_centered_requires_two_hydrogens_per_oxygen(self) -> None:
        lattice = np.diag([20.0, 20.0, 20.0])
        with self.assertRaisesRegex(ValueError, "two H atoms per O"):
            find_hydrogen_bonds_acceptor_centered(
                np.zeros((3, 3)), lattice, np.asarray([0, 1]), np.asarray([2])
            )


if __name__ == "__main__":
    unittest.main()
