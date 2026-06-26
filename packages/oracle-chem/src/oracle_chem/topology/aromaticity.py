"""
Aromaticity detection based on:
- ring topology
- planarity
- atomic eligibility derived from synthons

This module does NOT infer electronic structure.
It classifies rings and atoms using available descriptors.
"""

import numpy as np


class Aromaticity:
    """
    Aromaticity detector.
    """

    def __init__(
        self,
        graph,
        discrete_graph,
        ring_set,
        synthons=None,
        force_aromatic=False,
        planarity_tol=0.1,
    ):
        self.graph = graph
        self.dgraph = discrete_graph
        self.ringset = ring_set
        self.synthons = synthons
        self.force_aromatic = force_aromatic
        self.planarity_tol = planarity_tol

        self.Z = graph.Z
        self.coords = graph.coords

        self.aromatic_atoms = set()
        self.aromatic_bonds = set()

        self._analyze()

    # --------------------------------------------------------
    # Core logic
    # --------------------------------------------------------

    def _analyze(self):
        for ring in self.ringset:
            if not self._is_planar(ring):
                continue

            if not self._is_ring_aromatic(ring):
                continue

            for i in ring.atoms:
                self.aromatic_atoms.add(i)

            for (i, j) in ring.bonds:
                self.aromatic_bonds.add((i, j))

    # --------------------------------------------------------
    # Planarity
    # --------------------------------------------------------

    def _is_planar(self, ring):
        coords = np.array([self.coords[i] for i in ring.atoms])
        centroid = coords.mean(axis=0)
        coords -= centroid

        _, _, vh = np.linalg.svd(coords, full_matrices=False)
        normal = vh[-1]

        distances = np.abs(coords @ normal)
        return distances.max() < self.planarity_tol

    # --------------------------------------------------------
    # Aromaticity decision
    # --------------------------------------------------------

    def _is_ring_aromatic(self, ring):
        if self.force_aromatic:
            return True

        for i in ring.atoms:
            if not self._is_atom_aromatic(i):
                return False

        return True

    def _is_atom_aromatic(self, i):
        """
        Decide whether atom i can participate in aromaticity.

        Aromaticity is treated as a discrete attribute derived from:
        - ring membership
        - planarity (checked at ring level)
        - atomic ability to sustain pi delocalization
        """

        # Preferred path: derive from synthons
        if self.synthons is not None:
            # Use continuous pi-delocalization as eligibility proxy
            d_pi = self.synthons.delocalization(i)

            # Minimal threshold: must be non-zero and not purely localized
            return d_pi > 1e-3

        # Fallback: conservative element/coordination rules
        Z = self.Z[i]
        deg = len(self.dgraph.neighbors(i))

        if Z == 6:
            return deg == 3
        if Z in (7, 15):
            return deg in (2, 3)
        if Z in (8, 16):
            return deg == 2

        return False

