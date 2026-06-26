class Ring:
    """
    Topological + geometrical ring object.

    Ring coordinates are intended to be built from
    internal angles (valence and dihedral angles),
    not from Cartesian puckering coordinates.
    """

    def __init__(self, index, atoms, coords):
        self.index = index
        self.atoms = list(atoms)            # ordered cyclic list
        self.coords = np.asarray(coords)    # (N,3)

        # Topology
        self.bonds = self._build_bonds()
        self.adjacent_rings_atoms = set()
        self.adjacent_rings_bonds = set()

    # ----------------------------------------------------------
    # Cyclic topology helpers
    # ----------------------------------------------------------

    def cyclic_triplets(self):
        """
        Triplets (i-1, i, i+1) for valence angles.
        """
        n = len(self.atoms)
        return [
            (self.atoms[(i - 1) % n],
             self.atoms[i],
             self.atoms[(i + 1) % n])
            for i in range(n)
        ]

    def cyclic_quartets(self):
        """
        Quartets (i-1, i, i+1, i+2) for dihedral angles.
        """
        n = len(self.atoms)
        return [
            (self.atoms[(i - 1) % n],
             self.atoms[i],
             self.atoms[(i + 1) % n],
             self.atoms[(i + 2) % n])
            for i in range(n)
        ]

    # ----------------------------------------------------------
    # Numerical values (optional, but useful)
    # ----------------------------------------------------------

    def valence_angles(self):
        """
        Return list of valence angles (radians) along the ring.
        """
        angles = []
        for (i, j, k) in self.cyclic_triplets():
            angles.append(self._angle(i, j, k))
        return angles

    def dihedral_angles(self):
        """
        Return list of dihedral angles (radians) along the ring.
        """
        dihedrals = []
        for (i, j, k, l) in self.cyclic_quartets():
            dihedrals.append(self._dihedral(i, j, k, l))
        return dihedrals

    # ----------------------------------------------------------
    # Geometry primitives
    # ----------------------------------------------------------

    def _angle(self, i, j, k):
        rji = self.coords[self.atoms.index(i)] - self.coords[self.atoms.index(j)]
        rjk = self.coords[self.atoms.index(k)] - self.coords[self.atoms.index(j)]
        return _vector_angle(rji, rjk)

    def _dihedral(self, i, j, k, l):
        p = [self.coords[self.atoms.index(x)] for x in (i, j, k, l)]
        return _dihedral_angle(*p)
