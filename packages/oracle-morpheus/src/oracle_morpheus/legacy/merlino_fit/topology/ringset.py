# ============================================================
# RingSet class for Merlino (ROBUST VERSION)
#
# Responsibilities:
#   - detect all simple cycles in the molecular graph
#   - optionally limit ring size (ring_max)
#   - build Ring objects
#   - manage atom-ring and bond-ring mappings
#   - detect fused / connected rings
#
# Notes:
#   - speed is not a priority
#   - accuracy and completeness are preferred
# ============================================================

from collections import defaultdict

try:
    from .rings import Ring
except ImportError:
    from rings import Ring


class RingSet:
    """
    Container and manager for all rings in a molecular graph.
    """

    def __init__(self, graph, coords=None, ring_max=None):
        """
        Parameters
        ----------
        graph : object
            Molecular graph object.
            Must provide:
              - graph.natoms or graph.n_atoms
              - graph.neighbors(atom)
        coords : ndarray, optional
            Cartesian coordinates of the system.
        ring_max : int, optional
            Maximum ring size to consider.
            If None, no size limit is applied.
        """
        self.graph = graph
        self.coords = coords
        self.ring_max = ring_max

        # --- storage ---
        self.rings = []
        self.atom_to_rings = defaultdict(list)
        self.bond_to_rings = defaultdict(list)

        # --- build ---
        self._detect_rings()
        self._build_connectivity()

    # ============================================================
    # Ring detection (ROBUST DFS)
    # ============================================================

    def _detect_rings(self):
        """
        Detect all simple cycles up to ring_max using DFS.
        """
        seen_cycles = set()
        ring_index = 0

        nat = getattr(self.graph, "natoms",
                      getattr(self.graph, "n_atoms", None))
        if nat is None:
            raise AttributeError("Graph must define natoms or n_atoms")

        for start in range(nat):
            self._dfs_cycles(
                start=start,
                current=start,
                visited=[start],
                seen_cycles=seen_cycles,
                ring_index_ref=[ring_index],
            )
            ring_index = len(self.rings)

    def _dfs_cycles(self, start, current, visited, seen_cycles, ring_index_ref):
        """
        Depth-first search for cycles.
        """
        if self.ring_max is not None and len(visited) > self.ring_max:
            return

        for nbr in self.graph.neighbors(current):
            if nbr == start and len(visited) >= 3:
                cycle = visited[:]
                canonical = self._canonical_cycle(cycle)

                if canonical in seen_cycles:
                    continue

                seen_cycles.add(canonical)

                atoms = list(canonical)
                ring = Ring(
                    index=ring_index_ref[0],
                    atoms=atoms,
                    coords=self.coords,
                )

                self.rings.append(ring)

                for a in atoms:
                    self.atom_to_rings[a].append(ring_index_ref[0])

                for b in ring.bonds:
                    self.bond_to_rings[b].append(ring_index_ref[0])

                ring_index_ref[0] += 1
                continue

            if nbr in visited:
                continue

            # enforce ordering to avoid mirrored duplicates
            if nbr < start:
                continue

            self._dfs_cycles(
                start=start,
                current=nbr,
                visited=visited + [nbr],
                seen_cycles=seen_cycles,
                ring_index_ref=ring_index_ref,
            )

    # ============================================================
    # Canonicalization
    # ============================================================

    def _canonical_cycle(self, cycle):
        """
        Return a canonical tuple representation of a cycle.
        """
        n = len(cycle)
        rotations = []

        for i in range(n):
            r1 = cycle[i:] + cycle[:i]
            r2 = list(reversed(r1))
            rotations.append(tuple(r1))
            rotations.append(tuple(r2))

        return min(rotations)

    # ============================================================
    # Ring connectivity
    # ============================================================

    def _build_connectivity(self):
        """
        Detect fused / connected rings.
        """
        n = len(self.rings)

        for i in range(n):
            ri = self.rings[i]
            for j in range(i + 1, n):
                rj = self.rings[j]

                if ri.shares_atoms_with(rj):
                    ri.add_connected_ring(j)
                    rj.add_connected_ring(i)

    # ============================================================
    # Public API
    # ============================================================

    def __len__(self):
        return len(self.rings)

    def __iter__(self):
        return iter(self.rings)

    def get_ring(self, index):
        return self.rings[index]

    def rings_of_atom(self, atom):
        return self.atom_to_rings.get(atom, [])

    def rings_of_bond(self, bond):
        return self.bond_to_rings.get(tuple(sorted(bond)), [])

    def fused_rings(self, index):
        """
        Return indices of rings fused to the given ring.
        """
        return list(self.rings[index].connected_rings)

    def __repr__(self):
        lim = self.ring_max if self.ring_max is not None else "∞"
        return f"<RingSet: {len(self.rings)} rings (max={lim})>"
