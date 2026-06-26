from rdkit import Chem
from rdkit.Chem import SanitizeFlags
import numpy as np

from .continuous_graph import bond_order
from .discrete_graph import DiscreteGraph
from .ringset import RingSet
from .atomic_synthons import AtomicSynthons
from .aromaticity import Aromaticity
from .topology_validator import validate_topology


class ContinuousGraphAdapter:
    def __init__(self, coords, Z):
        self.coords = coords
        self.Z = Z
        self.natoms = len(Z)

        neighbors = [list(range(self.natoms)) for _ in range(self.natoms)]
        for i in range(self.natoms):
            neighbors[i].remove(i)

        self.BO = np.zeros((self.natoms, self.natoms))
        cache = {}

        for i in range(self.natoms):
            for j in range(i + 1, self.natoms):
                bo = bond_order(i, j, Z, coords, neighbors, cache)
                self.BO[i, j] = self.BO[j, i] = bo


def build_rdkit_mol(coords, atomic_numbers, force_aromatic=True, **_):
    cg = ContinuousGraphAdapter(coords, atomic_numbers)
    dg = DiscreteGraph(cg)
    validate_topology(dg)

    rings = RingSet(dg, coords=coords)

    neighbors = [list(dg.adjacency[i]) for i in range(dg.natoms)]

    synthons = AtomicSynthons(
        Z=atomic_numbers,
        coords=coords,
        neighbors=neighbors,
    )

    arom = Aromaticity(
        cg,
        dg,
        rings,
        synthons=synthons,
        force_aromatic=force_aromatic,
    )

    mol = Chem.RWMol()
    atom_map = {}

    # --------------------------------------------------------
    # Identify H atoms bound to carbon (C–H) to be removed
    # --------------------------------------------------------
    skip_H = set()
    for i, Zi in enumerate(atomic_numbers):
        if int(Zi) != 1:
            continue
        for j in dg.adjacency[i]:
            if int(atomic_numbers[j]) == 6:
                skip_H.add(i)
                break

    # --------------------------------------------------------
    # Add atoms (unchanged, except skipping C–H hydrogens)
    # --------------------------------------------------------
    for i, Z in enumerate(atomic_numbers):
        if i in skip_H:
            continue

        atom = Chem.Atom(int(Z))
        atom.SetNoImplicit(True)

        # Canonical discrete signature
        atom.SetProp("synthon_signature", synthons.canonical_signature_str(i))

        Zc, A, NED, D = synthons.canonical_signature(i)

        atom.SetIntProp("synthon_Z", Zc)
        atom.SetBoolProp("synthon_aromatic", bool(A))
        atom.SetIntProp("synthon_NED", int(NED))
        atom.SetIntProp("synthon_degree", int(D))
        atom.SetDoubleProp("Zeff", synthons.Zeff(i))

        rd_idx = mol.AddAtom(atom)
        atom_map[i] = rd_idx

    # --------------------------------------------------------
    # Add bonds (skip those involving removed C–H hydrogens)
    # --------------------------------------------------------
    for (i, j) in dg.bonds:
        if i in skip_H or j in skip_H:
            continue

        rd_i, rd_j = atom_map[i], atom_map[j]

        if (i, j) in arom.aromatic_bonds:
            mol.AddBond(rd_i, rd_j, Chem.BondType.AROMATIC)
            mol.GetAtomWithIdx(rd_i).SetIsAromatic(True)
            mol.GetAtomWithIdx(rd_j).SetIsAromatic(True)
            continue

        bo = cg.BO[i, j]
        bond = (
            Chem.BondType.TRIPLE if bo >= 2.1 else
            Chem.BondType.DOUBLE if bo >= 1.3 else
            Chem.BondType.SINGLE
        )
        mol.AddBond(rd_i, rd_j, bond)

    Chem.SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_SETAROMATICITY)
    return mol.GetMol()


def xyz_to_smiles(coords, atomic_numbers, **kwargs):
    mol = build_rdkit_mol(coords, atomic_numbers, **kwargs)
    return Chem.MolToSmiles(
        mol,
        canonical=True,
        kekuleSmiles=False,
        allHsExplicit=False,
        allBondsExplicit=False,
        # originale
        # kekuleSmiles=True,
        # allHsExplicit=True,
        # allBondsExplicit=True,
    )
