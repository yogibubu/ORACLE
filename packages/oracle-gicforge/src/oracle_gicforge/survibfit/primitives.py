from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

import numpy as np

from .geometry import angle, dihedral, oop, linear_components
from .bmat import bond_grad, angle_grad, dihedral_grad, oop_grad, linear_grad, finite_diff_grad


@dataclass(frozen=True)
class Primitive:
    kind: str  # bond, angle, dihedral, linear_bend, out_of_plane, frag_trans, frag_rot
    atoms: Tuple[int, ...]
    mode: int = 0  # for linear_bend: -1 or -2 to distinguish components; for frag: axis index
    ref: Tuple[int, ...] = ()


def build_primitives(discrete_graph, coords, linear_threshold=np.deg2rad(170.0)) -> List[Primitive]:
    """Generate primitive internal coordinates from discrete topology.

    Uses bonds/angles/dihedrals, plus linear-bend and out-of-plane primitives.
    """
    n = discrete_graph.natoms
    bonds = [Primitive("bond", (i, j)) for (i, j) in discrete_graph.bonds]

    # angles (including linear detection)
    angles = []
    linears = []
    for j in range(n):
        neigh = sorted(discrete_graph.adjacency[j])
        for a in range(len(neigh)):
            for b in range(a + 1, len(neigh)):
                i = neigh[a]
                k = neigh[b]
                ang = angle(i, j, k, coords)
                if ang >= linear_threshold:
                    linears.append(Primitive("linear_bend", (i, j, k), mode=-1))
                    linears.append(Primitive("linear_bend", (i, j, k), mode=-2))
                else:
                    angles.append(Primitive("angle", (i, j, k)))

    # dihedrals
    diheds = []
    for j in range(n):
        for k in discrete_graph.adjacency[j]:
            if k == j:
                continue
            # neighbors of j excluding k
            neigh_j = [i for i in discrete_graph.adjacency[j] if i != k]
            # neighbors of k excluding j
            neigh_k = [l for l in discrete_graph.adjacency[k] if l != j]
            for i in neigh_j:
                for l in neigh_k:
                    diheds.append(Primitive("dihedral", (i, j, k, l)))

    # out-of-plane
    oops = []
    for j in range(n):
        neigh = list(discrete_graph.adjacency[j])
        if len(neigh) < 3:
            continue
        for i in neigh:
            others = [x for x in neigh if x != i]
            for a in range(len(others)):
                for b in range(a + 1, len(others)):
                    k = others[a]
                    l = others[b]
                    oops.append(Primitive("out_of_plane", (i, j, k, l)))

    return bonds + angles + linears + diheds + oops


def eval_primitive(p: Primitive, coords: np.ndarray):
    if p.kind == "bond":
        i, j = p.atoms
        return np.linalg.norm(coords[i] - coords[j])
    if p.kind == "angle":
        i, j, k = p.atoms
        return angle(i, j, k, coords)
    if p.kind == "dihedral":
        i, j, k, l = p.atoms
        return dihedral(i, j, k, l, coords)
    if p.kind == "out_of_plane":
        i, j, k, l = p.atoms
        return oop(i, j, k, l, coords)
    if p.kind == "linear_bend":
        i, j, k = p.atoms
        c1, c2 = linear_components(i, j, k, coords)
        return c1 if p.mode == -1 else c2
    if p.kind == "frag_trans":
        frag = np.array(p.atoms, dtype=int)
        ref = np.array(p.ref, dtype=int)
        c_frag = coords[frag].mean(axis=0)
        c_ref = coords[ref].mean(axis=0)
        return (c_frag - c_ref)[p.mode]
    if p.kind == "frag_rot":
        frag = np.array(p.atoms, dtype=int)
        ref = np.array(p.ref, dtype=int)

        def _frame(atoms_idx):
            # Geometric centroid + principal axes of inertia (TRIC convention)
            x = coords[atoms_idx] - coords[atoms_idx].mean(axis=0)
            I = np.zeros((3, 3))
            for v in x:
                I += (np.dot(v, v) * np.eye(3) - np.outer(v, v))
            evals, evecs = np.linalg.eigh(I)
            order = np.argsort(evals)
            F = evecs[:, order]
            # right-handed frame
            if np.linalg.det(F) < 0.0:
                F[:, -1] *= -1.0
            return F

        F_ref = _frame(ref)
        F_frag = _frame(frag)
        R = F_ref.T @ F_frag
        tr = np.trace(R)
        if tr > 0.0:
            s = 0.5 / np.sqrt(tr + 1.0)
            qw = 0.25 / s
            qx = (R[2, 1] - R[1, 2]) * s
            qy = (R[0, 2] - R[2, 0]) * s
            qz = (R[1, 0] - R[0, 1]) * s
        else:
            if R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
                s = 2.0 * np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2])
                qw = (R[2, 1] - R[1, 2]) / s
                qx = 0.25 * s
                qy = (R[0, 1] + R[1, 0]) / s
                qz = (R[0, 2] + R[2, 0]) / s
            elif R[1, 1] > R[2, 2]:
                s = 2.0 * np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2])
                qw = (R[0, 2] - R[2, 0]) / s
                qx = (R[0, 1] + R[1, 0]) / s
                qy = 0.25 * s
                qz = (R[1, 2] + R[2, 1]) / s
            else:
                s = 2.0 * np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1])
                qw = (R[1, 0] - R[0, 1]) / s
                qx = (R[0, 2] + R[2, 0]) / s
                qy = (R[1, 2] + R[2, 1]) / s
                qz = 0.25 * s
        quat = np.array([qw, qx, qy, qz])
        if quat[0] < 0.0:
            quat = -quat
        # rotation vector from quaternion
        qw, qx, qy, qz = quat
        v = np.array([qx, qy, qz])
        vn = np.linalg.norm(v)
        if vn < 1e-12:
            rotvec = np.zeros(3)
        else:
            ang = 2.0 * np.arctan2(vn, qw)
            rotvec = v / vn * ang
        return rotvec[p.mode]
    raise ValueError(f"Unknown primitive kind: {p.kind}")


def eval_primitives(prims, coords: np.ndarray):
    return np.array([eval_primitive(p, coords) for p in prims], dtype=float)


def grad_primitive(p: Primitive, coords: np.ndarray, fd_step=1e-4):
    if p.kind == "bond":
        i, j = p.atoms
        return bond_grad(i, j, coords)
    if p.kind == "angle":
        i, j, k = p.atoms
        return angle_grad(i, j, k, coords)
    if p.kind == "dihedral":
        i, j, k, l = p.atoms
        return dihedral_grad(i, j, k, l, coords)
    if p.kind == "out_of_plane":
        i, j, k, l = p.atoms
        return oop_grad(i, j, k, l, coords)
    if p.kind == "linear_bend":
        i, j, k = p.atoms
        return linear_grad(i, j, k, coords, mode=p.mode)
    if p.kind in ("frag_trans", "frag_rot"):
        return finite_diff_grad(lambda x: eval_primitive(p, x), coords, h=fd_step)
    raise ValueError(f"Unknown primitive kind: {p.kind}")
