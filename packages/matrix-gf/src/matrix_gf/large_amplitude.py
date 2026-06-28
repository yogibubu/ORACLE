from __future__ import annotations

from dataclasses import dataclass, field
import math
import re

import numpy as np

from .harmonic import CM_PER_HARTREE, solve_wilson_gf


DEFAULT_LARGE_AMPLITUDE_FREQUENCY_CUTOFF_CM = 250.0
DEFAULT_LARGE_AMPLITUDE_FG_COUPLING_TOLERANCE = 0.05
DEFAULT_TORSION_DOUBLE_BOND_ORDER_THRESHOLD = 1.45
DEFAULT_G_INVERSE_RCOND = 1.0e-12
DEFAULT_LARGE_AMPLITUDE_FAMILIES = (
    "torsion",
    "ring_puckering",
    "butterfly",
    "oop",
    "linear",
    "special",
)


@dataclass(frozen=True)
class LargeAmplitudeTopologyContext:
    """Topology/synthon context fixed by preprocessing and consumed by GF only."""

    atomic_numbers: tuple[int, ...] = ()
    bonds: tuple[tuple[int, int], ...] = ()
    bond_orders: dict[tuple[int, int], float] = field(default_factory=dict)
    synthon_signatures: dict[int, tuple[object, ...]] = field(default_factory=dict)
    coordinates_angstrom: tuple[tuple[float, float, float], ...] = ()
    bond_order_source: str = ""

    def graph(self) -> dict[int, set[int]]:
        graph: dict[int, set[int]] = {}
        for left, right in self.bonds:
            graph.setdefault(left, set()).add(right)
            graph.setdefault(right, set()).add(left)
        return graph

    def bond_order(self, left: int, right: int) -> float | None:
        key = _pair_key(left, right)
        return self.bond_orders.get(key)

    def synthon_signature(self, atom: int) -> tuple[object, ...]:
        signature = self.synthon_signatures.get(atom)
        if signature is not None:
            return signature
        if 1 <= atom <= len(self.atomic_numbers):
            return (self.atomic_numbers[atom - 1],)
        return (atom,)


@dataclass(frozen=True)
class LargeAmplitudeCoordinate:
    index: int
    name: str
    irrep: str
    family: str
    label: str
    local_frequency_cm: float | None = None
    active: bool = True
    status: str = "ACTIVE"


@dataclass(frozen=True)
class LargeAmplitudeBlock:
    label: str
    family: str
    indices: tuple[int, ...]
    frequencies_cm: tuple[float, ...]
    g_inverse_block: tuple[tuple[float, ...], ...]
    g_inverse_source: str
    max_f_coupling_to_rest: float
    max_g_coupling_to_rest: float
    relative_f_coupling_to_rest: float
    relative_g_coupling_to_rest: float

    @property
    def relative_fg_coupling_to_rest(self) -> float:
        """Dimensionless diagnostic that treats F and G couplings symmetrically."""
        return max(self.relative_f_coupling_to_rest, self.relative_g_coupling_to_rest)


@dataclass(frozen=True)
class LargeAmplitudeModeContribution:
    mode: int
    frequency_cm: float
    ped_percent: float


@dataclass(frozen=True)
class LargeAmplitudeDVRCandidate:
    index: int
    name: str
    irrep: str
    family: str
    status: str
    frequency_cm: float | None
    fg_coupling_to_rest: float
    central_bond: tuple[int, int] | None = None
    periodicity: int | None = None
    minimum_rad: float | None = None
    force_constant_hartree: float | None = None
    fourier_amplitude_cm: float | None = None
    barrier_cm: float | None = None
    g_inverse_diagonal: float | None = None
    g_inverse_source: str = ""
    reason: str = ""


@dataclass(frozen=True)
class GFLargeAmplitudeAnalysis:
    coordinates: tuple[LargeAmplitudeCoordinate, ...]
    blocks: tuple[LargeAmplitudeBlock, ...]
    mode_contributions: tuple[LargeAmplitudeModeContribution, ...]
    dvr_candidates: tuple[LargeAmplitudeDVRCandidate, ...] = ()
    g_inverse: tuple[tuple[float, ...], ...] = ()
    g_inverse_source: str = ""
    families: tuple[str, ...] = DEFAULT_LARGE_AMPLITUDE_FAMILIES
    frequency_cutoff_cm: float | None = DEFAULT_LARGE_AMPLITUDE_FREQUENCY_CUTOFF_CM
    fg_coupling_tolerance: float = DEFAULT_LARGE_AMPLITUDE_FG_COUPLING_TOLERANCE

    @property
    def coordinate_count(self) -> int:
        return len(self.coordinates)

    @property
    def active_coordinate_count(self) -> int:
        return sum(1 for coordinate in self.coordinates if coordinate.active)


def gic_coordinate_family(name: str, label: str) -> str:
    """Return the MATRIX coordinate family encoded in a frozen GIC name/label."""
    text = f"{name} {label}".lower()
    if any(token in text for token in ("rpck", "qpck", "phip", "pck")):
        return "ring_puckering"
    if any(token in text for token in ("btfl", "butterfly")):
        return "butterfly"
    if any(token in text for token in ("cybe", "cyclic_bend", "ring_bend")):
        return "ring_bend"
    if any(token in text for token in ("fragment", "frag", "centroid", "center", "centre")):
        return "special"
    if re.search(r"\br\s*\(", text) or "str" in text or "stretch" in text or "bond(" in text:
        return "stretch"
    if re.search(r"\ba\s*\(", text) or "bend" in text or "angle(" in text:
        return "bend"
    if re.search(r"\bd\s*\(", text) or "tors" in text or "dih" in text or "dihedral(" in text:
        return "torsion"
    if (
        re.search(r"\bu\s*\(", text)
        or "oupl" in text
        or "oop" in text
        or "improper" in text
        or "out_of_plane(" in text
    ):
        return "oop"
    if re.search(r"\bl\s*\(", text) or "linear_bend" in text or "lin" in text:
        return "linear"
    return ""


def large_amplitude_analysis_from_gf_matrices(
    *,
    force_constants: np.ndarray,
    g_matrix: np.ndarray,
    frequencies_cm: np.ndarray,
    ped: np.ndarray,
    gic_labels: tuple[str, ...],
    gic_names: tuple[str, ...] = (),
    gic_irreps: tuple[str, ...] = (),
    families: tuple[str, ...] = DEFAULT_LARGE_AMPLITUDE_FAMILIES,
    frequency_cutoff_cm: float | None = DEFAULT_LARGE_AMPLITUDE_FREQUENCY_CUTOFF_CM,
    fg_coupling_tolerance: float = DEFAULT_LARGE_AMPLITUDE_FG_COUPLING_TOLERANCE,
    topology_context: LargeAmplitudeTopologyContext | None = None,
) -> GFLargeAmplitudeAnalysis:
    """Analyze large-amplitude GIC subspaces directly in the non-redundant GF basis."""
    f_mat = np.asarray(force_constants, dtype=float)
    g_mat = np.asarray(g_matrix, dtype=float)
    if f_mat.shape != g_mat.shape or f_mat.ndim != 2 or f_mat.shape[0] != f_mat.shape[1]:
        raise ValueError("large-amplitude analysis needs square F/G matrices of the same shape")
    ncoord = f_mat.shape[0]
    g_inverse, g_inverse_source = _invert_g_matrix(g_mat)
    names = _padded(gic_names, ncoord)
    irreps = _padded(gic_irreps, ncoord, fill="UNK")
    labels = _padded(gic_labels, ncoord)
    selected_families = tuple(dict.fromkeys(family for family in families if family))
    selected = set(selected_families)
    coordinates = []
    by_family: dict[str, list[int]] = {family: [] for family in selected_families}
    for index in range(ncoord):
        family = gic_coordinate_family(names[index], labels[index])
        if family not in selected:
            continue
        local_frequency = _local_coordinate_frequency_cm(f_mat, g_mat, index)
        active, status = _coordinate_activity_status(
            local_frequency,
            frequency_cutoff_cm=frequency_cutoff_cm,
        )
        coordinates.append(
            LargeAmplitudeCoordinate(
                index=index + 1,
                name=names[index] or f"GIC{index + 1:03d}",
                irrep=irreps[index] or "UNK",
                family=family,
                label=labels[index],
                local_frequency_cm=local_frequency,
                active=active,
                status=status,
            )
        )
        by_family.setdefault(family, []).append(index)

    blocks = [
        _large_amplitude_block(
            label=family,
            family=family,
            indices=tuple(indices),
            force_constants=f_mat,
            g_matrix=g_mat,
            g_inverse=g_inverse,
            g_inverse_source=g_inverse_source,
        )
        for family, indices in by_family.items()
        if indices
    ]
    all_indices = tuple(coordinate.index - 1 for coordinate in coordinates)
    if len(blocks) > 1 and all_indices:
        blocks.append(
            _large_amplitude_block(
                label="all_large_amplitude",
                family="mixed_large_amplitude",
                indices=all_indices,
                force_constants=f_mat,
                g_matrix=g_mat,
                g_inverse=g_inverse,
                g_inverse_source=g_inverse_source,
            )
        )

    ped_values = np.asarray(ped, dtype=float)
    mode_contributions: list[LargeAmplitudeModeContribution] = []
    if all_indices and ped_values.ndim == 2 and ped_values.shape[0] == ncoord:
        freqs = np.asarray(frequencies_cm, dtype=float).reshape(-1)
        for mode in range(min(ped_values.shape[1], freqs.size)):
            mode_contributions.append(
                LargeAmplitudeModeContribution(
                    mode=mode + 1,
                    frequency_cm=float(freqs[mode]),
                    ped_percent=float(np.sum(ped_values[list(all_indices), mode])),
                )
            )
    dvr_candidates = tuple(
        _dvr_candidate(
            coordinate,
            force_constants=f_mat,
            g_matrix=g_mat,
            g_inverse=g_inverse,
            g_inverse_source=g_inverse_source,
            topology_context=topology_context,
            fg_coupling_tolerance=fg_coupling_tolerance,
        )
        for coordinate in coordinates
    )
    return GFLargeAmplitudeAnalysis(
        coordinates=tuple(coordinates),
        blocks=tuple(blocks),
        mode_contributions=tuple(mode_contributions),
        dvr_candidates=dvr_candidates,
        g_inverse=tuple(tuple(float(value) for value in row) for row in g_inverse),
        g_inverse_source=g_inverse_source,
        families=selected_families,
        frequency_cutoff_cm=frequency_cutoff_cm,
        fg_coupling_tolerance=float(fg_coupling_tolerance),
    )


def large_amplitude_topology_context_from_arrays(
    *,
    atomic_numbers: np.ndarray,
    coordinates_angstrom: np.ndarray,
    bond_order_overrides: dict[tuple[int, int], float] | None = None,
    bond_order_source: str = "",
) -> LargeAmplitudeTopologyContext:
    """Build the large-amplitude context with the shared MATRIX topology engine."""
    from matrix_chem import build_topology_objects

    numbers = tuple(int(value) for value in np.asarray(atomic_numbers, dtype=int).reshape(-1))
    coords = np.asarray(coordinates_angstrom, dtype=float)
    _continuous, discrete, _ringset, synthons, _aromaticity = build_topology_objects(
        coords,
        numbers,
        bond_order_overrides=bond_order_overrides,
        bond_order_source=bond_order_source or "Topology Pauling continuous model",
    )
    bonds = tuple((int(i) + 1, int(j) + 1) for i, j in getattr(discrete, "bonds", ()))
    bond_orders: dict[tuple[int, int], float] = {}
    for i in range(len(numbers)):
        for j in range(i + 1, len(numbers)):
            value = float(synthons.bond_order(i, j))
            if value > 0.05:
                bond_orders[(i + 1, j + 1)] = value
    signatures = {
        idx + 1: tuple(synthons.canonical_signature(idx)) for idx in range(len(numbers))
    }
    return LargeAmplitudeTopologyContext(
        atomic_numbers=numbers,
        bonds=bonds,
        bond_orders=bond_orders,
        synthon_signatures=signatures,
        coordinates_angstrom=tuple(tuple(float(x) for x in row) for row in coords),
        bond_order_source=bond_order_source,
    )


def _large_amplitude_block(
    *,
    label: str,
    family: str,
    indices: tuple[int, ...],
    force_constants: np.ndarray,
    g_matrix: np.ndarray,
    g_inverse: np.ndarray,
    g_inverse_source: str,
) -> LargeAmplitudeBlock:
    index = np.asarray(indices, dtype=int)
    f_sub = force_constants[np.ix_(index, index)]
    g_sub = g_matrix[np.ix_(index, index)]
    g_inv_sub = g_inverse[np.ix_(index, index)]
    gf = solve_wilson_gf(f_sub, g_sub, scale_to_cm=True)
    f_coupling, f_relative = _coupling_to_rest(force_constants, indices)
    g_coupling, g_relative = _coupling_to_rest(g_matrix, indices)
    return LargeAmplitudeBlock(
        label=label,
        family=family,
        indices=tuple(int(value) + 1 for value in indices),
        frequencies_cm=tuple(float(value) for value in gf.frequencies_cm),
        g_inverse_block=tuple(tuple(float(value) for value in row) for row in g_inv_sub),
        g_inverse_source=g_inverse_source,
        max_f_coupling_to_rest=f_coupling,
        max_g_coupling_to_rest=g_coupling,
        relative_f_coupling_to_rest=f_relative,
        relative_g_coupling_to_rest=g_relative,
    )


def _dvr_candidate(
    coordinate: LargeAmplitudeCoordinate,
    *,
    force_constants: np.ndarray,
    g_matrix: np.ndarray,
    g_inverse: np.ndarray,
    g_inverse_source: str,
    topology_context: LargeAmplitudeTopologyContext | None,
    fg_coupling_tolerance: float,
) -> LargeAmplitudeDVRCandidate:
    index0 = coordinate.index - 1
    _f_coupling, f_relative = _coupling_to_rest(force_constants, (index0,))
    _g_coupling, g_relative = _coupling_to_rest(g_matrix, (index0,))
    fg_relative = max(f_relative, g_relative)
    base = {
        "index": coordinate.index,
        "name": coordinate.name,
        "irrep": coordinate.irrep,
        "family": coordinate.family,
        "frequency_cm": coordinate.local_frequency_cm,
        "fg_coupling_to_rest": fg_relative,
        "g_inverse_diagonal": _g_inverse_diagonal(g_inverse, index0),
        "g_inverse_source": g_inverse_source,
    }
    if not coordinate.active:
        return LargeAmplitudeDVRCandidate(
            **base,
            status=coordinate.status,
            reason=coordinate.status,
        )
    if coordinate.family != "torsion":
        return LargeAmplitudeDVRCandidate(
            **base,
            status="ACTIVE_BLOCK_DVR",
            reason="NON_TORSIONAL_LARGE_AMPLITUDE",
        )
    torsions = _torsion_terms(f"{coordinate.name} {coordinate.label}")
    if not torsions:
        return LargeAmplitudeDVRCandidate(
            **base,
            status="ACTIVE_TORSION_NO_PRIMITIVE",
            reason="NO_D_PRIMITIVE_LABEL",
        )
    central_bonds = tuple(dict.fromkeys(_pair_key(term[1], term[2]) for term in torsions))
    if len(central_bonds) != 1:
        return LargeAmplitudeDVRCandidate(
            **base,
            status="ACTIVE_COUPLED_TORSIONS",
            reason="MULTIPLE_CENTRAL_BONDS",
        )
    central_bond = central_bonds[0]
    if topology_context is None:
        return LargeAmplitudeDVRCandidate(
            **base,
            status="PENDING_TOPOLOGY",
            central_bond=central_bond,
            reason="NO_TOPOLOGY_CONTEXT",
        )
    bond_order = topology_context.bond_order(*central_bond)
    if bond_order is not None and bond_order >= DEFAULT_TORSION_DOUBLE_BOND_ORDER_THRESHOLD:
        return LargeAmplitudeDVRCandidate(
            **base,
            status="EXCLUDED_HIGH_BOND_ORDER",
            central_bond=central_bond,
            reason=f"BOND_ORDER_GE_{DEFAULT_TORSION_DOUBLE_BOND_ORDER_THRESHOLD:g}",
        )
    periodicity = _torsion_periodicity(topology_context, central_bond)
    minimum_rad = _torsion_minimum_rad(topology_context, torsions[0])
    if fg_relative > float(fg_coupling_tolerance):
        return LargeAmplitudeDVRCandidate(
            **base,
            status="ACTIVE_COUPLED_DVR",
            central_bond=central_bond,
            periodicity=periodicity,
            minimum_rad=minimum_rad,
            reason="FG_COUPLING_ABOVE_TOLERANCE",
        )
    force_constant = float(force_constants[index0, index0])
    if force_constant <= 0.0:
        return LargeAmplitudeDVRCandidate(
            **base,
            status="ACTIVE_TORSION_UNSTABLE",
            central_bond=central_bond,
            periodicity=periodicity,
            minimum_rad=minimum_rad,
            force_constant_hartree=force_constant,
            reason="NON_POSITIVE_CURVATURE",
        )
    amplitude_cm = force_constant * CM_PER_HARTREE / float(periodicity * periodicity)
    return LargeAmplitudeDVRCandidate(
        **base,
        status="ACTIVE_TORSION_1D",
        central_bond=central_bond,
        periodicity=periodicity,
        minimum_rad=minimum_rad,
        force_constant_hartree=force_constant,
        fourier_amplitude_cm=amplitude_cm,
        barrier_cm=2.0 * amplitude_cm,
        reason="ONE_TERM_FOURIER",
    )


def _invert_g_matrix(g_matrix: np.ndarray) -> tuple[np.ndarray, str]:
    values = np.asarray(g_matrix, dtype=float)
    symmetric = 0.5 * (values + values.T)
    try:
        inverse = np.linalg.inv(symmetric)
        source = "G_INVERSE"
        if not np.all(np.isfinite(inverse)):
            raise np.linalg.LinAlgError("non-finite inverse")
    except np.linalg.LinAlgError:
        inverse = np.linalg.pinv(symmetric, rcond=DEFAULT_G_INVERSE_RCOND)
        source = "G_PSEUDOINVERSE"
    return 0.5 * (inverse + inverse.T), source


def _g_inverse_diagonal(g_inverse: np.ndarray, index: int) -> float | None:
    if index < 0 or index >= g_inverse.shape[0]:
        return None
    value = float(g_inverse[index, index])
    if not math.isfinite(value):
        return None
    return value


def _local_coordinate_frequency_cm(
    force_constants: np.ndarray,
    g_matrix: np.ndarray,
    index: int,
) -> float | None:
    try:
        gf = solve_wilson_gf(
            np.asarray([[force_constants[index, index]]], dtype=float),
            np.asarray([[g_matrix[index, index]]], dtype=float),
            scale_to_cm=True,
        )
    except ValueError:
        return None
    if gf.frequencies_cm.size == 0:
        return None
    return float(gf.frequencies_cm[0])


def _coordinate_activity_status(
    frequency_cm: float | None,
    *,
    frequency_cutoff_cm: float | None,
) -> tuple[bool, str]:
    if frequency_cm is None:
        return False, "EXCLUDED_INVALID_LOCAL_GF"
    if frequency_cutoff_cm is None:
        return True, "ACTIVE_NO_FREQUENCY_CUTOFF"
    cutoff = float(frequency_cutoff_cm)
    if abs(float(frequency_cm)) <= cutoff:
        return True, "ACTIVE"
    return False, "EXCLUDED_FREQUENCY"


def _coupling_to_rest(matrix: np.ndarray, indices: tuple[int, ...]) -> tuple[float, float]:
    values = np.asarray(matrix, dtype=float)
    selected = np.zeros(values.shape[0], dtype=bool)
    selected[list(indices)] = True
    rest = ~selected
    if not np.any(rest):
        return 0.0, 0.0
    coupling = values[np.ix_(selected, rest)]
    max_abs = float(np.max(np.abs(coupling))) if coupling.size else 0.0
    diag_scale = float(np.max(np.abs(np.diag(values)))) if values.size else 0.0
    return max_abs, max_abs / max(diag_scale, 1.0e-30)


def _torsion_terms(text: str) -> tuple[tuple[int, int, int, int], ...]:
    matches = re.findall(
        r"\bD\s*\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)",
        text,
        flags=re.IGNORECASE,
    )
    return tuple(tuple(int(value) for value in match) for match in matches)


def _torsion_periodicity(
    context: LargeAmplitudeTopologyContext,
    central_bond: tuple[int, int],
) -> int:
    left, right = central_bond
    graph = context.graph()
    left_period = _substituent_equivalence_period(
        sorted(atom for atom in graph.get(left, set()) if atom != right),
        context,
    )
    right_period = _substituent_equivalence_period(
        sorted(atom for atom in graph.get(right, set()) if atom != left),
        context,
    )
    return max(1, math.lcm(max(1, left_period), max(1, right_period)))


def _substituent_equivalence_period(
    atoms: list[int],
    context: LargeAmplitudeTopologyContext,
) -> int:
    if not atoms:
        return 1
    classes: dict[tuple[object, ...], int] = {}
    for atom in atoms:
        signature = context.synthon_signature(atom)
        classes[signature] = classes.get(signature, 0) + 1
    return max(classes.values(), default=1)


def _torsion_minimum_rad(
    context: LargeAmplitudeTopologyContext,
    atoms: tuple[int, int, int, int],
) -> float | None:
    coords = np.asarray(context.coordinates_angstrom, dtype=float)
    if coords.shape != (len(context.atomic_numbers), 3):
        return None
    if any(atom < 1 or atom > len(context.atomic_numbers) for atom in atoms):
        return None
    p0, p1, p2, p3 = (coords[atom - 1] for atom in atoms)
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    norm = float(np.linalg.norm(b1))
    if norm <= 1.0e-12:
        return None
    b1 = b1 / norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    v_norm = float(np.linalg.norm(v))
    w_norm = float(np.linalg.norm(w))
    if v_norm <= 1.0e-12 or w_norm <= 1.0e-12:
        return None
    v = v / v_norm
    w = w / w_norm
    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1, v), w))
    return math.atan2(y, x)


def _pair_key(left: int, right: int) -> tuple[int, int]:
    return (int(left), int(right)) if int(left) < int(right) else (int(right), int(left))


def _padded(values: tuple[str, ...], size: int, *, fill: str = "") -> tuple[str, ...]:
    if len(values) >= size:
        return tuple(values[:size])
    return tuple(values) + tuple(fill for _ in range(size - len(values)))
