from __future__ import annotations

from dataclasses import dataclass
import re

import numpy as np

from .harmonic import solve_wilson_gf


DEFAULT_LARGE_AMPLITUDE_FAMILIES = (
    "torsion",
    "ring_puckering",
    "butterfly",
    "oop",
    "linear",
    "special",
)


@dataclass(frozen=True)
class LargeAmplitudeCoordinate:
    index: int
    name: str
    irrep: str
    family: str
    label: str


@dataclass(frozen=True)
class LargeAmplitudeBlock:
    label: str
    family: str
    indices: tuple[int, ...]
    frequencies_cm: tuple[float, ...]
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
class GFLargeAmplitudeAnalysis:
    coordinates: tuple[LargeAmplitudeCoordinate, ...]
    blocks: tuple[LargeAmplitudeBlock, ...]
    mode_contributions: tuple[LargeAmplitudeModeContribution, ...]
    families: tuple[str, ...] = DEFAULT_LARGE_AMPLITUDE_FAMILIES

    @property
    def coordinate_count(self) -> int:
        return len(self.coordinates)


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
) -> GFLargeAmplitudeAnalysis:
    """Analyze large-amplitude GIC subspaces directly in the non-redundant GF basis."""
    f_mat = np.asarray(force_constants, dtype=float)
    g_mat = np.asarray(g_matrix, dtype=float)
    if f_mat.shape != g_mat.shape or f_mat.ndim != 2 or f_mat.shape[0] != f_mat.shape[1]:
        raise ValueError("large-amplitude analysis needs square F/G matrices of the same shape")
    ncoord = f_mat.shape[0]
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
        coordinates.append(
            LargeAmplitudeCoordinate(
                index=index + 1,
                name=names[index] or f"GIC{index + 1:03d}",
                irrep=irreps[index] or "UNK",
                family=family,
                label=labels[index],
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
    return GFLargeAmplitudeAnalysis(
        coordinates=tuple(coordinates),
        blocks=tuple(blocks),
        mode_contributions=tuple(mode_contributions),
        families=selected_families,
    )


def _large_amplitude_block(
    *,
    label: str,
    family: str,
    indices: tuple[int, ...],
    force_constants: np.ndarray,
    g_matrix: np.ndarray,
) -> LargeAmplitudeBlock:
    index = np.asarray(indices, dtype=int)
    f_sub = force_constants[np.ix_(index, index)]
    g_sub = g_matrix[np.ix_(index, index)]
    gf = solve_wilson_gf(f_sub, g_sub, scale_to_cm=True)
    f_coupling, f_relative = _coupling_to_rest(force_constants, indices)
    g_coupling, g_relative = _coupling_to_rest(g_matrix, indices)
    return LargeAmplitudeBlock(
        label=label,
        family=family,
        indices=tuple(int(value) + 1 for value in indices),
        frequencies_cm=tuple(float(value) for value in gf.frequencies_cm),
        max_f_coupling_to_rest=f_coupling,
        max_g_coupling_to_rest=g_coupling,
        relative_f_coupling_to_rest=f_relative,
        relative_g_coupling_to_rest=g_relative,
    )


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


def _padded(values: tuple[str, ...], size: int, *, fill: str = "") -> tuple[str, ...]:
    if len(values) >= size:
        return tuple(values[:size])
    return tuple(values) + tuple(fill for _ in range(size - len(values)))
