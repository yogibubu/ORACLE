from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np

from .models import AnharmonicInput
from .vci import (
    QuarticForceField,
    VCIOptions,
    VCIResult,
    _mode_powers,
    _selected_force_field,
    _selected_quanta_limits,
    _term_element,
    _x_matrix_power,
    force_field_from_anharmonic_input,
    generate_vibrational_basis,
    solve_vci,
)


@dataclass(frozen=True)
class VPT2State:
    basis_state: tuple[int, ...]
    harmonic_cm: float
    first_order_cm: float
    second_order_cm: float

    @property
    def energy_cm(self) -> float:
        return self.harmonic_cm + self.first_order_cm + self.second_order_cm


@dataclass(frozen=True)
class VPT2Result:
    states: tuple[VPT2State, ...]
    basis: tuple[tuple[int, ...], ...]
    options: VCIOptions

    @property
    def energies_cm(self) -> np.ndarray:
        return np.array([state.energy_cm for state in self.states], dtype=float)

    @property
    def excitation_energies_cm(self) -> np.ndarray:
        energies = self.energies_cm
        return energies - energies[0]


@dataclass(frozen=True)
class VPT2VCIComparison:
    vpt2: VPT2Result
    vci: VCIResult
    energy_differences_cm: np.ndarray
    excitation_differences_cm: np.ndarray


def _harmonic_energy(freqs: np.ndarray, state: tuple[int, ...]) -> float:
    return float(np.dot(freqs, np.array(state, dtype=float) + 0.5))


def _term_matrix_element(
    indices: tuple[int, ...],
    left: tuple[int, ...],
    right: tuple[int, ...],
    powers: dict[int, np.ndarray],
) -> float:
    mode_powers = _mode_powers(indices, len(left))
    return _term_element(left, right, tuple(powers[p] for p in mode_powers))


def _intermediate_basis(n_modes: int, max_quanta_per_mode: int) -> tuple[tuple[int, ...], ...]:
    return tuple(product(range(max_quanta_per_mode + 1), repeat=n_modes))


def solve_vpt2(
    force_field: QuarticForceField,
    max_quanta: int,
    n_roots: int | None = None,
    *,
    options: VCIOptions | None = None,
    intermediate_extra_quanta: int = 4,
    denominator_tolerance_cm: float = 1.0e-8,
) -> VPT2Result:
    """Compute VPT2 energies on the canonical ORACLE QFF model.

    The implementation is deliberately solver-format independent. Quartic terms
    contribute at first order and cubic terms at second order.
    """
    opts = options or VCIOptions()
    opts.validate()
    mode_min, mode_max = _selected_quanta_limits(force_field, opts)
    qff = _selected_force_field(force_field, opts)
    freqs = np.asarray(qff.harmonic_frequencies_cm, dtype=float)
    basis = generate_vibrational_basis(
        len(freqs),
        max_quanta,
        freqs,
        opts.basis_energy_cutoff_cm,
        opts.max_basis_states,
        mode_min,
        mode_max,
        opts.excitation_class_limits,
    )
    if n_roots is not None:
        basis = basis[:n_roots]
    max_n = max(max(state) for state in basis) + intermediate_extra_quanta
    powers = {power: _x_matrix_power(max_n, power) for power in range(5)}
    intermediates = _intermediate_basis(len(freqs), max_n)

    quartic_terms = tuple((tuple(sorted(k)), float(v)) for k, v in qff.quartic_cm.items())
    cubic_terms = tuple((tuple(sorted(k)), float(v)) for k, v in qff.cubic_cm.items())
    states: list[VPT2State] = []
    for state in basis:
        harmonic = _harmonic_energy(freqs, state)
        first = 0.0
        for indices, coeff in quartic_terms:
            first += coeff * _term_matrix_element(indices, state, state, powers)

        second = 0.0
        for other in intermediates:
            if other == state:
                continue
            coupling = 0.0
            for indices, coeff in cubic_terms:
                coupling += coeff * _term_matrix_element(indices, other, state, powers)
            if coupling == 0.0:
                continue
            denom = harmonic - _harmonic_energy(freqs, other)
            if abs(denom) <= denominator_tolerance_cm:
                continue
            second += coupling * coupling / denom
        states.append(VPT2State(state, harmonic, first, second))
    return VPT2Result(states=tuple(states), basis=basis, options=opts)


def solve_vpt2_from_anharmonic_input(
    input_data: AnharmonicInput,
    max_quanta: int,
    n_roots: int | None = None,
    **kwargs: object,
) -> VPT2Result:
    return solve_vpt2(
        force_field_from_anharmonic_input(input_data),
        max_quanta=max_quanta,
        n_roots=n_roots,
        **kwargs,
    )


def compare_vpt2_vci(
    force_field: QuarticForceField,
    max_quanta: int,
    n_roots: int,
    *,
    options: VCIOptions | None = None,
    vci_method: str = "dense",
) -> VPT2VCIComparison:
    vpt2 = solve_vpt2(force_field, max_quanta=max_quanta, n_roots=n_roots, options=options)
    vci = solve_vci(force_field, max_quanta=max_quanta, n_roots=n_roots, options=options, method=vci_method)
    n = min(len(vpt2.energies_cm), len(vci.energies_cm))
    return VPT2VCIComparison(
        vpt2=vpt2,
        vci=vci,
        energy_differences_cm=vci.energies_cm[:n] - vpt2.energies_cm[:n],
        excitation_differences_cm=vci.excitation_energies_cm[:n] - vpt2.excitation_energies_cm[:n],
    )


def compare_vpt2_vci_from_anharmonic_input(
    input_data: AnharmonicInput,
    max_quanta: int,
    n_roots: int,
    **kwargs: object,
) -> VPT2VCIComparison:
    return compare_vpt2_vci(
        force_field_from_anharmonic_input(input_data),
        max_quanta=max_quanta,
        n_roots=n_roots,
        **kwargs,
    )
