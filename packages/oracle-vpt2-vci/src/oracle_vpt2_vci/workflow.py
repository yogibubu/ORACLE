from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from oracle_gf import GFResult, mass_weighted_cartesian_hessian, solve_wilson_gf

from .gaussian_qff import FCHKData, anharmonic_input_from_gaussian_fchk, hessian_input_from_gaussian_fchk, read_gaussian_fchk_qff, read_indexed_qff_text
from .vci import VCIResult, solve_vci, zero_anharmonic_force_field


@dataclass(frozen=True)
class VPT2VCIRun:
    gaussian_data: FCHKData
    gf: GFResult
    vci: VCIResult


def run_python_vci_from_gaussian_fchk(
    path: Path,
    *,
    qff_path: Path | None = None,
    max_quanta: int = 2,
    n_roots: int = 10,
) -> VPT2VCIRun:
    """Run the independent Python GF/VCI path from a Gaussian FCHK file.

    If `qff_path` is provided, indexed cubic/quartic normal-coordinate terms
    are read from that file and used by VCI. Otherwise a harmonic VCI reference
    is run from the Gaussian frequencies.
    """
    data = read_gaussian_fchk_qff(Path(path))
    hessian_input = hessian_input_from_gaussian_fchk(Path(path))
    anharmonic_input = anharmonic_input_from_gaussian_fchk(Path(path))
    mw_hessian = mass_weighted_cartesian_hessian(hessian_input.cartesian_hessian, hessian_input.masses_amu)
    gf = solve_wilson_gf(mw_hessian, np.eye(mw_hessian.shape[0]))
    frequencies = (
        anharmonic_input.anharmonic_frequencies_cm
        if anharmonic_input.anharmonic_frequencies_cm.size
        else anharmonic_input.harmonic_frequencies_cm
    )
    qff = read_indexed_qff_text(qff_path, frequencies) if qff_path is not None else zero_anharmonic_force_field(frequencies)
    vci = solve_vci(qff, max_quanta=max_quanta, n_roots=n_roots)
    return VPT2VCIRun(gaussian_data=data, gf=gf, vci=vci)
