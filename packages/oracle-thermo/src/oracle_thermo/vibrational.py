from __future__ import annotations

from collections.abc import Iterable
from math import exp

import numpy as np

from oracle_chem import Phy, get_physical_constants

from .models import ThermoContribution


def vibrational_thermo(
    frequencies_cm1: Iterable[float],
    *,
    T_K: float,
    cutoff_cm1: float = 10.0,
    keep_low_positive: bool = False,
) -> ThermoContribution:
    """Harmonic vibrational thermochemistry, ported from Merlino thermo_vib."""
    if T_K <= 0.0:
        raise ValueError("vibrational thermochemistry requires T_K > 0")

    freq = np.asarray(tuple(float(value) for value in frequencies_cm1), dtype=float)
    nfreq_total = int(freq.size)
    if keep_low_positive:
        freq_used = freq[freq > 0.0]
    else:
        freq_used = freq[freq > float(cutoff_cm1)]
    nfreq_used = int(freq_used.size)

    if nfreq_used == 0:
        reason = "no frequencies above cutoff" if nfreq_total else "no #VIBRATIONAL frequencies"
        return ThermoContribution(
            Q_dimless=1.0,
            U_kJmol=0.0,
            H_kJmol=0.0,
            S_JmolK=0.0,
            Cv_JmolK=0.0,
            Cp_JmolK=0.0,
            available=False,
            reason=reason,
            diagnostics={
                "T_K": float(T_K),
                "nfreq_total": nfreq_total,
                "nfreq_used": nfreq_used,
                "cutoff_cm1": float(cutoff_cm1),
                "keep_low_positive": bool(keep_low_positive),
            },
        )

    constants = get_physical_constants()
    h = constants[Phy.PLANCK]
    c = constants[Phy.C_LIGHT]
    kB = constants[Phy.BOLTZMANN]
    NA = constants[Phy.AVOGADRO]
    R = kB * NA

    nu_hz = freq_used * c
    x = (h * nu_hz) / (kB * T_K)
    x = np.clip(x, 0.0, 700.0)
    emx = np.exp(-x)
    ex = np.exp(x)

    ln1m = np.log1p(-emx)
    lnQ = float(np.sum(-0.5 * x - ln1m))
    Qvib = exp(lnQ)

    U = float(np.sum(h * nu_hz * (0.5 + 1.0 / (ex - 1.0))) * NA)
    H = U
    S = float(R * np.sum((x * emx) / (1.0 - emx) - ln1m))
    Cv = float(R * np.sum((x * x * ex) / ((ex - 1.0) ** 2)))

    return ThermoContribution(
        Q_dimless=Qvib,
        U_kJmol=U / 1000.0,
        H_kJmol=H / 1000.0,
        S_JmolK=S,
        Cv_JmolK=Cv,
        Cp_JmolK=Cv,
        diagnostics={
            "T_K": float(T_K),
            "lnQ": lnQ,
            "nfreq_total": nfreq_total,
            "nfreq_used": nfreq_used,
            "cutoff_cm1": float(cutoff_cm1),
            "keep_low_positive": bool(keep_low_positive),
        },
    )
