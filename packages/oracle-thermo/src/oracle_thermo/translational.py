from __future__ import annotations

from math import exp, log, pi

from oracle_chem import Phy, Structure, get_physical_constants
from oracle_core import BasicSection

from .models import ThermoContribution


def translational_thermo(structure: Structure, basic: BasicSection) -> ThermoContribution:
    """Ideal-gas translational thermochemistry, ported from Merlino thermo_trasl."""
    T = float(basic.temperature_K)
    P_atm = float(basic.pressure_atm)
    if T <= 0.0:
        raise ValueError("translational thermochemistry requires T_K > 0")
    if P_atm <= 0.0:
        raise ValueError("translational thermochemistry requires P_ATM > 0")

    constants = get_physical_constants()
    kB = constants[Phy.BOLTZMANN]
    h = constants[Phy.PLANCK]
    NA = constants[Phy.AVOGADRO]
    R = kB * NA
    mass_kg = float(structure.total_mass_isotope) * constants[Phy.TO_KG]
    pressure_pa = P_atm * 101325.0

    lnQ = (
        1.5 * log((2.0 * pi * mass_kg * kB * T) / (h * h))
        + log(kB * T / pressure_pa)
    )
    U = 1.5 * R * T
    H = 2.5 * R * T
    Cv = 1.5 * R
    Cp = 2.5 * R
    S = R * (lnQ + 2.5)

    return ThermoContribution(
        Q_dimless=exp(lnQ),
        U_kJmol=U / 1000.0,
        H_kJmol=H / 1000.0,
        S_JmolK=S,
        Cv_JmolK=Cv,
        Cp_JmolK=Cp,
        diagnostics={
            "T_K": T,
            "P_ATM": P_atm,
            "mass_amu": float(structure.total_mass_isotope),
            "lnQ": lnQ,
        },
    )
