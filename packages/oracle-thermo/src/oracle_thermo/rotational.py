from __future__ import annotations

from math import erf, exp, log, pi, sqrt

from oracle_chem import Phy, get_physical_constants

from .models import ThermoContribution


def qrot_quantum_linear(Beff_Hz: float, T_K: float, Jmax: int = 200) -> float:
    beta = _h() / (_kB() * T_K)
    return sum(
        (2 * J + 1) * exp(-beta * Beff_Hz * J * (J + 1)) for J in range(Jmax + 1)
    )


def qrot_quantum_spherical(B_Hz: float, T_K: float, Jmax: int = 200) -> float:
    beta = _h() / (_kB() * T_K)
    return sum(
        (2 * J + 1) ** 2 * exp(-beta * B_Hz * J * (J + 1)) for J in range(Jmax + 1)
    )


def qrot_quantum_symmetric(
    A_Hz: float,
    B_Hz: float,
    T_K: float,
    Jmax: int = 150,
) -> float:
    beta = _h() / (_kB() * T_K)
    q = 0.0
    for J in range(Jmax + 1):
        prefactor = 2 * J + 1
        JJ1 = J * (J + 1)
        for K in range(-J, J + 1):
            q += prefactor * exp(-beta * (B_Hz * JJ1 + (A_Hz - B_Hz) * K * K))
    return q


def rotational_thermo(
    A_MHz: float | None,
    B_MHz: float | None,
    C_MHz: float | None,
    rotor_type: str,
    *,
    T_K: float,
    sigma: int | None = None,
    Tq: float = 200.0,
    Tc: float = 250.0,
    Jmax_linear: int = 200,
    Jmax_spherical: int = 200,
    Jmax_symmetric: int = 150,
    dT_num: float = 0.05,
) -> ThermoContribution:
    """Rotational thermochemistry, preserving Merlino's quantum/classical crossover."""
    if T_K <= 0.0:
        raise ValueError("rotational thermochemistry requires T_K > 0")
    constants = tuple(_finite_constant(value) for value in (A_MHz, B_MHz, C_MHz))
    A_MHz, B_MHz, C_MHz = _complete_rotational_constants(*constants)
    sigma_i = max(int(sigma or 1), 1)
    rotor_label = _canonical_rotor_type(rotor_type)

    if not _rotational_constants_are_usable(A_MHz, B_MHz, C_MHz, rotor_label):
        return ThermoContribution(
            Q_dimless=1.0,
            U_kJmol=0.0,
            H_kJmol=0.0,
            S_JmolK=0.0,
            Cv_JmolK=0.0,
            Cp_JmolK=0.0,
            available=False,
            reason="no usable rotational constants",
            diagnostics={
                "T_K": float(T_K),
                "sigma": sigma_i,
                "A_MHz": A_MHz,
                "B_MHz": B_MHz,
                "C_MHz": C_MHz,
                "rotor_type": rotor_label,
            },
        )

    A_Hz = A_MHz * 1.0e6
    B_Hz = B_MHz * 1.0e6
    C_Hz = C_MHz * 1.0e6
    Beff_Hz = max(B_Hz, C_Hz)
    kind = _rotor_kind(rotor_label)

    def lnQ_quantum(Tloc: float) -> float:
        if kind == "linear":
            Q = qrot_quantum_linear(Beff_Hz, Tloc, Jmax_linear)
        elif kind == "spherical":
            Q = qrot_quantum_spherical(B_Hz, Tloc, Jmax_spherical)
        else:
            Q = qrot_quantum_symmetric(A_Hz, Beff_Hz, Tloc, Jmax_symmetric)
        return log(Q / sigma_i)

    lnQq = lnQ_quantum(T_K)
    lnQc = _lnq_classical(A_Hz, B_Hz, Beff_Hz, T_K, sigma_i, kind)

    if T_K <= Tq:
        lnQ = lnQq
        dlnQ, d2lnQ = _dlnQ_num(lnQ_quantum, T_K, dT=dT_num)
        regime = "quantum"
    elif T_K >= Tc:
        lnQ = lnQc
        dlnQ, d2lnQ = _dlnq_classical(T_K, kind)
        regime = "classical"
    else:
        w, wp, wpp = _mix_weights(T_K, Tq, Tc)
        dlnQq, d2lnQq = _dlnQ_num(lnQ_quantum, T_K, dT=dT_num)
        dlnQc, d2lnQc = _dlnq_classical(T_K, kind)
        lnQ = (1.0 - w) * lnQq + w * lnQc
        dlnQ = (1.0 - w) * dlnQq + w * dlnQc + wp * (lnQc - lnQq)
        d2lnQ = (
            (1.0 - w) * d2lnQq
            + w * d2lnQc
            + 2.0 * wp * (dlnQc - dlnQq)
            + wpp * (lnQc - lnQq)
        )
        regime = "mixed"

    R = _R()
    U = R * T_K * T_K * dlnQ
    S = R * (lnQ + T_K * dlnQ)
    Cv = R * (2.0 * T_K * dlnQ + T_K * T_K * d2lnQ)
    Qrot = exp(lnQ)

    return ThermoContribution(
        Q_dimless=Qrot,
        U_kJmol=U / 1000.0,
        H_kJmol=U / 1000.0,
        S_JmolK=S,
        Cv_JmolK=Cv,
        Cp_JmolK=Cv,
        diagnostics={
            "Qrot": Qrot,
            "lnQ": lnQ,
            "T_K": float(T_K),
            "sigma": sigma_i,
            "rotor_type": rotor_label,
            "rotor_kind": kind,
            "regime": regime,
            "A_MHz": A_MHz,
            "B_MHz": B_MHz,
            "C_MHz": C_MHz,
        },
    )


def _finite_constant(value: float | None) -> float | None:
    if value is None:
        return None
    value = float(value)
    return value if value > 0.0 else None


def _complete_rotational_constants(
    A_MHz: float | None,
    B_MHz: float | None,
    C_MHz: float | None,
) -> tuple[float, float, float]:
    if A_MHz is None and B_MHz is None and C_MHz is None:
        return 0.0, 0.0, 0.0
    if B_MHz is None:
        B_MHz = A_MHz if A_MHz is not None else C_MHz
    if A_MHz is None:
        A_MHz = B_MHz
    if C_MHz is None:
        C_MHz = B_MHz
    return float(A_MHz or 0.0), float(B_MHz or 0.0), float(C_MHz or 0.0)


def _canonical_rotor_type(rotor_type: str) -> str:
    label = (rotor_type or "").strip().lower()
    mapping = {
        "linear_top": "linear",
        "linear rotor": "linear",
        "spherical_top": "spherical",
        "spherical rotor": "spherical",
        "symmetric_top_prolate": "symmetric_prolate",
        "symmetric_top_oblate": "symmetric_oblate",
        "asymmetric_top_quasi_prolate": "asymmetric_prolate",
        "asymmetric_top_quasi_oblate": "asymmetric_oblate",
    }
    return mapping.get(label, label or "unknown")


def _rotor_kind(rotor_type: str) -> str:
    if rotor_type == "linear":
        return "linear"
    if rotor_type == "spherical":
        return "spherical"
    return "symmetric"


def _rotational_constants_are_usable(
    A_MHz: float,
    B_MHz: float,
    C_MHz: float,
    rotor_type: str,
) -> bool:
    if rotor_type == "linear":
        return max(B_MHz, C_MHz) > 0.0
    if rotor_type == "spherical":
        return B_MHz > 0.0
    return A_MHz > 0.0 and max(B_MHz, C_MHz) > 0.0


def _lnq_classical(
    A_Hz: float,
    B_Hz: float,
    Beff_Hz: float,
    T_K: float,
    sigma: int,
    kind: str,
) -> float:
    h = _h()
    kB = _kB()
    if kind == "linear":
        theta = h * Beff_Hz / kB
        return log(T_K / (sigma * theta))
    if kind == "spherical":
        theta = h * B_Hz / kB
        return 0.5 * log(pi) + 1.5 * log(T_K / theta) - log(float(sigma))
    thetaA = h * A_Hz / kB
    thetaB = h * Beff_Hz / kB
    return (
        0.5 * log(pi)
        + 1.5 * log(T_K)
        - log(float(sigma))
        - 0.5 * log(thetaA * thetaB * thetaB)
    )


def _dlnq_classical(T_K: float, kind: str) -> tuple[float, float]:
    if kind == "linear":
        return 1.0 / T_K, -1.0 / (T_K * T_K)
    return 1.5 / T_K, -1.5 / (T_K * T_K)


def _dlnQ_num(lnQ_func, T_K: float, *, dT: float = 0.05) -> tuple[float, float]:
    Tp = T_K + dT
    Tm = max(T_K - dT, 1.0e-6)
    lnQp = lnQ_func(Tp)
    lnQ0 = lnQ_func(T_K)
    lnQm = lnQ_func(Tm)
    d1 = (lnQp - lnQm) / (Tp - Tm)
    d2 = 2.0 * (
        (lnQp - lnQ0) / (Tp - T_K) - (lnQ0 - lnQm) / (T_K - Tm)
    ) / (Tp - Tm)
    return d1, d2


def _mix_weights(T_K: float, Tq: float, Tc: float) -> tuple[float, float, float]:
    T0 = 0.5 * (Tq + Tc)
    dTmix = 0.2 * (Tc - Tq)
    x = (T_K - T0) / dTmix
    ex = exp(-x * x)
    w = 0.5 * (1.0 + erf(x))
    wp = ex / (sqrt(pi) * dTmix)
    wpp = (-2.0 * x) * ex / (sqrt(pi) * dTmix**2)
    return w, wp, wpp


def _constants():
    return get_physical_constants()


def _h() -> float:
    return _constants()[Phy.PLANCK]


def _kB() -> float:
    return _constants()[Phy.BOLTZMANN]


def _R() -> float:
    constants = _constants()
    return constants[Phy.BOLTZMANN] * constants[Phy.AVOGADRO]
