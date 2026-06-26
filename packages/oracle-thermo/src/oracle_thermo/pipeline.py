from __future__ import annotations

from pathlib import Path

from oracle_chem import Structure, read_enriched_xyz
from oracle_chem.rotational import rotational_info
from oracle_core import read_basic_section
from oracle_rovib import read_rotational_section, read_vibrational_section

from .contracts import write_thermo_section
from .models import THERMO_KEYS, ThermoContribution, ThermoResult, ThermoSection
from .rotational import rotational_thermo
from .translational import translational_thermo
from .vibrational import vibrational_thermo


def thermo_pipeline(
    xyzin: Path | str,
    *,
    report: bool = True,
    report_path: Path | None = None,
    write_section: bool = True,
    cutoff_cm1: float = 10.0,
    keep_low_positive: bool = False,
) -> ThermoResult:
    target = Path(xyzin)
    if target.is_dir():
        raise IsADirectoryError(f"thermo_pipeline expects an xyzin file, got directory: {target}")

    basic = read_basic_section(target)
    geometry = read_enriched_xyz(target)
    structure = Structure(
        symbols=list(geometry.atoms),
        coords=[tuple(row) for row in geometry.coordinates_angstrom],
        isotopes=None,
    )
    rotational_section = read_rotational_section(target)
    vibrational_section = read_vibrational_section(target)

    translational = translational_thermo(structure, basic)
    rotational = _rotational_contribution(structure, basic, rotational_section)
    vibrational = vibrational_thermo(
        vibrational_section.frequencies_cm1,
        T_K=float(basic.temperature_K),
        cutoff_cm1=cutoff_cm1,
        keep_low_positive=keep_low_positive,
    )
    total = _total_contribution(translational, rotational, vibrational)

    section = ThermoSection(
        translational=translational,
        rotational=rotational,
        vibrational=vibrational,
        total=total,
    )
    result = ThermoResult(
        basic=basic,
        rotational_section=rotational_section,
        vibrational_section=vibrational_section,
        section=section,
    )
    if write_section:
        write_thermo_section(target, section)
    if report:
        output = report_path or target.parent / "thermo.report"
        output.write_text("\n".join(thermo_report_lines(result)).rstrip() + "\n", encoding="utf-8")
    return result


def run_thermo_on_xyzin(
    xyzin: Path | str,
    *,
    report: bool = True,
    report_path: Path | None = None,
    write_section: bool = True,
    cutoff_cm1: float = 10.0,
    keep_low_positive: bool = False,
) -> ThermoResult:
    return thermo_pipeline(
        xyzin,
        report=report,
        report_path=report_path,
        write_section=write_section,
        cutoff_cm1=cutoff_cm1,
        keep_low_positive=keep_low_positive,
    )


def thermo_report_lines(result: ThermoResult) -> list[str]:
    lines = [
        "THERMO PIPELINE REPORT",
        "keys: Q_dimless U_kJmol H_kJmol S_JmolK Cv_JmolK Cp_JmolK",
        f"T_K = {float(result.basic.temperature_K):.8g}",
        f"P_ATM = {float(result.basic.pressure_atm):.8g}",
    ]
    for label, contribution in (
        ("TRASL", result.translational),
        ("ROT", result.rotational),
        ("VIB", result.vibrational),
        ("TOT", result.total),
    ):
        if contribution is None:
            continue
        lines.append("")
        lines.append(label)
        for key in THERMO_KEYS:
            value = getattr(contribution, key)
            if value is not None:
                lines.append(f"{key} = {float(value):.12g}")
        if not contribution.available:
            lines.append(f"available = 0")
            lines.append(f"reason = {contribution.reason}")
    return lines


def _rotational_contribution(
    structure: Structure,
    basic,
    section,
) -> ThermoContribution:
    constants = _resolved_rotational_constants(structure, section)
    rotor_type = section.rotor_type
    if not rotor_type:
        rotor_type = rotational_info(structure, isotopic=True)["rotor_type"]
    return rotational_thermo(
        constants[0],
        constants[1],
        constants[2],
        rotor_type,
        T_K=float(basic.temperature_K),
        sigma=section.symmetry_number,
    )


def _resolved_rotational_constants(structure: Structure, section) -> tuple[float, float, float]:
    if section.A_MHz is None and section.B_MHz is None and section.C_MHz is None:
        info = rotational_info(structure, isotopic=True)
        return float(info["A"]), float(info["B"]), float(info["C"])
    B = section.B_MHz
    if B is None:
        B = section.A_MHz if section.A_MHz is not None else section.C_MHz
    A = section.A_MHz if section.A_MHz is not None else B
    C = section.C_MHz if section.C_MHz is not None else B
    return float(A or 0.0), float(B or 0.0), float(C or 0.0)


def _total_contribution(*contributions: ThermoContribution | None) -> ThermoContribution:
    totals = {}
    for key in THERMO_KEYS:
        if key == "Q_dimless":
            q = 1.0
            present = False
            for contribution in contributions:
                value = None if contribution is None else contribution.Q_dimless
                if value is None:
                    continue
                q *= float(value)
                present = True
            totals[key] = q if present else None
            continue

        value_sum = 0.0
        present = False
        for contribution in contributions:
            value = None if contribution is None else getattr(contribution, key)
            if value is None:
                continue
            value_sum += float(value)
            present = True
        totals[key] = value_sum if present else None

    return ThermoContribution(
        Q_dimless=totals["Q_dimless"],
        U_kJmol=totals["U_kJmol"],
        H_kJmol=totals["H_kJmol"],
        S_JmolK=totals["S_JmolK"],
        Cv_JmolK=totals["Cv_JmolK"],
        Cp_JmolK=totals["Cp_JmolK"],
        diagnostics={"components": ("trasl", "rot", "vib")},
    )
