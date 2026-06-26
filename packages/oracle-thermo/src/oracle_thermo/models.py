from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

from oracle_core import BasicSection
from oracle_rovib import RotationalSection, VibrationalSection


THERMO_KEYS = ("Q_dimless", "U_kJmol", "H_kJmol", "S_JmolK", "Cv_JmolK", "Cp_JmolK")
THERMO_LABELS = ("trasl", "rot", "vib", "tot")


@dataclass(frozen=True)
class ThermoContribution:
    Q_dimless: float | None = None
    U_kJmol: float | None = None
    H_kJmol: float | None = None
    S_JmolK: float | None = None
    Cv_JmolK: float | None = None
    Cp_JmolK: float | None = None
    available: bool = True
    reason: str = "ok"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    def as_dict(self, *, include_diagnostics: bool = False) -> dict[str, object]:
        values: dict[str, object] = {
            "Q_dimless": self.Q_dimless,
            "U_kJmol": self.U_kJmol,
            "H_kJmol": self.H_kJmol,
            "S_JmolK": self.S_JmolK,
            "Cv_JmolK": self.Cv_JmolK,
            "Cp_JmolK": self.Cp_JmolK,
            "available": self.available,
            "reason": self.reason,
        }
        if include_diagnostics:
            values.update(self.diagnostics)
        return values


@dataclass(frozen=True)
class ThermoSection:
    translational: ThermoContribution | None = None
    rotational: ThermoContribution | None = None
    vibrational: ThermoContribution | None = None
    total: ThermoContribution | None = None
    schema: str = "oracle.xyz.thermo.v1"

    def contribution(self, label: str) -> ThermoContribution | None:
        normalized = label.strip().lower()
        if normalized in {"trasl", "trans", "translational"}:
            return self.translational
        if normalized in {"rot", "rotational"}:
            return self.rotational
        if normalized in {"vib", "vibrational"}:
            return self.vibrational
        if normalized in {"tot", "total"}:
            return self.total
        raise KeyError(label)


@dataclass(frozen=True)
class ThermoResult:
    basic: BasicSection
    rotational_section: RotationalSection
    vibrational_section: VibrationalSection
    section: ThermoSection

    @property
    def translational(self) -> ThermoContribution | None:
        return self.section.translational

    @property
    def rotational(self) -> ThermoContribution | None:
        return self.section.rotational

    @property
    def vibrational(self) -> ThermoContribution | None:
        return self.section.vibrational

    @property
    def total(self) -> ThermoContribution | None:
        return self.section.total
