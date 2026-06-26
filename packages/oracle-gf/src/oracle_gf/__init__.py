"""Frozen-GIC GF/PED analysis workflows."""

from importlib import import_module

from .harmonic import GFResult, mass_weighted_cartesian_hessian, solve_wilson_gf
from .internal import (
    BOHR_TO_ANGSTROM,
    GFLocalOptions,
    InternalGFResult,
    PEDTable,
    gf_from_cartesian_hessian_and_gic_b_matrix,
    gf_from_cartesian_hessian_and_oracle_gics,
    gf_from_gaussian_fchk_with_oracle_gics,
    gf_from_hessian_input_and_gic_definition,
    gf_from_hessian_input_and_xyzin,
    gf_from_hessian_input_with_oracle_gics,
    gic_labels_from_u,
    local_force_constant_mask,
    primitive_label,
    pulay_scale_internal_hessian,
    topology_bonds_from_xyzin,
)
from .models import HessianInput
from .nonbonded import (
    nonbonded_cartesian_hessian_correction,
    synthon_charges_from_xyzin,
)

_SERVICE_EXPORTS = {
    "GFReport",
    "GFLocalOptions",
    "format_gf_report",
    "gf_csv_tables",
    "pulay_scaling_factors",
    "run_gf_report_from_fchk",
    "run_xyzin_gf_report_from_fchk",
    "write_csv_tables",
}


def __getattr__(name: str):
    if name in _SERVICE_EXPORTS:
        return getattr(import_module(".service", __name__), name)
    raise AttributeError(name)


__all__ = [
    "BOHR_TO_ANGSTROM",
    "GFReport",
    "GFResult",
    "HessianInput",
    "InternalGFResult",
    "PEDTable",
    "format_gf_report",
    "gf_csv_tables",
    "gf_from_cartesian_hessian_and_gic_b_matrix",
    "gf_from_cartesian_hessian_and_oracle_gics",
    "gf_from_gaussian_fchk_with_oracle_gics",
    "gf_from_hessian_input_and_gic_definition",
    "gf_from_hessian_input_and_xyzin",
    "gf_from_hessian_input_with_oracle_gics",
    "gic_labels_from_u",
    "local_force_constant_mask",
    "mass_weighted_cartesian_hessian",
    "nonbonded_cartesian_hessian_correction",
    "primitive_label",
    "pulay_scale_internal_hessian",
    "pulay_scaling_factors",
    "run_gf_report_from_fchk",
    "run_xyzin_gf_report_from_fchk",
    "solve_wilson_gf",
    "synthon_charges_from_xyzin",
    "topology_bonds_from_xyzin",
    "write_csv_tables",
]
