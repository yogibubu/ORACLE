"""VPT2/VCI contracts for ORACLE."""

from oracle_core import ScientificValidationError

from .contracts import (
    DavidsonSettings,
    ForceFieldSource,
    VCIRequest,
    VPT2VCIInventory,
    inventory_vpt2_vci_backends,
)
from .gaussian_qff import (
    FCHKData,
    anharmonic_input_from_gaussian_fchk,
    lower_to_symmetric,
    read_gaussian_fchk_qff,
    read_indexed_qff_text,
)
from .davidson import DavidsonResult, davidson_lowest
from .gui_service import (
    VPT2VCIReport,
    format_vpt2_vci_report,
    load_force_field,
    run_vpt2_vci_report,
    vpt2_vci_csv_tables,
    write_csv_tables,
)
from .models import AnharmonicInput
from .validation import validate_force_field
from .vci import (
    QuarticForceField,
    VCIBlockInfo,
    VCIOptions,
    VCIResult,
    VCIStateContribution,
    build_vci_hamiltonian,
    force_field_from_anharmonic_input,
    generate_vibrational_basis,
    solve_vci,
    solve_vci_from_anharmonic_input,
    zero_anharmonic_force_field,
)
from .vpt2 import (
    VPT2Result,
    VPT2State,
    VPT2VCIComparison,
    compare_vpt2_vci,
    compare_vpt2_vci_from_anharmonic_input,
    solve_vpt2,
    solve_vpt2_from_anharmonic_input,
)
from .workflow import VPT2VCIRun, run_python_vci_from_gaussian_fchk

__all__ = [
    "DavidsonSettings",
    "DavidsonResult",
    "AnharmonicInput",
    "FCHKData",
    "ForceFieldSource",
    "QuarticForceField",
    "ScientificValidationError",
    "VCIBlockInfo",
    "VCIOptions",
    "VCIRequest",
    "VCIResult",
    "VCIStateContribution",
    "VPT2VCIInventory",
    "VPT2VCIRun",
    "VPT2Result",
    "VPT2State",
    "VPT2VCIComparison",
    "VPT2VCIReport",
    "anharmonic_input_from_gaussian_fchk",
    "build_vci_hamiltonian",
    "compare_vpt2_vci",
    "compare_vpt2_vci_from_anharmonic_input",
    "davidson_lowest",
    "force_field_from_anharmonic_input",
    "format_vpt2_vci_report",
    "generate_vibrational_basis",
    "inventory_vpt2_vci_backends",
    "load_force_field",
    "lower_to_symmetric",
    "read_gaussian_fchk_qff",
    "read_indexed_qff_text",
    "run_python_vci_from_gaussian_fchk",
    "run_vpt2_vci_report",
    "solve_vci",
    "solve_vci_from_anharmonic_input",
    "solve_vpt2",
    "solve_vpt2_from_anharmonic_input",
    "validate_force_field",
    "vpt2_vci_csv_tables",
    "write_csv_tables",
    "zero_anharmonic_force_field",
]
