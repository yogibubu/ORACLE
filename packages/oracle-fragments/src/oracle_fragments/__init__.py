"""Topology-backed fragment contracts for ORACLE."""

from .contracts import (
    ORACLE_XYZ_ASSEMBLY_SCHEMA,
    ORACLE_XYZ_FRAGMENT_LIBRARY_SCHEMA,
    ORACLE_XYZ_FRAGMENTS_SCHEMA,
    FragmentDefinition,
    FragmentContractError,
    FragmentRecord,
    build_fragment_definition_from_xyzin,
    fragment_build_section_lines,
    fragment_plan_section_lines,
    read_fragment_records,
    validate_fragment_prerequisites,
    write_fragment_build_section,
    write_fragment_plan_section,
)

__all__ = [
    "FragmentDefinition",
    "FragmentContractError",
    "FragmentRecord",
    "ORACLE_XYZ_ASSEMBLY_SCHEMA",
    "ORACLE_XYZ_FRAGMENT_LIBRARY_SCHEMA",
    "ORACLE_XYZ_FRAGMENTS_SCHEMA",
    "build_fragment_definition_from_xyzin",
    "fragment_build_section_lines",
    "fragment_plan_section_lines",
    "read_fragment_records",
    "validate_fragment_prerequisites",
    "write_fragment_build_section",
    "write_fragment_plan_section",
]
