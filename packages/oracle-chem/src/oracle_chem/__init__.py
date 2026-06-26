"""Atoms, geometry, topology, rings, symmetry and synthons for ORACLE."""

from .geometry import MolecularGeometry
from .geometry_io import (
    GeometryParseError,
    normalize_atom_symbol,
    parse_xyz_lines,
    read_enriched_xyz,
    read_geometry,
    read_xyz,
)
from .zmatrix import (
    ZMatrix,
    ZMatrixAtom,
    parse_zmatrix_text,
    read_zmatrix,
    zmatrix_to_geometry,
)
from .babel import (
    BabelPreprocessResult,
    SymmetryThresholds,
    preprocess_to_enriched_xyz,
)
from .validation import (
    ORACLE_XYZ_VALIDATION_SCHEMA,
    ValidationMessage,
    ValidationResult,
    validate_enriched_molecule,
    validation_section_lines,
    write_validation_section,
)
from .topology.atomic_synthons import AtomicSynthons
from .topology.pipeline import build_topology_objects

__all__ = [
    "AtomicSynthons",
    "GeometryParseError",
    "MolecularGeometry",
    "ORACLE_XYZ_VALIDATION_SCHEMA",
    "BabelPreprocessResult",
    "SymmetryThresholds",
    "ValidationMessage",
    "ValidationResult",
    "ZMatrix",
    "ZMatrixAtom",
    "build_topology_objects",
    "normalize_atom_symbol",
    "parse_xyz_lines",
    "preprocess_to_enriched_xyz",
    "read_enriched_xyz",
    "read_geometry",
    "read_xyz",
    "read_zmatrix",
    "parse_zmatrix_text",
    "validate_enriched_molecule",
    "validation_section_lines",
    "write_validation_section",
    "zmatrix_to_geometry",
]
