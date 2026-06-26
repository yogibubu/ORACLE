"""Atoms, geometry, topology, rings, symmetry and synthons for ORACLE."""

from .average_atomic_masses import atomic_mass
from .geometry import MolecularGeometry
from .geometry_io import (
    GeometryParseError,
    normalize_atom_symbol,
    parse_xyz_lines,
    read_enriched_xyz,
    read_geometry,
    read_xyz,
    read_xyz_atoms_coords,
    write_xyz,
)
from .inertia import center_of_mass, inertia_tensor, principal_moments
from .isotopes_table import Isotope, get_default_isotope, get_isotope, get_isotopes
from .physical_constants import Phy, get_physical_constants
from .rotational import rotational_constants_MHz
from .structure import Structure
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
from .symmetry import (
    MolecularSymmetry,
    SymmetryOperation,
    analyze_molecular_symmetry,
    symmetry_section_lines,
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
    "Isotope",
    "MolecularGeometry",
    "MolecularSymmetry",
    "ORACLE_XYZ_VALIDATION_SCHEMA",
    "Phy",
    "BabelPreprocessResult",
    "Structure",
    "SymmetryThresholds",
    "SymmetryOperation",
    "ValidationMessage",
    "ValidationResult",
    "ZMatrix",
    "ZMatrixAtom",
    "atomic_mass",
    "build_topology_objects",
    "center_of_mass",
    "get_default_isotope",
    "get_isotope",
    "get_isotopes",
    "get_physical_constants",
    "analyze_molecular_symmetry",
    "inertia_tensor",
    "normalize_atom_symbol",
    "parse_xyz_lines",
    "principal_moments",
    "preprocess_to_enriched_xyz",
    "read_enriched_xyz",
    "read_geometry",
    "read_xyz",
    "read_xyz_atoms_coords",
    "read_zmatrix",
    "rotational_constants_MHz",
    "parse_zmatrix_text",
    "validate_enriched_molecule",
    "validation_section_lines",
    "write_validation_section",
    "write_xyz",
    "symmetry_section_lines",
    "zmatrix_to_geometry",
]
