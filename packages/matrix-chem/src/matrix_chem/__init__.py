"""Atoms, geometry, topology, rings, symmetry and synthons for MATRIX."""

from .average_atomic_masses import atomic_mass
from .geometry import MolecularGeometry
from .geometry_io import (
    GeometryParseError,
    detect_qm_output_format,
    normalize_atom_symbol,
    parse_xyz_lines,
    read_enriched_xyz,
    read_geometry,
    read_geometry_with_kind,
    read_xyz,
    read_xyz_atoms_coords,
    write_xyz,
)
from .structure_files import read_mol2, read_molfile
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
from .link import (
    BabelPreprocessResult,
    LinkPreprocessResult,
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
    MATRIX_XYZ_VALIDATION_SCHEMA,
    ORACLE_XYZ_VALIDATION_SCHEMA,
    ValidationMessage,
    ValidationResult,
    validate_enriched_molecule,
    validation_section_lines,
    write_validation_section,
)
from .topology.atomic_synthons import AtomicSynthons
from .topology.contracts import (
    MATRIX_XYZ_FRAGMENTS_SCHEMA,
    MATRIX_XYZ_SYNTHONS_SCHEMA,
    MATRIX_XYZ_TOPOLOGY_SCHEMA,
    ORACLE_XYZ_FRAGMENTS_SCHEMA,
    ORACLE_XYZ_SYNTHONS_SCHEMA,
    ORACLE_XYZ_TOPOLOGY_SCHEMA,
)
from .topology.pipeline import build_topology_objects
from .topology.snapshot import (
    TOPOLOGY_SNAPSHOT_SCHEMA,
    TopologySnapshotComparison,
    compare_topology_snapshot_entry,
    topology_report_lines,
    topology_snapshot_document,
    topology_snapshot_from_xyzin,
    write_topology_report,
    write_topology_snapshot,
)

__all__ = [
    "AtomicSynthons",
    "GeometryParseError",
    "Isotope",
    "MolecularGeometry",
    "MolecularSymmetry",
    "MATRIX_XYZ_FRAGMENTS_SCHEMA",
    "MATRIX_XYZ_SYNTHONS_SCHEMA",
    "MATRIX_XYZ_TOPOLOGY_SCHEMA",
    "MATRIX_XYZ_VALIDATION_SCHEMA",
    "ORACLE_XYZ_FRAGMENTS_SCHEMA",
    "ORACLE_XYZ_SYNTHONS_SCHEMA",
    "ORACLE_XYZ_TOPOLOGY_SCHEMA",
    "ORACLE_XYZ_VALIDATION_SCHEMA",
    "Phy",
    "BabelPreprocessResult",
    "LinkPreprocessResult",
    "Structure",
    "SymmetryThresholds",
    "SymmetryOperation",
    "TOPOLOGY_SNAPSHOT_SCHEMA",
    "TopologySnapshotComparison",
    "ValidationMessage",
    "ValidationResult",
    "ZMatrix",
    "ZMatrixAtom",
    "atomic_mass",
    "build_topology_objects",
    "center_of_mass",
    "compare_topology_snapshot_entry",
    "get_default_isotope",
    "get_isotope",
    "get_isotopes",
    "get_physical_constants",
    "analyze_molecular_symmetry",
    "inertia_tensor",
    "detect_qm_output_format",
    "normalize_atom_symbol",
    "parse_xyz_lines",
    "principal_moments",
    "preprocess_to_enriched_xyz",
    "read_enriched_xyz",
    "read_geometry",
    "read_geometry_with_kind",
    "read_mol2",
    "read_molfile",
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
    "topology_report_lines",
    "topology_snapshot_document",
    "topology_snapshot_from_xyzin",
    "write_topology_report",
    "write_topology_snapshot",
    "zmatrix_to_geometry",
]
