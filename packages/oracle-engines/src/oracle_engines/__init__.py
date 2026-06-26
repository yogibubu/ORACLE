"""Fortran and vendored backend discovery wrappers."""

from .gicforge import (
    LEGACY_GICFORGE_FILES,
    FortranBackendLayout,
    gicforge_fortran_layout,
    legacy_gicforge_source_paths,
    validate_legacy_gicforge_sources,
)

__all__ = [
    "FortranBackendLayout",
    "LEGACY_GICFORGE_FILES",
    "gicforge_fortran_layout",
    "legacy_gicforge_source_paths",
    "validate_legacy_gicforge_sources",
]
