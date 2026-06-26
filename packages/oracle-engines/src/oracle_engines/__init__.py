"""Fortran and vendored backend discovery wrappers."""

from .gicforge import (
    LEGACY_GICFORGE_FILES,
    FortranBackendLayout,
    gicforge_fortran_layout,
    legacy_gicforge_source_paths,
    validate_legacy_gicforge_sources,
)
from .legacy_gicforge import (
    LegacyGICForgeRun,
    legacy_gicforge_executable,
    read_legacy_gicforge_run,
    run_legacy_gicforge,
)

__all__ = [
    "FortranBackendLayout",
    "LEGACY_GICFORGE_FILES",
    "LegacyGICForgeRun",
    "gicforge_fortran_layout",
    "legacy_gicforge_executable",
    "legacy_gicforge_source_paths",
    "read_legacy_gicforge_run",
    "run_legacy_gicforge",
    "validate_legacy_gicforge_sources",
]
