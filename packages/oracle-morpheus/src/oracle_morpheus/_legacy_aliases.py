from __future__ import annotations

import importlib
import sys


_ALIASES = (
    "geometry",
    "merlino_core",
    "merlino_fortran",
    "merlino_fit",
    "topology",
    "merlino_gic",
)


def install_legacy_aliases() -> None:
    """Resolve Merlino import names to the vendored ORACLE MORPHEUS copy."""
    for name in _ALIASES:
        module = importlib.import_module(f"oracle_morpheus.legacy.{name}")
        sys.modules.setdefault(name, module)
