# ADR-0009: MATRIX Framework Naming Transition

Date: 2026-06-27

## Status

Accepted. Compatibility aliases are active; physical package renames remain
deferred until the scientific contracts are fully stable.

## Context

The current repository and runtime packages still use ORACLE-oriented names
because the active refactory is stabilizing scientific contracts, test coverage
and compatibility aliases. Once the refactory is complete, the naming should
separate the framework, user-facing GUI and scientific tools more clearly.

## Decision

After the refactory is stable:

- **MATRIX** becomes the framework/package family name:
  **Molecular Analysis Toolkit for Reusable Integrated eXperiments**.
- **ORACLE** remains the user-facing GUI/application name:
  **Operator for Routing, Analysis, Control, Launch and Exploration**.
- **GICForge** is renamed to **NEO**:
  **Nonredundant Equivariant Orthogonalizer**.
- Other tools should use Matrix-saga character names when a rename is useful
  and does not obscure the scientific contract.

NEO is assigned to GICForge rather than to the GUI because that tool is the
coordinate engine that builds, projects, symmetrizes and reduces the molecular
internal-coordinate representation. ORACLE is a better GUI name because the GUI
routes users through project state, launches tools and exposes the scientific
state without owning parser or kernel logic.

No runtime package is renamed immediately. Current imports, CLIs, docs and
tests keep their existing names until the compatibility surface is stable.

The first compatibility layer is active:

- `matrix` is an installable console alias for the framework CLI.
- `oracle neo ...` is an alias for `oracle gicforge ...`.
- `neo ...` is an installable console alias for the GICForge/NEO coordinate
  tool.

`oracle_core.tool_contracts` records both current names and planned names. This
is the source of truth for migration planning until actual package renames are
performed.

## Consequences

- Public contracts remain stable during scientific porting.
- Renames happen as a planned compatibility migration, not as incidental churn.
- Future docs can introduce tools as, for example, "NEO, formerly GICForge",
  only after compatibility aliases are in place.
- The GUI may be branded ORACLE while running on the MATRIX framework.
- Internal Python packages keep their `oracle_*` import names for now, so the
  migration does not break downstream scripts or tests.
