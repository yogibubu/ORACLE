# ADR-0009: MATRIX Framework Naming Transition

Date: 2026-06-27

## Status

Accepted as a future refactory constraint.

## Context

The current repository and runtime packages still use ORACLE-oriented names
because the active refactory is stabilizing scientific contracts, test coverage
and compatibility aliases. Once the refactory is complete, the naming should
separate the framework, user-facing GUI and scientific tools more clearly.

## Decision

After the refactory is stable:

- **MATRIX** becomes the framework/package family name.
- **ORACLE** remains the user-facing GUI/application name.
- **GICForge** is renamed to **NEO**.
- Other tools should use Matrix-saga character names when a rename is useful
  and does not obscure the scientific contract.

No runtime package is renamed immediately. Current imports, CLIs, docs and
tests keep their existing names until the compatibility surface is stable.

`oracle_core.tool_contracts` records both current names and planned names. This
is the source of truth for migration planning until actual package renames are
performed.

## Consequences

- Public contracts remain stable during scientific porting.
- Renames happen as a planned compatibility migration, not as incidental churn.
- Future docs can introduce tools as, for example, "NEO, formerly GICForge",
  only after compatibility aliases are in place.
- The GUI may be branded ORACLE while running on the MATRIX framework.

