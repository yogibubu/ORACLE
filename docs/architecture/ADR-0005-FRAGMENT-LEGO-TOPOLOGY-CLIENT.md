# ADR-0005: Fragment/LEGO Workflows Are Topology Clients

Date: 2026-06-26

## Status

Accepted

## Context

ORACLE needs a Nano-LEGO-like utility that can:

- fragment a molecule into reusable chemistry units;
- search a fragment library such as LCB25-derived fragments;
- assemble candidate structures from compatible fragments.

This tool is scientifically useful only if all ORACLE modules agree on the same
connectivity, rings, aromaticity, synthons and atom identifiers. Recomputing
topology inside the fragment utility would recreate the duplication that ORACLE
is meant to remove.

## Decision

Create `oracle-fragments` as the package boundary for future fragmentation,
fragment search and fragment assembly. It is not a topology engine.

`oracle-fragments` must consume enriched XYZ files that already contain:

- `#TOPOLOGY` with schema `oracle.xyz.topology.v1`;
- `#SYNTHONS` with schema `oracle.xyz.synthons.v1`.

It may then write:

- `#FRAGMENTS` with schema `oracle.xyz.fragments.v1`;
- `#FRAGMENT_LIBRARY` with schema `oracle.xyz.fragment_library.v1`;
- `#ASSEMBLY` with schema `oracle.xyz.assembly.v1`.

The initial implementation only records a planned `#FRAGMENTS` section after
validating the topology/synthon prerequisites. Fragment extraction, library
indexing and assembly algorithms are deferred until the shared topology contract
is robust across ORACLE cases and LCB25.

The initial command is:

```bash
python -m oracle fragments plan molecule.xyzin
```

## Consequences

- LCB25 can be cached and preprocessed now, but fragment indexing waits for the
  robust shared topology pass.
- GUI, MORPHEUS, GICForge and future data workflows will consume the same
  fragment sections instead of building their own molecular graph logic.
- Any future Python or Fortran implementation of fragment algorithms must be a
  backend behind this package contract and use the same input sections.
