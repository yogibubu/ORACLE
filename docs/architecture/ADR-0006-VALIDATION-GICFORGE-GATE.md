# ADR-0006: Validation Gate Before GICForge

Date: 2026-06-26

## Status

Accepted

## Context

After ORACLE-Babel has imported a molecule, determined symmetry, built
topology, saved synthons and optionally planned fragments, the molecule should
be validated before downstream scientific tools examine it.

GICForge is the first major post-validation tool. It constructs non-redundant
GICs, optionally symmetrizes them and writes frozen coordinate sections for
MORPHEUS, GF/PED, Gaussian writers and later workflows.

## Decision

`oracle-chem` owns `#VALIDATION` with schema `oracle.xyz.validation.v1`.

Validation consumes the existing enriched XYZ state:

- ordinary XYZ block;
- `#SYMMETRY`;
- `#TOPOLOGY`;
- `#SYNTHONS`;
- optional `#FRAGMENTS`.

It writes `STATUS PASS`, `WARN` or `FAIL` plus diagnostic messages. It does not
recompute topology or synthons.

`oracle-gicforge` must require `#VALIDATION` with `STATUS PASS` before writing
`#GIC` or `#SYCART`. Planned GICForge sections are allowed before the full
GICForge algorithm is migrated, but they must record dependencies and
symmetrization intent.

## Initial Commands

```bash
python -m oracle validate molecule.xyzin
python -m oracle gicforge plan molecule.xyzin --symmetrize --sycart
```

## Consequences

- Downstream tools do not inspect raw geometry or topology directly before the
  molecule passes validation.
- If Avogadro edits the XYZ block, dependent sections must be regenerated:
  `#TOPOLOGY`, `#SYNTHONS`, `#FRAGMENTS`, `#VALIDATION`, `#GIC` and `#SYCART`.
- Future robust validation can add geometry hashes and more chemical
  diagnostics without changing the gate shape.
