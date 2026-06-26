# ORACLE

**Operational Recognition of Atomistic Connectivity and Local Environments**

ORACLE is the suite-level refactor of the current Merlino/MORPHEUS codebase.
The first releases keep Merlino compatibility aliases while the scientific
packages are separated behind stable service, CLI and manifest contracts.

## Initial Scope

- `oracle-core`: workspace layout, configuration, manifests, checksums and
  common errors.
- `oracle-chem`: atoms, masses, geometry, topology, rings and symmetry.
  Descriptor topology and atomic synthons live here as shared chemistry
  primitives. It owns the canonical `MolecularGeometry` model and plain/enriched
  XYZ parsers.
- `oracle-gicforge`: non-redundant GIC construction, frozen schemas, SYCART and
  B-matrix evaluation.
- `oracle-morpheus`: semiexperimental geometry refinement, constraints,
  predicates, parameter classes, diagnostics and ensemble refinement.
- `oracle-gf`: frozen-GIC GF/PED analysis.
- `oracle-gaussian`: Gaussian input/output adapters.
- `oracle-engines`: Fortran and vendored backend discovery/build wrappers.
- `oracle-dvr`: scan/grid to DVR workflows.
- `oracle-vpt2-vci`: QFF, VPT2/VCI and Davidson workflows.
- `oracle-gui`: GUI controllers and views only.

## Migration Rule

Scientific behavior is migrated package by package. Existing `merlino_*`
imports stay valid until ORACLE-native tests cover the new public APIs.

## Core Architecture Rule

ORACLE modules must not reinvent shared operations. Common tasks such as
sectioned XYZ I/O, atom and isotope data, topology, symmetry, GIC construction,
Gaussian parsing, backend execution and manifests belong to shared libraries.

Geometry and QM-program parsers are unified: they return shared ORACLE models
and live in parser packages, not GUI or workflow modules.

The preprocessing layer is ORACLE-Babel: it imports external sources, writes an
Avogadro-compatible enriched XYZ, determines symmetry with explicit thresholds,
builds topology once, and saves descriptors/synthons for downstream tools.

Python and strict Fortran77 implementations may intentionally coexist for the
same scientific kernel. In that case they are backends behind one ORACLE service
contract, with shared schemas, shared enriched XYZ sections and shared
regression or identity tests.

The enriched XYZ file is the canonical communication object between modules:
tools append or replace their own uppercase sections and preserve all unrelated
sections. See
`docs/architecture/ADR-0001-SHARED-LIBRARIES-AND-XYZ-CONTAINER.md`.

## Workspace Contract

Project workspaces use:

```text
inputs/
runs/
outputs/
reports/
cache/
logs/
```

Every workflow run writes an `oracle.run.v1` manifest containing input/output
paths, SHA256 hashes, parameters, backend metadata and status.
