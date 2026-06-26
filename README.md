# ORACLE

**Operational Recognition of Atomistic Connectivity and Local Environments**

ORACLE is the suite-level refactor of the legacy MORPHEUS/structural-analysis
codebase. The first releases keep compatibility aliases while the scientific
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
- `oracle-babel`: import adapters for external sources and molecular databases
  such as LCB25.
- `oracle-fragments`: topology-backed fragmentation, fragment-library search
  and future assembly contracts.
- `oracle-rovib`: rotational, vibrational and rovibrational compatibility
  sections.
- `oracle-thermo`: thermochemistry from enriched XYZ state.
- `oracle-engines`: Fortran and vendored backend discovery/build wrappers.
- `oracle-dvr`: scan/grid to DVR workflows.
- `oracle-vpt2-vci`: QFF, VPT2/VCI and Davidson workflows.
- `oracle-gui`: GUI controllers and views only.

LCB25 geometries are managed as a reproducible local cache rather than
committed source files:

```bash
python -m oracle lcb25 fetch
```

ORACLE shell helpers are available without editing personal startup files
automatically. `oracle-set` creates `ORACLE_VENV` when missing and installs the
core runtime dependencies, including RDKit for SMILES imports:

```bash
source /Users/vincenzobarone/ORACLE/scripts/oracle_env.sh
oracle-set
oracle-run --help
oracle-test-all
```

The demanding `test_molecules` corpus is versioned under
`tests/fixtures/test_molecules` for GICForge/parser regressions.

```bash
python -m oracle gicforge corpus
python -m oracle gicforge corpus-audit
python -m oracle gicforge corpus --format paths --suffix .inp
```

Fragment workflows are planned against existing topology/synthon sections:

```bash
python -m oracle babel preprocess molecule.inp molecule.xyzin
python -m oracle fragments plan molecule.xyzin
python -m oracle fragments build molecule.xyzin
```

Validation is the gate before GICForge:

```bash
python -m oracle validate molecule.xyzin
python -m oracle gicforge plan molecule.xyzin --symmetrize --sycart
python -m oracle gicforge build molecule.xyzin --symmetrize --sycart
python -m oracle gicforge bmatrix molecule.xyzin bmat.out
python -m oracle gicforge gaussian-input molecule.xyzin job.gjf
```

## Migration Rule

Scientific behavior is migrated package by package. Existing `oracle_*`
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
SMILES imports are handled in ORACLE-Babel through RDKit when the active
environment provides it.

Python and strict Fortran77 implementations may intentionally coexist for the
same scientific kernel. In that case they are backends behind one ORACLE service
contract, with shared schemas, shared enriched XYZ sections and shared
regression or identity tests.

The enriched XYZ file is the canonical communication object between modules:
tools append or replace their own uppercase sections and preserve all unrelated
sections. See
`docs/architecture/ADR-0001-SHARED-LIBRARIES-AND-XYZ-CONTAINER.md`.

Scientific tools also support standalone `xyzin` mode: a sufficiently populated
enriched XYZ can be passed directly to SEFit/MORPHEUS, GF/PED, Thermo, DVR or
anharmonic workflows without rerunning the whole ORACLE preprocessing pipeline.

```bash
python -m oracle rovib summarize molecule.xyzin
```

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
