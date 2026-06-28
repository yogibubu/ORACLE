# MATRIX

**Molecular Analysis Toolkit for Reusable Integrated eXperiments**

MATRIX is the suite-level framework for the refactor of the legacy
MORPHEUS/structural-analysis codebase. The ORACLE name is kept for the
user-facing GUI/orchestrator. The first releases keep compatibility aliases
while the scientific packages are separated behind stable service, CLI and
manifest contracts.

## Initial Scope

- `matrix-core`: workspace layout, configuration, manifests, checksums and
  common errors.
- `matrix-chem`: atoms, masses, geometry, topology, rings and symmetry.
  Descriptor topology and atomic synthons live here as shared chemistry
  primitives. It owns the canonical `MolecularGeometry` model and plain/enriched
  XYZ parsers.
- `matrix-neo`: non-redundant GIC construction, frozen schemas, SYCART and
  B-matrix evaluation.
- `matrix-morpheus`: semiexperimental geometry refinement, constraints,
  predicates, parameter classes, diagnostics and ensemble refinement.
- `matrix-gf`: frozen-GIC GF/PED analysis.
- `matrix-qm`: shared QM data sections for Cartesian Hessians, normal modes,
  QFF/anharmonic force fields, electronic data and normalized QM properties.
- `matrix-gaussian`: Gaussian input/output adapters, job status/run helpers,
  FCHK/QFF promotion and Gaussian log rovibrational promotion.
- `matrix-molpro`: Molpro launch helpers and output adapters returning shared
  MATRIX geometry and enriched XYZ state.
- `matrix-orca`: ORCA launch/status helpers and output adapters for final
  geometry plus Cartesian Hessian data when printed by ORCA.
- `matrix-mrcc`: MRCC output adapters returning shared MATRIX geometry and
  enriched XYZ state.
- `matrix-link`: LINK preprocessing/import
  adapter for external sources and molecular databases such as LCB25.
- `matrix-fragments`: topology-backed fragmentation, fragment-library search
  and future assembly contracts.
- `matrix-rovib`: rotational, vibrational and rovibrational compatibility
  sections.
- `matrix-thermo`: thermochemistry from enriched XYZ state.
- `matrix-engines`: Fortran and vendored backend discovery/build wrappers.
- `matrix-dvr`: scan/grid to DVR workflows.
- `matrix-vpt2-vci`: QFF, VPT2/VCI and Davidson workflows.
- `matrix-oracle`: GUI controllers and views only.

LCB25 geometries are managed as a reproducible local cache rather than
committed source files:

```bash
python -m matrix lcb25 fetch
```

MATRIX shell helpers are available without editing personal startup files
automatically. `matrix-set` creates `MATRIX_VENV` when missing and installs the
core runtime dependencies, including RDKit for SMILES imports and
Pandas/SymPy for the vendored WMS-Rot rotational engine:

```bash
source /Users/vincenzobarone/MATRIX/scripts/matrix_env.sh
matrix-set
matrix-run --help
matrix-test-all
```

Full MATRIX installation notes, including visualization programs
Avogadro, Molden, MOrbVis browser fallback and XQuartz setup, are in
[`docs/INSTALL_MATRIX.md`](docs/INSTALL_MATRIX.md).

Tool contracts, manual links and missing-section guidance are available from
the same metadata used by the GUI:

```bash
python -m matrix help
python -m matrix help gf --xyzin molecule.xyzin
python -m matrix manuals NEO --format markdown
python -m matrix properties summary molecule.xyzin
```

The demanding `test_molecules` corpus is versioned under
`tests/fixtures/test_molecules` for GICForge/parser regressions.

```bash
python -m matrix gicforge corpus
python -m matrix gicforge corpus-audit
python -m matrix gicforge corpus --format paths --suffix .inp
```

Fragment workflows are planned against existing topology/synthon sections:

```bash
python -m matrix link preprocess molecule.inp molecule.xyzin
python -m matrix fragments plan molecule.xyzin
python -m matrix fragments build molecule.xyzin
```

Validation is the gate before GICForge:

```bash
python -m matrix validate molecule.xyzin
python -m matrix gicforge plan molecule.xyzin --symmetrize --sycart
python -m matrix gicforge build molecule.xyzin --symmetrize --sycart
python -m matrix gicforge bmatrix molecule.xyzin bmat.out
python -m matrix gicforge gaussian-input molecule.xyzin job.gjf
```

QM adapters and Gaussian job helpers are available directly:

```bash
python -m matrix gaussian summary gauin.log
python -m matrix gaussian status .
python -m matrix gaussian run . --executable gdv
python -m matrix gaussian formchk gicforge.chk gicforge.fchk
python -m matrix gaussian fchk-summary gicforge.fchk
python -m matrix gaussian promote-fchk gaussian.fchk molecule.xyzin
python -m matrix gaussian promote-rovib gauin.log molecule.xyzin
python -m matrix molpro status molpro_workdir
python -m matrix molpro run molpro_workdir --input molecule.com
python -m matrix molpro promote molpro.out molecule.xyzin
python -m matrix orca status orca_workdir
python -m matrix orca run orca_workdir --input molecule.inp
python -m matrix orca summary orca.out
python -m matrix orca promote orca.out molecule.xyzin
python -m matrix mrcc promote mrcc.out molecule.xyzin
```

Remote QM jobs on `oracle` use the same launcher contract, with SSH/SCP used
only as transport and no private parser on the remote machine:

```bash
python -m matrix qm remote-submit job.gjf --engine gdv32 --host enzo@oracle
python -m matrix qm remote-status --host enzo@oracle
python -m matrix qm remote-fetch JOB_NAME --host enzo@oracle --dest runs
python -m matrix qm remote-fetch JOB_NAME --host enzo@oracle \
  --dest runs --promote molpro --xyzin molecule.xyzin
python -m matrix qm remote-fetch JOB_NAME --host enzo@oracle \
  --dest runs --promote orca --xyzin molecule.xyzin
```

Thermo and rovibrational utilities run from the same enriched `xyzin` state:

```bash
python -m matrix thermo molecule.xyzin
python -m matrix rovib vibin molecule.xyzin --fchk gaussian.fchk
python -m matrix rovib coriolis molecule.xyzin --out coriolis.report
python -m matrix rovib qcent molecule.xyzin --out qcent.report
python -m matrix rovib dos molecule.xyzin
python -m matrix rovib dos-rovib molecule.xyzin
python -m matrix gf --xyzin molecule.xyzin
python -m matrix vpt2-vci --xyzin molecule.xyzin --run-dir runs/vpt2_vci
python -m matrix vpt2-vci --collect molecule.xyzin
python -m matrix dvr prepare scan.log --outdir runs/dvr --xyzin molecule.xyzin
python -m matrix dvr run scan.log --outdir runs/dvr --xyzin molecule.xyzin
python -m matrix dvr run --xyzin molecule.xyzin
python -m matrix dvr collect molecule.xyzin
```

## Migration Rule

Scientific behavior is migrated package by package. Existing `oracle_*`
imports stay valid until MATRIX-native tests cover the new public APIs.

## Core Architecture Rule

MATRIX modules must not reinvent shared operations. Common tasks such as
sectioned XYZ I/O, atom and isotope data, topology, symmetry, GIC construction,
QM-program parsing, backend execution and manifests belong to shared libraries.

Geometry and QM-program parsers are unified: they return shared MATRIX models
and live in parser packages, not GUI or workflow modules.

The preprocessing layer is LINK: it imports external sources, writes an
Avogadro-compatible enriched XYZ, determines symmetry with explicit thresholds,
builds topology once, and saves descriptors/synthons for downstream tools.
SMILES imports are handled in LINK through RDKit when the active environment
provides it. The legacy `babel` subcommand remains a compatibility alias for
`matrix link ...`; `python -m oracle ...` remains a deprecated CLI alias for
`python -m matrix ...`.

Python and strict Fortran77 implementations may intentionally coexist for the
same scientific kernel. In that case they are backends behind one MATRIX service
contract, with shared schemas, shared enriched XYZ sections and shared
regression or identity tests.

The enriched XYZ file is the canonical communication object between modules:
tools append or replace their own uppercase sections and preserve all unrelated
sections. See
`docs/architecture/ADR-0001-SHARED-LIBRARIES-AND-XYZ-CONTAINER.md`.

Scientific tools also support standalone `xyzin` mode: a sufficiently populated
enriched XYZ can be passed directly to SEFit/MORPHEUS, GF/PED, Thermo, DVR or
anharmonic workflows without rerunning the whole LINK preprocessing pipeline.

```bash
python -m matrix rovib summarize molecule.xyzin
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
