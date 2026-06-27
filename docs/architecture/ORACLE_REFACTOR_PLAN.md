# ORACLE Refactor Plan

Date: 2026-06-26

## Goal

Create a clean ORACLE repository from the legacy workspace without losing the
working scientific code, manuals, regression fixtures or GUI workflows.

ORACLE historically meant **Operational Recognition of Atomistic Connectivity
and Local Environments**. During the active refactory it remains the repository
and GUI/application name. After the refactory is stable, the framework/package
family should be renamed **MATRIX** (**Molecular Analysis Toolkit for Reusable
Integrated eXperiments**), GICForge should become **NEO** (**Nonredundant
Equivariant Orthogonalizer**), and ORACLE should remain the GUI/user-facing
application (**Operator for Routing, Analysis, Control, Launch and
Exploration**). Existing module names remain compatibility aliases until the
package boundaries are stable and covered by tests. See
`ADR-0009-MATRIX-NAMING-TRANSITION.md`.

## Current Situation

The repository already contains useful ORACLE package boundaries:

- `oracle_core`: configuration, workspace layout, manifests and shared errors.
- `oracle_gicforge`: GICForge service, frozen GIC schemas and B-matrix evaluation.
- `oracle_morpheus`: MORPHEUS/SEfit single- and multi-structure refinement.
- `oracle_gf`: frozen-GIC GF/PED analysis.
- `oracle_qm`: shared QM tensor sections for Hessians, normal modes and QFF.
- `oracle_gaussian`: Gaussian adapters.
- `oracle_engines`: backend discovery and wrappers.
- `oracle_dvr`: DVR workflow wrapper.
- `oracle_vpt2_vci`: QFF, VPT2/VCI and Davidson prototypes.
- `oracle_gui`: newer workflow dashboard.

The same tree also still contains legacy or mixed-responsibility areas:

- `gui` and `advanced` still orchestrate scientific services directly.
- `geometry`, `topology` and the legacy fitting stack overlap and create
  cross-package dependencies.
- `working`, `tmp`, generated paper artifacts and LaTeX build files pollute
  the checkout.
- `doc/papers` mixes source manuscripts, generated tables, analysis results and
  build products.
- `puckering_dvr` is a large vendored backend and should be treated as an
  engine, not as ordinary ORACLE application code.

The current worktree is dirty. Before moving files, preserve it as a historical
state and separate code changes from generated outputs.

## Manuals Read

The `newmsr_overleaf` manuals define the target scientific contracts:

- MORPHEUS Manual: single-structure semiexperimental refinement, constraints,
  predicates, parameter classes, GIC and symmetry-Cartesian coordinate models,
  diagnostics and CLI/GUI workflow.
- GICForge Manual: deterministic non-redundant GIC construction, ring and
  butterfly coordinates, symmetry adaptation, SYCART and the Python/Fortran
  identity contract.
- GF/PED Manual: Cartesian Hessian plus frozen GIC definition, Pulay-style
  scaling, Wilson GF/PED reports and CSV outputs.
- Multi-Structure MORPHEUS Manual: shared class-correction refinement across
  related molecules or conformers, priors, hard constraints, synthon `Zeff`
  typing, ensemble diagnostics and paper artifact generation.

These manuals should become ORACLE contract documentation, not side artifacts.

## Proposed ORACLE Repository Layout

Use a monorepo first. Split into multiple repositories only after interfaces,
tests and release boundaries are stable.

```text
ORACLE/
  pyproject.toml
  README.md
  scripts/
    oracle_env.sh
  docs/
    architecture/
    manuals/
    papers/
    archive/oracle3/
  packages/
    oracle-core/
      src/oracle_core/
    oracle-chem/
      src/oracle_chem/
    oracle-gicforge/
      src/oracle_gicforge/
    oracle-morpheus/
      src/oracle_morpheus/
    oracle-gf/
      src/oracle_gf/
    oracle-qm/
      src/oracle_qm/
    oracle-gaussian/
      src/oracle_gaussian/
    oracle-molpro/
      src/oracle_molpro/
    oracle-mrcc/
      src/oracle_mrcc/
    oracle-fragments/
      src/oracle_fragments/
    oracle-rovib/
      src/oracle_rovib/
    oracle-thermo/
      src/oracle_thermo/
    oracle-engines/
      src/oracle_engines/
    oracle-dvr/
      src/oracle_dvr/
    oracle-vpt2-vci/
      src/oracle_vpt2_vci/
    oracle-gui/
      src/oracle_gui/
  engines/
    fortran/gicforge/
    fortran/dvr/
    fortran/vpt2_vci/
    vendored/puckering_dvr/
  examples/
    morpheus/
    gicforge/
    gf/
    dvr/
    vpt2_vci/
  benchmarks/
    semiexp_msr/
  tests/
    fixtures/test_molecules/
    integration/
    regression/
  tools/
    oracle_run.py
    migrate_oracle.py
```

## Central Architecture Constraint

All ORACLE modules must reuse the same libraries for the same tasks. The suite
must not grow duplicate XYZ parsers, topology builders, isotope tables, GIC
builders, Gaussian parsers, manifest writers or backend launchers.

This does not forbid intentional dual scientific kernels. Some modules may have
both Python and strict Fortran77 implementations. In ORACLE these are treated as
backends behind one service contract, not as separate tools with separate data
models. They must consume the same enriched XYZ sections, use the same shared
libraries around the kernel, record backend metadata in the manifest and share
identity/regression tests.

The canonical communication file is an enriched XYZ container. It starts with a
plain XYZ block and is progressively enriched by named uppercase sections. Each
tool owns only its own section, replaces only that section and preserves all
others. External formats are import/export adapters; downstream ORACLE modules
consume the enriched XYZ container.

The implemented standalone tool contracts are recorded in
`oracle_core.tool_contracts`. That registry is intentionally small and
machine-readable: it lists current package names, future names where decided,
required/optional sections, produced sections and the canonical standalone CLI
entry point for each tool.

This constraint is formalized in
`ADR-0001-SHARED-LIBRARIES-AND-XYZ-CONTAINER.md` and
`ORACLE_XYZ_CONTAINER.md`.

The preprocessing layer is ORACLE-Babel. It normalizes external sources into an
Avogadro-compatible enriched XYZ, determines symmetry with explicit thresholds,
builds topology once, saves synthons/descriptors and leaves downstream tools to
reuse the saved sections.

## Subproject Responsibilities

`oracle-core`

- Workspace layout: `inputs/`, `runs/`, `outputs/`, `reports/`, `cache`,
  `logs`.
- Run manifests, checksums, config, logging and shared typed errors.
- Shared sectioned enriched-XYZ utilities.
- No chemistry, GUI, Gaussian or Fortran-specific logic.

`oracle-chem`

- Atoms, masses, isotopes, XYZ, inertia, topology, rings, symmetry and primitive
  geometry utilities.
- Owns canonical geometry data models and plain/enriched XYZ parsers.
- Owns the unified Z-matrix parser. GUI, Gaussian and legacy adapters must call
  this parser rather than keeping private Z-matrix readers.
- This should absorb the stable parts of `geometry`, `topology` and topology
  pieces now under `oracle_morpheus`.
- Descriptor topology and atomic synthons are first-class chemistry primitives,
  not side tools. `AtomicSynthons`, descriptor parameters, aromaticity and
  fragment/synthon signatures belong here so MORPHEUS, GICForge, fragment
  search, GUI diagnostics and future ML/data workflows reuse the same
  descriptors.
- Owns `#VALIDATION`, the post-preprocessing gate that checks the enriched XYZ
  state before GICForge and later tools consume the molecule.
- `topology_reporting` is migrated as source material but should not become a
  public API until the remaining `survibfit` reporting dependencies are moved
  or replaced by ORACLE services.

`oracle-gicforge`

- GIC construction, ring numbering, symmetry labels, SYCART, frozen schemas,
  B-matrix evaluation and Python/Fortran comparison.
- Owns the GICForge public service, normalized schema and Python/Fortran77
  backend contract.
- Starts only from an enriched XYZ with `#VALIDATION STATUS PASS`; it then
  writes frozen `#GIC` and optional `#SYCART` sections for downstream tools.
- May create Gaussian input on request, but the file-format serialization is
  delegated to `oracle-gaussian`.

`oracle-morpheus`

- Single-structure and multi-structure semiexperimental refinement.
- Consumes frozen coordinate models from `oracle-gicforge` or
  symmetry-Cartesian data.
- Owns constraints, predicates, parameter classes, robust least squares,
  diagnostics, reports and reference-library search.

`oracle-gf`

- Frozen-GIC Hessian transformation, Wilson GF/PED, Pulay scaling, reports and
  CSV tables.
- Must stay independent from VPT2/VCI.

`oracle-qm`

- Shared enriched-XYZ contracts for Cartesian Hessians, normal modes and QFF
  data promoted by external QM adapters.
- Owns `#CARTESIAN_HESSIAN`, `#NORMAL_MODES` and `#QFF`.
- Scientific tools consume these sections instead of reparsing QM output files.

`oracle-gaussian`

- Gaussian input writing, log/FCHK/QFF parsing and normalized correction tables.
- Gaussian is an adapter, not an internal model.
- Gaussian geometry/log parsers must return shared ORACLE models such as
  `oracle_chem.MolecularGeometry`; downstream tools must not parse Gaussian text
  directly.
- Gaussian Z-matrix input delegates to the shared `oracle-chem` Z-matrix parser.
- Gaussian input generation from GICForge consumes enriched XYZ plus `#GIC`.
- Gaussian rovibrational log data is promoted to shared `#VIBRATIONAL`,
  `#ROTATIONAL` and `#DELTABVIB` sections before GF, Thermo, SEfit/MORPHEUS or
  anharmonic workflows consume it.

`oracle-molpro` / `oracle-mrcc`

- Program-specific output adapters for Molpro and MRCC geometry import.
- Return shared `oracle_chem.MolecularGeometry` objects and feed ORACLE-Babel
  preprocessing, which writes `#BASIC`, `#SOURCE`, `#SYMMETRY`, `#TOPOLOGY` and
  `#SYNTHONS`.
- Downstream GF, SEfit/MORPHEUS and anharmonic workflows must consume enriched
  XYZ sections, not Molpro/MRCC text.

`oracle-fragments`

- Nano-LEGO-like fragmentation, fragment-library search and future assembly
  workflows.
- It is a client of `#TOPOLOGY` and `#SYNTHONS`, not a separate topology
  engine.
- LCB25 molecules are first imported/preprocessed by ORACLE-Babel, then indexed
  here as full references or topology/synthon fragments.
- Query molecule fragmentation and LCB25 fragment-library lookup are inverse
  uses of the same service contract.

`oracle-rovib`

- Rotational and vibrational section contracts, DeltaVib/alpha bridge values,
  Coriolis and Q-cent compatibility.
- Provides standalone `xyzin` readers/writers for `#ROTATIONAL` and
  `#VIBRATIONAL`, plus summary tooling for existing ORACLE containers.
- Migrates ORACLE `geometry/rotational_pipeline.py`, `vibrational.py`,
  `vib_anh.py`, `rovib_pipeline.py`, `coriolis.py` and `qcent.py`.
- Treats external CeDiTT/alpha-resonances payloads as imported data that enrich
  `#ROTATIONAL`, not as private parser logic in downstream tools.

`oracle-thermo`

- Thermochemistry from enriched XYZ state.
- Reads `#BASIC`, `#ROTATIONAL` and optional `#VIBRATIONAL`.
- Owns `#THERMO` and migrates ORACLE `geometry/thermo_*` modules.

`oracle-engines`

- Discovery, build checks and subprocess wrappers for Fortran and vendored
  engines.
- GUI and scientific packages should call this layer instead of executing
  binaries directly.
- Backend metadata and executable/source checksums are recorded in manifests.

`oracle-dvr`

- Scan/grid to DVR workflow, Cremer-Pople labeling and output readers.
- Treat `puckering_dvr` as a vendored engine or adapter backend.

`oracle-vpt2-vci`

- Normal-mode QFF, VPT2/VCI basis management and Davidson diagonalization.
- Consumes normalized force-field data from `oracle-gaussian`.

`oracle-gui`

- PySide6 controllers and views only.
- No scientific algorithms, no direct Fortran calls, no duplicate parsers.

## Pipeline Contract

Every project workflow should follow the same shape:

```text
inputs/
  user-provided geometries, job files, Gaussian logs, observation tables
runs/
  timestamped immutable run directories
outputs/
  canonical machine-readable outputs
reports/
  human-readable text, HTML, LaTeX or PDF summaries
cache/
  reusable derived data
logs/
  application and backend logs
```

Every run writes an `oracle.run.v1` manifest with:

- workflow name and schema version
- input paths and SHA256 hashes
- output paths and SHA256 hashes
- parameters
- backend executable/source metadata
- git commit and dirty flag when available
- status and user-facing messages

During compatibility migration, accept and emit `oracle.run.v1` where needed,
but the new ORACLE manifest should be the forward contract.

Every workflow should also accept or produce an enriched XYZ file as the
canonical state handoff. Workflow-specific files may be generated in `runs/`,
`outputs/` or `reports/`, but the reusable state must be written back to the
appropriate XYZ section through the shared section API.

Standalone mode is mandatory for scientific packages. SEFit/MORPHEUS, GF/PED,
Thermo, DVR and VPT2/VCI must accept a sufficiently populated `xyzin` directly
and validate only the sections they consume; a full ORACLE project workspace is
recommended but not required for those direct runs.

## Migration Phases

### Phase 0: Freeze and Inventory

1. Create a clean archive or branch from the current dirty `oracle-legacy` state.
2. Split pending changes into groups: code, examples, benchmarks, papers,
   generated artifacts and local runtime outputs.
3. Add missing ignore rules for LaTeX build files and runtime outputs that are
   currently showing up in `git status`.
4. Copy the manuals from `newmsr_overleaf` into `doc/manuals` or confirm that
   the existing `doc/manuals` copy is canonical.

Exit criteria:

- Historical state preserved.
- No generated runtime output mixed with source changes.
- Existing `./freeze_check.sh` behavior documented, even if not yet green.

### Phase 1: ORACLE Skeleton

1. Create the ORACLE monorepo skeleton with workspace-aware `pyproject.toml`.
2. Add package placeholders and compatibility aliases from legacy imports to
   ORACLE packages.
3. Move only docs and metadata first; do not move scientific code in the first
   commit.
4. Add `python -m oracle` as the initial ORACLE CLI.

Exit criteria:

- `python -m oracle --help` works.
- Existing compatibility imports still work.
- No scientific behavior changed.

### Phase 2: Core and Engines

1. Extract core services into `oracle-core`.
2. Extract engine discovery and backend build checks into `oracle-engines`.
3. Normalize manifest schema to `oracle.run.v1` while preserving legacy
   manifest readers.
4. Replace direct subprocess calls in services with engine wrappers.

Exit criteria:

- Core tests pass.
- GICForge/DVR/Fortran build discovery tests pass.
- Manifest compatibility test covers both schemas.

### Phase 3: Chemistry Foundation

1. Merge stable `geometry`, `topology` and selected `oracle_morpheus.topology`
   services into `oracle-chem`.
2. Remove circular dependencies between topology and `oracle_morpheus`.
3. Define one public molecular model for atoms, coordinates, masses, graph,
   rings and symmetry metadata.

Exit criteria:

- `oracle-gicforge` and `oracle-morpheus` can consume `oracle-chem`.
- No package imports legacy `geometry` or `topology` directly except
  compatibility wrappers.

### Phase 4: GICForge and MORPHEUS

1. Move GICForge services and frozen schema models into `oracle-gicforge`.
2. Move MORPHEUS single-structure and ensemble workflows into
   `oracle-morpheus`.
3. Keep the manuals as contract docs and add examples that match each manual.
4. Add regression fixtures for GICForge, SYCART, single MORPHEUS and ensemble
   MORPHEUS. The imported ORACLE `test_molecules` corpus is the first
   demanding GICForge/parser fixture set.

Exit criteria:

- `oracle gicforge define`, `oracle morpheus fit` and
  `oracle morpheus ensemble` run without importing GUI modules.
- Existing semiexp tests pass under ORACLE names and ORACLE aliases.

### Phase 5: Analysis Engines

1. Move GF/PED into `oracle-gf`.
2. Move Gaussian adapters into `oracle-gaussian`.
3. Move DVR orchestration into `oracle-dvr`.
4. Move VPT2/VCI into `oracle-vpt2-vci`.
5. Keep GF independent from VPT2/VCI and keep Gaussian as an adapter.

Exit criteria:

- `oracle gf`, `oracle gaussian summary`, `oracle dvr`, `oracle vpt2-vci`
  expose non-GUI CLI workflows.
- Each writes an ORACLE manifest and has small fixture tests.

### Phase 6: GUI Rewire

1. Replace legacy `gui` and `advanced` direct logic with calls into
   `oracle-*` services.
2. Keep GUI state, view models and file selection in `oracle-gui`.
3. Remove scientific parsing, fitting and direct Fortran execution from GUI
   classes.

Exit criteria:

- GUI smoke tests pass.
- Service/CLI tests cover the scientific behavior behind each GUI action.

### Phase 7: Cleanup and Release

1. Move old manuscripts and historical files into `docs/archive/oracle3`.
2. Remove duplicated generated artifacts from source control.
3. Keep benchmark inputs and golden outputs, but regenerate reports through
   explicit commands.
4. Publish ORACLE 0.1 as a compatibility release with ORACLE aliases.

Exit criteria:

- Clean `git status` after a documented build/test run.
- New clone can run tests without local `working/` state.
- Documentation explains legacy-to-ORACLE mapping.

## Immediate First Tasks

1. Commit or archive the current dirty state before any moves.
2. Add a repository hygiene patch:
   ignore LaTeX build products, keep generated paper outputs out of ordinary
   status, and untrack runtime-only files after review.
3. Create ORACLE skeleton in a separate directory or branch.
4. Copy transition plans, package architecture notes, repository layout notes
   and manuals into ORACLE docs.
5. Implement only `oracle-core` and compatibility aliases first.
6. Run the smallest validation set before each migration step.

## Risk Register

- Dirty worktree: high risk of losing active scientific edits if cleanup and
  migration are mixed.
- GUI coupling: medium-high risk; GUI imports many scientific packages and must
  be rewired last.
- Topology/GIC coupling: high risk; `oracle_gicforge` and `oracle_morpheus` still
  depend on `oracle_morpheus` and legacy topology code.
- Generated paper artifacts: medium risk; useful for publications but noisy for
  source control.
- Fortran/vendored engines: medium risk; treat as stable engines with narrow
  wrappers.
- Rename risk: high if modules are renamed too early. Keep ORACLE aliases until
  ORACLE tests are green.

## Validation Matrix

- Core: workspace layout, config, manifest checksum tests.
- Chemistry: XYZ, masses, topology, ring numbering and symmetry tests.
- GICForge: Python/Fortran identity, SYCART, frozen schema and B-matrix tests.
- MORPHEUS: single-structure fit, constraints, predicates, class constraints,
  robust loss, leave-one-out and ensemble tests.
- GF/PED: frozen-GIC FCHK workflow, scaling and CSV export tests.
- Gaussian: input writer, log/FCHK/QFF parser tests.
- DVR: Gaussian log/grid to DVR args, manifest and output reader tests.
- VPT2/VCI: QFF load, basis construction, VPT2, VCI and Davidson tests.
- GUI: smoke tests only after service/CLI tests pass.
