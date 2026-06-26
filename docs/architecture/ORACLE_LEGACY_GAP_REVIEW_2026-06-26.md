# ORACLE legacy Gap Review

Date: 2026-06-26

Source reviewed: `/Users/vincenzobarone/oracle-legacy`

## Summary

The main package families are represented in ORACLE, but the current ORACLE
plan under-described several ORACLE responsibilities:

- standalone scientific runs from an existing `xyzin`;
- rovibrational/DeltaVib compatibility utilities;
- thermochemistry;
- local semiexperimental reference libraries;
- paper/benchmark artifact generation;
- Gaussian job launching/status helpers.

These are now tracked as explicit architecture responsibilities. They should be
migrated after the shared parser, topology and section contracts are stable.

## Reviewed ORACLE Areas

Core/package areas:

- `oracle_core`: CLI, config, manifests, workspace, `xyzin` sections and
  geometry helpers.
- `oracle_gicforge`: GICForge Python/Fortran service, frozen GIC schema, symmetry,
  B-matrix evaluation.
- `oracle_morpheus`: SEFit/MORPHEUS, ensemble fitting, MSR import, Kraitchman,
  reference library search, paper benchmarks.
- `oracle_gf`: HessianInput, Wilson GF/PED, internal-coordinate transform,
  reports and CSV export.
- `oracle_vpt2_vci`: QFF/FCHK adapters, VPT2, VCI, Davidson and workflow.
- `oracle_dvr`: DVR request and manifest preparation.
- `oracle_gaussian`: Gaussian job discovery/status plus log summaries.
- `oracle_engines`: backend discovery and build checks.
- `geometry`: rotational, vibrational, rovib, Coriolis, Q-cent and
  thermochemistry utilities.
- `topology` and `oracle_morpheus/topology`: graph/rings/synthons/descriptors.
- `advanced` and `gui`: GUI orchestration, launchers and viewers.

Data/docs:

- `data/se_geometries`: local SE reference geometry library.
- `benchmarks/semiexp_msr`: regression and paper benchmark inputs/goldens.
- `doc/XYZIN_FORMAT.md`: historical `xyzin` section contract.
- `doc/ROVIB_COMPATIBILITY_INTERFACE_2026-03-24.md`: DeltaVib/alpha bridge
  contract.
- `doc/VPT2_VCI_CORE.md`: GF/VPT2/VCI separation and adapter policy.
- `doc/FORTRAN_BACKENDS.md`: active Fortran backends and Python ownership.

## Gaps Added To ORACLE Plan

`oracle-rovib`

- Owns normalized `#ROTATIONAL`, `#VIBRATIONAL`, `#DELTABVIB`, `#CORIOLIS`
  and `#QCENT` contracts.
- Migrates ORACLE `geometry/rotational_pipeline.py`, `vibrational.py`,
  `vib_anh.py`, `rovib_pipeline.py`, `coriolis.py` and `qcent.py`.
- Treats CeDiTT/alpha-resonances output as an external compatibility payload,
  not as ORACLE-owned theory unless explicitly migrated later.

`oracle-thermo`

- Owns `#THERMO`.
- Migrates `geometry/thermo_trasl.py`, `thermo_rot.py`, `thermo_vib.py`,
  `thermo_pipeline.py` and `thermo_writer.py`.
- Reads `#BASIC`, `#ROTATIONAL` and optional `#VIBRATIONAL` from `xyzin`.

Standalone `xyzin` mode

- SEFit/MORPHEUS, GF/PED, Thermo, DVR and VPT2/VCI must be runnable from a
  sufficiently populated `xyzin` file.
- ORACLE pipeline commands create good state, but package CLIs must not require
  the user to rerun import/preprocess when the file already contains the needed
  sections.

Reference libraries

- LCB25 is already cached in ORACLE.
- ORACLE `data/se_geometries` remains a separate source to import or merge
  into the ORACLE library/index layer after topology/synthon validation.

Gaussian job helpers

- ORACLE has parsers and an input writer scaffold.
- ORACLE's `oracle_gaussian.jobs` behavior, including `gauin`/`gauin.gjf`
  selection and normal-termination status, remains to migrate into
  `oracle-gaussian`.

Paper/benchmark generation

- `semiexp-benchmark`, `semiexp-ensemble-paper` and generated publication
  artifacts are not core runtime services.
- They should migrate as benchmark/report commands after MORPHEUS core APIs are
  stable.

## Not Missing, But Still Pending

- Python/Fortran77 dual backends are already allowed by ADR-0002.
- GICForge, GF/PED, DVR and VPT2/VCI package placeholders already exist.
- Fragment/Nano-LEGO direction is covered by `oracle-fragments`.
- Unified geometry/Gaussian/Z-matrix parsing is started but not complete for
  all QM programs (`Molpro`, `MRCC`, FCHK/QFF details still need migration).
- GIC symmetry is an ORACLE strength and must preserve Merlino3.0 behavior.
  The current ORACLE projector covers tested low-symmetry, `C/D/Dnh/Dnd`
  families plus matrix-classified `Td`, `Oh` and `Ih`. Remaining parity work is
  strict Fortran77 projector diagnostics equivalent to Python and regression
  fixtures for ring/butterfly/special coordinates under symmetry.

## Legacy Namespace Cleanup 2026-06-26

Closed in the first MORPHEUS cleanup pass:

- `merlino_core` service imports in MORPHEUS were replaced by `oracle_core`
  (`manifest`, `numerics`, `xyzin_geometry`, `isotopologues`, `paths`,
  scientific errors and BDPCS3 parameters).
- Generic geometry utilities used by SEfit/Kraitchman were moved to
  `oracle_chem` (`average_atomic_masses`, `isotopes_table`,
  `physical_constants`, `structure`, `inertia`, `rotational`).
- XYZ tuple I/O and XYZ writing now come from `oracle_chem.geometry_io` instead
  of `survibfit.modify_geom`.
- Topology/radius/element imports in MORPHEUS and the GIC kernel now use
  `oracle_chem.topology`.
- GICForge executable discovery now uses `oracle_engines` instead of
  `merlino_fortran`.
- Runtime aliases for `merlino_core`, `merlino_fortran`, `geometry` and
  `topology` were removed.

Remaining intentional boundary:

- `merlino_fit.survibfit` and `merlino_gic` are still the active SEfit/GIC
  kernel interfaces. They should be ported as real ORACLE kernel modules or
  replaced by adapters to `oracle_gicforge`; they should not be renamed
  cosmetically.
- Only `merlino_fit` and `merlino_gic` remain as runtime aliases, because they
  are still the active SEfit/GIC kernel boundary.

Validation for this cleanup pass: `PYTHONPATH=. pytest -q` returned
`125 passed`.
