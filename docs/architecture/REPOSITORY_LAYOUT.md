# Repository Layout

This file is the operational map for the ORACLE repository.

## Root Policy

The repository root contains only launchers, freeze scripts, git metadata, and
the top-level README. New implementation files should not be added at root.

## Runtime Code

- `packages/oracle-core/src/oracle_core/`: CLI, manifests, workspace handling,
  common errors and shared sectioned-xyz helpers.
- `packages/oracle-chem/src/oracle_chem/`: chemistry primitives, topology and
  preprocessing services used by all tools.
- `packages/oracle-gicforge/src/oracle_gicforge/`: frozen GIC definition,
  symmetry adaptation and reusable B-matrix evaluation.
- `packages/oracle-gf/src/oracle_gf/`: harmonic Cartesian-Hessian/GIC-B/GF/PED
  package, physically separated from VPT2/VCI.
- `packages/oracle-vpt2-vci/src/oracle_vpt2_vci/`: normal-mode anharmonic QFF,
  VPT2, VCI and Davidson services.
- `packages/oracle-dvr/src/oracle_dvr/`: scan/grid to DVR workflow services.
- `packages/oracle-morpheus/src/oracle_morpheus/`: semiexperimental refinement
  and SEfit workflows.
- `packages/oracle-gui/src/oracle_gui/`: GUI controllers and views.
- `engines/fortran/`: active Fortran backends. GICForge lives in
  `engines/fortran/gicforge`; DVR lives in `engines/fortran/dvr`; independent
  source kernels for harmonic GF, normal-mode VPT2/VCI and Davidson live in
  `engines/fortran/vpt2_vci`.
- `engines/puckering_dvr/`: vendored Python backend for mass-weighted path DVR
  analysis from completed Gaussian scan/path logs.

## Generated And Local Files

- `working/`: GUI/runtime work area. Ignored by git.
- `projects/`: local data libraries. Ignored by git.
- `archives/freezes/`: local freeze archives and checksums. Ignored by git.
- `engines/fortran/*/build/`: compiler logs and build intermediates. Ignored by git.

## Documentation

- `docs/architecture/ORACLE_LEGACY_GAP_REVIEW_2026-06-26.md`: current
  migration/gap status snapshot.
- `docs/architecture/PACKAGE_ARCHITECTURE.md`: package boundaries and workflow
  data flow.
- `docs/architecture/ORACLE_GIC_METHOD.md`: GIC construction, reduction,
  symmetrization and Python/Fortran parity notes.
- `docs/archive/`: cleanup, freeze and triage history.
- `docs/reports/`: larger reports and TeX artifacts.
- Historical ORACLE legacy notes are not duplicated in ORACLE; use the frozen
  `oracle-legacy` archive for that material.

## Fortran Build Policy

Fortran sources are kept in fixed-form Fortran-compatible style. Build scripts
should use `-std=legacy` and suppress compiler deprecation noise from historical
constructs, while still failing on real compilation/link errors.

`engines/fortran/gicforge/compile_legacy` is the canonical GICForge build command.
It updates:

- `engines/fortran/gicforge/build/gicforge_legacy`

`engines/fortran/dvr/compile_MAC` is the canonical Fortran DVR build command. It updates:

- `engines/fortran/dvr/build/path_dvr`
- `bin/path_dvr.x`

`engines/fortran/vpt2_vci/compile_check` verifies the independent GF/VCI/Davidson
source kernels. It updates:

- `engines/fortran/vpt2_vci/build/gf_core.o`
- `engines/fortran/vpt2_vci/build/vci_core.o`
- `engines/fortran/vpt2_vci/build/vpt2_core.o`
- `engines/fortran/vpt2_vci/build/davidson_core.o`
