# ORACLE Refactor Plan

## Baseline

ORACLE starts from the final ORACLE legacy working baseline:

- Code baseline: `c704017 feat(dvr): add fortran77 path solver`
- Closure docs: `836a418 docs: clarify oracle3 closure point`
- Final ORACLE legacy freeze: `archives/freezes/oracle-legacy_freeze_20260620_141839.tar.gz`

ORACLE legacy is frozen. Do not continue structural work there.

## Refactor Goal

ORACLE must make the main functional areas independently maintainable. The GUI
should orchestrate workflows, not contain scientific or backend-specific logic.
Fortran programs should be called through narrow Python wrappers. Gaussian and
DVR parsing/writing should have one source of truth.

The planned new scientific capabilities are limited and explicit:

- Harmonic GF analysis from the independent `oracle_gf` Python package and
  `fortran/vpt2_vci/gf_core.f`.
- VPT2/VCI from Gaussian quartic force fields, with a Davidson diagonalizer for
  large VCI spaces.
- Semiexperimental equilibrium-geometry determination from least-squares fits of
  experimental rotational constants for isotopologues, using QM vibrational
  corrections.

## Target Packages

- `oracle_core`: configuration, paths, logging, manifests, job status and common
  exceptions.
- `oracle_chem`: atoms, XYZ, masses, isotopes, topology, rings and symmetry
  data structures.
- `oracle_gicforge`: ring numbering, GIC construction, Gaussian GIC coordinate
  definitions, optional symmetry adaptation, frozen GIC schemas and reusable
  B-matrix evaluation.
- `oracle_gf`: Cartesian-Hessian to GIC-B/GF/PED transformations, Pulay
  scaling and harmonic report/CSV services. This package is physically
  separated from VPT2/VCI.
- `oracle_gaussian`: Gaussian input writing, log parsing, scan extraction,
  GIC-value extraction and job metadata.
- `oracle_engines`: executable/source discovery, build checks, subprocess
  wrappers and normalized error reporting for GICForge, DVR and source-library
  Fortran kernels.
- `oracle_dvr`: Gaussian-log/grid to DVR workflows, Cremer-Pople mapping,
  Fortran bridge integration and output readers.
- `oracle_vpt2_vci`: Gaussian quartic force-field extraction, VPT2/VCI input
  preparation, VCI basis control, Fortran backend orchestration and Davidson
  diagonalization outputs.
- `oracle_morpheus`: semiexperimental equilibrium-geometry fits from
  isotopologue rotational constants, QM vibrational corrections and
  least-squares diagnostics.
- `oracle_gui`: PySide6 windows and controllers only; all scientific work goes
  through service interfaces.
- `oracle_data`: local libraries, catalog indexes and curated chemical data.

The current ORACLE legacy folders stay in place until each area has been migrated
behind tests.

## Interface Contracts

### GICForge

Input:

- Cartesian XYZ file.
- Optional backend options: symmetry, B matrix, Gaussian route/resources.

Output:

- Gaussian `.gjf`.
- Readable GIC report.
- Optional B matrix.
- Manifest JSON with executable, input checksum, options and output paths.

### DVR

Input:

- Completed Gaussian log or prepared grid CSV.
- Boundary/solver settings.
- Optional ring/Cremer-Pople labeling settings.

Output:

- Levels CSV.
- Vectors/profile CSV.
- Summary text.
- Figures where requested.
- Manifest JSON with Gaussian-log checksum, solver, backend and output paths.

### Gaussian

Input:

- Molecular model plus route/resources/workflow options.

Output:

- `.gjf` input.
- Completed `.log` where run locally.
- Parsed scan/path records.
- Manifest JSON.

### VPT2/VCI

Input:

- Gaussian output containing a quadratic/cubic/quartic force field, or an
  equivalent normalized force-field file produced by `oracle_gaussian`.
- Normal-mode metadata, harmonic frequencies and transformation data.
- Basis/cutoff settings for VCI.
- Requested number of roots and Davidson convergence settings.

Output:

- VPT2 energies and spectroscopic constants where available.
- VCI levels and dominant basis-state coefficients.
- Davidson convergence report.
- Manifest JSON with Gaussian input checksum, basis/cutoff settings, executable
  versions and output paths.

Implementation rule:

- Reuse the existing Fortran force-field/VPT2/VCI code first.
- Add Davidson as an isolated numerical component with a matrix-vector product
  interface; avoid dense VCI diagonalization for large spaces.
- Python handles parsing, input preparation, job orchestration and output
  normalization.

### Semiexperimental Equilibrium Geometry

Input:

- Experimental rotational constants for multiple isotopologues.
- Isotopic compositions and atomic masses.
- QM vibrational corrections, normally `Delta B_vib`, from Gaussian or a
  normalized ORACLE correction table.
- Initial equilibrium geometry and fit constraints.

Output:

- Fitted semiexperimental equilibrium geometry.
- Parameter covariance/correlation diagnostics.
- Residuals by isotopologue and rotational constant.
- Corrected equilibrium rotational constants.
- Manifest JSON with experimental data checksums, QM correction source,
  weighting scheme and constraints.

Implementation rule:

- Keep least-squares fitting in Python unless a Fortran backend is demonstrably
  needed.
- Keep vibrational-correction parsing in `oracle_gaussian`; `oracle_morpheus`
  should consume normalized correction data.
- Support constrained fits and fixed parameters from the beginning, because
  isotopologue data are often insufficient for a fully free structure.

## Migration Sequence

1. Add the new package skeletons and compatibility imports.
2. Move path/config/logging helpers into `oracle_core`.
3. Consolidate geometry/topology/ring primitives in `oracle_chem`.
4. Move GIC generation and ring numbering into `oracle_gicforge`.
5. Introduce `oracle_engines` wrappers for `gicforge.x`, `path_dvr.x` and
   source-only Fortran kernels such as `fortran/vpt2_vci/gf_core.f`.
6. Move Gaussian parsing/writing into `oracle_gaussian`.
7. Move DVR workflow orchestration into `oracle_dvr`.
8. Inventory existing harmonic/anharmonic Fortran code and define normalized
   input and output files.
9. Add the semiexperimental geometry data model and least-squares interface.
10. Rewire GUI controllers to call service interfaces.
11. Finalize ORACLE branding as a compatibility-managed migration: suite name,
    documentation, manifests, GUI labels and public CLI aliases first;
    package/module renames only after stable `oracle_*` compatibility imports
    and regression tests exist.
    Working expansion:
    **ORACLE = Operational Recognition of Atomistic Connectivity and Local
    Environments**.
    ORACLE is the full computational-chemistry framework: perception,
    topology, synthons/local environments, coordinate representations,
    symmetry-aware descriptors, computational spectroscopy, data/ML-assisted
    workflows, and downstream engines. Within ORACLE, the MORPHEUS block
    should identify the GICForge+SEfit layer for non-redundant internal
    coordinate generation, symmetry handling, constraints/classes/predicates,
    and semiexperimental refinement/model generation.
    Include a separate branding task to design a new polished Oracle logo
    suitable for the GUI, documentation, repository, presentations and future
    MORPHEUS/ERC material.
12. Remove compatibility wrappers only after tests cover the new imports.

Each step should end with a small commit and a green validation run.

## Rules

- Do not add new scientific logic directly inside GUI classes.
- Do not duplicate Gaussian log parsers.
- Do not call Fortran executables directly from GUI code.
- Do not hardcode absolute paths except in local environment setup.
- Keep all file-based backend contracts documented and tested.
- Prefer integration tests that exercise real input/output files.
- Keep GICForge and DVR Fortran77 backends independent compiled tools.
- Keep VPT2/VCI numerical kernels independent from Gaussian parsing and GUI
  orchestration.
- Keep semiexperimental geometry fitting independent from any specific QM
  package by consuming normalized vibrational-correction tables.
- New workflows must expose a non-Qt service or CLI path before GUI wiring.
- Every workflow run should write a `oracle.run.v1` manifest with input/output
  checksums, parameters and backend metadata.
- Project workspaces should use `inputs/`, `runs/`, `outputs/`, `reports/`,
  `cache/` and `logs/` rather than writing new files directly into `working/`.

## First Milestone

The first real milestone is a no-behavior-change extraction of `oracle_core`
and `oracle_engines`, with the existing GUI and command-line workflows still
passing `./freeze_check.sh`.
