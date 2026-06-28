# MATRIX TODO

## Current Release Maintenance

- The first release is documented in `docs/architecture/MATRIX_FIRST_RELEASE_CLOSURE.md`.
- Naming cleanup is the first maintenance task before new scientific features:
  public documentation, manuals and user-facing commands must consistently use
  MATRIX for the framework/container, ORACLE for the GUI/orchestrator and NEO
  for the GICForge coordinate engine. Runtime compatibility aliases stay active
  until package renames are planned and tested.
- Keep the user-facing manuals aligned with implemented command-line and GUI
  workflows: LINK/preprocessing, NEO/GICForge, GF/PED, MORPHEUS/SEFit,
  rotational spectroscopy/WMS-Rot, vibrational spectroscopy, thermo tables,
  VPT2/VCI, DVR, QM adapters and the ORACLE GUI.
- When code changes, update the relevant manual in the same commit as tests.
- Do not move parsing, topology, diagonalization or plotting logic into manual
  examples. Examples must use the public MATRIX commands and shared libraries.

## Next Release Scope

- Fragment search, robust fragmentation/assembly and the nano-LEGO workflow are
  next-release work. The current release keeps only the existing fragment
  sections and fragment-aware consumers needed by implemented tools.
- TRINITY external energy/gradient geometry optimization remains a next-release
  production feature. The current release keeps its schema, command skeleton and
  GUI entry point so workflows can be designed without blocking the implemented
  spectroscopy/refinement tools.
- LCB25-backed fragment libraries should stay in the next-release fragment
  track, including one-time database caching, fragment indexing and assembly
  heuristics.
- Any next-release fragment/nano-LEGO/TRINITY implementation must consume the
  existing topology, synthon, GIC and xyzin contracts instead of introducing
  private parsers or molecular graph logic.

## Post-First-Release Scientific Roadmap

- Add PySCF as an optional electronic-structure backend with MATRIX-managed
  installation notes, job launcher, result normalization and analysis adapters.
  PySCF outputs must be promoted into the same shared QM/electronic sections
  used by Gaussian, Molpro and MRCC.
- Add VSCF for harmonically coupled anharmonic oscillators in local modes. The
  implementation should reuse existing GF/GIC force constants, QFF/DVR
  contracts and the shared diagonalizer instead of introducing private
  vibrational data structures.
- Extend the implemented GF scaling-class frontend with automatic class
  suggestions from `#SYNTHONS`. MATRIX already supports class records in
  `--scale-file`, `--scale-class` and `--scale-preview`; future work should add
  default classes for XH stretchings, XY stretchings, HXH/HXY/XYZ bendings,
  out-of-plane coordinates, Q and phi ring coordinates and dihedrals. Store
  method-specific factors in a small library when calibrated values are known;
  otherwise require user-supplied factors. A finer optional split by synthon
  descriptors can refine these default classes later. Do not port `vibint`
  Gaussian-log parsing or private GF reconstruction.
- Add vibronic spectroscopy workflows, including electronic-transition,
  normal-mode and Franck-Condon/Herzberg-Teller data contracts, with plotting
  and publication exports in the electronic/vibrational workbenches.
- Add ionization-potential and electron-affinity workflows. Supported modes
  should include vertical and adiabatic IP/EA, spin/state bookkeeping,
  provenance of charge states and connection to electronic spectra.
- Verify where Watson A/S rotational-constant conversions are implemented,
  including quartic and sextic terms, and document the public entry points.
- Add run manifests for every multi-step workflow so each publication figure or
  table can be traced back to input sections, command-line options, external
  executables and MATRIX version.
- Import or redesign the validated one-mode Gaussian rovibrational probe from
  `/Users/vincenzobarone/merlino_fit/scripts/gaussian_one_mode_probe.py` as a
  MATRIX diagnostic workflow.  It should reuse the shared Gaussian/QFF adapters,
  normal-mode contracts and diagonalizer, and should keep the validated
  `dI/dQ`, `I(Q)` and `1/I(Q)` checks.
- Before deleting `/Users/vincenzobarone/newmerlin`, extract the remaining
  useful historical material listed in
  `docs/architecture/MATRIX_LEGACY_REPO_AUDIT_2026-06-28.md`: `vci1d.f`,
  `allvib.f`, `rates.f`, `source/newrot/newkra.f`, `der12.F`, old notes and
  the `.inp` fixtures not yet in the MATRIX corpus.
- Add a release-quality benchmark corpus for end-to-end workflows, including
  small fast fixtures for CI and larger demanding molecules for periodic
  numerical audits.
- Add packaging/versioning policy for the MATRIX transition, including final
  command aliases, deprecation messages for ORACLE/Merlino names and a
  reproducible environment lock for optional heavy dependencies.
- Execute the MATRIX naming transition in stages:
  1. freeze the canonical naming table in `ADR-0009`;
  2. update top-level docs, manuals and GUI labels;
  3. make public commands prefer `matrix ...`, `link ...`, `neo ...` and
     `matrix-oracle` style entry points while keeping old aliases;
  4. rename distributions/packages only after tests and downstream scripts are
     covered by compatibility shims;
  5. remove legacy ORACLE/Merlino wording only after one deprecation cycle.
- Keep `ruff check packages tests` as a CI correctness gate. The active ruff
  profile gates syntax/undefined-name failures now; broader style cleanup for
  legacy ports, scientific index names and vendored WMS-Rot code remains a
  separate maintenance track and must not be mixed with scientific feature
  changes.

## WMS-Rot Homologation

- Remove the hard `pandas` dependency from MATRIX-facing WMS-Rot services after
  the line-list, assignment and CSV layers have typed ORACLE records and golden
  regression tests.
- Keep `pandas` available only as a temporary compatibility dependency for the
  imported first-party WMS-Rot engine.
- Move WMS-Rot adapter boundaries toward ordinary MATRIX services: shared
  xyzin sections, run manifests, publication exports, shared diagonalizer and
  shared GUI workflow state.
- Do not add Gaussian, QM-output or topology parsing inside WMS-Rot. All such
  data must enter through shared MATRIX adapters and enriched xyzin sections.

## DVR Diagonalization

- Keep large DVR Hamiltonian diagonalizations behind
  `matrix_core.diagonalizer`.
- Do not add private DVR diagonalizer wrappers in new Python DVR workflows.
- Add GPU/performance regression cases when realistic large DVR fixtures are
  available.
- Implement the large-amplitude metric policy before treating GF-generated
  DVR plans as production Hamiltonians.  The equilibrium `G_INV_BLOCK` from
  GF is only a reference; production DVR must either rebuild B/G and
  \(G^{-1}_{SS}\) at each grid geometry (`GRID_RECOMPUTE`) or include
  derivative corrections to the metric (`METRIC_DERIVATIVES`).  A constant
  reference metric is allowed only as an explicitly declared screening/test
  approximation.

## Experimental Spectrum Databases

- Add a shared experimental-spectrum adapter layer for external databases
  instead of hardcoding database queries inside spectroscopy GUIs.
- Vibrational: NIST IR data may be fetched automatically only for gas-phase
  JCAMP records; condensed-phase and missing records require user instruction.
- Rotational: add database-backed comparison for predicted line lists against
  gas-phase microwave/submillimeter catalogs such as CDMS/JPL/Splatalogue when
  a reliable molecule identifier and frequency/unit contract are available.
- Keep database provenance, query URL, phase/state and selection filters in the
  export metadata so publication plots remain traceable.
