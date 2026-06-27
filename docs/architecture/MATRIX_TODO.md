# MATRIX TODO

## Current Release Completion

- Finish user-facing manuals for the implemented tools before calling MATRIX
  complete at this release level.
- Keep the manuals focused on the tools that already have tested command-line
  and GUI workflows: ORACLE-Babel/preprocessing, NEO/GICForge, GF/PED,
  MORPHEUS/SEFit, rotational spectroscopy/WMS-Rot, vibrational spectroscopy,
  thermo/kinetics tables, VPT2/VCI, DVR, QM adapters and the ORACLE GUI.
- Each manual must document the shared xyzin sections consumed and produced by
  the tool, standalone command usage, GUI workflow usage, publication exports,
  required external programs and the regression fixtures that define expected
  behavior.
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
  `oracle_core.diagonalizer`.
- Do not add private DVR diagonalizer wrappers in new Python DVR workflows.
- Add GPU/performance regression cases when realistic large DVR fixtures are
  available.

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
