# MATRIX First Release Closure

Status: closed for the first release baseline.

This release establishes MATRIX as the framework layer around the implemented
ORACLE tools. The release is closed on the condition that the repository stays
green under the standard test suite and the manuals in `docs/manuals` remain
generated from their LaTeX sources.

## Included Scope

- ORACLE-Babel preprocessing from XYZ, QM adapters and SMILES/RDKit.
- Shared enriched `xyzin` section contracts.
- NEO/GICForge GIC construction, symmetrization, special coordinates and B
  matrices.
- GF/PED, including `#GF_PED`, local force-field mode and GIC force-constant
  scaling.
- MORPHEUS/SEFit and multi-structure MORPHEUS workflows.
- Rotational spectroscopy through the local WMS-Rot engine.
- Vibrational IR/Raman/VCD/ROA spectra, mirrored comparisons, NIST gas-phase
  IR import and hybrid `harmonic(level1)+[anharmonic-harmonic](level2)` spectra
  with normal-mode overlap checks.
- Thermochemistry tables and rovibrational DOS workflows.
- VPT2/VCI from normalized `#QFF`.
- DVR preparation/run state and shared diagonalizer routing.
- QM adapters for shared Hessian, normal-mode, QFF, electronic, transition and
  orbital sections.
- ORACLE GUI as the user-facing MATRIX application.

## Deferred Scope

The following items start the next release and are not blockers for this one:

- robust fragment search and fragmentation/assembly;
- nano-LEGO;
- LCB25 fragment indexing and assembly heuristics;
- production TRINITY external energy/gradient optimization loops;
- PySCF launcher and analysis adapters;
- VSCF for harmonically coupled anharmonic oscillators in local modes;
- vibronic spectroscopy;
- ionization-potential and electron-affinity workflows;
- online help/manual generation for the most important CLI and GUI tools;
- rotational database comparison adapters for CDMS/JPL/Splatalogue;
- deeper WMS-Rot homologation and removal of temporary `pandas` compatibility
  boundaries.

## Documentation Gate

The first-release manual set is:

- `matrix_first_release_manual.pdf`;
- `gicforge_manual.pdf`;
- `gf_manual.pdf`;
- `morpheus_manual.pdf`;
- `multistructure_manual.pdf`.

These manuals are the release contract for implemented tools. Future work must
update the relevant manual at the same time as code and tests.
