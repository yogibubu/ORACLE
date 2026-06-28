# MATRIX Manuals

This directory contains the user-facing manuals for the current MATRIX
workflows.

| Manual | Use it for | Main files |
| --- | --- | --- |
| MATRIX First Release Manual | Operational guide for all implemented first-release tools, shared xyzin contracts, release boundary and next-release exclusions. | `matrix_first_release_manual.tex`, `matrix_first_release_manual.pdf` |
| MORPHEUS Manual | Single-structure semiexperimental refinement, constraints, predicates, parameter classes, active coordinate models and fit diagnostics. | `morpheus_manual.tex`, `morpheus_manual.pdf` |
| NEO / GICForge Manual | Automatic construction of non-redundant GICs, ring coordinates, butterfly coordinates, out-of-plane policy, symmetry adaptation, SYCART and the Python/Fortran coordinate contract. | `gicforge_manual.tex`, `gicforge_manual.pdf` |
| GF/PED Manual | Wilson GF analysis from Cartesian Hessians, frozen-GIC Hessian transformation, Pulay-style scaling, frequencies, internal force constants and PED tables. | `gf_manual.tex`, `gf_manual.pdf` |
| Multi-Structure MORPHEUS Manual | Shared class-correction refinement across related molecules or conformers, priors, hard constraints, synthon `Zeff` typing and ensemble diagnostics. | `multistructure_manual.tex`, `multistructure_manual.pdf` |
| Remote QM Execution on oracle | Operational guide for launching Gaussian GDV32/G16, Molpro 2025.3 and ORCA jobs on the private `oracle` Linux host, copying files and promoting outputs through MATRIX adapters. | `oracle_qm_remote_manual.tex`, `oracle_qm_remote_manual.pdf` |

Recommended reading order for a new workflow:

1. Read `matrix_first_release_manual.pdf` for the first-release tool map,
   command overview and next-release boundary.
2. Read `morpheus_manual.pdf` for the general refinement model.
3. Read `gicforge_manual.pdf` if the coordinate basis, symmetry labels or ring
   coordinates need to be inspected or debugged.
4. Read `gf_manual.pdf` for harmonic GF/PED analysis using the same frozen GIC
   definition.
5. Read `multistructure_manual.pdf` for homologous-series or conformer-pair
   refinements with shared class corrections.
6. Read `oracle_qm_remote_manual.pdf` before running heavy Gaussian, Molpro or
   ORCA calculations on the remote `oracle` machine.

All PDF files are generated from the adjacent LaTeX sources with `latexmk`.
For in-program help, use `matrix help`, `matrix help TOOL`,
`matrix manuals TOOL --format markdown` or the ORACLE GUI `Help` tab. These
views are generated from the same tool-contract metadata used to validate xyzin
sections, so missing-section guidance stays aligned with the CLI and GUI.

## Manual Completion Plan

The current first-release manual pass is closed by `matrix_first_release_manual`
plus the specialist manuals listed above and the generated online help exposed
through the CLI and GUI. Fragment search, nano-LEGO assembly, production
TRINITY optimization, PySCF adapters, VSCF, vibronic spectroscopy and IP/EA
workflows belong to the next release and do not block this release.

| Manual | Status | Scope |
| --- | --- | --- |
| LINK / Preprocessing | Covered in release manual | Geometry/QM/SMILES import, RDKit use, symmetry/topology execution, Avogadro handoff and generated xyzin sections. |
| NEO / GICForge | Specialist manual plus release manual | Rename transition from GICForge to NEO, fragment-aware GICs already implemented, symmetry projector, B matrix, Python/Fortran parity and golden tests. |
| GF/PED | Specialist manual plus release manual | Standalone xyzin mode, Gaussian FCHK/log Hessian promotion, `#GF_PED`, geometry/frequency diagnostics, local force-field model, electrostatic/vdW subtraction, GIC force-constant scaling and publication tables. |
| MORPHEUS / SEFit | Specialist manual plus release manual | Standalone and GUI workflows, local `se_geometries`, multi-structure refinements and exported diagnostics. |
| Rotational Spectroscopy / WMS-Rot | Covered in release manual | Local WMS-Rot engine, line-list generation, shared diagonalizer policy, database-comparison roadmap and publication exports. |
| Vibrational Spectroscopy | Covered in release manual | IR/Raman/VCD/ROA plotting, NIST gas-phase import policy, mirrored comparisons and hybrid level1+level2 spectra with normal-mode overlap checks. |
| Thermochemistry / Kinetics | Covered in release manual | Thermo tables, rovibrational DOS, current `#THERMO` outputs and planned `#KINETICS` extension boundary. |
| VPT2/VCI | Covered in release manual | QFF input contract, VPT2/VCI run modes, CSV/report outputs and relation to the vibrational workbench. |
| DVR | Covered in release manual | Direct DVR runs, post-run normalization, shared diagonalizer use and `#DVR` GUI state. |
| QM Adapters | Covered in release manual | One-adapter-per-format policy, Gaussian/Molpro/MRCC/ORCA promotion commands, shared `#CARTESIAN_HESSIAN`, `#NORMAL_MODES`, `#QFF`, electronic, orbital and `#PROPERTIES` sections, including quadrupole coupling promotion. |
| Remote QM Execution on oracle | Specialist manual plus release manual | Remote directory layout, `matrix-submit`, `matrix-status`, GDV32/G16 policy, Molpro 2025.3, ORCA module use, file transfer and adapter promotion. |
| ORACLE GUI | Covered in release manual | Project dashboard, spectroscopy workbenches, structure/synthon view, external viewers, missing-section guidance and export workflows. |
