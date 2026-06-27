# ORACLE/MORPHEUS Manuals

This directory contains the user-facing manuals for the current MORPHEUS and
GIC-related workflows.

| Manual | Use it for | Main files |
| --- | --- | --- |
| MORPHEUS Manual | Single-structure semiexperimental refinement, constraints, predicates, parameter classes, active coordinate models and fit diagnostics. | `morpheus_manual.tex`, `morpheus_manual.pdf` |
| GICForge Manual | Automatic construction of non-redundant GICs, ring coordinates, butterfly coordinates, out-of-plane policy, symmetry adaptation, SYCART and the Python/Fortran coordinate contract. | `gicforge_manual.tex`, `gicforge_manual.pdf` |
| GF/PED Manual | Wilson GF analysis from Cartesian Hessians, frozen-GIC Hessian transformation, Pulay-style scaling, frequencies, internal force constants and PED tables. | `gf_manual.tex`, `gf_manual.pdf` |
| Multi-Structure MORPHEUS Manual | Shared class-correction refinement across related molecules or conformers, priors, hard constraints, synthon `Zeff` typing and ensemble diagnostics. | `multistructure_manual.tex`, `multistructure_manual.pdf` |

Recommended reading order for a new workflow:

1. Read `morpheus_manual.pdf` for the general refinement model.
2. Read `gicforge_manual.pdf` if the coordinate basis, symmetry labels or ring
   coordinates need to be inspected or debugged.
3. Read `gf_manual.pdf` for harmonic GF/PED analysis using the same frozen GIC
   definition.
4. Read `multistructure_manual.pdf` for homologous-series or conformer-pair
   refinements with shared class corrections.

All PDF files are generated from the adjacent LaTeX sources with `latexmk`.

## Manual Completion Plan

Before the current MATRIX release is considered complete, add or refresh manuals
for the implemented tools below. Fragment search, nano-LEGO assembly and
production TRINITY optimization belong to the next release and should not block
this manual pass.

| Manual | Status | Scope |
| --- | --- | --- |
| ORACLE-Babel / Preprocessing | Needed | Geometry/QM/SMILES import, RDKit use, symmetry/topology execution, Avogadro handoff and generated xyzin sections. |
| NEO / GICForge | Existing, refresh | Rename transition from GICForge to NEO, fragment-aware GICs already implemented, symmetry projector, B matrix, Python/Fortran parity and golden tests. |
| GF/PED | Existing, refresh | Standalone xyzin mode, `#GF_PED`, local force-field model, electrostatic/vdW subtraction, GIC force-constant scaling and publication tables. |
| MORPHEUS / SEFit | Existing, refresh | Standalone and GUI workflows, local `se_geometries`, multi-structure refinements and exported diagnostics. |
| Rotational Spectroscopy / WMS-Rot | Needed | Local WMS-Rot engine, line-list generation, shared diagonalizer policy, database-comparison roadmap and publication exports. |
| Vibrational Spectroscopy | Needed | IR/Raman/VCD/ROA plotting, NIST gas-phase import policy, mirrored comparisons and hybrid level1+level2 spectra with normal-mode overlap checks. |
| Thermochemistry / Kinetics | Needed | Thermo tables, rovibrational DOS, current `#THERMO` outputs and planned `#KINETICS` extension boundary. |
| VPT2/VCI | Needed | QFF input contract, VPT2/VCI run modes, CSV/report outputs and relation to the vibrational workbench. |
| DVR | Needed | Direct DVR runs, post-run normalization, shared diagonalizer use and `#DVR` GUI state. |
| QM Adapters | Needed | One-adapter-per-format policy, Gaussian/Molpro/MRCC promotion commands, shared `#CARTESIAN_HESSIAN`, `#NORMAL_MODES`, `#QFF`, electronic and orbital sections. |
| ORACLE GUI | Needed | Project dashboard, spectroscopy workbenches, structure/synthon view, external viewers, missing-section guidance and export workflows. |
