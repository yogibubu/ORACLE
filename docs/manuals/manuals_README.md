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
