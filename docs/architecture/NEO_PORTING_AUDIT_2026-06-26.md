# NEO/GICForge Porting Audit

Date: 2026-06-26

## Rule

ORACLE ports Merlino3.0 first. A behavior can change only when Merlino3.0 has a
documented fragility or bug, and the ORACLE change must be covered by a
regression test and noted in the method documentation.

## Sources Reviewed

- `/Users/vincenzobarone/merlino3.0/merlino_gic/model.py`
- `/Users/vincenzobarone/merlino3.0/merlino_gic/gic_symmetry.py`
- `/Users/vincenzobarone/merlino3.0/merlino_fit/survibfit/primitives.py`
- `/Users/vincenzobarone/merlino3.0/merlino_fit/survibfit/transforms.py`
- `/Users/vincenzobarone/merlino3.0/fortran/gicforge`
- `/Users/vincenzobarone/Desktop/newmsr_overleaf/gicforge_manual.tex`
- `/Users/vincenzobarone/Desktop/newmsr_overleaf/main.tex`
- `/Users/vincenzobarone/Desktop/newmsr_overleaf/manuals_README.md`

## Porting Matrix

| Area | Merlino3.0 source | ORACLE status | Action |
| --- | --- | --- | --- |
| Ordinary primitives | `survibfit/primitives.py`, `mkprim.f` | Present: stretch, bend, linear bend, torsion, out-of-plane | Keep parity tests against analytic/FD B rows. |
| Ring bend/torsion blocks | `survibfit/transforms.py`, `mkcyc.f`, `mksalc.f`, Overleaf GICForge manual | Present: `CYCLIC_BEND`, `RING_PUCKER_COMPONENT`/`RPck`, `CONDENSED_RING_TORSION`, `BUTTERFLY`; Gaussian export derives `QPck`/`PhiP` from selected `RPck` pairs with the Merlino `PrtPckQP` rule; projector golden tests cover `RPck` SALCs and derived `QPck`/`PhiP` labels | Next: add more chemically diverse corpus molecules, not new core machinery. |
| Butterfly coordinates | `mksalc.f:BtFly`, `survibfit/transforms.py:ring_butterfly_u` | Present: fused-ring shared-bond torsions are classified as `BUTTERFLY`; corpus tests keep fused ring rank with RPck enabled | Next: compare coefficient vectors with executable Merlino output if the legacy backend later emits machine-readable ring-block coefficients. |
| Non-redundant reduction | `gicprune.f`, `locsvd.f`, `transforms.py` | Present: analytic B-row MGS, protected special-first policy, ring families not mixed with ordinary blocks; fused corpus tests cover selected `RPck` value and B rows | Next: extend golden cases to bridged saturated rings where bend rank can saturate before puckering rows. |
| Symmetrization | `gic_symmetry.py`, `symmetry_global.py`, `gic_type_symmetry.f`, `symm.f` | Present: Merlino label-only parity plus matrix projectors; no type mixing; total-symmetric subset stored; Fortran audit summaries now expose projector status, symmetry-group counts, special-coordinate groups, mixed-family groups, total-symmetric GIC counts, SALC coefficient counts and coefficient-normalization residuals | Next: compare coefficient vectors with executable Merlino output if the legacy backend later emits machine-readable projector coefficients. |
| Fragment/TRIC coordinates | Merlino TRIC roadmap, Overleaf manual, geomeTRIC reference model | Present: fragment center distance, fragment center-atom distance, translations, orientations, analytic B rows, Gaussian symbolic export | Next: add atom-frame angle/torsion and center-frame tilt/orientation modes. |
| Ring/bond/interaction centers | Fragment roadmap, topology ring docs | Present: bond centers, ring centers, atom-center distance candidates, analytic chain-rule B row | Next: add center-angle, center-torsion and hapticity/coordination center scoring. |
| Python/Fortran parity | `doc/GIC_PYTHON_FORTRAN_COMPARISON.md`, Fortran GICForge sources | Present for core GICForge: legacy source vendored and compiled; `frag_tric_bmat.f` mirrors special-coordinate B rows; executable harness runs Merlino on fused corpus cases and compares final rank, GIC labels and B-row subspaces while reporting projector no-mixing diagnostics | Next: keep extending the corpus; do not add new core machinery unless a regression exposes a real gap. |
| `xyzin` standalone restart | Merlino frozen schema docs | Present: frozen `#GIC`, `#SYCART`, report, B-matrix from saved sections | Next: extend downstream GF/SEFit/Thermo/VPT2 commands to consume the saved GIC section directly. |

## New Regression Anchors

- `tests/test_matrix_neo.py::test_gicforge_label_only_characters_match_merlino3`
  freezes Merlino3.0 label-only character behavior.
- `tests/test_matrix_neo.py::test_gicforge_classifies_ring_and_butterfly_primitives_like_merlino`
  freezes the Merlino ring-family separation for cyclic bend, RPck ring
  puckering, condensed-ring torsion and butterfly torsion.
- `tests/test_matrix_neo.py::test_gicforge_ring_puckering_coefficients_match_merlino_six_ring`
  freezes the six-membered-ring RPck coefficient vectors against Merlino3.0.
- `tests/test_matrix_neo.py::test_gicforge_ring_puckering_numeric_corpus_fused`
  compares selected RPck values and analytic B rows numerically on fused
  naphthalene/phenanthrene/pyrene corpus cases.
- `tests/test_matrix_neo.py::test_gicforge_point_group_projector_symmetrizes_ring_puckering_components`
  freezes point-group projector behavior for selected `RPck` sources and checks
  that symmetrized Gaussian `QPck`/`PhiP` labels remain attached.
- `tests/test_matrix_fortran_neo.py::test_legacy_merlino_ring_and_butterfly_blocks_remain_reference`
  prevents accidental removal of the strict Fortran77 ring/butterfly and
  `PrtPckQP` reference routines.
- `tests/test_matrix_fortran_neo.py::test_legacy_merlino_executable_bmatrix_span_matches_oracle_corpus`
  runs the vendored Merlino executable on naphthalene, phenanthrene and pyrene,
  then compares final rank and the Wilson-B row space with ORACLE.
- `tests/test_matrix_fortran_neo.py::test_default_fortran_audit_covers_official_golden_parity_roles`
  binds the machine-readable golden registry to the default Python/Fortran audit
  set so fused, bridged, spiro and parity-designated cases cannot silently drop
  out of the periodic gate.
- GIC reports and Fortran audit summaries now expose selected families,
  protected coordinates, skipped singular/dependent primitive details, symmetry
  source blocks and nontrivial SALC coefficient normalizations.

## Known Merlino Fragilities To Correct In ORACLE

- Label-only symmetry depends on operation spelling. ORACLE preserves it for
  compatibility, but uses matrix-classified projectors when operation matrices
  are available.
- Ring primitive generation in Merlino is split across Python and Fortran paths.
  ORACLE must centralize the saved topology/ring source in `xyzin` and make all
  backends consume the same ring records.
- Fragment/center derivatives were not one shared service in Merlino. ORACLE
  centralizes these rows in `matrix-neo` and mirrors them in
  `frag_tric_bmat.f`.

## Remaining Work Outside The Closed Core

1. Add golden corpus cases for pyridine, coronene-like fused rings and
   norcamphor/testosterone-like bridged systems.
2. Add center-angle, center-torsion and atom-frame coordinates for
   ring/metal/H-bond cases after the topology-center records are finalized.
3. Expose full strict-Fortran projector diagnostics and coefficient traces if
   downstream debugging needs operation-by-operation comparison with executable
   Merlino output.
