# Migration Inventory

Date: 2026-06-26

Source repository: `/Users/vincenzobarone/merlino3.0`

Target skeleton: `/Users/vincenzobarone/ORACLE`

## Current Source State

The Merlino tree is a dirty working checkout on branch
`codex/provout-polycyclics`. The dirty state includes:

- active code changes in `merlino_core/cli.py`, `merlino_semiexp/__init__.py`,
  `merlino_semiexp/tests/test_ensemble.py` and new reference-library files;
- semiexperimental examples and benchmark input changes;
- generated analysis and paper artifacts under `doc/papers/ensemble_jpcl`;
- generated/new MSR paper artifacts under `doc/papers/newmsr`;
- a new `data/` tree;
- a new `doc/manuals/` tree;
- the new ORACLE planning document.

Do not bulk move or reset this tree. First separate source edits from generated
outputs.

## Observed Package Coupling

Highest coupling edges from Python imports:

- `gui -> merlino_semiexp`
- `gui -> merlino_fit`
- `gui -> geometry`
- `topology -> merlino_fit`
- `merlino_fit -> topology`
- `merlino_gic -> merlino_fit`
- `merlino_semiexp -> merlino_fit`
- `merlino_semiexp -> geometry`

This confirms that ORACLE should extract `core`, `chem`, `gicforge` and
`morpheus` before rewiring the GUI.

## Source Areas

Keep as migration sources:

- `merlino_core`
- `merlino_gic`
- `merlino_semiexp`
- `merlino_gf`
- `merlino_gaussian`
- `merlino_fortran`
- `merlino_dvr`
- `merlino_vpt2_vci`
- `merlino_gui`
- stable parts of `geometry`, `topology` and `merlino_fit`
- `merlino_fit/topology`, especially `AtomicSynthons`,
  `descriptor_parameters`, continuous/discrete topology, rings, aromaticity and
  topology writers. The top-level `topology/` package in Merlino is mostly a
  compatibility wrapper around this canonical implementation.
- `fortran`
- `puckering_dvr`
- `examples`
- `benchmarks`
- `doc/manuals`

Treat as runtime/generated unless explicitly reviewed:

- `working`
- `tmp`
- LaTeX build files
- generated paper fragments and analysis CSV/JSON outputs
- local GUI/runtime output directories

## Manual Contracts Imported

The following manual sources were copied from
`/Users/vincenzobarone/Desktop/newmsr_overleaf` into `docs/manuals`:

- `morpheus_manual.tex` and `.pdf`
- `gicforge_manual.tex` and `.pdf`
- `gf_manual.tex` and `.pdf`
- `multistructure_manual.tex` and `.pdf`
- `manuals_README.md`

## Immediate Hygiene Patch For Merlino

Recommended next patch in `merlino3.0`:

1. Add ignore coverage for paper build products and local analysis outputs.
2. Keep `doc/manuals/*.tex` and `doc/manuals/*.pdf`; ignore their auxiliary
   build products.
3. Decide whether generated `doc/papers/*/analysis` outputs are source
   artifacts, benchmark goldens or regenerated paper outputs.
4. Move local runtime outputs out of source control only after preserving the
   current dirty state.

## Parser Consolidation Targets

Known geometry/QM parser duplicates in Merlino should become compatibility
wrappers around ORACLE:

- `gui/xyz_reader.py`
- `gui/zmat_reader.py`
- `gui/readers.py`
- `gui/gaussian.py`
- `gui/molpro.py`
- `gui/mrcc.py`
- `merlino_core/xyzin_geometry.py`
- `merlino_gaussian/parsers.py`
- `merlino_semiexp/geometry_input.py`
- `merlino_fit/survibfit/gaussian_log.py`
- `merlino_vpt2_vci/gaussian_qff.py`

No new parser logic should be added to these paths.
