# Migration Inventory

Date: 2026-06-26

Source repository: `/Users/vincenzobarone/oracle-legacy`

Target skeleton: `/Users/vincenzobarone/ORACLE`

## Current Source State

The ORACLE tree is a dirty working checkout on branch
`codex/provout-polycyclics`. The dirty state includes:

- active code changes in `oracle_core/cli.py`, `oracle_morpheus/__init__.py`,
  `oracle_morpheus/tests/test_ensemble.py` and new reference-library files;
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

- `gui -> oracle_morpheus`
- `gui -> oracle_morpheus`
- `gui -> geometry`
- `topology -> oracle_morpheus`
- `oracle_morpheus -> topology`
- `oracle_gicforge -> oracle_morpheus`
- `oracle_morpheus -> oracle_morpheus`
- `oracle_morpheus -> geometry`

This confirms that ORACLE should extract `core`, `chem`, `gicforge` and
`morpheus` before rewiring the GUI.

## Source Areas

Keep as migration sources:

- `oracle_core`
- `oracle_gicforge`
- `oracle_morpheus`
- `oracle_gf`
- `oracle_gaussian`
- `oracle_engines`
- `oracle_dvr`
- `oracle_vpt2_vci`
- `oracle_gui`
- stable parts of `geometry`, `topology` and `oracle_morpheus`
- `geometry/rotational_pipeline.py`, `geometry/vibrational.py`,
  `geometry/vib_anh.py`, `geometry/rovib_pipeline.py`,
  `geometry/coriolis.py` and `geometry/qcent.py` as `oracle-rovib` sources
- `geometry/thermo_trasl.py`, `geometry/thermo_rot.py`,
  `geometry/thermo_vib.py`, `geometry/thermo_pipeline.py` and
  `geometry/thermo_writer.py` as `oracle-thermo` sources
- `oracle_morpheus/topology`, especially `AtomicSynthons`,
  `descriptor_parameters`, continuous/discrete topology, rings, aromaticity and
  topology writers. The top-level `topology/` package in ORACLE is mostly a
  compatibility wrapper around this canonical implementation.
- `fortran`
- `puckering_dvr`
- `examples`
- `benchmarks`
- `doc/manuals`
- `data/se_geometries`, as a separate local SE reference-library source to
  review alongside LCB25

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

## Immediate Hygiene Patch For ORACLE

Recommended next patch in `oracle-legacy`:

1. Add ignore coverage for paper build products and local analysis outputs.
2. Keep `doc/manuals/*.tex` and `doc/manuals/*.pdf`; ignore their auxiliary
   build products.
3. Decide whether generated `doc/papers/*/analysis` outputs are source
   artifacts, benchmark goldens or regenerated paper outputs.
4. Move local runtime outputs out of source control only after preserving the
   current dirty state.

## Parser Consolidation Targets

Known geometry/QM parser duplicates in ORACLE should become compatibility
wrappers around ORACLE:

- `gui/xyz_reader.py`
- `gui/zmat_reader.py`
- `gui/readers.py`
- `gui/gaussian.py`
- `gui/molpro.py`
- `gui/mrcc.py`
- `oracle_core/xyzin_geometry.py`
- `oracle_gaussian/parsers.py`
- `oracle_morpheus/geometry_input.py`
- `oracle_morpheus/survibfit/gaussian_log.py`
- `oracle_vpt2_vci/gaussian_qff.py`

No new parser logic should be added to these paths.

## Standalone XYzin Workflows

ORACLE scientific tools can be launched directly from a populated `xyzin`
file. ORACLE must preserve this mode. SEFit/MORPHEUS, GF/PED, Thermo, DVR and
VPT2/VCI should validate the sections they need and run without forcing a full
ORACLE preprocessing workspace in the same command.
