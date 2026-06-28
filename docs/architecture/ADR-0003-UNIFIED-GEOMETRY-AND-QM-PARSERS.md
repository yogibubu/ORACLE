# ADR 0003: Unified Geometry And QM Parsers

Date: 2026-06-26

## Decision

ORACLE must have one parser stack for geometries and QM-program files. There
must not be parallel XYZ readers, Gaussian readers, Z-matrix readers, GUI-only
readers or workflow-specific geometry parsers.

Every parser must return shared MATRIX data models, starting with
`matrix_chem.MolecularGeometry`. Workflow packages may enrich that model, but
they must not define their own incompatible geometry containers.

## Rules

- Geometry parsers live in shared libraries, not GUI classes or workflow
  scripts.
- QM-program adapters live in program-specific packages such as
  `matrix-gaussian`, `matrix-molpro`, `matrix-orca` or `matrix-mrcc`, but
  return shared MATRIX models.
- The enriched XYZ container remains the canonical handoff file.
- Parser diagnostics should be explicit and user-facing.
- A parser may preserve program-specific metadata, but coordinates, atom labels,
  charge, multiplicity and fixed-coordinate hints must use shared fields.
- New parser code must be covered by fixtures before any GUI or workflow uses it.
- Old ORACLE parser modules should become compatibility wrappers around the
  ORACLE parser stack during migration.

## Initial Parser Ownership

- `matrix-chem`: plain XYZ, enriched XYZ first block, canonical
  `MolecularGeometry`, element normalization, Z-matrix parsing and geometry
  serialization.
- `matrix-gaussian`: Gaussian `.com`/`.gjf` Cartesian input, Gaussian log/out
  summaries, Gaussian Z-matrix input through the shared Z-matrix parser and
  FCHK/QFF adapters.
- `matrix-molpro`: Molpro output geometry, charge and multiplicity adapters.
- `matrix-orca`: ORCA output geometry, charge, multiplicity, energy,
  frequencies and Cartesian Hessian adapters when the output prints those data.
- `matrix-mrcc`: MRCC output geometry, charge and multiplicity adapters.
- Future program packages must consume the same shared data model and must not
  add workflow-local parsers.

## Not Allowed

- GUI code parsing coordinate files directly.
- MORPHEUS, GICForge, GF/PED, DVR or VPT2/VCI each keeping private geometry
  parsers.
- Separate atom-symbol tables in different modules.
- Separate Gaussian Cartesian-input parsers in semiexperimental and GUI code.
- Separate Z-matrix parsers in GUI, Gaussian or workflow code.
- Direct downstream consumption of raw QM text when an ORACLE adapter exists.

## Migration Targets From ORACLE

Known ORACLE parser locations to replace or wrap:

- `gui/xyz_reader.py`
- `gui/zmat_reader.py`
- `gui/readers.py`
- `gui/gaussian.py`
- `gui/molpro.py`
- `gui/mrcc.py`
- `matrix_core/xyzin_geometry.py`
- `matrix_gaussian/parsers.py`
- `matrix_morpheus/geometry_input.py`
- `matrix_morpheus/survibfit/gaussian_log.py`
- `matrix_vpt2_vci/gaussian_qff.py`

The first ORACLE migration step is not to delete these files in ORACLE, but to
make the ORACLE parser stack the only new implementation path.
