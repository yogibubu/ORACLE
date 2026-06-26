# ADR 0003: Unified Geometry And QM Parsers

Date: 2026-06-26

## Decision

ORACLE must have one parser stack for geometries and QM-program files. There
must not be parallel XYZ readers, Gaussian readers, Z-matrix readers, GUI-only
readers or workflow-specific geometry parsers.

Every parser must return shared ORACLE data models, starting with
`oracle_chem.MolecularGeometry`. Workflow packages may enrich that model, but
they must not define their own incompatible geometry containers.

## Rules

- Geometry parsers live in shared libraries, not GUI classes or workflow
  scripts.
- QM-program adapters live in program-specific packages such as
  `oracle-gaussian`, `oracle-molpro` or `oracle-mrcc`, but return shared ORACLE
  models.
- The enriched XYZ container remains the canonical handoff file.
- Parser diagnostics should be explicit and user-facing.
- A parser may preserve program-specific metadata, but coordinates, atom labels,
  charge, multiplicity and fixed-coordinate hints must use shared fields.
- New parser code must be covered by fixtures before any GUI or workflow uses it.
- Old Merlino parser modules should become compatibility wrappers around the
  ORACLE parser stack during migration.

## Initial Parser Ownership

- `oracle-chem`: plain XYZ, enriched XYZ first block, canonical
  `MolecularGeometry`, element normalization and geometry serialization.
- `oracle-gaussian`: Gaussian `.com`/`.gjf` Cartesian input, Gaussian log/out
  summaries, FCHK/QFF adapters.
- Future program packages: Molpro, MRCC and other QM formats. They must consume
  the same shared data model.

## Not Allowed

- GUI code parsing coordinate files directly.
- MORPHEUS, GICForge, GF/PED, DVR or VPT2/VCI each keeping private geometry
  parsers.
- Separate atom-symbol tables in different modules.
- Separate Gaussian Cartesian-input parsers in semiexperimental and GUI code.
- Direct downstream consumption of raw QM text when an ORACLE adapter exists.

## Migration Targets From Merlino

Known Merlino parser locations to replace or wrap:

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

The first ORACLE migration step is not to delete these files in Merlino, but to
make the ORACLE parser stack the only new implementation path.

