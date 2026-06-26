# ADR-0007: Post-GICForge Analysis Pipeline

Date: 2026-06-26

## Status

Accepted

## Context

After validation, GICForge constructs the frozen coordinate model for the
molecule. Several workflows can then start:

- GF/PED and vibrational analysis from a Cartesian Hessian;
- rotational analysis and optional SEFit/MORPHEUS;
- thermochemistry;
- anharmonic workflows such as QFF, VPT2 and VCI.

Many of these workflows need data from quantum-mechanical program outputs.
Those outputs must be read through ORACLE adapters, not separately by each
tool.

## Decision

GICForge is the first post-validation producer. It writes frozen `#GIC` and
optional `#SYCART` sections. Downstream tools consume those sections and the
shared QM adapters.

`oracle-gaussian` owns Gaussian file-format I/O:

- Gaussian input writing from enriched XYZ plus `#GIC`;
- Gaussian log/FCHK/QFF/Hessian adapters as they are migrated;
- normalized data models for downstream tools.

The user-facing GICForge workflow may request Gaussian input, but the actual
file serialization is delegated to `oracle-gaussian` to avoid duplicate
Gaussian writers.

Initial command:

```bash
python tools/oracle_run.py gicforge gaussian-input molecule.xyzin job.gjf \
  --route "#p b3lyp/def2svp opt freq"
```

Downstream ownership:

- `oracle-gf`: consumes frozen `#GIC` plus Cartesian Hessian data.
- `oracle-morpheus`: consumes rotational data, corrections and frozen
  coordinate models for SEFit.
- future `oracle-thermo`: consumes normalized frequencies, rotational
  constants and electronic energies.
- `oracle-vpt2-vci`: consumes normalized QFF/anharmonic force-field data.

## Consequences

- No downstream package reads Gaussian text directly when an ORACLE adapter
  exists.
- GICForge can produce Gaussian input as a workflow operation without owning a
  private Gaussian writer.
- Cartesian Hessian, rotational, thermochemical and anharmonic data should be
  promoted to schema-controlled ORACLE models before use by scientific tools.
