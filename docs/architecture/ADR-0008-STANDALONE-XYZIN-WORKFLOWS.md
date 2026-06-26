# ADR-0008: Scientific Tools Support Standalone XYzin Workflows

Date: 2026-06-26

## Status

Accepted

## Context

ORACLE tools are not only steps inside a single GUI/pipeline. SEFit/MORPHEUS,
GF/PED, thermochemistry and anharmonic programs can also be run independently
when the user already has a valid `xyzin`/enriched XYZ file.

The ORACLE pipeline remains useful for creating and validating that file, but
it must not become a hard runtime dependency for every scientific tool.

## Decision

Every scientific ORACLE package must expose a standalone service and CLI path
that accepts an enriched XYZ/`xyzin` file directly.

Standalone tools must:

- use `oracle-core` section APIs;
- require only the sections that their workflow actually consumes;
- fail with explicit missing-section diagnostics when prerequisites are absent;
- preserve unrelated sections when writing results;
- optionally write manifests when a workspace/run directory is supplied, but
  not require a full ORACLE project workspace.

The pipeline view and the standalone view are the same data contract:

```text
xyzin/enriched XYZ
  -> package validates required sections
  -> package reads optional QM adapter outputs
  -> package writes its own section/report
```

Examples:

- `oracle-rovib summarize` reads `#BASIC`, `#ROTATIONAL` and `#VIBRATIONAL`
  directly from an existing `xyzin`.
- MORPHEUS/SEFit reads geometry plus `#ISOTOPOLOGUES`, and consumes `#GIC` or
  `#SYCART` depending on coordinate model.
- GF/PED reads frozen `#GIC` plus `#CARTESIAN_HESSIAN`.
- Thermochemistry reads `#BASIC`, `#ROTATIONAL` and optionally
  `#VIBRATIONAL`.
- VPT2/VCI reads `#QFF` in standalone mode and does not rebuild GICs. Gaussian
  FCHK and indexed QFF text remain adapter/compatibility entry points.
- DVR writes and reads `#DVR` request/manifest/output state. Gaussian scan logs
  remain adapter inputs for `oracle dvr prepare`, not private parsers inside
  downstream DVR clients.

## Consequences

- `xyzin` is a public interchange file, not an internal cache.
- Existing ORACLE fixtures that already contain enough sections remain useful
  regression tests.
- ORACLE-Babel and validation are recommended producers of clean state, but
  standalone tools can consume externally prepared compatible state.
- Package CLIs should prefer `--xyzin` for direct standalone mode.
- GF/PED standalone mode consumes `#CARTESIAN_HESSIAN` plus frozen `#GIC`.
- VPT2/VCI standalone loaders consume `#QFF` when present.
- DVR standalone orchestration stores discoverable run state in `#DVR`.
