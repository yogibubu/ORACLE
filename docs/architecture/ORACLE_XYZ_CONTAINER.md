# ORACLE Enriched XYZ Container

The enriched XYZ file is the canonical communication file between ORACLE
modules.

It starts as a normal XYZ file and is enriched step by step by tools that append
or replace named sections. This is the pipeline spine: each tool can be run
independently, but it consumes the same container and preserves sections owned
by other tools.

## Non-Negotiable Rules

- Use shared ORACLE APIs to read, write and replace sections.
- Preserve unrelated sections exactly when updating a file.
- Do not duplicate parsers in GUI, scripts or individual scientific modules.
- Treat external file formats as adapters. Convert them into enriched XYZ
  sections before downstream consumption.
- Give every section a `SCHEMA` line.
- Keep section names uppercase.

## Minimal Flow

```text
plain XYZ
  -> oracle-chem adds #TOPOLOGY and #SYMMETRY
  -> oracle-fragments adds #FRAGMENTS, #FRAGMENT_LIBRARY or #ASSEMBLY
  -> oracle-chem adds #VALIDATION
  -> oracle-gicforge adds #GIC and optionally #SYCART
  -> branches:
       oracle-gaussian writes Gaussian inputs or imports QM outputs
       oracle-gf adds #GF_PED from Cartesian Hessian plus frozen #GIC
       oracle-rovib adds #ROTATIONAL, #VIBRATIONAL, #DELTABVIB, #CORIOLIS or #QCENT
       oracle-thermo adds #THERMO from BASIC/ROTATIONAL/VIBRATIONAL state
       oracle-morpheus adds #ISOTOPOLOGUES and #MORPHEUS
       oracle-vpt2-vci adds #VPT2_VCI from normalized anharmonic data
       oracle-dvr adds #DVR
```

The GUI should display and orchestrate this state, not own a parallel data
model.

The same file is also the standalone input contract. A tool may be run directly
from an already prepared `xyzin` file when the sections it needs are present.

## GIC State

`#GIC` stores the frozen coordinate contract, not only a list of Gaussian input
lines. In built files it includes:

- selected primitive coordinates;
- final frozen GICs, including linear-combination coefficients;
- point-group/symmetry-group metadata;
- total-symmetric irrep and active total-symmetric GIC list;
- reduction diagnostics;
- symmetrization diagnostics;
- Gaussian ReadGIC text generated from the frozen state.

Downstream modules should consume the frozen GICs and diagnostics directly and
use `[GAUSSIAN_GIC]` only when writing Gaussian inputs. Optimizers and
least-squares refinements should use `TOTAL_SYMMETRIC_GICS` for
symmetry-preserving active variables and reevaluate the B matrix from the frozen
definition at each geometry step.

## Compatibility

The ORACLE `xyzin` format is the historical source:

- ordinary XYZ block first;
- uppercase appended sections;
- section replacement must preserve all unrelated sections.

ORACLE keeps this behavior and upgrades legacy section schemas from
`oracle.xyzin.*` to `oracle.xyz.*` as modules migrate.
