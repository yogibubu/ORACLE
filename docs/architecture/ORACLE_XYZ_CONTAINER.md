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

## Gaussian Topology Overrides

When ORACLE-Babel imports a Gaussian `.log` or `.out`, `oracle-gaussian` may
write `#GAUSSIAN_TOPOLOGY` before the topology/synthon pass. Only two Gaussian
quantities are accepted as topology overrides:

- CM5 atomic charges;
- Mayer bond orders.

If CM5 charges are present, `#SYNTHONS` uses them as the atomic charge column
and records `CHARGE_SOURCE Gaussian CM5`. If CM5 charges are absent, synthons
use the ORACLE electronegativity charge model and record
`CHARGE_SOURCE Synthons electronegativity model`.

If Mayer bond orders are present, topology and synthons use them and record
`BOND_ORDER_SOURCE Gaussian Mayer`. If Mayer bond orders are absent, topology
uses the ORACLE continuous Pauling bond-order model and records
`BOND_ORDER_SOURCE Topology Pauling continuous model`. Gaussian total bond
orders are not a fallback source.

## Gaussian Rovibrational Promotion

Gaussian log/out text is converted once by `oracle-gaussian` before downstream
tools run. `oracle gaussian promote-rovib` promotes harmonic frequencies, IR
intensities, anharmonic chi matrices, rotational constants and vibrational
rotational corrections into shared ORACLE sections:

- `#VIBRATIONAL` stores frequencies, intensities and optional chi values.
- `#ROTATIONAL` stores rotational constants, point group, temperature and the
  DeltaBvib bridge values consumed by rovibrational utilities.
- `#DELTABVIB` stores the normalized DeltaBvib values and, when available, the
  Gaussian alpha rows used to compute them.

GF, Thermo, SEfit/MORPHEUS and anharmonic workflows consume these sections
rather than reparsing Gaussian output.

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
