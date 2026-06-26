Fortran77 Path DVR
==================

`path_dvr.f` is the active Fortran77 DVR kernel for ORACLE.

It deliberately does not parse Gaussian output and does not generate paths.
Python prepares equally spaced grids from supported Gaussian scan/path output,
then this solver reads only `dvrin`.

The diagonalizer is `DVRHQRII` in `dvr_hqrii.f`, a local copy of the GICForge
`HQRII1` Householder/QR inverse-iteration routine. Jacobi diagonalization is
not used in this backend.

Input format:

```text
MODE
...
```

`MODE=1`: one-dimensional grid DVR.

```text
NPTS NLEVELS IBOUND
q_1  V_1
...
q_N  V_N
```

- `q`: mass-weighted coordinate in atomic units.
- `V`: relative potential in cm-1.
- `IBOUND=0`: non-periodic sinc DVR.
- `IBOUND=1`: periodic Fourier DVR.

`MODE=2`: one-dimensional distributed Gaussian basis.

```text
NGRID NBASIS NLEVELS
q_1   V_1
...
q_N   V_N
```

`MODE=3`: two-dimensional product grid DVR.

```text
N1 N2 NLEVELS IBOUND1 IBOUND2 G11 G22 G12
q1_1 ... q1_N1
q2_1 ... q2_N2
V_11 ... V_N1N2
```

Build:

```bash
cd engines/fortran/dvr
./compile_MAC
```

Outputs:

- `dvrout`
- `dvr_levels.csv`
- `dvr_vectors.csv`
