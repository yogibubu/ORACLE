# ORACLE Shared Diagonalization

ORACLE keeps dense Hermitian diagonalization behind one shared API:

```python
from oracle_core import diagonalize_hermitian, eigh_arrays, eigvalsh_array
```

Scientific tools should call this API instead of `numpy.linalg.eigh` when the
matrix is part of the main numerical workload. This currently covers:

- Wilson GF and symmetry-block GF diagonalizations.
- VCI dense blocks and Davidson projected subspaces.
- Vibro-rotational Hessian normal-mode extraction.
- The vendored WMS-Rot rotational Hamiltonian engine.

Small geometric diagonalizations, SVDs and inertia tensor utilities can stay
local when GPU transfer overhead would dominate the calculation.

## Backend Selection

The default backend policy is `auto`:

1. Use a GPU backend for sufficiently large matrices when one is already
   available in the active environment.
2. Otherwise use SciPy LAPACK.
3. Fall back to NumPy LAPACK when SciPy is not available.

Supported backend names are:

- `auto`
- `cpu`
- `gpu`
- `scipy`
- `numpy`
- `cupy`
- `torch-cuda`
- `torch-mps`

Environment controls:

```bash
export ORACLE_DIAGONALIZER_BACKEND=auto
export ORACLE_DIAGONALIZER_GPU_MIN_SIZE=128
export ORACLE_DIAGONALIZER_STRICT_GPU=0
```

`ORACLE_DIAGONALIZER_STRICT_GPU=1` makes GPU failures fatal. With the default
value, a failed GPU call falls back to SciPy CPU so production workflows keep
running.

`oracle-set` does not install CuPy or PyTorch automatically. GPU packages are
hardware- and driver-specific, so they must be installed explicitly in the same
environment when acceleration is wanted. On NVIDIA/CUDA systems install the
matching CuPy or PyTorch CUDA build. On Apple Silicon, PyTorch with MPS support
can be used when the operation is supported by the local PyTorch/MPS stack.

## Pandas And SymPy

Pandas and SymPy are runtime dependencies because the imported first-party
WMS-Rot engine uses them:

- Pandas stores, filters and exports rotational line-list tables.
- SymPy provides symbolic/angular-momentum algebra helpers used by the
  rotational Hamiltonian code.

They are not mathematically required by GF, VCI or the shared diagonalizer. If
we later replace the WMS-Rot table layer with native ORACLE arrays/records,
Pandas can become optional for non-rotational workflows. For now it stays in
`oracle-set` to preserve the validated WMS-Rot behavior.

## Extending The Rotational Hamiltonian

Additional rotational Hamiltonian terms can be added without changing the GUI
contract. The right extension point is the local WMS-Rot engine and the shared
`#ROTATIONAL`/WMS-Rot input sections:

1. Define the physical term, units and parameter names.
2. Add the operator/matrix element implementation in the WMS-Rot engine.
3. Add parser/export support in the ORACLE rotational adapter.
4. Store the normalized parameters in the xyzin rotational section.
5. Add golden tests against known spectra or the original WMS-Rot web workflow.

Ordinary centrifugal-distortion-like terms fit this model directly. Hindered
internal rotors are possible too, but they require an explicit rotor basis and
coupling model; they should be treated as a larger Hamiltonian extension rather
than as only one scalar correction.
