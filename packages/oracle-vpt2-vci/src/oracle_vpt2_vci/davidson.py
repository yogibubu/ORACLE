from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


MatVec = Callable[[np.ndarray], np.ndarray]


@dataclass(frozen=True)
class DavidsonResult:
    """Lowest eigenpairs from an independent symmetric Davidson iteration."""

    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    residual_norms: np.ndarray
    iterations: int
    converged: bool


def _orthonormalize(candidates: np.ndarray, reference: np.ndarray | None = None, tol: float = 1.0e-12) -> np.ndarray:
    vecs: list[np.ndarray] = []
    refs = [] if reference is None else [reference[:, i] for i in range(reference.shape[1])]
    for col in range(candidates.shape[1]):
        v = candidates[:, col].astype(float, copy=True)
        for q in refs:
            v -= q * float(np.dot(q, v))
        for q in vecs:
            v -= q * float(np.dot(q, v))
        norm = float(np.linalg.norm(v))
        if norm > tol:
            vecs.append(v / norm)
    if not vecs:
        return np.zeros((candidates.shape[0], 0), dtype=float)
    return np.column_stack(vecs)


def davidson_lowest(
    matvec: MatVec,
    diagonal: np.ndarray,
    *,
    n_roots: int,
    max_subspace: int = 80,
    max_iter: int = 200,
    convergence: float = 1.0e-8,
) -> DavidsonResult:
    """Compute the lowest symmetric eigenpairs using a standalone Davidson solver.

    The solver knows nothing about Gaussian, GDV, VCI or force-field storage. It
    only requires a matrix-vector product and an approximate diagonal.
    """
    diag = np.asarray(diagonal, dtype=float)
    n = diag.size
    if n_roots < 1 or n_roots > n:
        raise ValueError("n_roots must be between 1 and the matrix size")
    if max_subspace < n_roots:
        raise ValueError("max_subspace must be >= n_roots")
    if convergence <= 0.0:
        raise ValueError("convergence must be positive")

    guess_indices = np.argsort(diag)[:n_roots]
    guesses = np.zeros((n, n_roots), dtype=float)
    for col, idx in enumerate(guess_indices):
        guesses[idx, col] = 1.0
    v = _orthonormalize(guesses)
    av = np.column_stack([matvec(v[:, i]) for i in range(v.shape[1])])

    theta = np.zeros(n_roots, dtype=float)
    ritz = np.zeros((n, n_roots), dtype=float)
    residual_norms = np.full(n_roots, np.inf, dtype=float)

    for iteration in range(1, max_iter + 1):
        projected = v.T @ av
        evals, evecs = np.linalg.eigh((projected + projected.T) * 0.5)
        theta = evals[:n_roots]
        coeff = evecs[:, :n_roots]
        ritz = v @ coeff
        aritz = av @ coeff
        residuals = aritz - ritz * theta[None, :]
        residual_norms = np.linalg.norm(residuals, axis=0)
        if bool(np.all(residual_norms <= convergence)):
            return DavidsonResult(theta, ritz, residual_norms, iteration, True)

        corrections: list[np.ndarray] = []
        for root in range(n_roots):
            if residual_norms[root] <= convergence:
                continue
            denom = theta[root] - diag
            small = np.abs(denom) < 1.0e-10
            denom = denom.copy()
            denom[small] = np.where(denom[small] < 0.0, -1.0e-10, 1.0e-10)
            corrections.append(residuals[:, root] / denom)
        if not corrections:
            break
        add = _orthonormalize(np.column_stack(corrections), v)
        if add.shape[1] == 0:
            break

        if v.shape[1] + add.shape[1] > max_subspace:
            v = _orthonormalize(ritz)
            av = np.column_stack([matvec(v[:, i]) for i in range(v.shape[1])])
            add = _orthonormalize(np.column_stack(corrections), v)
            if add.shape[1] == 0:
                continue

        v = np.column_stack((v, add))
        av_add = np.column_stack([matvec(add[:, i]) for i in range(add.shape[1])])
        av = np.column_stack((av, av_add))

    return DavidsonResult(theta, ritz, residual_norms, max_iter, False)
