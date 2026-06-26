#!/usr/bin/env python3
"""Analyze the tetrahydrofuran GIC pseudorotation scan.

The script extracts the relaxed scan from a Gaussian log, builds a
symmetrized one-dimensional potential along Qtan, and diagonalizes a
four-well periodic Hamiltonian using the microwave-paper kinetic
coefficient B_phi = 3.19 cm^-1.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh


ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "gaussian_outputs" / "scan.log"
OUT = ROOT / "gaussian_outputs"
FIGS = ROOT / "figs"

HARTREE_TO_CM = 219474.63137
B_PHI = 3.19
STEP_RAD = 0.087266463


def parse_hf_arrays(text: str) -> list[list[float]]:
    arrays: list[list[float]] = []
    pattern = re.compile(r"HF=([^\\]+(?:\\n?[^\\]+)*?)\\RMSD=")
    for match in pattern.finditer(text):
        raw = match.group(1).replace("\n", "").replace(" ", "")
        values = []
        for item in raw.split(","):
            if item:
                values.append(float(item))
        arrays.append(values)
    return arrays


def parse_qtan_blocks(text: str) -> list[list[float]]:
    blocks: list[list[float]] = []
    for chunk in text.split("Summary of Optimized Potential Surface Scan")[1:]:
        qtan: list[float] = []
        for line in chunk.splitlines():
            if line.strip().startswith("Qtan"):
                qtan.extend(float(x) for x in re.findall(r"[-+]?\d+\.\d+", line))
            elif qtan and "Largest change from initial coordinates" in line:
                break
        if qtan:
            blocks.append(qtan)
    return blocks


def parse_final_rotational_constants(text: str) -> dict[str, list[list[float]]]:
    data = {"positive": [], "negative": []}
    job: str | None = None
    waiting = False

    for line in text.splitlines():
        if "positive t" in line:
            job = "positive"
        elif "negative t" in line:
            job = "negative"

        if "Optimization completed." in line and job:
            waiting = True
            continue

        if waiting and "Rotational constants (GHZ):" in line:
            values = [float(x) for x in re.findall(r"[-+]?\d+\.\d+", line)]
            if len(values) >= 3 and job:
                data[job].append(values[-3:])
            waiting = False

    return data


def finite_difference_periodic(
    x_abs: np.ndarray,
    v_abs: np.ndarray,
    a_ghz: np.ndarray,
    b_ghz: np.ndarray,
    c_ghz: np.ndarray,
    n_grid: int = 2000,
    n_levels: int = 24,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Four equivalent wells obtained by periodic repetition of one branch."""
    xmax = float(x_abs[-1])
    period = 2.0 * xmax
    length = 4.0 * period
    dx = length / n_grid
    grid = np.linspace(-length / 2.0, length / 2.0 - dx, n_grid)

    v_interp = PchipInterpolator(x_abs, v_abs)
    a_interp = PchipInterpolator(x_abs, a_ghz)
    b_interp = PchipInterpolator(x_abs, b_ghz)
    c_interp = PchipInterpolator(x_abs, c_ghz)

    reduced = ((grid + xmax) % period) - xmax
    q_abs_grid = np.abs(reduced)
    potential = v_interp(q_abs_grid)

    main = 2.0 * B_PHI / dx**2 + potential
    off = -B_PHI / dx**2 * np.ones(n_grid - 1)
    hamiltonian = diags(
        [main, off, off, [-B_PHI / dx**2], [-B_PHI / dx**2]],
        [0, -1, 1, n_grid - 1, -(n_grid - 1)],
        format="csr",
    )

    values, vectors = eigsh(hamiltonian, k=n_levels, which="SA", tol=1.0e-10)
    order = np.argsort(values)
    values = values[order]
    vectors = vectors[:, order]

    # eigsh returns normalized Euclidean vectors; expectation values are sums
    # over grid points because the same normalization cancels in numerator and
    # denominator.
    weights = vectors**2
    a_eff = weights.T @ a_interp(q_abs_grid)
    b_eff = weights.T @ b_interp(q_abs_grid)
    c_eff = weights.T @ c_interp(q_abs_grid)
    abc_eff = np.vstack([a_eff, b_eff, c_eff]).T

    return grid, potential, values, vectors, abc_eff


def local_levels(x_abs: np.ndarray, v_abs: np.ndarray, n_grid: int = 1600) -> np.ndarray:
    xmax = float(x_abs[-1])
    interp = PchipInterpolator(x_abs, v_abs)
    grid = np.linspace(-xmax, xmax, n_grid)
    dx = grid[1] - grid[0]
    potential = interp(np.abs(grid))[1:-1]
    n = potential.size
    main = 2.0 * B_PHI / dx**2 + potential
    off = -B_PHI / dx**2 * np.ones(n - 1)
    hamiltonian = diags([main, off, off], [0, -1, 1], format="csr")
    return np.sort(eigsh(hamiltonian, k=8, which="SA", return_eigenvectors=False, tol=1.0e-10))


def main() -> None:
    text = LOG.read_text(errors="ignore")
    hf_arrays = parse_hf_arrays(text)
    qtan_blocks = parse_qtan_blocks(text)
    rot_blocks = parse_final_rotational_constants(text)

    if not hf_arrays:
        raise RuntimeError(f"No Gaussian HF archive arrays found in {LOG}")

    branch_energy = max(hf_arrays, key=len)
    n_points = len(branch_energy)

    matching_qtan = [block for block in qtan_blocks if len(block) == n_points]
    if matching_qtan:
        qtan = np.abs(np.array(matching_qtan[0], dtype=float))
    else:
        qtan = np.arange(n_points, dtype=float) * STEP_RAD

    energies = np.array(branch_energy, dtype=float)
    rel_cm = (energies - energies.min()) * HARTREE_TO_CM

    rot_candidates = [vals for vals in rot_blocks.values() if len(vals) == n_points]
    if rot_candidates:
        abc = np.array(rot_candidates[0], dtype=float)
    else:
        abc = np.full((n_points, 3), np.nan)

    order = np.argsort(qtan)
    x_abs = qtan[order]
    v_abs = rel_cm[order]
    abc = abc[order]

    # The endpoint is a plateau; the numerical maximum of the interpolant is
    # the barrier estimate for this symmetrized one-dimensional model.
    dense_x = np.linspace(0.0, x_abs[-1], 1201)
    interp = PchipInterpolator(x_abs, v_abs)
    dense_v = interp(dense_x)
    imax = int(np.argmax(dense_v))
    barrier_x = dense_x[imax]
    barrier_v = dense_v[imax]

    fit_matrix = np.vstack([x_abs**2, x_abs**4, x_abs**6, x_abs**8]).T
    poly_coeffs = np.linalg.lstsq(fit_matrix[1:], v_abs[1:], rcond=None)[0]
    harmonic_spacing = 2.0 * math.sqrt(poly_coeffs[0] * B_PHI)

    grid, potential, periodic_levels, _, abc_eff = finite_difference_periodic(
        x_abs, v_abs, abc[:, 0], abc[:, 1], abc[:, 2]
    )
    local = local_levels(x_abs, v_abs)

    profile_csv = OUT / "thf_pseudorotation_profile.csv"
    with profile_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "point",
                "abs_Qtan_rad",
                "abs_Qtan_deg",
                "energy_hartree",
                "relative_energy_cm-1",
                "A_GHz",
                "B_GHz",
                "C_GHz",
            ]
        )
        for i, (x, e, v, row) in enumerate(zip(x_abs, energies[order], v_abs, abc), start=1):
            writer.writerow([i, x, math.degrees(x), e, v, row[0], row[1], row[2]])

    levels_csv = OUT / "thf_pseudorotation_levels.csv"
    with levels_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "model",
                "state",
                "energy_above_minimum_cm-1",
                "energy_above_ground_cm-1",
                "A_eff_GHz",
                "B_eff_GHz",
                "C_eff_GHz",
            ]
        )
        for i, (energy, row) in enumerate(zip(periodic_levels, abc_eff)):
            writer.writerow(
                [
                    "four_well_periodic_symmetric",
                    i,
                    energy,
                    energy - periodic_levels[0],
                    row[0],
                    row[1],
                    row[2],
                ]
            )
        for i, energy in enumerate(local):
            writer.writerow(["single_well_dirichlet", i, energy, energy - local[0], "", "", ""])

    transitions_csv = OUT / "thf_pseudorotation_transitions.csv"
    with transitions_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["model", "transition", "delta_cm-1"])
        for i in range(len(periodic_levels) - 1):
            writer.writerow(
                [
                    "four_well_periodic_symmetric",
                    f"{i + 1}<-{i}",
                    periodic_levels[i + 1] - periodic_levels[i],
                ]
            )
        for i in range(len(local) - 1):
            writer.writerow(["single_well_dirichlet", f"{i + 1}<-{i}", local[i + 1] - local[i]])

    summary = OUT / "thf_pseudorotation_summary.txt"
    with summary.open("w") as handle:
        handle.write("Tetrahydrofuran pseudorotation scan analysis\n")
        handle.write("============================================\n\n")
        handle.write(f"Input log: {LOG}\n")
        handle.write(f"Unique scan points used: {n_points}\n")
        handle.write(f"Kinetic coefficient B_phi: {B_PHI:.4f} cm^-1\n")
        handle.write(f"Barrier estimate: {barrier_v:.3f} cm^-1 at |Qtan| = {math.degrees(barrier_x):.2f} deg\n")
        handle.write(f"Harmonic spacing from even polynomial curvature: {harmonic_spacing:.3f} cm^-1\n")
        handle.write("\nEven polynomial fit V = c2*x^2 + c4*x^4 + c6*x^6 + c8*x^8, x in radians:\n")
        for name, value in zip(["c2", "c4", "c6", "c8"], poly_coeffs):
            handle.write(f"  {name} = {value:.10f} cm^-1\n")
        handle.write("\nFour-well periodic symmetric levels, relative to ground:\n")
        for i, energy in enumerate(periodic_levels[:12]):
            handle.write(f"  {i:2d}  {energy - periodic_levels[0]:10.4f} cm^-1\n")
        handle.write("\nSingle-well local levels, relative to ground:\n")
        for i, energy in enumerate(local[:8]):
            handle.write(f"  {i:2d}  {energy - local[0]:10.4f} cm^-1\n")

    full_x_deg = np.linspace(-90.0, 90.0, 2200)
    distance_from_minimum_deg = 90.0 - np.abs(full_x_deg)
    qtan_from_folded = np.radians(distance_from_minimum_deg / 90.0 * math.degrees(x_abs[-1]))
    full_v = interp(qtan_from_folded)

    extended_csv = OUT / "thf_pseudorotation_profile_minus90_90.csv"
    with extended_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["symmetry_coordinate_deg", "relative_energy_cm-1"])
        for x, v in zip(full_x_deg, full_v):
            writer.writerow([x, v])

    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    ax.plot(full_x_deg, full_v, color="#0a6470", lw=2.3, label="symmetry-extended interpolation")
    ax.axhline(barrier_v, color="#6f7780", lw=0.8, ls="--")
    ax.text(
        12.0,
        barrier_v + 1.4,
        f"barrier {barrier_v:.1f} cm$^{{-1}}$",
        fontsize=9,
        color="#31373d",
    )

    for energy in periodic_levels[:8]:
        ax.hlines(
            energy,
            -90.0,
            90.0,
            color="#2a77be",
            lw=0.9,
            alpha=0.70,
        )

    level_labels = [
        ("0-3", periodic_levels[3]),
        ("4", periodic_levels[4]),
        ("5,6", periodic_levels[5]),
        ("7", periodic_levels[7]),
    ]
    for label, energy in level_labels:
        ax.text(
            85.0,
            energy,
            label,
            va="center",
            fontsize=8,
            color="#2a77be",
        )

    for marker in (-90.0, 0.0, 90.0):
        ax.axvline(marker, color="#0a6470", lw=0.55, alpha=0.18)

    ax.set_xlabel(r"Symmetry coordinate / degrees")
    ax.set_ylabel(r"Energy / cm$^{-1}$")
    ax.set_title("Tetrahydrofuran relaxed pseudorotation scan", fontsize=14)
    ax.set_xlim(-90.0, 90.0)
    ax.set_ylim(-2.0, max(62.0, barrier_v + 7.0))
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    fig.savefig(FIGS / "thf_pseudorotation_potential.pdf")
    fig.savefig(FIGS / "thf_pseudorotation_potential.png", dpi=240)

    print(summary)
    print(profile_csv)
    print(levels_csv)
    print(transitions_csv)
    print(extended_csv)
    print(FIGS / "thf_pseudorotation_potential.pdf")


if __name__ == "__main__":
    main()
