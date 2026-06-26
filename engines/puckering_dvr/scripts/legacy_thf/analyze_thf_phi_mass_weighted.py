#!/usr/bin/env python3
"""Analyze a THF pseudorotation scan using a mass-weighted path coordinate.

The Gaussian scan is assumed to contain one optimized structure for each fixed
pseudorotational phase phi.  Phi is only the label used to generate the path:
the one-dimensional Hamiltonian is built on the arc length in mass-weighted
Cartesian coordinates, so the reduced mass of the coordinate is one.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "gaussian_outputs"
FIGS = ROOT / "figs"

HARTREE_TO_CM = 219474.6313705
AMU_TO_ME = 1822.888486209
BOHR_TO_ANGSTROM = 0.529177210903
MASS_WEIGHTED_A_TO_AU = math.sqrt(AMU_TO_ME) / BOHR_TO_ANGSTROM

MASSES_AMU = {
    1: 1.00782503223,
    6: 12.0,
    8: 15.99491461957,
}

MW_REFERENCE = {
    "delta_e_01_cm-1": 21307.71 / 29979.2458,
    "delta_e_23_cm-1": 61205.69 / 29979.2458,
    "barrier_high_flexible_cm-1": 45.03,
    "barrier_low_flexible_cm-1": 15.61,
    "barrier_high_simple_cm-1": 43.1,
    "barrier_low_simple_cm-1": 18.2,
}


@dataclass
class OptPoint:
    atoms: np.ndarray
    coords_angstrom: np.ndarray
    energy_hartree: float


def parse_orientation_block(lines: list[str], start: int) -> tuple[np.ndarray, np.ndarray, int]:
    """Parse a Gaussian Input/Standard orientation block."""
    i = start + 1
    dash_count = 0
    while i < len(lines):
        if set(lines[i].strip()) == {"-"}:
            dash_count += 1
            if dash_count == 2:
                i += 1
                break
        i += 1

    atoms: list[int] = []
    coords: list[list[float]] = []
    while i < len(lines):
        line = lines[i]
        if set(line.strip()) == {"-"}:
            break
        fields = line.split()
        if len(fields) >= 6:
            atoms.append(int(fields[1]))
            coords.append([float(fields[3]), float(fields[4]), float(fields[5])])
        i += 1

    if not atoms:
        raise ValueError(f"Empty orientation block near line {start + 1}")
    return np.array(atoms, dtype=int), np.array(coords, dtype=float), i


def parse_gaussian_optimized_points(log_path: Path) -> list[OptPoint]:
    lines = log_path.read_text(errors="ignore").splitlines()
    points: list[OptPoint] = []
    last_atoms: np.ndarray | None = None
    last_coords: np.ndarray | None = None
    last_energy: float | None = None

    i = 0
    while i < len(lines):
        line = lines[i]
        if "SCF Done:" in line:
            match = re.search(r"=\s*([-+]?\d+\.\d+)", line)
            if match:
                last_energy = float(match.group(1))

        if "Input orientation:" in line or "Standard orientation:" in line:
            last_atoms, last_coords, i = parse_orientation_block(lines, i)

        if "Optimization completed." in line:
            if last_atoms is not None and last_coords is not None and last_energy is not None:
                points.append(
                    OptPoint(
                        atoms=last_atoms.copy(),
                        coords_angstrom=last_coords.copy(),
                        energy_hartree=last_energy,
                    )
                )
        i += 1

    return points


def masses_for_atoms(atoms: np.ndarray) -> np.ndarray:
    missing = sorted({int(atom) for atom in atoms if int(atom) not in MASSES_AMU})
    if missing:
        raise ValueError(f"No mass is defined for atomic numbers: {missing}")
    return np.array([MASSES_AMU[int(atom)] for atom in atoms], dtype=float)


def mass_weighted_align(
    mobile: np.ndarray, reference: np.ndarray, masses: np.ndarray
) -> np.ndarray:
    """Align mobile onto reference by a mass-weighted Kabsch rotation."""
    weights = masses / masses.sum()
    mobile_centroid = np.sum(mobile * weights[:, None], axis=0)
    reference_centroid = np.sum(reference * weights[:, None], axis=0)
    x = mobile - mobile_centroid
    y = reference - reference_centroid
    covariance = x.T @ (masses[:, None] * y)
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    if np.linalg.det(u @ vt) < 0.0:
        correction[-1, -1] = -1.0
    rotation = u @ correction @ vt
    return x @ rotation + reference_centroid


def mass_weighted_arc_lengths(
    geometries: list[np.ndarray], masses: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    s_sqrtamu_angstrom = np.zeros(len(geometries), dtype=float)
    step_lengths = np.zeros(len(geometries), dtype=float)
    for i in range(1, len(geometries)):
        aligned = mass_weighted_align(geometries[i], geometries[i - 1], masses)
        diff = aligned - geometries[i - 1]
        ds = math.sqrt(float(np.sum(masses[:, None] * diff * diff)))
        step_lengths[i] = ds
        s_sqrtamu_angstrom[i] = s_sqrtamu_angstrom[i - 1] + ds
    return s_sqrtamu_angstrom, step_lengths


def dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.atan2(y, x)


def ring_puckering_diagnostics(
    coords: np.ndarray, ring_indices: list[int]
) -> tuple[np.ndarray, float, float, float, float]:
    r = coords[ring_indices]
    torsions = np.array(
        [
            dihedral(r[0], r[1], r[2], r[3]),
            dihedral(r[1], r[2], r[3], r[4]),
            dihedral(r[2], r[3], r[4], r[0]),
            dihedral(r[3], r[4], r[0], r[1]),
            dihedral(r[4], r[0], r[1], r[2]),
        ],
        dtype=float,
    )
    t0, t1, t2, t3, t4 = torsions
    scale = math.sqrt(2.0 / 5.0)
    qc = scale * (
        t0
        - 0.8090169944 * t1
        + 0.3090169944 * t2
        + 0.3090169944 * t3
        - 0.8090169944 * t4
    )
    qs = scale * (
        0.5877852523 * t1
        - 0.9510565163 * t2
        + 0.9510565163 * t3
        - 0.5877852523 * t4
    )
    q = math.hypot(qc, qs)
    phi = math.degrees(math.atan2(qs, qc)) % 360.0
    return torsions, qc, qs, q, phi


def solve_periodic_hamiltonian(
    s_au: np.ndarray,
    potential_cm: np.ndarray,
    *,
    cells: int,
    n_grid: int,
    n_levels: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = s_au - s_au[0]
    period = float(x[-1])
    if period <= 0.0:
        raise ValueError("The mass-weighted path length is zero")
    if len(np.unique(x)) != len(x):
        raise ValueError("The mass-weighted coordinate contains duplicate points")

    interpolator = PchipInterpolator(x, potential_cm)
    length = period * cells
    grid = np.linspace(0.0, length, n_grid, endpoint=False)
    dx = grid[1] - grid[0]
    phase = np.mod(grid, period)
    potential_grid = interpolator(phase)

    kinetic_main = HARTREE_TO_CM / dx**2
    kinetic_off = -0.5 * HARTREE_TO_CM / dx**2
    main = kinetic_main + potential_grid
    off = kinetic_off * np.ones(n_grid - 1)
    hamiltonian = diags(
        [main, off, off, [kinetic_off], [kinetic_off]],
        [0, -1, 1, n_grid - 1, -(n_grid - 1)],
        format="csr",
    )

    values = eigsh(
        hamiltonian,
        k=min(n_levels, n_grid - 2),
        which="SA",
        return_eigenvectors=False,
        tol=1.0e-10,
    )
    values = np.sort(values)
    return grid, potential_grid, values


def write_profile(
    path: Path,
    points: list[OptPoint],
    target_phi: np.ndarray,
    s_sqrtamu_angstrom: np.ndarray,
    step_lengths: np.ndarray,
    s_au: np.ndarray,
    rel_cm: np.ndarray,
    ring_indices: list[int],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "point",
                "target_phi_deg",
                "s_sqrtamu_angstrom",
                "step_sqrtamu_angstrom",
                "s_au",
                "energy_hartree",
                "relative_energy_cm-1",
                "T0_deg",
                "T1_deg",
                "T2_deg",
                "T3_deg",
                "T4_deg",
                "QC_rad",
                "QS_rad",
                "q_rad",
                "diagnostic_phi_deg",
            ]
        )
        for i, point in enumerate(points):
            torsions, qc, qs, q, phi = ring_puckering_diagnostics(
                point.coords_angstrom, ring_indices
            )
            writer.writerow(
                [
                    i + 1,
                    target_phi[i],
                    s_sqrtamu_angstrom[i],
                    step_lengths[i],
                    s_au[i],
                    point.energy_hartree,
                    rel_cm[i],
                    *[math.degrees(value) for value in torsions],
                    qc,
                    qs,
                    q,
                    phi,
                ]
            )


def write_levels(path: Path, levels_by_cells: dict[int, np.ndarray]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cells", "level", "energy_cm-1", "relative_to_ground_cm-1"])
        for cells, levels in levels_by_cells.items():
            ground = float(levels[0])
            for i, level in enumerate(levels):
                writer.writerow([cells, i, level, level - ground])


def write_comparison(
    path: Path,
    barriers: dict[str, float],
    levels_by_cells: dict[int, np.ndarray],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "calculated_cm-1", "mw_reference_cm-1", "difference_cm-1"])
        writer.writerow(
            [
                "barrier_endpoint_average",
                barriers["endpoint_average_cm-1"],
                MW_REFERENCE["barrier_high_flexible_cm-1"],
                barriers["endpoint_average_cm-1"] - MW_REFERENCE["barrier_high_flexible_cm-1"],
            ]
        )
        writer.writerow(
            [
                "barrier_phi90",
                barriers["phi90_cm-1"],
                MW_REFERENCE["barrier_low_flexible_cm-1"],
                barriers["phi90_cm-1"] - MW_REFERENCE["barrier_low_flexible_cm-1"],
            ]
        )
        for cells, levels in levels_by_cells.items():
            if len(levels) >= 4:
                ground = levels[0]
                d01 = levels[1] - levels[0]
                d23 = levels[3] - levels[2]
                writer.writerow(
                    [
                        f"delta_e_01_cells_{cells}",
                        d01,
                        MW_REFERENCE["delta_e_01_cm-1"],
                        d01 - MW_REFERENCE["delta_e_01_cm-1"],
                    ]
                )
                writer.writerow(
                    [
                        f"delta_e_23_cells_{cells}",
                        d23,
                        MW_REFERENCE["delta_e_23_cm-1"],
                        d23 - MW_REFERENCE["delta_e_23_cm-1"],
                    ]
                )
                writer.writerow([f"zero_point_cells_{cells}", ground, "", ""])


def write_summary(
    path: Path,
    log_path: Path,
    points: list[OptPoint],
    target_phi: np.ndarray,
    s_sqrtamu_angstrom: np.ndarray,
    s_au: np.ndarray,
    barriers: dict[str, float],
    levels_by_cells: dict[int, np.ndarray],
) -> None:
    lines: list[str] = []
    lines.append("THF pseudorotation analysis with mass-weighted coordinate")
    lines.append(f"Gaussian log: {log_path}")
    lines.append(f"Optimized points: {len(points)}")
    lines.append(
        f"Target phi range: {target_phi[0]:.3f} to {target_phi[-1]:.3f} degrees"
    )
    lines.append(
        "Path length: "
        f"{s_sqrtamu_angstrom[-1]:.8f} sqrt(amu) Angstrom = {s_au[-1]:.8f} a.u."
    )
    lines.append("Reduced mass of coordinate s: 1 by construction")
    lines.append("")
    lines.append("Barriers relative to the computed minimum:")
    lines.append(f"  phi = 0 endpoint:   {barriers['phi0_cm-1']:.4f} cm^-1")
    lines.append(f"  phi = 90:           {barriers['phi90_cm-1']:.4f} cm^-1")
    lines.append(f"  phi = 180 endpoint: {barriers['phi180_cm-1']:.4f} cm^-1")
    lines.append(
        f"  endpoint average:   {barriers['endpoint_average_cm-1']:.4f} cm^-1"
    )
    lines.append("")
    lines.append("Microwave-paper reference values:")
    lines.append("  flexible model barriers: 45.03 and 15.61 cm^-1")
    lines.append("  measured splittings: 0.710749 and 2.041602 cm^-1")
    lines.append("")
    for cells, levels in levels_by_cells.items():
        lines.append(f"Periodic Hamiltonian using {cells} cell(s) of the 0-180 path:")
        ground = float(levels[0])
        for i, level in enumerate(levels[:8]):
            lines.append(f"  level {i:2d}: {level - ground:10.6f} cm^-1 above ground")
        if len(levels) >= 4:
            lines.append(f"  Delta E01: {levels[1] - levels[0]:.6f} cm^-1")
            lines.append(f"  Delta E23: {levels[3] - levels[2]:.6f} cm^-1")
        lines.append("")
    path.write_text("\n".join(lines) + "\n")


def plot_potential(
    path_pdf: Path,
    path_png: Path,
    target_phi: np.ndarray,
    s_sqrtamu_angstrom: np.ndarray,
    rel_cm: np.ndarray,
    barriers: dict[str, float],
) -> None:
    dense_phi = np.linspace(float(target_phi[0]), float(target_phi[-1]), 1200)
    phi_interp = PchipInterpolator(target_phi, rel_cm)
    dense_v_phi = phi_interp(dense_phi)

    dense_s = np.linspace(float(s_sqrtamu_angstrom[0]), float(s_sqrtamu_angstrom[-1]), 1200)
    s_interp = PchipInterpolator(s_sqrtamu_angstrom, rel_cm)
    dense_v_s = s_interp(dense_s)

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6), constrained_layout=True)
    axes[0].plot(dense_phi, dense_v_phi, color="#1f5a7a", lw=2.2)
    axes[0].axvline(90.0, color="#8f2d2d", lw=1.0, ls="--")
    axes[0].set_xlabel(r"$\phi$ / degrees")
    axes[0].set_ylabel(r"$V$ / cm$^{-1}$")
    axes[0].set_title("Barrier positions")
    axes[0].text(
        0.03,
        0.93,
        (
            f"endpoints {barriers['endpoint_average_cm-1']:.1f}\n"
            f"90 deg {barriers['phi90_cm-1']:.1f}"
        ),
        transform=axes[0].transAxes,
        va="top",
        ha="left",
        fontsize=9,
    )

    axes[1].plot(dense_s, dense_v_s, color="#3f6f3b", lw=2.2)
    axes[1].set_xlabel(r"$s$ / $\sqrt{\mathrm{amu}}$ Angstrom")
    axes[1].set_ylabel(r"$V$ / cm$^{-1}$")
    axes[1].set_title("Hamiltonian coordinate")
    axes[1].text(
        0.04,
        0.93,
        r"$\mu_s = 1$",
        transform=axes[1].transAxes,
        va="top",
        ha="left",
        fontsize=10,
    )

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(color="#d8d8d8", lw=0.5, alpha=0.7)

    fig.savefig(path_pdf)
    fig.savefig(path_png, dpi=240)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze a THF phi scan with a mass-weighted arc-length coordinate."
    )
    parser.add_argument(
        "--log",
        type=Path,
        default=OUT / "thf_phi_scan.log",
        help="Gaussian log containing constrained optimizations along phi.",
    )
    parser.add_argument(
        "--prefix",
        default="thf_phi_mass_weighted",
        help="Prefix for CSV, summary, and figure outputs.",
    )
    parser.add_argument(
        "--ring",
        default="1,2,3,4,5",
        help="One-based atom indices for O-C-C-C-C ring atoms.",
    )
    parser.add_argument("--phi-start", type=float, default=0.0)
    parser.add_argument("--phi-step", type=float, default=5.0)
    parser.add_argument("--grid", type=int, default=2400)
    parser.add_argument("--levels", type=int, default=16)
    parser.add_argument(
        "--cells",
        default="1,2",
        help="Comma-separated number of repeated 0-180 cells for the Hamiltonian.",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Parse the log and report the number of optimized points without writing outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    log_path = args.log.expanduser().resolve()
    if not log_path.exists():
        raise FileNotFoundError(
            f"{log_path} does not exist. Put the completed phi scan there or pass --log."
        )

    points = parse_gaussian_optimized_points(log_path)
    if len(points) < 5:
        raise RuntimeError(
            f"Only {len(points)} optimized points were found in {log_path}; "
            "a 0-180 phi scan at 5 degrees should contain 37 points."
        )

    atoms = points[0].atoms
    for point in points[1:]:
        if not np.array_equal(point.atoms, atoms):
            raise RuntimeError("Atom order changes between optimized structures.")

    ring_indices = [int(item.strip()) - 1 for item in args.ring.split(",")]
    if len(ring_indices) != 5 or min(ring_indices) < 0 or max(ring_indices) >= len(atoms):
        raise ValueError("--ring must contain five valid one-based atom indices")

    if args.check_only:
        print(f"optimized_points={len(points)}")
        print(f"first_energy_hartree={points[0].energy_hartree:.12f}")
        print(f"last_energy_hartree={points[-1].energy_hartree:.12f}")
        return

    masses = masses_for_atoms(atoms)
    geometries = [point.coords_angstrom for point in points]
    s_sqrtamu_angstrom, step_lengths = mass_weighted_arc_lengths(geometries, masses)
    s_au = s_sqrtamu_angstrom * MASS_WEIGHTED_A_TO_AU

    energies = np.array([point.energy_hartree for point in points], dtype=float)
    rel_cm = (energies - energies.min()) * HARTREE_TO_CM
    target_phi = args.phi_start + args.phi_step * np.arange(len(points), dtype=float)
    if target_phi[0] > 1.0e-8 or target_phi[-1] < 180.0 - 1.0e-8:
        raise RuntimeError(
            "The target phi labels must cover 0, 90, and 180 degrees. "
            f"With {len(points)} points, phi-start={args.phi_start} and "
            f"phi-step={args.phi_step} give {target_phi[0]:.3f} to "
            f"{target_phi[-1]:.3f} degrees."
        )

    phi_interp = PchipInterpolator(target_phi, rel_cm)
    phi0 = float(phi_interp(0.0))
    phi90 = float(phi_interp(90.0))
    phi180 = float(phi_interp(180.0))
    barriers = {
        "phi0_cm-1": phi0,
        "phi90_cm-1": phi90,
        "phi180_cm-1": phi180,
        "endpoint_average_cm-1": 0.5 * (phi0 + phi180),
    }

    cells_values = [int(item.strip()) for item in args.cells.split(",") if item.strip()]
    levels_by_cells: dict[int, np.ndarray] = {}
    for cells in cells_values:
        _, _, levels = solve_periodic_hamiltonian(
            s_au,
            rel_cm,
            cells=cells,
            n_grid=args.grid * cells,
            n_levels=args.levels,
        )
        levels_by_cells[cells] = levels

    OUT.mkdir(exist_ok=True)
    FIGS.mkdir(exist_ok=True)
    write_profile(
        OUT / f"{args.prefix}_profile.csv",
        points,
        target_phi,
        s_sqrtamu_angstrom,
        step_lengths,
        s_au,
        rel_cm,
        ring_indices,
    )
    write_levels(OUT / f"{args.prefix}_levels.csv", levels_by_cells)
    write_comparison(OUT / f"{args.prefix}_mw_comparison.csv", barriers, levels_by_cells)
    write_summary(
        OUT / f"{args.prefix}_summary.txt",
        log_path,
        points,
        target_phi,
        s_sqrtamu_angstrom,
        s_au,
        barriers,
        levels_by_cells,
    )
    plot_potential(
        FIGS / f"{args.prefix}_potential.pdf",
        FIGS / f"{args.prefix}_potential.png",
        target_phi,
        s_sqrtamu_angstrom,
        rel_cm,
        barriers,
    )

    print(f"Wrote outputs with prefix {args.prefix}")
    print(f"Reduced mass of s: 1")
    print(
        "Path length: "
        f"{s_sqrtamu_angstrom[-1]:.8f} sqrt(amu) Angstrom = {s_au[-1]:.8f} a.u."
    )


if __name__ == "__main__":
    main()
