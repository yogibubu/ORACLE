#!/usr/bin/env python3
"""Run the Fortran77 DVR kernel on a grid prepared by the Python workflow."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


def _normalized(name: str) -> str:
    return name.strip().lower().replace("/", "_").replace(" ", "_")


def read_grid_csv(path: Path) -> tuple[list[float], list[float]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Empty grid CSV: {path}")
        by_name = {_normalized(name): name for name in reader.fieldnames}
        q_key = by_name.get("q_au") or by_name.get("grid_au") or by_name.get("s_au")
        v_key = (
            by_name.get("potential_cm-1")
            or by_name.get("potential_cm1")
            or by_name.get("v_cm-1")
            or by_name.get("v_cm1")
        )
        if q_key is None or v_key is None:
            raise ValueError(
                f"Grid CSV must contain q_au/grid_au/s_au and potential_cm-1 columns: {path}"
            )
        q: list[float] = []
        v: list[float] = []
        for row in reader:
            q.append(float(row[q_key]))
            v.append(float(row[v_key]))
    if len(q) < 2:
        raise ValueError("Fortran DVR requires at least two grid points")
    return q, v


def read_grid2d_csv(
    path: Path,
    q1_key: str | None,
    q2_key: str | None,
    potential_key: str | None,
) -> tuple[list[float], list[float], list[float]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Empty 2D grid CSV: {path}")
        by_name = {_normalized(name): name for name in reader.fieldnames}
        q1_col = q1_key or by_name.get("q1") or by_name.get("q1_au")
        q2_col = q2_key or by_name.get("q2") or by_name.get("q2_au")
        v_col = (
            potential_key
            or by_name.get("potential_cm-1")
            or by_name.get("potential_cm1")
            or by_name.get("energy_cm-1")
            or by_name.get("energy_cm1")
        )
        if q1_col is None or q2_col is None or v_col is None:
            raise ValueError("2D CSV needs q1, q2 and potential/energy columns")
        rows = [(float(row[q1_col]), float(row[q2_col]), float(row[v_col])) for row in reader]
    q1 = sorted({row[0] for row in rows})
    q2 = sorted({row[1] for row in rows})
    values = {(a, b): v for a, b, v in rows}
    if len(values) != len(q1) * len(q2):
        raise ValueError("2D CSV is not a complete rectangular grid")
    vgrid = [values[(a, b)] for a in q1 for b in q2]
    return q1, q2, vgrid


def write_dvrin(
    path: Path,
    q: list[float],
    v: list[float],
    levels: int,
    boundary: str,
    mode: str,
    basis_size: int,
) -> None:
    ibound = 1 if boundary == "periodic" else 0
    with path.open("w") as handle:
        if mode == "gaussian":
            handle.write("2\n")
            handle.write(f"{len(q)} {basis_size} {levels}\n")
        else:
            handle.write("1\n")
            handle.write(f"{len(q)} {levels} {ibound}\n")
        for qi, vi in zip(q, v):
            handle.write(f"{qi:.16e} {vi:.16e}\n")


def write_dvrin_2d(
    path: Path,
    q1: list[float],
    q2: list[float],
    v: list[float],
    levels: int,
    boundary1: str,
    boundary2: str,
    g11: float,
    g22: float,
    g12: float,
) -> None:
    ib1 = 1 if boundary1 == "periodic" else 0
    ib2 = 1 if boundary2 == "periodic" else 0
    with path.open("w") as handle:
        handle.write("3\n")
        handle.write(f"{len(q1)} {len(q2)} {levels} {ib1} {ib2} {g11:.16e} {g22:.16e} {g12:.16e}\n")
        for value in q1:
            handle.write(f"{value:.16e}\n")
        for value in q2:
            handle.write(f"{value:.16e}\n")
        for value in v:
            handle.write(f"{value:.16e}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run ORACLE Fortran77 path DVR on a prepared grid CSV.")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--grid-csv", type=Path)
    source.add_argument("--grid2d-csv", type=Path)
    parser.add_argument("--exe", type=Path, default=Path(__file__).resolve().parents[3] / "bin" / "path_dvr.x")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", default="fortran_dvr")
    parser.add_argument("--levels", type=int, default=12)
    parser.add_argument("--boundary", choices=("periodic", "nonperiodic"), default="nonperiodic")
    parser.add_argument("--mode", choices=("grid", "gaussian"), default="grid")
    parser.add_argument("--basis-size", type=int, default=24)
    parser.add_argument("--q1-key")
    parser.add_argument("--q2-key")
    parser.add_argument("--potential-key")
    parser.add_argument("--boundary1", choices=("periodic", "nonperiodic"), default="nonperiodic")
    parser.add_argument("--boundary2", choices=("periodic", "nonperiodic"), default="nonperiodic")
    parser.add_argument("--g11", type=float, default=1.0)
    parser.add_argument("--g22", type=float, default=1.0)
    parser.add_argument("--g12", type=float, default=0.0)
    args = parser.parse_args()

    exe = args.exe.expanduser().resolve()
    if not exe.exists():
        raise FileNotFoundError(f"Fortran DVR executable not found: {exe}")

    run_dir = args.outdir.expanduser().resolve() / f"{args.prefix}_fortran_run"
    run_dir.mkdir(parents=True, exist_ok=True)
    for stale_name in ("dvrout", "dvr_levels.csv", "dvr_vectors.csv", "stdout.txt", "stderr.txt"):
        stale = run_dir / stale_name
        if stale.exists():
            stale.unlink()
    if args.grid2d_csv:
        q1, q2, v2d = read_grid2d_csv(
            args.grid2d_csv.expanduser().resolve(),
            args.q1_key,
            args.q2_key,
            args.potential_key,
        )
        write_dvrin_2d(
            run_dir / "dvrin",
            q1,
            q2,
            v2d,
            args.levels,
            args.boundary1,
            args.boundary2,
            args.g11,
            args.g22,
            args.g12,
        )
    else:
        q, v = read_grid_csv(args.grid_csv.expanduser().resolve())
        write_dvrin(run_dir / "dvrin", q, v, args.levels, args.boundary, args.mode, args.basis_size)

    result = subprocess.run(
        [str(exe)],
        cwd=run_dir,
        text=True,
        capture_output=True,
        check=False,
    )
    (run_dir / "stdout.txt").write_text(result.stdout)
    (run_dir / "stderr.txt").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"Fortran DVR failed with exit code {result.returncode}; see {run_dir}"
        )

    for name in ("dvrout", "dvr_levels.csv", "dvr_vectors.csv"):
        src = run_dir / name
        if src.exists():
            dst = args.outdir / f"{args.prefix}_{name}"
            shutil.copy2(src, dst)
            print(dst)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
