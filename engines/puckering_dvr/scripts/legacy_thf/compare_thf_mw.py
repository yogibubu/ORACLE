#!/usr/bin/env python3
"""Compare the THF scan model with the microwave analysis of Meyer et al."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "gaussian_outputs"
FIGS = ROOT / "figs"


MW_BARRIERS = {
    "high_barrier": 45.03,
    "low_barrier": 15.61,
}

MW_SPLITTINGS = {
    "DeltaE01": 0.710749,
    "DeltaE23": 2.041602,
}

MW_ROT_MHZ = {
    0: (7096.770, 6976.371, 4008.219),
    1: (7092.859, 6982.939, 4008.377),
    2: (7131.432, 6920.466, 3999.042),
    3: (7127.863, 6927.011, 3997.073),
}

MW_FIR_CM = {
    (5, 3): 24.5,
    (7, 5): 23.0,
    (8, 6): 19.9,
    (9, 7): 28.0,
    (10, 8): 27.3,
    (11, 9): 34.9,
    (12, 10): 33.6,
    (13, 12): 41.6,
    (15, 14): 48.2,
    (17, 16): 54.4,
}


def read_levels() -> list[dict[str, float | int | str]]:
    rows = []
    with (OUT / "thf_pseudorotation_levels.csv").open() as handle:
        for row in csv.DictReader(handle):
            if row["model"] != "four_well_periodic_symmetric":
                continue
            rows.append(
                {
                    "state": int(row["state"]),
                    "E0": float(row["energy_above_ground_cm-1"]),
                    "Emin": float(row["energy_above_minimum_cm-1"]),
                    "A": float(row["A_eff_GHz"]) * 1000.0,
                    "B": float(row["B_eff_GHz"]) * 1000.0,
                    "C": float(row["C_eff_GHz"]) * 1000.0,
                }
            )
    rows.sort(key=lambda item: int(item["state"]))
    return rows


def read_barrier() -> float:
    summary = (OUT / "thf_pseudorotation_summary.txt").read_text()
    for line in summary.splitlines():
        if line.startswith("Barrier estimate:"):
            return float(line.split()[2])
    raise RuntimeError("Barrier estimate not found in summary")


def write_csv(path: Path, rows: list[list[object]], header: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    levels = read_levels()
    calc_barrier = read_barrier()

    calc_splittings = {
        "DeltaE01": levels[1]["E0"] - levels[0]["E0"],
        "DeltaE23": levels[3]["E0"] - levels[2]["E0"],
    }

    barrier_rows = [
        [
            "high_barrier",
            MW_BARRIERS["high_barrier"],
            calc_barrier,
            calc_barrier - MW_BARRIERS["high_barrier"],
            100.0 * (calc_barrier - MW_BARRIERS["high_barrier"]) / MW_BARRIERS["high_barrier"],
            "Calculated branch corresponds to the high-barrier side of the MW flexible model.",
        ],
        [
            "low_barrier",
            MW_BARRIERS["low_barrier"],
            "",
            "",
            "",
            "Not resolved by the present one-branch symmetric scan.",
        ],
    ]
    write_csv(
        OUT / "thf_mw_barrier_comparison.csv",
        barrier_rows,
        ["quantity", "MW_cm-1", "calculated_cm-1", "calc_minus_MW_cm-1", "percent_error", "comment"],
    )

    splitting_rows = []
    for key, observed in MW_SPLITTINGS.items():
        calculated = calc_splittings[key]
        splitting_rows.append(
            [
                key,
                observed,
                calculated,
                calculated - observed,
                100.0 * (calculated - observed) / observed,
            ]
        )
    write_csv(
        OUT / "thf_mw_splitting_comparison.csv",
        splitting_rows,
        ["quantity", "MW_cm-1", "calculated_cm-1", "calc_minus_MW_cm-1", "percent_error"],
    )

    rotational_rows = []
    for state in range(4):
        calc = levels[state]
        obs = MW_ROT_MHZ[state]
        for label, observed, calculated in zip(["A", "B", "C"], obs, [calc["A"], calc["B"], calc["C"]]):
            rotational_rows.append(
                [
                    state,
                    label,
                    observed,
                    calculated,
                    calculated - observed,
                    100.0 * (calculated - observed) / observed,
                ]
            )
    write_csv(
        OUT / "thf_mw_rotational_comparison.csv",
        rotational_rows,
        ["state", "constant", "MW_MHz", "calculated_MHz", "calc_minus_MW_MHz", "percent_error"],
    )

    fir_rows = []
    for (upper, lower), observed in MW_FIR_CM.items():
        if upper >= len(levels) or lower >= len(levels):
            calculated = ""
            diff = ""
            percent = ""
        else:
            calculated = levels[upper]["E0"] - levels[lower]["E0"]
            diff = calculated - observed
            percent = 100.0 * diff / observed
        fir_rows.append([f"{upper}<-{lower}", observed, calculated, diff, percent])
    write_csv(
        OUT / "thf_mw_far_ir_comparison.csv",
        fir_rows,
        ["transition", "MW_cm-1", "calculated_cm-1", "calc_minus_MW_cm-1", "percent_error"],
    )

    with (OUT / "thf_mw_comparison_summary.txt").open("w") as handle:
        handle.write("Comparison with tetrahydrofuran microwave pseudorotation data\n")
        handle.write("=============================================================\n\n")
        handle.write("Reference data: Meyer, Lopez, Alonso, Melandri, Favero, and Caminati, J. Chem. Phys. 111, 7871 (1999).\n")
        handle.write("Calculated data: B3LYP-D3BJ/def2-TZVP relaxed GIC scan plus a one-dimensional Hamiltonian with B_phi = 3.19 cm^-1.\n\n")
        handle.write(f"High barrier: calculated {calc_barrier:.2f} cm^-1 vs MW flexible-model {MW_BARRIERS['high_barrier']:.2f} cm^-1.\n")
        handle.write("Low barrier: not available from the present one-branch symmetric scan; MW flexible-model value is 15.61 cm^-1.\n")
        handle.write(
            f"DeltaE01: calculated {calc_splittings['DeltaE01']:.3f} cm^-1 vs MW {MW_SPLITTINGS['DeltaE01']:.3f} cm^-1.\n"
        )
        handle.write(
            f"DeltaE23: calculated {calc_splittings['DeltaE23']:.3f} cm^-1 vs MW {MW_SPLITTINGS['DeltaE23']:.3f} cm^-1.\n\n"
        )

        mue_rot = sum(abs(float(row[4])) for row in rotational_rows) / len(rotational_rows)
        handle.write(f"Mean unsigned error on A/B/C for states 0-3: {mue_rot:.2f} MHz.\n")
        handle.write("Interpretation: the local scan captures the high-barrier scale and the first tunneling splitting, but the missing low-barrier branch limits higher splittings, far-infrared intervals, and state-dependent rotational constants.\n")

    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.4))

    axes[0].bar(["MW high", "calc"], [MW_BARRIERS["high_barrier"], calc_barrier], color=["#8ca8c4", "#0a6470"])
    axes[0].axhline(MW_BARRIERS["low_barrier"], color="#d45a2c", lw=1.4, ls="--", label="MW low")
    axes[0].set_ylabel(r"Barrier / cm$^{-1}$")
    axes[0].set_title("Barrier")
    axes[0].legend(frameon=False, fontsize=8)

    split_labels = list(MW_SPLITTINGS.keys())
    x = range(len(split_labels))
    axes[1].bar([i - 0.18 for i in x], [MW_SPLITTINGS[k] for k in split_labels], width=0.34, color="#8ca8c4", label="MW")
    axes[1].bar([i + 0.18 for i in x], [calc_splittings[k] for k in split_labels], width=0.34, color="#0a6470", label="calc")
    axes[1].set_xticks(list(x), split_labels)
    axes[1].set_ylabel(r"Splitting / cm$^{-1}$")
    axes[1].set_title("Tunnel Splittings")
    axes[1].legend(frameon=False, fontsize=8)

    states = list(range(4))
    obs_b = [MW_ROT_MHZ[s][1] for s in states]
    calc_b = [levels[s]["B"] for s in states]
    axes[2].plot(states, obs_b, "o-", color="#8ca8c4", label="MW B")
    axes[2].plot(states, calc_b, "s-", color="#0a6470", label="calc B")
    axes[2].set_xlabel("Pseudorotational state")
    axes[2].set_ylabel("B / MHz")
    axes[2].set_title("Rotational Constant")
    axes[2].legend(frameon=False, fontsize=8)

    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGS / "thf_mw_comparison.pdf")
    fig.savefig(FIGS / "thf_mw_comparison.png", dpi=240)

    print(OUT / "thf_mw_comparison_summary.txt")
    print(OUT / "thf_mw_barrier_comparison.csv")
    print(OUT / "thf_mw_splitting_comparison.csv")
    print(OUT / "thf_mw_rotational_comparison.csv")
    print(OUT / "thf_mw_far_ir_comparison.csv")
    print(FIGS / "thf_mw_comparison.pdf")


if __name__ == "__main__":
    main()
