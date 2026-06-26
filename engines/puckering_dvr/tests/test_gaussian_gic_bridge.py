from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


def load_mw_path_dvr():
    script = Path(__file__).resolve().parents[1] / "scripts" / "mw_path_dvr.py"
    spec = importlib.util.spec_from_file_location("mw_path_dvr", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def orientation_block(coords):
    lines = [
        " Input orientation:",
        " ---------------------------------------------------------------------",
        " Center     Atomic      Atomic             Coordinates (Angstroms)",
        " Number     Number       Type             X           Y           Z",
        " ---------------------------------------------------------------------",
    ]
    for i, (x, y, z) in enumerate(coords, start=1):
        lines.append(f" {i:6d} {6:10d} {0:11d} {x:15.6f} {y:11.6f} {z:11.6f}")
    lines.append(" ---------------------------------------------------------------------")
    return "\n".join(lines)


def test_gaussian_gic_values_are_profiled_and_mapped_to_cp(tmp_path):
    dvr = load_mw_path_dvr()
    assert dvr.canonical_ring_indices([2, 3, 4, 0, 1]) == [0, 1, 2, 3, 4]
    assert dvr.canonical_ring_indices([3, 2, 1, 0, 4]) == [0, 1, 2, 3, 4]

    theta = np.linspace(0.0, 2.0 * np.pi, 5, endpoint=False)
    base = np.column_stack([np.cos(theta), np.sin(theta), np.zeros(5)])
    coords1 = base.copy()
    coords1[:, 2] = [0.15, -0.10, 0.12, -0.08, -0.09]
    coords2 = base.copy()
    coords2[:, 2] = [-0.12, 0.16, -0.07, -0.11, 0.14]

    log = tmp_path / "scan.log"
    log.write_text(
        "\n".join(
            [
                " Charge = 0 Multiplicity = 1",
                orientation_block(coords1),
                " SCF Done:  E(RB3LYP) =  -100.0000000000     A.U. after    1 cycles",
                " ! Name     Definition           Value          Derivative Info.                !",
                " ! QPck0001 GIC-477              0.8000         estimate D2E/DX2                !",
                " ! PhiP0001 GIC-478              0.1000         estimate D2E/DX2                !",
                " Step number     1 out of a maximum of  20 on scan point   1 out of   2",
                orientation_block(coords2),
                " SCF Done:  E(RB3LYP) =  -100.0010000000     A.U. after    1 cycles",
                " QPck0001     0.80000   0.00000   0.00000   0.00000   0.00000   0.90000",
                " PhiP0001     0.10000   0.00000   0.00000   0.00000   0.00000   0.30000",
                " Step number     1 out of a maximum of  20 on scan point   2 out of   2",
                " Normal termination of Gaussian",
            ]
        )
        + "\n"
    )

    structures = dvr.read_gaussian_log(log, selection="all")
    assert len(structures) == 2
    assert structures[0].props["gic_QPck0001"] == 0.8
    assert structures[1].props["gic_PhiP0001"] == 0.3

    ring = [0, 1, 2, 3, 4]
    raw_energy, rel_energy = dvr.potential_cm(structures, "energy_hartree", "hartree")
    bridges = dvr.fit_gic_to_cremer_pople_bridges(structures, ring)
    assert 1 in bridges

    profile = tmp_path / "profile.csv"
    zeros = np.zeros(len(structures))
    dvr.write_profile(
        profile,
        structures,
        zeros,
        zeros,
        zeros,
        zeros,
        raw_energy,
        rel_energy,
        "energy_hartree",
        [],
        ring,
        bridges,
    )
    text = profile.read_text()
    assert "gic_QPck0001" in text
    assert "CP_m2_q_angstrom" in text
    assert "CP_from_GIC_m2_q_angstrom" in text
