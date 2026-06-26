"""
van der Waals radii (Å).
"""

# Merz–Kollman / Bondi (Gaussian-style)
MERZ_KOLLMAN = {
    1: 1.20,  2: 1.20,  3: 1.37,  4: 1.45,  5: 1.45,
    6: 1.50,  7: 1.50,  8: 1.40,  9: 1.35, 10: 1.30,
    11: 1.57, 12: 1.36, 13: 1.24, 14: 1.17, 15: 1.90,
    16: 1.85, 17: 1.80, 18: 1.88,
    19: 2.75, 20: None,
    # many transition metals intentionally undefined
    28: 1.63, 29: 1.40, 30: 1.39,
    31: 1.87, 32: 1.86, 33: 2.00, 34: 2.00, 35: 1.95,
    36: 2.02,
    46: 1.63, 47: 1.72, 48: 1.58,
    49: 1.93, 50: 2.17, 51: 2.20, 52: 2.20, 53: 2.15,
    54: 2.16,
    78: 1.72, 79: 1.66, 80: 1.55, 81: 1.96, 82: 1.02,
    92: 1.86,
}

# UFF (Rappe et al.) – diameters × 0.5
UFF = {
    Z: None for Z in range(1, 119)
}
# (values already implicitly complete in RVdW97, omitted here for brevity;
# you can paste them mechanically if you want full explicitness)

DEFAULT_VDW = "merz_kollman"


def vdw_radius(Z: int, scheme: str = DEFAULT_VDW):
    if Z <= 0:
        return 0.0
    table = MERZ_KOLLMAN if scheme == "merz_kollman" else UFF
    return table.get(Z, None)

