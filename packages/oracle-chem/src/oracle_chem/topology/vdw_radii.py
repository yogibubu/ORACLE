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

# UFF (Rappe et al.) diameters converted to radii; ORACLE RVdW97 table.
_UFF_DIAMETERS = [
    2.886, 2.362, 2.451, 2.745, 4.083, 3.851, 3.660, 3.500, 3.364, 3.243,
    2.983, 3.021, 4.499, 4.295, 4.147, 4.035, 3.947, 3.868, 3.812, 3.399,
    3.295, 3.175, 3.144, 3.023, 2.961, 2.912, 2.872, 2.834, 3.495, 2.763,
    4.383, 4.280, 4.230, 4.205, 4.189, 4.141, 4.114, 3.641, 3.345, 3.124,
    3.165, 3.052, 2.998, 2.963, 2.929, 2.899, 3.148, 2.848, 4.463, 4.392,
    4.420, 4.470, 4.500, 4.404, 4.517, 3.703, 3.522, 3.556, 3.606, 3.575,
    3.547, 3.520, 3.493, 3.368, 3.451, 3.428, 3.409, 3.391, 3.374, 3.355,
    3.640, 3.141, 3.170, 3.069, 2.954, 3.120, 2.840, 2.754, 3.293, 2.705,
    4.347, 4.297, 4.370, 4.709, 4.750, 4.765, 4.900, 3.677, 3.478, 3.396,
    3.424, 3.395, 3.424, 3.424, 3.381, 3.326, 3.339, 3.313, 3.299, 3.286,
    3.274, 3.248, 3.236,
]
UFF = {Z: None for Z in range(1, 119)}
UFF.update({Z: 0.5 * value for Z, value in enumerate(_UFF_DIAMETERS, start=1)})

# UFF vdW well depths D_i in kcal/mol. Values are loaded where currently needed
# by GF non-bonded corrections; unsupported elements raise a clear error.
UFF_WELL_DEPTH_KCAL = {
    1: 0.044,
    2: 0.056,
    3: 0.025,
    4: 0.085,
    5: 0.180,
    6: 0.105,
    7: 0.069,
    8: 0.060,
    9: 0.050,
    10: 0.042,
    11: 0.030,
    12: 0.111,
    13: 0.505,
    14: 0.402,
    15: 0.305,
    16: 0.274,
    17: 0.227,
    18: 0.185,
}

DEFAULT_VDW = "merz_kollman"


def vdw_radius(Z: int, scheme: str = DEFAULT_VDW):
    if Z <= 0:
        return 0.0
    table = MERZ_KOLLMAN if scheme == "merz_kollman" else UFF
    return table.get(Z, None)


def uff_well_depth_kcal(Z: int):
    return UFF_WELL_DEPTH_KCAL.get(int(Z))
