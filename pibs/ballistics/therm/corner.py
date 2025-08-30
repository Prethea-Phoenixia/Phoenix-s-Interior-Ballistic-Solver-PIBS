import logging
from math import exp
from os import path

from ..num import quadratic
from .nasa7 import Reaction, Specie

# create logger
logger = logging.getLogger("CORNER")
logger.setLevel(logging.INFO)

# create console handler and set level to info
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)

Specie.read(path.join(path.dirname(__file__), "nasa7.dat"))

CO2 = Specie.get("CO2")
H2O = Specie.get("H2O")
CO = Specie.get("CO")
H2 = Specie.get("H2")
H = Specie.get("H")
N2 = Specie.get("N2")
OH = Specie.get("OH")
NO = Specie.get("NO")
N = Specie.get("N")
O = Specie.get("O")
O2 = Specie.get("O2")
CH4 = Specie.get("CH4")
NH3 = Specie.get("NH3")

# TABLE 3.9 Corrections To Equilibrium Constants from Corner.
# T in Kelvin, -DeltaB in cc/(gm.mol), -DeltaC/2 in (cc/(gm.mol))^2

negDeltaTable = [
    [1600, 34.2, 490],
    [1700, 33.8, 460],
    [1800, 33.5, 435],
    [1900, 33.2, 410],
    [2000, 33.0, 390],
    [2100, 32.8, 370],
    [2200, 32.6, 355],
    [2300, 32.5, 340],
    [2400, 32.4, 325],
    [2500, 32.2, 310],
    [2600, 32.1, 300],
    [2700, 32.0, 290],
    [2800, 31.9, 280],
    [2900, 31.8, 270],
    [3000, 31.7, 260],
    [3100, 31.6, 255],
    [3200, 31.5, 245],
    [3300, 31.4, 235],
    [3400, 31.4, 230],
    [3500, 31.3, 225],
    [3600, 31.2, 215],
    [3700, 31.1, 210],
    [3800, 31.1, 205],
    [3900, 31.0, 200],
    [4000, 30.9, 195],
]


"""
# TABLE 2.04 from Hunt
# T in Kelvin, K0(T)..... K8(T)
# PDF page 111 Book page 91
#
# P = V/RT
#
# K0         = [CO] [H2O] / [CO2] [H2]            CO2 + H2 <-> CO + H2O
# K1 sqrt(P) = [H] / [H2]^0.5                       0.5 H2 <->  H
# K2 sqrt(P) = [OH] [H2]^0.5 / [H2O]                  H2O  <-> OH + 0.5 H2
# K3 sqrt(P) = [NO] [H2] / [H2O] [N2]^0.5     H2O + 0.5 N2 <-> NO + H2
# K4 sqrt(P) = [N] / [N2]^0.5                       0.5 N2 <-> N
# K5   P     = [O] [H2] / [H2O]                       H2O  <-> O + H2
# K6   P     = [O2] [H2]^2 / [H2O]^2                 2H2O  <-> O2 + 2H2
# K7  P^-2   = [CH4][CO2] / [CO]^2 [H2]^2     2 CO + 2 H2  <-> CH4 + CO2
# K8  P^-1   = [NH3] / [N2]^0.5 [H2]^1.5   0.5 N2 + 1.5 H2 <-> NH3
#
#
# fmt: off
equilibriumKT = [
 [800, 0.2478,      0,      0,      0,      0,      0,     0,   31.25, 2.91e-3],
 [1000,0.7286,      0,      0,      0,      0,      0,     0, 3.72e-2, 5.64e-4],
 [1200, 1.435,      0,      0,      0,      0,      0,     0, 3.99e-4, 1.85e-4],
 [1400, 2.270,      0,      0,      0,      0,      0,     0, 1.55e-5, 8.2e-5],
 [1500, 2.704,      0,      0,      0,      0,      0,     0,  4.2e-6, 5.9e-5],
 [1600, 3.132,   6e-5,      0,      0,      0,      0,     0, 1.35e-6, 4.5e-5],
 [1700, 3.555,   2e-4,      0,      0,      0,      0,     0,  4.9e-7, 3.5e-5],
 [1800, 3.975,   4e-4,   2e-5,      0,      0,      0,     0,  2.0e-7, 2.8e-5],
 [1900, 4.385,   9e-4,   6e-5,      0,      0,      0,     0,  9.2e-8, 2.3e-5],
 [2000, 4.782, 0.0016,   2e-4,      0,      0,      0,     0,  4.5e-8, 1.9e-5],
 [2100, 5.161, 0.0031,   4e-4,   2e-5,      0,      0,     0,  2.4e-8, 1.6e-5],
 [2200, 5.520, 0.0057,   8e-4,   5e-5,      0,      0,     0,  1.3e-8, 1.4e-5],
 [2300, 5.852, 0.0097, 0.0016,   1e-4,   1e-5,   1e-5,     0,       0, 1.2e-5],
 [2400, 6.155, 0.0159, 0.0029,   2e-4,   3e-5,   3e-5,   1e-5,      0, 1.1e-5],
 [2500, 6.433, 0.0251, 0.0052,   4e-4,   7e-5,   9e-5,   3e-5,      0, 1.0e-5],
 [2600, 6.694, 0.0383, 0.0090,   7e-4,   2e-4,   2e-4,   1e-4,      0, 0.9e-5],
 [2700, 6.939, 0.0566, 0.0146, 0.0012,   3e-4,   5e-4,   2e-4,      0, 0.8e-5],
 [2800, 7.167, 0.0814, 0.0231, 0.0020,   5e-4, 0.0012,   5e-4,      0, 0.8e-5],
 [2900, 7.379, 0.1143, 0.0355, 0.0034,   8e-4, 0.0026, 0.0010,      0, 0.7e-5],
 [3000, 7.574, 0.1574, 0.0529, 0.0055, 0.0014, 0.0053, 0.0020,      0, 0.6e-5],
 [3100, 7.753, 0.2125, 0.0768, 0.0086, 0.0022, 0.0103, 0.0038,      0,      0],
 [3200, 7.917, 0.2813, 0.1089, 0.0129, 0.0035, 0.0190, 0.0071,      0,      0],
 [3300, 8.068, 0.3658, 0.1513, 0.0188, 0.0053, 0.0339, 0.0126,      0,      0],
 [3400, 8.205, 0.4682, 0.2064, 0.0271, 0.0078, 0.0586, 0.0216,      0,      0],
 [3500, 8.330, 0.5910, 0.2760, 0.0387, 0.0113, 0.0978, 0.0358,      0,      0],
 [3600, 8.443, 0.7367, 0.3626, 0.0539, 0.0161, 0.1591, 0.0580,      0,      0],
 [3700, 8.544, 0.9079, 0.4693, 0.0734, 0.0225, 0.2516, 0.0917,      0,      0],
 [3800, 8.633, 1.1070, 0.6001, 0.0982, 0.0309, 0.3886, 0.1412,      0,      0],
 [3900, 8.710, 1.3360, 0.7589, 0.1300, 0.0417, 0.5878, 0.2132,      0,      0],
 [4000, 8.775, 1.5970, 0.9495, 0.1693, 0.0554, 0.8711, 0.3157,      0,      0],
]

# fmt: on
"""

k0 = Reaction("water-gas", LHS={CO2: 1, H2: 1}, RHS={CO: 1, H2O: 1})
k1 = Reaction("hydrogen decomposition", LHS={H2: 0.5}, RHS={H: 1})
k2 = Reaction("water-hydroxyl", LHS={H2O: 1}, RHS={H2: 0.5, OH: 1})
k3 = Reaction("water-nitroxide", LHS={H2O: 1, N2: 0.5}, RHS={NO: 1, H2: 1})
k4 = Reaction("nitrogen decomposiiton", LHS={N2: 0.5}, RHS={N: 1})
k5 = Reaction("water-oxygen radical", LHS={H2O: 1}, RHS={O: 1, H2: 1})
k6 = Reaction("water-decomposition", LHS={H2O: 1}, RHS={O2: 0.5, H2: 1})
k7 = Reaction(
    "methane synthesis w/ carbon dioxide",
    LHS={CO: 2, H2: 2},
    RHS={CH4: 1, CO2: 1},
)
k8 = Reaction("ammonia synthesis", LHS={N2: 0.5, H2: 1.5}, RHS={NH3: 1})

# TABLE 3.8 Corrections to Pressure from Corner
#            B                    C
# T in K, H2, N2/CO, CO2, H2O, H2, N2/CO, CO2, H2O

BCTable = [
    [1600, 16.4, 32.1, 45.7, -4.2, 20, 210, 1385, 220],
    [1700, 16.3, 32.3, 47.3, -2.5, 20, 200, 1305, 210],
    [1800, 16.2, 32.4, 48.7, -1.1, 20, 190, 1235, 195],
    [1900, 16.1, 32.6, 49.9, 0.2, 20, 180, 1170, 185],
    [2000, 16.0, 32.6, 50.9, 1.2, 15, 170, 1110, 175],
    [2100, 15.9, 32.7, 51.8, 2.2, 15, 160, 1055, 170],
    [2200, 15.8, 32.7, 52.6, 3.0, 15, 155, 1010, 160],
    [2300, 15.7, 32.8, 53.2, 3.7, 15, 150, 965, 155],
    [2400, 15.6, 32.8, 53.8, 4.4, 15, 140, 925, 145],
    [2500, 15.6, 32.8, 54.4, 5.0, 15, 135, 885, 140],
    [2600, 15.5, 32.7, 54.8, 5.5, 15, 130, 855, 135],
    [2700, 15.4, 32.7, 55.3, 6.0, 15, 125, 825, 130],
    [2800, 15.3, 32.7, 55.6, 6.4, 10, 120, 795, 125],
    [2900, 15.3, 32.6, 56.0, 6.8, 10, 120, 765, 120],
    [3000, 15.2, 32.6, 56.2, 7.1, 10, 115, 740, 120],
    [3100, 15.1, 32.6, 56.5, 7.5, 10, 110, 720, 115],
    [3200, 15.0, 32.5, 56.7, 7.7, 10, 105, 695, 110],
    [3300, 15.0, 32.4, 56.9, 8.0, 10, 105, 675, 105],
    [3400, 14.9, 32.4, 57.1, 8.3, 10, 100, 650, 105],
    [3500, 14.8, 32.3, 57.3, 8.5, 10, 95, 635, 100],
    [3600, 14.8, 32.3, 57.4, 8.7, 10, 95, 615, 100],
    [3700, 14.7, 32.2, 57.5, 8.9, 10, 90, 600, 95],
    [3800, 14.7, 32.2, 57.6, 9.1, 10, 90, 585, 95],
    [3900, 14.6, 32.1, 57.7, 9.3, 10, 85, 570, 90],
    [4000, 14.5, 32.0, 57.8, 9.4, 10, 85, 555, 90],
]


# Table 3.2 from Corner
# Add 1.987 for constant pressure values.
# Mean molecular heats over the temperature range of 300 deg K to T deg K"
# "at zero density from Corner
# in calories per gram-mole per degree at constant volume
# T (K), CO2, H2O, CO, H2, N2, OH, NO, O2, CH4, NH3
# for monoatomic gas, take 2.980

"""
# fmt:off
MMHTable = [
 [800,  8.896,  6.599, 5.244, 5.019, 5.179, 5.092, 5.432, 5.570,  9.870, 8.322],
 [1000, 9.409,  6.883, 5.403, 5.059, 5.326, 5.136, 5.597, 5.759, 11.129, 9.021],
 [1200, 9.824,  7.176, 5.553, 5.117, 5.468, 5.197, 5.744, 5.915, 12.243, 9.672],
 [1400, 10.165, 7.467, 5.684, 5.191, 5.598, 5.288, 5.871, 6.044, 13.235, 10.270],
 [1500, 10.313, 7.612, 5.743, 5.232, 5.658, 5.337, 5.927, 6.100, 13.683, 10.548],
 [1600, 10.449, 7.752, 5.799, 5.275, 5.714, 5.385, 5.979, 6.152, 14.098, 10.811],
 [1700, 10.574, 7.887, 5.851, 5.318, 5.767, 5.432, 6.026, 6.202, 14.486, 11.061],
 [1800, 10.690, 8.017, 5.899, 5.362, 5.816, 5.477, 6.069, 6.248, 14.849, 11.298],
 [1900, 10.797, 8.142, 5.944, 5.407, 5.862, 5.520, 6.109, 6.292, 15.188, 11.523],
 [2000, 10.896, 8.263, 5.986, 5.451, 5.906, 5.563, 6.147, 6.335, 15.506, 11.735],
 [2100, 10.989, 8.380, 6.026, 5.495, 5.946, 5.606, 6.183, 6.375, 15.798, 11.936],
 [2200, 11.075, 8.492, 6.063, 5.538, 5.984, 5.648, 6.217, 6.413, 16.083, 12.126],
 [2300, 11.156, 8.599, 6.098, 5.580, 6.020, 5.691, 6.250, 6.450, 16.347, 12.307],
 [2400, 11.233, 8.702, 6.131, 5.621, 6.054, 5.734, 6.281, 6.486, 16.595, 12.478],
 [2500, 11.305, 8.801, 6.162, 5.662, 6.086, 5.774, 6.309, 6.522, 16.828, 12.640],
 [2600, 11.373, 8.896, 6.191, 5.702, 6.117, 5.812, 6.335, 6.557],
 [2700, 11.437, 8.988, 6.218, 5.740, 6.146, 5.849, 6.360, 6.591],
 [2800, 11.498, 9.076, 6.244, 5.778, 6.174, 5.885, 6.384, 6.625],
 [2900, 11.556, 9.160, 6.269, 5.815, 6.200, 5.922, 6.406, 6.658],
 [3000, 11.611, 9.241, 6.293, 5.851, 6.225, 5.957, 6.427, 6.690],
 [3100, 11.664, 9.319, 6.315, 5.886, 6.249, 5.990, 6.448, 6.721],
 [3200, 11.714, 9.395, 6.336, 5.920, 6.271, 6.023, 6.468, 6.752],
 [3300, 11.762, 9.468, 6.356, 5.953, 6.293, 6.054, 6.487, 6.782],
 [3400, 11.808, 9.538, 6.376, 5.985, 6.314, 6.085, 6.505, 6.812],
 [3500, 11.852, 9.605, 6.394, 6.017, 6.333, 6.114, 6.522, 6.841],
 [3600, 11.895, 9.670, 6.412, 6.047, 6.352, 6.143, 6.539, 6.869],
 [3700, 11.936, 9.733, 6.429, 6.077, 6.370, 6.170, 6.555, 6.896],
 [3800, 11.975, 9.794, 6.446, 6.106, 6.385, 6.197, 6.570, 6.923],
 [3900, 12.013, 9.852, 6.462, 6.135, 6.404, 6.223, 6.585, 6.949],
 [4000, 12.050, 9.908, 6.477, 6.163, 6.420, 6.248, 6.599, 6.974],
]
# fmt:on
"""


"""
Table 2.08 from Hunt
E1 /100
H2, N2 & CO, CO2, H2O

for E2 * 1e^-4:
H2:     3
N2, CO: 34,
CO2:    220
H2O:    35
"""

E1Table = [
    [1600, 45, -110, -870, -935],
    [1700, 50, -95, -840, -895],
    [1800, 55, -80, -810, -855],
    [1900, 65, -70, -785, -825],
    [2000, 70, -55, -755, -795],
    [2100, 75, -40, -730, -770],
    [2200, 80, -25, -700, -745],
    [2300, 90, -15, -675, -725],
    [2400, 95, -0, -645, -705],
    [2500, 100, 15, -620, -685],
    [2600, 105, 25, -595, -665],
    [2700, 110, 40, -565, -650],
    [2800, 120, 55, -540, -635],
    [2900, 125, 65, -515, -620],
    [3000, 130, 80, -490, -605],
    [3100, 135, 95, -465, -590],
    [3200, 140, 105, -440, -580],
    [3300, 145, 120, -415, -565],
    [3400, 150, 130, -390, -555],
    [3500, 160, 145, -365, -540],
    [3600, 165, 155, -340, -530],
    [3700, 170, 170, -315, -520],
    [3800, 175, 185, -290, -510],
    [3900, 180, 195, -265, -500],
    [4000, 185, 210, -240, -490],
]


"""
(X): nbr. molecule per unit mass of propellant
{Y}: nbr. atom per unit mass of propellant

unless otherwise stated, mol/gram is assumed for
propellant work derived form Hunt.

Eq. 2.04
(CO) * (H2O) / [(CO2) * (H2)] = K0

{C} = (CO) + (CO2)                 | + (CH4)
{N} = 2 (N2)                       | + (NO) + (N) + (NH3)
{H} = 2 (H2O) + 2 (H2) +           | + (H) + (OH) + 3(NH3) + 4(CH4)
{O} = (CO) + (H2O) + 2 (CO2)       | + (OH) + (NO) + (O) + 2 (O2)

Major                              Minor
(CO) = {C} - (CO2)                 | - (CH4)

dissociation considered, in descending order of significance:

[H] = [H2]^0.5                         (V/RT)^0.5 K1(T)
[OH] = [H2O][H2]^-0.5                  (V/RT)^0.5 K2(T)
[NO] = [H2O][N2]^0.5[H2]^-1            (V/RT)^0.5 K3(T)
[N] = [N2]^0.5                         (V/RT)^0.5 K4(T)
[O] = [H2O][H2]^-1                     (V/RT)     K5(T)
[O2] = ([H2O]/[H2])^2                  (V/RT)     K6(T)
[CH4] = [CO]^2[H2]^2/[CO2]             (V/RT)^-2  K7(T)
[NH3] = [N2]^0.5[H2]^1.5               (V/RT)^-1  K8(T)
"""


def balance(Hf, T, Ci, Hi, Oi, Ni, V=1 / 0.1, its=200, tol=1e-9):
    """
    Ci, Hi, Oi, Ni are in mol/g.
    Consequently,
    CO2j, OHj, Hj, NOj, O2j, Oj, Nj are also in mol/g

    T: temperature in Kelvin
    V: volume per unit mass of gas, in cc/g
    1/V is also commonly known as Delta, or load density, if the propellant
    solid covolume can be ignored.
    in literature, usually 1/V = 0.2 g/cc = 200 kg/m^3 is adopted for IB
    purposes. However, very often covolume values are given for a specific
    load condition.

        > a value of 0.089(0.084) g/cc is adopted for propellant referenced
        from ADA043778. When this value is adopted, prediciton is good to
        within 1% of tabulated.

        >It is possible that for modern software, even if the user sets a 1/V
        of 0.2 g/cc, the program will apply its own averaging routine resulting
        in different result than is calculated here, given the apparent accuracy
        of all other parameter estimation performed here, and the closeness the
        composition can be solved to published data.

    Covolume is suitably invariant over typical internal ballistic
    temperature, range, however it varies not-insignificantly by load density
    in test chambers, and thus must be chosen with reference to gas density
    achieved in actual gun system.

    n: number of gram-molecules per unit mass of gas, in mol/g
    """
    if T > 4000 or T < max(v.fitTlow for v in (CO2, H2O, CO, H2, N2, OH, NO, O2, CH4, NH3)):
        raise ValueError("T Not supported")

    negDeltaB, neghalfDeltaC = 0, 0
    for i in range(len(negDeltaTable) - 1):
        Tlow, negDeltaBlow, neghalfDeltaClow = negDeltaTable[i]
        Thigh, negDeltaBhigh, neghalfDeltaChigh = negDeltaTable[i + 1]
        if Tlow <= T <= Thigh:
            k = (T - Tlow) / (Thigh - Tlow)
            negDeltaB = negDeltaBlow * (1 - k) + negDeltaBhigh * k
            neghalfDeltaC = neghalfDeltaClow * (1 - k) + neghalfDeltaChigh * k
            break

    BC = [0 for _ in range(8)]  # this disables covolume compensation for major products for T < 1600K
    for i in range(len(BCTable) - 1):
        Tlow, *BClow = BCTable[i]
        Thigh, *BChigh = BCTable[i + 1]
        if Tlow <= T <= Thigh:
            k = (T - Tlow) / (Thigh - Tlow)
            BC = [vi * (1 - k) + vj * k for vi, vj in zip(BClow, BChigh)]
            break

    BH2, BN2CO, BCO2, BH2O, CH2, CN2CO, CCO2, CH2O = BC

    DeltaT = T - 300

    HH = HO = HN = 2.980 * DeltaT  # monoatomic
    HCO2, HH2O, HCO, HH2, HN2, HOH, HNO, HO2, HCH4, HNH3 = (
        v.getMMH(T) * DeltaT for v in (CO2, H2O, CO, H2, N2, OH, NO, O2, CH4, NH3)
    )

    E1 = (0 for _ in range(4))
    for i in range(len(E1Table) - 1):
        Tlow, *E1low = E1Table[i]
        Thigh, *E1high = E1Table[i + 1]
        if Tlow <= T <= Thigh:
            k = (T - Tlow) / (Thigh - Tlow)
            E1 = [(vi * (1 - k) + vj * k) * 1e2 for vi, vj in zip(E1low, E1high)]
            break

    E1H2, E1N2, E1CO2, E1H2O = E1
    E1CO = E1N2
    """
    Second part of Table 3.10 from Corner, Page 130 PDF/110 Book
    "Corrections to Internal Energies"
    """
    if T < 1600:
        E2H2 = E2N2 = E2CO = E2CO2 = E2H2O = 0
    else:
        E2H2 = 3e4
        E2N2 = E2CO = 34e4
        E2CO2 = 220e4
        E2H2O = 35e4

    R = 82.06  # in cc atm /(K mol)
    sqrtVdivRT = (V / (R * T)) ** 0.5

    # initialize minor species.
    Hj, OHj, NOj, Nj, Oj, O2j, CH4j, NH3j = 0, 0, 0, 0, 0, 0, 0, 0

    # set up initial value for iteration: only major products
    G = Ci
    H = Oi - Ci
    I = 0.5 * Hi + Ci - Oi

    n = Ci + 0.5 * Hi + 0.5 * Ni  # mol/g n = None
    if negDeltaB == 0 and neghalfDeltaC == 0:
        K0 = k0(T)
    else:
        K0 = k0(T) * exp(n / V * negDeltaB + (n / V) ** 2 * neghalfDeltaC)

    # print(quadratic((1 - K0), -(G + H + K0 * I), G * H), Ci)
    soln_1, soln_2 = quadratic((1 - K0), -(G + H + K0 * I), G * H)  # two roots to the quadratic
    CO2j = None
    for soln in (soln_1, soln_2):
        if Ci >= soln >= 0:
            CO2j = soln if CO2j is None else min(soln, CO2j)
    if CO2j is None:
        raise ValueError(
            "Cannot get a first-estimate for CO2 molar concentration. "
            + "This is most likely a result of the mixture specified leaning towards "
            + "fuel poor and oxygen rich."
        )

    """get a first-order estimate of CO2 concentration assuming
    only major products (CO2, H2O, H2, CO, N2) by solving the
    equilibrium quadratic of the water-gas balance. This equilibrium
    is conveniently not dependent on pressure (since gas count
    remains constant throughout), as specified by Corner.

    Note: as a result of these assumptions, the method is only good for when the oxygen
    balance of the proepllant skews in fuel-rich (oxygen poor).

    N2j = 0.5 * Ni
    COj = G - CO2j
    H2Oj = H - CO2j
    H2j = I + CO2j
    """
    CO2j_nmp = CO2j  # no minor product

    i = its
    CO2j_0 = None
    CO2j_1 = None
    epsilon_0 = None
    epsilon = None
    kappa = 0.1
    delta = None
    olds = None

    firstRun = True
    while i > 0:
        N2j = 0.5 * (Ni - Nj - NOj - NH3j)

        G = Ci - CH4j
        COj = G - CO2j

        H = Oi - Ci - OHj - NOj - Oj - 2 * O2j + CH4j
        H2Oj = H - CO2j

        # fmt:off
        I = (0.5 * Hi - Oi + Ci - 0.5 * Hj + 0.5 * OHj - 1.5 * NH3j
             - 3 * CH4j + NOj + Oj + 2 * O2j)
        # fmt: on
        H2j = I + CO2j

        epsilon = COj * H2Oj - K0 * (CO2j * H2j)  # error

        if any((COj < 0, H2Oj < 0, CO2j < 0, H2j < 0)):
            logger.warning("Kappa-")
            kappa = kappa**2
            if kappa < 1e-16:
                raise ValueError("kappa has degenerated to zero.")
            CO2j = CO2j_nmp
            CO2j_0 = None
            epsilon = None
            epsilon_0 = None
            Hj, OHj, NOj, Nj, Oj, O2j, CH4j, NH3j = 0, 0, 0, 0, 0, 0, 0, 0
            i = its
            firstRun = True
            continue

        # Minor, Dissociation Products
        Hj = H2j**0.5 * sqrtVdivRT * k1(T)
        OHj = H2Oj / H2j**0.5 * sqrtVdivRT * k2(T) * exp(-20 * n / V)
        NOj = H2Oj * N2j**0.5 / H2j * sqrtVdivRT * k3(T) * exp(-20 * n / V)
        # NO + H2 -> H2O + 0.5 * N2
        Nj = N2j**0.5 * sqrtVdivRT * k4(T)
        Oj = (H2Oj / H2j) * sqrtVdivRT**2 * k5(T)
        O2j = (H2Oj / H2j) ** 2 * sqrtVdivRT**2 * k6(T)

        CH4j_1 = COj**2 * H2j**2 / CO2j * sqrtVdivRT**-4 * k7(T)

        # fmt: off
        CH4max = (0.5 * Hi - Oi + Ci - 0.5 * Hj + 0.5 * OHj - 1.5 * NH3j
                  + NOj + Oj + 2 * O2j + CO2j) / 3  # ensure H2j > 0
        # fmt: on

        CH4max = min(CH4max, Ci - CO2j)  # ensure COj > 0
        CH4min = -Oi + Ci + OHj + NOj + Oj + 2 * O2j + CO2j  # ensure H2Oj > 0
        CH4min = max(CH4min, 0)
        CH4j_1 = max(min(CH4max, CH4j_1), CH4min)  # enforce

        CH4j = CH4j + kappa * (CH4j_1 - CH4j)
        # damp out Methane oscillation

        NH3j_1 = N2j**0.5 * H2j**1.5 * sqrtVdivRT**-2 * k8(T)
        # fmt: off
        NH3max = (0.5 * Hi - Oi + Ci - 0.5 * Hj + 0.5 * OHj
                  - 3 * CH4j + NOj + Oj + 2 * O2j + CO2j) / 1.5
        # fmt: on
        NH3max = min(Ni - Nj - NOj, NH3max)
        NH3j_1 = min(NH3max, NH3j_1)
        NH3j = NH3j + kappa * (NH3j_1 - NH3j)

        if firstRun:
            # fmt:off
            olds = (CO2j, COj, H2Oj, H2j, Hj, OHj, NOj, Nj, Oj, O2j, CH4j, NH3j, N2j)
            # fmt: on

        else:
            delta = sum(
                abs(v - v0)
                for v, v0 in zip(
                    # fmt:off
                    (CO2j_1 if CO2j_1 is not None else CO2j, COj, H2Oj, H2j, Hj, OHj, NOj, Nj, Oj, O2j, CH4j_1, NH3j_1, N2j),
                    # fmt: on
                    olds,
                )
            )
            # fmt:off
            olds = (CO2j, COj, H2Oj, H2j, Hj, OHj, NOj, Nj, Oj, O2j, CH4j, NH3j, N2j)
            # fmt: on

        # fmt: off
        n = (CO2j + COj + H2Oj + H2j + Hj + OHj + NOj + Nj + Oj + O2j + CH4j
             + NH3j + N2j)  # mol/g
        # fmt: on

        # print(CO2j, COj, H2Oj, H2j)
        if negDeltaB == 0 and neghalfDeltaC == 0:
            pass
        else:
            K0 = k0(T) * exp(n / V * negDeltaB + (n / V) ** 2 * neghalfDeltaC)

        if firstRun:
            epsilon_0 = epsilon
            CO2j_0 = CO2j

            CO2j *= 1

            # we forfeit updating the first time around

        else:
            if epsilon_0 == epsilon:
                logger.warning("Solving terminated due to oscillating conditions.Δ={:}".format(delta))
                break

            CO2j_1 = CO2j - epsilon * (CO2j - CO2j_0) / (epsilon - epsilon_0)

            epsilon_0 = epsilon
            CO2j_0 = CO2j
            CO2j += kappa * (CO2j_1 - CO2j)

            if delta < tol and i < 50:
                break

        if firstRun:
            firstRun = False

        i -= 1

    if i == 0:
        logger.warning("Equilibrium finding terminated on maximum iterations. Result may not be correct.")

    logger.info("Solved for equilibrium condition in {:} iteration(s).".format(its - i) + "Δ = {:.2e}.".format(delta))
    speciesList = [
        ("CO", COj * 28.01, COj),
        ("CO2", CO2j * 44.01, CO2j),
        ("H2O", H2Oj * 18.02, H2Oj),
        ("NO", NOj * 30.01, NOj),
        ("OH", OHj * 17.01, OHj),
        ("N2", N2j * (14.01 * 2), N2j),
        ("N", Nj * 14.01, Nj),
        ("O2", O2j * (16.0 * 2), O2j),
        ("O", Oj * 16.0, Oj),
        ("H2", H2j * 1.008 * 2, H2j),
        ("H", Hj * 1.008, Hj),
        ("CH4", CH4j * 16.03, CH4j),
        ("NH3", NH3j * 17.03, NH3j),
    ]

    speciesList.sort(key=lambda x: x[1], reverse=True)

    """
    Find the internal energy of the gaseous products for major constitutents
    E = MMH * (T-300K) + E1 * n/V + E2 * (n/V)**2
    Correct for real gas effects of internal heat and energy content.
    """
    HH2 += E1H2 * (n / V) + E2H2 * (n / V) ** 2
    HN2 += E1N2 * (n / V) + E2N2 * (n / V) ** 2
    HCO2 += E1CO2 * (n / V) + E2CO2 * (n / V) ** 2
    HH2O += E1H2O * (n / V) + E2H2O * (n / V) ** 2
    HCO += E1CO * (n / V) + E2CO * (n / V) ** 2

    # fmt:off
    E = sum((HCO2 * CO2j, HH2O * H2Oj, HCO * COj, HH2 * H2j, HN2 * N2j,
            HOH * OHj, HNO * NOj, HO2 * O2j, HH * Hj, HO * Oj,
            HN * Nj, HCH4 * CH4j, HNH3 * NH3j))
    # fmt: on

    """
    Covolume calculations:
    """

    B = BCO2 * CO2j + BN2CO * (N2j + COj) + BH2O * H2Oj + BH2 * H2j
    C = CCO2 * CO2j + CN2CO * (N2j + COj) + CH2O * H2Oj + CH2 * H2j

    """
    p = nRT/(V-b) = nRT/V * (1+ B/V + nC/V^2)
    1/(V-b) = (1 + B/V + nC/V^2) / V
    V - b = V / (1+ B/V + nC/V^2)
    b = V - V / (1+ B/V + nC/V^2)
      = (B + nC/V) / (1+ B/V + nC/V^2)
      = (BV^2 + nCV) / (V^2 + BV + nC)
    """
    b = (B * V**2 + n * C * V) / (V**2 + B * V + n * C)
    p = n * R * T / (V - b) / 9.869

    """add the heat of formation of the products
    according to their proportions
    Table 3.4 Energies of Formation of Products from Corner
    From Graphite, Hydrogen Gas, Nitrogen Gas, and Oxygen at
    300K and Zero Pressure, Constant Volume Condition, Energies
    in Kilocalories per Mole.
    """
    Hf_p = -(
        COj * 26.69e3
        + CO2j * 94.03e3
        + H2Oj * 57.50e3  # this is in the gaseous form!
        + CH4j * 17.27
        + NH3j * 10.38
        + NOj * -21.53e3
        + OHj * -9.31e3
        + Nj * -85.25e3
        + Oj * -58.75e3
        + Hj * -51.79e3
    )

    DeltaEF = Hf_p - Hf
    # print("DeltaEF", DeltaEF)
    zeta = E + DeltaEF
    # print("Zeta", zeta)

    f = n * T * 8.314  # gas constant: 8.314 j/ mol K

    """ n: 1/M, mol/g, gas moles per gram
        b: covolume, cc/g
        p: pressure implied by load density, MPa
        f: propellant force, J/g
        E: internal energy of products.
    """
    """
    Cv = E / (T - 300)
    gamma = (n * 1.987 / Cv) + 1
    """
    return zeta, speciesList, n, E, b, p, f


if __name__ == "__main__":

    def f(T, ld):
        _, speciesList, _, _, b, p, f = balance(Hf=-470, T=T, Ci=0.02238, Hi=0.03014, Ni=0.01044, Oi=0.03468, V=1 / ld)
        print(T, ld)

        print(" @ Product  %mass  1e-6 mol/g")
        print(
            *[
                "{:>2} : {:^6} {:<6.1%}".format(i, name, mass) + "{:>10,.2f}".format(num * 1e6)
                for i, (name, mass, num) in enumerate(speciesList)
            ],
            sep="\n",
        )
        print("covol:", b, "g/cc")

        print("press:", p, "MPa")
        print("force:", f, "J/g")

    f(1200, 0.2)
