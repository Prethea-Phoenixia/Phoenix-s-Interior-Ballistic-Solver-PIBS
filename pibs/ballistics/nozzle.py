from __future__ import annotations

from math import inf

from .num import dekker


def get_cf(gamma: float, sr: float, tol: float = 1e-5) -> float:
    """
    takes the adiabatic index and area ration between constriction throat
    and the exit, calculate the thrust factor Cf
    See Hunt (1953) Interior Ballistics, appendix I.A01-A03
    Sr = S/St
    Vr = V/Vt
    """
    vr_old = 0
    vr = 1
    while abs(vr - vr_old) / vr > tol:
        vr_old = vr
        vr = ((gamma + 1) / (gamma - 1) * (1 - 2 / (gamma + 1) * (vr * sr) ** (1 - gamma))) ** 0.5

    cf = (2 / (gamma + 1)) ** (gamma / (gamma - 1)) * (gamma * vr + sr ** (1 - gamma) * vr ** (-gamma))

    return cf


def get_ar(gamma: float, pr: float) -> float:
    if pr <= 0 or pr >= 1:
        return inf
    return (
        (0.5 * (gamma + 1)) ** (1 / (gamma - 1))
        * pr ** (1 / gamma)
        * ((gamma + 1) / (gamma - 1) * (1 - pr ** ((gamma - 1) / gamma))) ** 0.5
    ) ** -1


def get_pr(gamma: float, ar: float, tol: float) -> tuple[float, float]:
    """
    Given an area ratio Ar, calculate the pressure ratio Pr for a
    converging-diverging nozzle. One area ratio corresponds to two
    pressure ratio, for subsonic and supersonic, respectively.

    Ar: Nozzle cs area at probe point x over throat area
        Ar = Ax / At
    Pr: pressure ratio at probe point x over upstream chamber pressure
        Pr = Px / Pc

    returns:
        Pr_sub: subsonic solution
        Pr_sup: supersonic solution
    """
    pr_c = (2 / (gamma + 1)) ** (gamma / (gamma - 1))

    if ar == 1:
        return pr_c, pr_c
    else:
        pr_sup = dekker(lambda pr: get_ar(gamma, pr), 0, pr_c, y=ar, tol=tol)
        pr_sub = dekker(lambda pr: get_ar(gamma, pr), pr_c, 1, y=ar, tol=tol)
        return pr_sub, pr_sup
