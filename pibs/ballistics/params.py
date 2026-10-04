from __future__ import annotations

import math
from dataclasses import dataclass
from math import log

from . import GunType, SolutionMethod
from .config import SimConfig
from .nozzle import get_cf
from .num import dekker, intg
from .prop import Propellant


def pidduck(wpm: float, k: float, tol: float) -> tuple[float, float]:
    """
    Pidduck's limiting solution to the Lagrange problem.
    wpm : w/(phi_1 * m), charge mass to equivalent corrected (fictitious) shot
          weight
    k   : adiabatic index of the gas, in practice this is not a great influence
    tol : numerical tolerance

    Pidduck's solution is reduced to that of M.A.Mamontov's solution at k -> 1,
    however numerical difficulty necessitates taking the limit.
    """
    if k < 1:
        raise ValueError("Invalid adiabatic index passed", k)

    def f(om: float, x: float) -> float:
        if k == 1:
            return math.exp(-om * x**2)
        else:
            return (1 - om * x**2) ** (1 / (k - 1))

    def g(om: float, x: float) -> float:
        return f(om, x) * x**2

    def f_omega(om: float) -> float:
        """
        Solve Ω by finding the root of:
        1
        ∫ (1-Ωξ²)^[1/(k-1)] dξ = (w/m) (k-1)/(2 k) (1-Ω)^[k/(k-1)]/Ω
        0
        金（2014）《枪炮内弹道学》(3-114) pp.160

        for the case of k -> 1:
        鲍廷钰，邱文坚（1995）《内弹道学》pp.196
        """
        if om == 0:
            return -math.inf

        i = intg(lambda x: f(om, x), 0, 1, tol)

        if k == 1:
            return i - 0.5 * wpm * math.exp(-om) / om
        else:
            return i - 0.5 * ((k - 1) / k) * wpm * ((1 - om) ** (k / (k - 1)) / om)

    omega = dekker(f_omega, 0, 1, tol=tol)

    if k == 1:
        labda_1 = (math.exp(omega) - 1) / wpm
    else:
        labda_1 = ((1 - omega) ** (k / (1 - k)) - 1) / wpm

    i_u = intg(lambda x: g(omega, x), 0, 1, tol)
    i_l = intg(lambda x: f(omega, x), 0, 1, tol)
    labda_2 = i_u / i_l

    return labda_1, labda_2


@dataclass(frozen=True)
class Params:
    """Derived ballistic parameters computed once from a SimConfig.

    load density            Δ := ω / V₀
    asymptotic velocity:    vⱼ := [2 f ω / (θ φ m)]^0.5

    the reduced system:
    reduced length:         l̄  := l / l₀ where l₀ := V₀ / S
    reduced pressure:       p̄ := p / (f Δ)
    reduced time:           t̄ := t vⱼ / l₀
    reduced velocity:       v̄ := v / vⱼ

    auxiliary definitions:
    l̄φ = 1 - Δ/ρₚ - Δ(α - 1/ρₚ) · ψ
    B    = [S² · eᵢ² / (f · ω · φ · m · uᵢ²)] · (fΔ)^(2(1-n))

    Note: for typographical purposes, ω has been laid out as w in code, and sometimes l̄ is rendered as λ, spelled as
        "labda" due to Python keyword lambda.

    Nozzle fields (chi_0 .. s_j) are None for conventional guns.
    labda_1/labda_2 (Lagrange-problem distribution coefficients) are None
    for recoilless guns, which assume the simplified phi = φ₁ + ω/(3m).
    """

    propellant: Propellant
    gun_type: GunType

    caliber: float
    e_1: float
    s: float
    m: float
    w: float
    vol_0: float
    p_0: float
    l_g: float
    chi_k: float
    l_0: float
    l_c: float
    delta: float
    tol: float
    phi_1: float
    psi_0: float
    z_0: float
    phi: float
    v_j: float
    b: float

    labda_1: float | None = None
    labda_2: float | None = None

    chi_0: float | None = None
    a_bar: float | None = None
    c_a: float | None = None
    c_f: float | None = None
    s_j_bar: float | None = None
    s_j: float | None = None

    @property
    def is_recoilless(self) -> bool:
        return self.gun_type == GunType.RECOILLESS


def derive_params(cfg: SimConfig) -> Params:
    prop = cfg.propellant

    caliber = cfg.caliber
    e_1 = 0.5 * cfg.web_thickness
    s = (0.5 * caliber) ** 2 * math.pi
    m = cfg.shot_mass
    w = cfg.charge_mass
    vol_0 = cfg.chamber_volume
    p_0 = cfg.start_pressure
    l_g = cfg.barrel_length
    chi_k = cfg.chambrage
    l_0 = vol_0 / s
    l_c = l_0 / chi_k
    delta = w / vol_0
    tol = cfg.tolerance
    phi_1 = 1 / (1 - cfg.drag_coefficient)

    psi_0 = (1 / delta - 1 / prop.rho_p) / (prop.f / p_0 + prop.alpha - 1 / prop.rho_p)
    if psi_0 <= 0:
        raise ValueError(
            "Initial burnup fraction is solved to be negative. This indicate an excessively high load density for start-pressure."
        )
    elif psi_0 >= 1:
        raise ValueError(
            "Initial burnup fraction is solved to be greater than unity. This indicate an excessively low loading density for start-pressure."
        )
    z_0 = dekker(prop.f_psi_z, 0, prop.z_b, y=psi_0, tol=tol)

    labda_1 = labda_2 = None
    if cfg.gun_type == GunType.RECOILLESS:
        phi = phi_1 + w / (3 * m)
    elif cfg.gun_type == GunType.CONVENTIONAL:
        if cfg.solution_method == SolutionMethod.LAGRANGE:
            labda_1, labda_2 = 1 / 2, 1 / 3
        elif cfg.solution_method == SolutionMethod.PIDDUCK:
            labda_1, labda_2 = pidduck(w / (phi_1 * m), prop.theta + 1, tol)
        elif cfg.solution_method == SolutionMethod.MAMONTOV:
            labda_1, labda_2 = pidduck(w / (phi_1 * m), 1, tol)
        else:
            raise ValueError("Unknown Solution")

        labda = l_g / l_0
        # 见《枪炮内弹道学》（金，2014）p.70 式
        cc = 1 - (1 - 1 / chi_k) * log(labda + 1) / labda
        phi = phi_1 + labda_2 * cc * w / m
    else:
        raise ValueError("unknown gun type")

    b = s**2 * e_1**2 / (prop.f * phi * w * m * prop.u_1**2) * (prop.f * delta) ** (2 * (1 - prop.n))

    v_j = (2 * prop.f * w / (prop.theta * phi * m)) ** 0.5

    chi_0 = a_bar = c_a = c_f = s_j_bar = s_j = None
    if cfg.gun_type == GunType.RECOILLESS:
        chi_0 = cfg.efficiency
        a_bar = cfg.expansion_ratio

        gamma = prop.theta + 1
        """
        Enforce the "recoilless condition" by setting the size of the
        throat.
        """
        phi_2 = 1
        c_a = (
            (0.5 * prop.theta * phi * m / w) ** 0.5
            * gamma**0.5
            * (2 / (gamma + 1)) ** (0.5 * (gamma + 1) / prop.theta)
            * phi_2
        )
        c_f = get_cf(gamma, a_bar, tol)
        s_j_bar = 1 / (c_f * chi_0)
        if s_j_bar > chi_k:
            raise ValueError(
                "Achieving recoilless condition necessitates a larger throat area than could be fit into breech face."
            )
        s_j = s_j_bar * s

    return Params(
        propellant=prop,
        gun_type=cfg.gun_type,
        caliber=caliber,
        e_1=e_1,
        s=s,
        m=m,
        w=w,
        vol_0=vol_0,
        p_0=p_0,
        l_g=l_g,
        chi_k=chi_k,
        l_0=l_0,
        l_c=l_c,
        delta=delta,
        tol=tol,
        phi_1=phi_1,
        psi_0=psi_0,
        z_0=z_0,
        phi=phi,
        v_j=v_j,
        b=b,
        labda_1=labda_1,
        labda_2=labda_2,
        chi_0=chi_0,
        a_bar=a_bar,
        c_a=c_a,
        c_f=c_f,
        s_j_bar=s_j_bar,
        s_j=s_j,
    )
