"""ODE right-hand sides of interior ballistics, in reduced (dimensionless) form.

Canonical state order: (t_bar, z, l_bar, v_bar, eta, tau)
    t_bar : reduced time
    z     : relative burnt web coordinate
    l_bar : reduced projectile travel
    v_bar : reduced velocity
    eta   : relative mass discharged through nozzle (recoilless only)
    tau   : reduced temperature (recoilless only)

Conventional guns carry dummy eta=0, tau=1 with zero derivatives so that both
gun types share one integration pipeline.

The master right-hand sides below are defined in the time domain. Integration
against another independent variable is obtained by the chain rule:

    dY/dx = (dY/dt) / (dx/dt)

via in_domain() - which is what the historical per-domain ODEs did anyway, in
algebraically pre-cancelled form. Note that the length-domain ODE has a 1/v_bar
factor that poses a starting problem, preventing its use from initial
conditions.
"""

from __future__ import annotations

from typing import Callable

# canonical indices of valid independent variables within the full state
T, Z, L, V = 0, 1, 2, 3

MasterRhs = Callable[[float, float, float, float, float, float], tuple[float, float, float, float, float]]
PBarGun = Callable[[float, float, float, float, float], float]
PBarRecoilless = Callable[[float, float, float, float, float, float | None], float]


def make_gun_rhs(*, theta: float, n: float, b: float, f_p_bar: PBarGun) -> MasterRhs:
    """
    the ODE system (2-72):

    dψ/dt =
        χ(1 + 2λZ + 3μZ²) √(θ/(2B)) · p̄^n,          when Z < 1
        (χ_s/Z_k)(1 + 2λ_s·Z/Z_k) √(θ/(2B)) · p̄^n,   when 1 ≤ Z < Z_k
        0,                                           when Z ≥ Z_k

    dZ/dt =
        √(θ/(2B)) · p̄^n,    when Z < Z_k
        0,                   when Z ≥ Z_k

    dl̄/dt̄ = v̄

    dv̄/dt̄ = (θ/2) · p̄

    dp̄/dt = [l₀ / ((l̄ + l̄_φ) · vⱼ)] · [1 + Δ(α - 1/ρ) · p̄] · (dψ/dt)
             - [(1 + θ) / (l̄ + l̄_φ)] · p̄ · v̄

    (dψ/dt is folded into the propellant form functions behind f_p_bar.)
    """

    def rhs(_t: float, z: float, l_bar: float, v_bar: float, _eta: float, _tau: float):
        p_bar = f_p_bar(z, l_bar, v_bar, 0.0, 1.0)
        dz = (0.5 * theta / b) ** 0.5 * p_bar**n
        dl_bar = v_bar
        dv_bar = theta * 0.5 * p_bar
        return dz, dl_bar, dv_bar, 0.0, 0.0

    return rhs


def make_recoilless_rhs(
    *,
    theta: float,
    n: float,
    b: float,
    c_a: float,
    s_j_bar: float,
    f_psi_z: Callable[[float], float],
    f_sigma_z: Callable[[float], float],
    f_p_bar: PBarRecoilless,
) -> MasterRhs:
    """Recoilless master rhs: the gun system extended with outflow fraction η
    and reduced temperature τ (adiabatic cooling of the chamber charge)."""

    def rhs(_t: float, z: float, l_bar: float, v_bar: float, eta: float, tau: float):
        psi = f_psi_z(z)
        d_psi = f_sigma_z(z)
        p_bar = f_p_bar(z, l_bar, v_bar, eta, tau, psi)

        dz = (0.5 * theta / b) ** 0.5 * p_bar**n
        dl_bar = v_bar
        dv_bar = theta * 0.5 * p_bar

        d_eta = c_a * s_j_bar * p_bar * tau**-0.5
        d_tau = ((1 - tau) * (d_psi * dz) - 2 * v_bar * dv_bar - theta * tau * d_eta) / (psi - eta)

        return dz, dl_bar, dv_bar, d_eta, d_tau

    return rhs


def in_domain(rhs: MasterRhs, x_index: int) -> Callable:
    """Re-parametrize a time-domain master rhs against another independent
    variable by the chain rule, dY/dx = (dY/dt) / (dx/dt).

    x_index: canonical index of the new independent variable (T, Z, L or V).
    The returned ODE takes (x, dep_state), where dep_state holds the five
    remaining canonical components in canonical order, and returns their
    x-derivatives in that same order.

    The L domain divides by dl/dt = v_bar; the 1/v_bar factor poses a starting
    problem that prevents its use from initial conditions.
    """

    def ode(x: float, dep: tuple[float, ...]) -> tuple[float, ...]:
        full = list(dep)
        full.insert(x_index, x)
        t_bar, z, l_bar, v_bar, eta, tau = full
        dz, dl_bar, dv_bar, d_eta, d_tau = rhs(t_bar, z, l_bar, v_bar, eta, tau)
        derivs = (1.0, dz, dl_bar, dv_bar, d_eta, d_tau)
        q = 1.0 / derivs[x_index]
        return tuple(d * q for i, d in enumerate(derivs) if i != x_index)

    return ode
