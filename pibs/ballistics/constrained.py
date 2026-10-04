from __future__ import annotations

import logging
import math
from math import log
from typing import Callable

from . import GunType, OptimizationTarget, Point, SolutionMethod
from .config import SimConfig
from .nozzle import get_cf
from .num import dekker, find_last_le, gss, merge_sorted_records, rkf
from .odes import T, V, Z, in_domain, make_gun_rhs, make_recoilless_rhs
from .params import pidduck


def validate_constrained_cfg(cfg: SimConfig) -> None:
    if any(
        (
            cfg.caliber <= 0,
            cfg.shot_mass <= 0,
            cfg.start_pressure <= 0,
            cfg.drag_coefficient < 0,
            cfg.drag_coefficient >= 1,
            cfg.chambrage < 1,
        )
    ):
        raise ValueError("Invalid parameters for constrained design")

    if any((cfg.design_pressure <= 0, cfg.design_velocity <= 0)):
        raise ValueError("Invalid design constraint")

    if cfg.gun_type == GunType.RECOILLESS and any((cfg.expansion_ratio < 1, cfg.efficiency > 1, cfg.efficiency <= 0)):
        raise ValueError("Invalid parameters for constrained design")


def minimum_load_fraction(cfg: SimConfig) -> float:
    """
    p = fΔ / (1 - αΔ)
    Δ = 1 / (f/p + α)
    """
    prop = cfg.propellant
    return (1 / (prop.f / cfg.design_pressure + prop.alpha)) / prop.rho_p * (1 + cfg.tolerance)


def solve_constrained(
    cfg: SimConfig,
    load_fraction: float,
    charge_mass_ratio: float,
    known_bore: bool,
    length_gun: float | None = None,
    logger: logging.Logger | None = None,
) -> tuple[float, float]:
    logger = logger or cfg.logger or logging.getLogger(__name__)
    validate_constrained_cfg(cfg)

    if load_fraction < minimum_load_fraction(cfg):
        raise ValueError(
            "Design pressure cannot be achieved, in the limit of closed bomb operation, at the current load fraction."
        )
    if load_fraction >= 1:
        raise ValueError("Chamber is overfull.")

    if known_bore:
        if length_gun is None:
            raise ValueError("known_bore option requires length_gun")

    if cfg.gun_type == GunType.CONVENTIONAL:
        return _solve_gun(
            cfg,
            load_fraction=load_fraction,
            charge_mass_ratio=charge_mass_ratio,
            known_bore=known_bore,
            length_gun=length_gun,
            logger=logger,
        )
    elif cfg.gun_type == GunType.RECOILLESS:
        return _solve_recoilless(
            cfg,
            load_fraction=load_fraction,
            charge_mass_ratio=charge_mass_ratio,
            known_bore=known_bore,
            length_gun=length_gun,
            logger=logger,
        )
    else:
        raise ValueError("unknown gun type")


def _solve_gun(
    cfg: SimConfig,
    load_fraction: float,
    charge_mass_ratio: float,
    known_bore: bool = False,
    length_gun: float | None = None,
    logger: logging.Logger | None = None,
) -> tuple[float, float]:
    logger = logger or logging.getLogger(__name__)
    prop = cfg.propellant
    m = cfg.shot_mass
    s = (cfg.caliber / 2) ** 2 * math.pi
    phi_1 = 1 / (1 - cfg.drag_coefficient)
    p_d, v_d = cfg.design_pressure, cfg.design_velocity
    chi_k = cfg.chambrage
    tol = cfg.tolerance
    theta, f, alpha, rho_p, u_1, n = prop.theta, prop.f, prop.alpha, prop.rho_p, prop.u_1, prop.n

    w = m * charge_mass_ratio
    v_0 = w / (rho_p * load_fraction)
    delta = w / v_0

    l_0 = v_0 / s

    if length_gun is None:
        l_bar_g_0 = cfg.max_length / l_0
    else:
        l_bar_g_0 = length_gun / l_0

    if cfg.solution_method == SolutionMethod.LAGRANGE:
        labda_1, labda_2 = 1 / 2, 1 / 3
    elif cfg.solution_method == SolutionMethod.PIDDUCK:
        labda_1, labda_2 = pidduck(w / (phi_1 * m), theta + 1, tol)
    elif cfg.solution_method == SolutionMethod.MAMONTOV:
        labda_1, labda_2 = pidduck(w / (phi_1 * m), 1, tol)
    else:
        raise ValueError("Unknown Solution")

    psi_0 = (1 / delta - 1 / rho_p) / (f / cfg.start_pressure + alpha - 1 / rho_p)
    z_0 = dekker(prop.f_psi_z, 0, prop.z_b, y=psi_0, tol=tol)

    def func_p_control_bar(z: float, l_bar: float, v_bar: float) -> float:
        psi = prop.f_psi_z(z)
        l_psi_bar = 1 - delta / rho_p - delta * (alpha - 1 / rho_p) * psi
        p_bar = max((psi - v_bar**2), 0) / max((l_bar + l_psi_bar), tol)
        if cfg.pressure_control_point == Point.PEAK_AVG:
            return p_bar
        else:
            cc_prime = (1 / chi_k + l_bar) / (1 + l_bar)
            labda_1_prime = labda_1 * cc_prime
            labda_2_prime = labda_2 * cc_prime

            factor_s = 1 + labda_2_prime * (w / (phi_1 * m))
            factor_b = (phi_1 * m + labda_2_prime * w) / (phi_1 * m + labda_1_prime * w)

            if cfg.pressure_control_point == Point.PEAK_SHOT:
                return p_bar / factor_s
            elif cfg.pressure_control_point == Point.PEAK_BREECH:
                return p_bar / factor_b
            else:
                raise ValueError(f"Unknown control {cfg.pressure_control_point}")

    p_bar_s = func_p_control_bar(z_0, 0, 0)
    p_bar_d = p_d / (f * delta)

    if p_bar_d < p_bar_s:
        raise ValueError("Design pressure cannot be lower than shot start pressure.")
    l_bar_d = cfg.max_length / l_0

    def abort_z(x: float, ys: tuple[float, ...], record: list[tuple[float, tuple[float, ...]]]) -> bool:
        z, (t_bar, l_bar, v_bar, _, _) = x, ys
        p_bar = func_p_control_bar(z, l_bar, v_bar)
        o_x, o_ys = record[-1]
        o_z = o_x
        ot_bar, ol_bar, ov_bar, _, _ = o_ys
        op_bar = func_p_control_bar(o_z, ol_bar, ov_bar)

        return (p_bar < op_bar) or (p_bar > 2 * p_bar_d)

    # iterate the chambrage correction cc against the solved gun length;
    # only cc and the length guess change between iterations
    cc: float | None = None
    it = 0
    while True:
        if cc is None:
            """
            张小兵，金志明（2014），《枪炮内弹道学》，北京理工大学出版社，pp 70 (1-128)
            """
            cc = 1 - (1 - 1 / chi_k) * log(l_bar_g_0 + 1) / l_bar_g_0

        phi = phi_1 + labda_2 * w / m * cc
        v_j: float = (2 * f * w / (theta * phi * m)) ** 0.5

        if v_j < v_d and not known_bore:
            raise ValueError(
                f"Design velocity cannot be achieved with the current load conditions. The 2nd ballistic limit for this loading \
            conditions is {v_j:.4g} m/s."
            )

        v_bar_d = v_d / v_j

        def f_p_guarded(z: float, l_bar: float, v_bar: float, _eta: float, _tau: float) -> float:
            psi = prop.f_psi_z(z)
            l_psi_bar = 1 - delta / rho_p - delta * (alpha - 1 / rho_p) * psi
            return max((psi - v_bar**2), 0) / max((l_bar + l_psi_bar), tol)

        def func_p_e_1(e_1: float) -> tuple[float, float, float, float, float]:
            b_e_1 = s**2 * e_1**2 / (f * phi * w * m * u_1**2) * (f * delta) ** (2 * (1 - n))
            ode_z = in_domain(make_gun_rhs(theta=theta, n=n, b=b_e_1, f_p_bar=f_p_guarded), Z)

            record = [(z_0, (0.0, 0.0, 0.0, 0.0, 1.0))]

            z_k, state_k, _ = rkf(
                d_func=ode_z,
                ini_val=(0.0, 0.0, 0.0, 0.0, 1.0),
                x_0=z_0,
                x_1=prop.z_b,
                rel_tol=tol,
                abort_func=abort_z,
                record=record,
            )
            t_bar_k, l_bar_k, v_bar_k = state_k[0], state_k[1], state_k[2]

            p_bar_j = func_p_control_bar(z_k, l_bar_k, v_bar_k)

            if p_bar_j >= 2 * p_bar_d:
                return p_bar_j - p_bar_d, z_k, t_bar_k, l_bar_k, v_bar_k

            if len(record) > 1:
                z_j = record[-2][0]
            else:
                z_j = z_0

            def func_p_z(z: float) -> tuple[float, float, float, float, float]:
                i = find_last_le(record, z)
                x = record[i][0]
                ys = record[i][1]

                r = []
                state = rkf(d_func=ode_z, ini_val=ys, x_0=x, x_1=z, rel_tol=tol, record=r)[1]
                merge_sorted_records(record, r)
                t_bar, l_bar, v_bar = state[0], state[1], state[2]
                return func_p_control_bar(z, l_bar, v_bar) - p_bar_d, z, t_bar, l_bar, v_bar

            z_p = gss(lambda z: func_p_z(z)[0], z_j, z_k, x_tol=tol, find_min=False)
            return func_p_z(z_p)

        probe_web = 0.5 * cfg.min_web
        dp_bar_probe = next_dp_bar_probe = func_p_e_1(probe_web)[0]
        next_probe_web = probe_web

        while dp_bar_probe * next_dp_bar_probe > 0:
            probe_web = next_probe_web
            dp_bar_probe = next_dp_bar_probe

            next_probe_web = probe_web * 2 if dp_bar_probe > 0 else probe_web * 0.5
            next_dp_bar_probe = func_p_e_1(next_probe_web)[0]

        # terminate on the pressure residual: the peak pressure should match
        # the design pressure to tol * p_d
        e_1_solved = dekker(lambda _e_1: func_p_e_1(_e_1)[0] + p_bar_d, probe_web, next_probe_web, y=p_bar_d, tol=tol)
        p_bar_dev, z_i, t_bar_i, l_bar_i, v_bar_i = func_p_e_1(e_1_solved)

        if known_bore:
            return e_1_solved, length_gun

        if v_j * v_bar_i > v_d:
            raise ValueError(f"Design velocity exceeded before peak pressure point (V = {v_bar_i * v_j:.4g} m/s).")

        b = s**2 * e_1_solved**2 / (f * phi * w * m * u_1**2) * (f * delta) ** (2 * (1 - n))

        def f_p_free(z: float, l_bar: float, v_bar: float, _eta: float, _tau: float) -> float:
            psi = prop.f_psi_z(z)
            l_psi_bar = 1 - delta / rho_p - delta * (alpha - 1 / rho_p) * psi
            return (psi - v_bar**2) / (l_bar + l_psi_bar)

        ode_v = in_domain(make_gun_rhs(theta=theta, n=n, b=b, f_p_bar=f_p_free), V)

        def abort_v(_: float, ys: tuple[float, ...], record: list[tuple[float, tuple[float, ...]]]) -> bool:
            t_bar, _, l_bar, _, _ = ys
            ot_bar = record[-1][1][0]
            return l_bar > l_bar_d or t_bar < ot_bar

        v_t_z_l_record = [(v_bar_i, (t_bar_i, z_i, l_bar_i, 0.0, 1.0))]
        try:
            v_bar_g, state_g, _ = rkf(
                d_func=ode_v,
                ini_val=v_t_z_l_record[-1][1],
                x_0=v_t_z_l_record[-1][0],
                x_1=v_bar_d,
                rel_tol=tol,
                abort_func=abort_v,
                record=v_t_z_l_record,
            )
            t_bar_g, z_g, l_bar_g = state_g[0], state_g[1], state_g[2]
        except ValueError:
            v_bar_m, state_m = v_t_z_l_record[-1]
            t_bar_m, z_m, l_bar_m = state_m[0], state_m[1], state_m[2]
            p_max = func_p_control_bar(z_m, l_bar_m, v_bar_m) * f * delta
            v_max = v_bar_m * v_j
            l_max: float = l_bar_m * l_0
            raise ValueError(
                "Integration appears to be approaching asymptote, "
                + f"last calculated to v = {v_max:.4g} m/s, x = {l_max:.4g} m, p = {p_max * 1e-6:.4g} MPa. "
                + "This indicates an excessive velocity target relative to pressure developed."
            )

        p_bar_g = func_p_control_bar(z_g, l_bar_g, v_bar_g)

        v_g = v_bar_g * v_j
        l_g: float = l_bar_g * l_0

        p_g = p_bar_g * f * delta
        if l_bar_g > l_bar_d:
            raise ValueError(
                f"Solution requires excessive tube length, last calculated to v = {v_g:.4g} m/s, x = {l_g:.4g} m, p = {p_g * 1e-6:.4g} MPa."
            )
        cc_n = 1 - (1 - 1 / chi_k) * log(l_bar_g + 1) / l_bar_g
        l_g = l_bar_g * l_0
        if abs((l_bar_g - l_bar_g_0) / min(l_bar_g, l_bar_g_0)) > tol:
            if it < cfg.max_iterations:
                cc = cc_n
                l_bar_g_0 = l_bar_g
                it += 1
                continue
            logger.warning(
                f"chambrage correction did not converge within {cfg.max_iterations} iterations; returning last iterate"
            )

        logger.info(
            f"ω/m={charge_mass_ratio:.2f}, Δ/ρ={load_fraction:.2f} -> e_1={e_1_solved * 1e3:.2f} mm, l_g={l_g * 1e3:.0f} mm ({it + 1} iters)"
        )
        return e_1_solved, l_g


def _solve_recoilless(
    cfg: SimConfig,
    load_fraction: float,
    charge_mass_ratio: float,
    known_bore: bool = False,
    length_gun: float | None = None,
    logger: logging.Logger | None = None,
) -> tuple[float, float]:
    logger = logger or logging.getLogger(__name__)
    prop = cfg.propellant
    m = cfg.shot_mass
    s = (cfg.caliber / 2) ** 2 * math.pi
    phi_1 = 1 / (1 - cfg.drag_coefficient)
    p_d, v_d = cfg.design_pressure, cfg.design_velocity
    chi_k = cfg.chambrage
    tol = cfg.tolerance
    theta, f, alpha, rho_p, u_1, n = prop.theta, prop.f, prop.alpha, prop.rho_p, prop.u_1, prop.n
    chi_0, a_bar = cfg.efficiency, cfg.expansion_ratio

    w = m * charge_mass_ratio
    vol_0 = w / (rho_p * load_fraction)

    delta = w / vol_0
    l_0 = vol_0 / s

    gamma = theta + 1
    phi = phi_1 + w / (3 * m)

    s_j_bar = 1 / (get_cf(gamma, a_bar, tol) * chi_0)
    if s_j_bar > chi_k:
        raise ValueError(
            "Achieving recoilless condition necessitates a larger throat area than could be fit into breech face."
        )
    s_j = s_j_bar * s

    k_0 = (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1))) * gamma**0.5

    phi_2 = 1
    c_a = (0.5 * theta * phi * m / w) ** 0.5 * k_0 * phi_2
    v_j: float = (2 * f * w / (theta * phi * m)) ** 0.5

    if v_j < v_d and not known_bore:
        raise ValueError(
            "propellant load too low to achieve design velocity, "
            + f"the 2nd ballistic limit for this loading conditions is {v_j:.4g} m/s,"
            + " and recoilless guns only achieve a part of that as well.",
        )

    psi_0 = (1 / delta - 1 / rho_p) / (f / cfg.start_pressure + alpha - 1 / rho_p)
    z_0 = dekker(lambda _z: prop.f_psi_z(_z), 0, prop.z_b, y=psi_0, tol=tol)

    def f_p_avg(z: float, l_bar: float, v_bar: float, eta: float, tau: float, psi: float | None = None) -> float:
        psi = psi if psi is not None else prop.f_psi_z(z)
        l_psi_bar = 1 - delta * ((1 - psi) / rho_p + alpha * (psi - eta))
        return tau / (l_bar + l_psi_bar) * (psi - eta)

    def f_p_bar_control(z: float, l_bar: float, v_bar: float, eta: float, tau: float) -> float:
        p_bar = f_p_avg(z, l_bar, v_bar, eta, tau)

        m_dot = c_a * v_j * s_j * p_bar * delta * tau**-0.5
        vb = m_dot * (vol_0 + s * l_bar * l_0) / (s * chi_k * w * (1 - eta))

        h_lim = 2 * phi_1 * m / (w * (1 - eta)) + 1
        h = h_lim if v_bar == 0 else min(h_lim, vb / (v_j * v_bar))

        if cfg.pressure_control_point == Point.PEAK_AVG:
            return p_bar
        else:
            p_s_bar = p_bar / (1 + w * (1 - eta) / (3 * phi_1 * m) * (1 - 0.5 * h))
            if cfg.pressure_control_point == Point.PEAK_SHOT:
                return p_s_bar
            elif cfg.pressure_control_point == Point.PEAK_STAG:
                return p_s_bar * (1 + w * (1 - eta) / (2 * phi_1 * m) * (1 + h) ** -1)
            elif cfg.pressure_control_point == Point.PEAK_BREECH:
                if h == h_lim:
                    return 0.0
                else:
                    return p_s_bar * (1 + w * (1 - eta) / (2 * phi_1 * m) * (1 - h))
        raise ValueError(f"unknown control {cfg.pressure_control_point}")

    p_bar_d = p_d / (f * delta)
    l_bar_d = cfg.max_length / l_0

    def abort_t(
        x: float,
        ys: tuple[float, ...],
        record: list[tuple[float, tuple[float, ...]]],
    ) -> bool:
        t_bar, (z, l_bar, v_bar, eta, tau) = x, ys
        p_bar_control = f_p_bar_control(z=z, l_bar=l_bar, v_bar=v_bar, eta=eta, tau=tau)
        p_avg_bar = f_p_avg(z, l_bar, v_bar, eta, tau)

        o_t_bar, (o_z, o_l_bar, o_v_bar, o_eta, o_tau) = record[-1]
        op_bar_control = f_p_bar_control(z=o_z, l_bar=o_l_bar, v_bar=o_v_bar, eta=o_eta, tau=o_tau)
        op_avg_bar = f_p_avg(o_z, o_l_bar, o_v_bar, o_eta, o_tau)

        return (p_bar_control > p_bar_d * (1 + 2 * tol)) or (
            (p_bar_control < op_bar_control) and (p_avg_bar < op_avg_bar)
        )

    def func_p_e1(e_1: float) -> tuple[float, float, tuple[float, ...]]:
        b_e_1 = (s**2 * e_1**2) / (f * phi * w * m * u_1**2) * (f * delta) ** (2 * (1 - n))
        ode_t = in_domain(
            make_recoilless_rhs(
                theta=theta,
                n=n,
                b=b_e_1,
                c_a=c_a,
                s_j_bar=s_j_bar,
                f_psi_z=prop.f_psi_z,
                f_sigma_z=prop.f_sigma_z,
                f_p_bar=f_p_avg,
            ),
            T,
        )

        record = [(0.0, (z_0, 0.0, 0.0, 0.0, 1.0))]

        aborted = False
        while not aborted:
            # march in chunks of 10 reduced-time units; abort is checked per
            # committed rkf step, so chunk size only sets record granularity
            t_bar_next = record[-1][0] + 10
            t_bar_k, state_k, aborted = rkf(
                d_func=ode_t,
                ini_val=record[-1][1],
                x_0=record[-1][0],
                x_1=t_bar_next,
                rel_tol=tol,
                abort_func=abort_t,
                record=record,
            )
            record.append((t_bar_k, tuple(state_k)))

        t_bar_k, (z_k, l_bar_k, v_bar_k, eta_k, tau_k) = record[-1]
        p_bar_k = f_p_bar_control(z_k, l_bar_k, v_bar_k, eta_k, tau_k)

        if p_bar_k > (1 + 2 * tol) * p_bar_d:
            return p_bar_k - p_bar_d, t_bar_k, (z_k, l_bar_k, v_bar_k, eta_k, tau_k)

        t_bar_j = 0.0

        def func_p_t_bar(t_bar: float) -> tuple[float, float, tuple[float, ...]]:
            i = find_last_le(record, t_bar)
            x = record[i][0]
            ys = record[i][1]

            r: list[tuple[float, tuple[float, ...]]] = []
            state = rkf(d_func=ode_t, ini_val=ys, x_0=x, x_1=t_bar, rel_tol=tol, record=r)[1]
            merge_sorted_records(record, r)
            z, l_bar, v_bar, eta, tau = state
            return f_p_bar_control(z, l_bar, v_bar, eta, tau) - p_bar_d, t_bar, tuple(state)

        t_p = gss(lambda t_bar: func_p_t_bar(t_bar)[0], t_bar_j, t_bar_k, x_tol=tol, find_min=False)
        return func_p_t_bar(t_p)

    probe_web = next_probe_web = 0.5 * cfg.min_web
    dp_bar_probe = next_dp_bar_probe = func_p_e1(probe_web)[0]

    while dp_bar_probe * next_dp_bar_probe > 0:
        probe_web = next_probe_web
        dp_bar_probe = next_dp_bar_probe
        next_probe_web = (probe_web * 2) if dp_bar_probe > 0 else (probe_web * 0.5)
        next_dp_bar_probe = func_p_e1(next_probe_web)[0]

    # terminate on the pressure residual: the peak pressure should match
    # the design pressure to tol * p_d
    e_1_solved = dekker(lambda _e_1: func_p_e1(_e_1)[0] + p_bar_d, probe_web, next_probe_web, y=p_bar_d, tol=tol)

    if known_bore:
        return e_1_solved, length_gun

    v_bar_d = v_d / v_j
    dp_bar_i, t_bar_i, (z_i, l_bar_i, v_bar_i, eta_i, tau_i) = func_p_e1(e_1_solved)

    if v_bar_i > v_bar_d:
        v = v_bar_i * v_j
        raise ValueError(f"Design velocity exceeded before peak pressure point (V = {v:.4g} m/s).")

    b = (s**2 * e_1_solved**2) / (f * phi * w * m * u_1**2) * (f * delta) ** (2 * (1 - n))

    ode_v = in_domain(
        make_recoilless_rhs(
            theta=theta,
            n=n,
            b=b,
            c_a=c_a,
            s_j_bar=s_j_bar,
            f_psi_z=prop.f_psi_z,
            f_sigma_z=prop.f_sigma_z,
            f_p_bar=f_p_avg,
        ),
        V,
    )

    def abort_v(
        x: float,
        ys: tuple[float, ...],
        record: list[tuple[float, tuple[float, ...]]],
    ) -> bool:
        t_bar, z, l_bar, eta, tau = ys
        ot_bar, *_ = record[-1][-1]
        return l_bar > l_bar_d or t_bar < ot_bar

    v_t_z_l_eta_tau_record = [(v_bar_i, (t_bar_i, z_i, l_bar_i, eta_i, tau_i))]
    try:
        v_bar_g, state_g, _ = rkf(
            d_func=ode_v,
            ini_val=v_t_z_l_eta_tau_record[-1][1],
            x_0=v_t_z_l_eta_tau_record[-1][0],
            x_1=v_bar_d,
            rel_tol=tol,
            abort_func=abort_v,
            record=v_t_z_l_eta_tau_record,
        )
        t_bar_g, z_g, l_bar_g, eta_g, tau_g = state_g

    except ValueError:
        v_bar_m, (t_bar_m, z_m, l_bar_m, eta_m, tau_m) = v_t_z_l_eta_tau_record[-1]
        p_max = f_p_bar_control(z_m, l_bar_m, v_bar_m, eta_m, tau_m) * f * delta
        v_max = v_bar_m * v_j
        l_max = l_bar_m * l_0
        raise ValueError(
            "Integration appears to be approaching asymptote, "
            + "last calculated to v = {:.4g} m/s, ".format(v_max)
            + "x = {:.4g} m, p = {:.4g} MPa. ".format(l_max, p_max * 1e-6)
            + "This indicates an excessive velocity target relative to pressure developed."
        )

    p_bar_g = f_p_bar_control(z_g, l_bar_g, v_bar_g, eta_g, tau_g)

    v_g = v_bar_g * v_j
    l_g: float = l_bar_g * l_0
    p_g = p_bar_g * f * delta

    if l_bar_g > l_bar_d:
        raise ValueError(
            "Solution requires excessive tube length, last calculated to "
            + "v = {:.4g} m/s, x = {:.4g} m, ".format(v_g, l_g)
            + "p = {:.4g} MPa.".format(p_g * 1e-6)
        )

    elif t_bar_g < t_bar_i:
        raise ValueError(
            "Target velocity is achieved before peak pressure point, "
            + "last calculated at v = {:.4g} m/s,".format(v_g)
            + "x = {:.4g} m, p = {:.4g} MPa. ".format(l_g, p_g * 1e-6)
        )

    logger.info(
        f"ω/m={charge_mass_ratio:.2f}, Δ/ρ={load_fraction:.2f} -> e_1={e_1_solved * 1e3:.2f} mm, l_g={l_g * 1e3:.0f} mm"
    )
    return e_1_solved, l_bar_g * l_0


def get_f(cfg: SimConfig, charge_mass_ratio: float, logger: logging.Logger | None = None) -> Callable:
    s = (cfg.caliber / 2) ** 2 * math.pi
    rho_p = cfg.propellant.rho_p

    def _f(load_fraction: float) -> tuple[float, float, float]:
        e_1_delta, l_g_delta = solve_constrained(
            cfg, load_fraction=load_fraction, charge_mass_ratio=charge_mass_ratio, known_bore=False, logger=logger
        )
        l_0 = (cfg.shot_mass * charge_mass_ratio / (rho_p * load_fraction)) / s
        return e_1_delta, l_g_delta, l_g_delta + l_0

    return _f


def maximum_load_fraction(cfg: SimConfig, charge_mass_ratio: float, logger: logging.Logger | None = None) -> float:
    f = get_f(cfg, charge_mass_ratio, logger=logger)
    x = minimum_load_fraction(cfg)
    delta = 1 - cfg.tolerance - x
    while delta > cfg.tolerance:
        try:
            f(x + delta)
            x += delta
        except ValueError:
            delta *= 0.5

    return x


def validate_load_fraction(
    cfg: SimConfig, charge_mass_ratio: float, proposed_lfs: list[float], logger: logging.Logger | None = None
) -> dict[float, tuple[float, float, float]]:
    # construct the load fraction ladder:
    f = get_f(cfg, charge_mass_ratio, logger=logger)
    proposed_lfs = [lf for lf in proposed_lfs if lf > minimum_load_fraction(cfg)]
    results = {}
    for lf in proposed_lfs:
        try:
            results[lf] = f(lf)
        except ValueError:
            break

    return results


def find_min_v(
    cfg: SimConfig,
    *,
    charge_mass_ratio: float,
    opt_target: OptimizationTarget = OptimizationTarget.MIN_BARR_VOLUME,
    logger: logging.Logger | None = None,
) -> tuple[float, float, float]:
    """
    find the minimum volume solution.
    """
    logger = logger or cfg.logger or logging.getLogger(__name__)
    logger.info("optimizing under constraints")
    low = minimum_load_fraction(cfg)
    logger.info(f"min Δ/ρ = {low:.3%}")

    high = maximum_load_fraction(cfg, charge_mass_ratio, logger=logger)
    logger.info(f"max Δ/ρ = {high:.3%}")

    if opt_target == OptimizationTarget.MIN_PROJ_TRAVEL:
        _f_index = 1
    elif opt_target == OptimizationTarget.MIN_BARR_VOLUME:
        _f_index = 2
    else:
        raise ValueError(f"Unknown target {opt_target}")

    _f = get_f(cfg, charge_mass_ratio, logger=logger)

    logger.info(f"Δ/ρ range: {low:.3%} - {high:.3%}")

    try:
        _f(low)
    except ValueError as e:
        raise ValueError(f"constraint not satisfied at minimum load fraction: {e}")

    lf = gss(lambda load_fraction: _f(load_fraction)[_f_index], low, high, x_tol=cfg.tolerance, find_min=True)

    e_1, l_g, _ = _f(lf)
    logger.info(f"optimal Δ/ρ = {lf:.2f}")
    return lf, e_1, l_g
