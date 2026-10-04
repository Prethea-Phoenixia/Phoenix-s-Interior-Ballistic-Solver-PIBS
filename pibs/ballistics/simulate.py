from __future__ import annotations

import logging
from dataclasses import dataclass
from math import exp, pi, tan
from typing import Callable

from . import Domain, Point
from .config import Material, SimConfig
from .generics import OutlineEntry, PressureProbePoint, PressureTraceEntry, SimResult, TableEntry
from .nozzle import get_pr
from .num import find_last_le, gss, merge_sorted_records, rkf
from .odes import L, MasterRhs, T, Z, in_domain, make_gun_rhs, make_recoilless_rhs
from .params import Params, derive_params

# Domain-specific ODEs are derived from the model's master rhs via
# in_domain(); see odes.py for the canonical state order.
# bar_data rows: (tag, t_bar, l_bar, z, v_bar, p_bar, eta, tau)

State = tuple[float, float, float, float, float]
BarRow = tuple[Point, float, float, float, float, float, float, float]


@dataclass
class _Model:
    """Per-gun-type closures and constants consumed by the shared pipeline."""

    f_p_bar: Callable[[float, float, float, float, float], float]  # (z, l_bar, v_bar, eta, tau) -> p_bar
    rhs: MasterRhs
    peak_tags: tuple[Point, ...]
    peak_bar: Callable[[Point, float, State, float], float]  # (tag, z, (t,l,v,eta,tau), p_bar) -> value to maximize
    build_row: Callable[[BarRow, int], tuple[TableEntry, PressureTraceEntry]]


def _gun_model(params: Params) -> _Model:
    prop = params.propellant
    theta, b, n = prop.theta, params.b, prop.n
    delta, rho_p, alpha = params.delta, prop.rho_p, prop.alpha
    p_scale = prop.f * delta

    def f_p_bar(z: float, l_bar: float, v_bar: float, eta: float, tau: float) -> float:
        psi = prop.f_psi_z(z)
        l_psi_bar = 1 - delta * ((1 - psi) / rho_p + (alpha * psi))
        p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)
        return p_bar

    def peak_bar(tag: Point, z: float, state: State, p_bar: float) -> float:
        if tag == Point.PEAK_AVG:
            return p_bar
        _, l_bar, _, _, _ = state
        ps_bar, pb_bar = _to_ps_pb(params, l_bar * params.l_0, p_bar)
        if tag == Point.PEAK_SHOT:
            return ps_bar
        elif tag == Point.PEAK_BREECH:
            return pb_bar
        raise ValueError(f"tag {tag} not handled.")

    def build_row(row: BarRow, trace_steps: int) -> tuple[TableEntry, PressureTraceEntry]:
        dtag, t_bar, l_bar, z, v_bar, p_bar, _, _ = row

        t = t_bar * (params.l_0 / params.v_j)
        l = l_bar * params.l_0
        psi = prop.f_psi_z(z)
        v = v_bar * params.v_j
        p = p_bar * p_scale
        ps, pb = _to_ps_pb(params, l, p)
        temp = _get_temperature(params, psi, l, p)

        p_line = []
        for i in range(trace_steps):
            x = i / trace_steps * (l + params.l_c)
            p_x, _ = _to_px_u(params, l, ps, pb, v, x)
            p_line.append(PressureProbePoint(x, p_x))
        p_line.append(PressureProbePoint(l + params.l_c, ps))

        return TableEntry(
            tag=dtag,
            time=t,
            travel=l,
            burnup=psi,
            velocity=v,
            breech_pressure=pb,
            avg_pressure=p,
            shot_pressure=ps,
            temperature=temp,
        ), PressureTraceEntry(dtag, temp, p_line)

    return _Model(
        f_p_bar=f_p_bar,
        rhs=make_gun_rhs(theta=theta, n=n, b=b, f_p_bar=f_p_bar),
        peak_tags=(Point.PEAK_AVG, Point.PEAK_SHOT, Point.PEAK_BREECH),
        peak_bar=peak_bar,
        build_row=build_row,
    )


def _recoilless_model(params: Params) -> _Model:
    prop = params.propellant
    theta, b, n = prop.theta, params.b, prop.n
    delta, rho_p, alpha = params.delta, prop.rho_p, prop.alpha
    c_a, s_j_bar = params.c_a, params.s_j_bar
    p_scale = prop.f * delta

    def f_p_bar(z: float, l_bar: float, v_bar: float, eta: float, tau: float, psi: float | None = None) -> float:
        psi = psi if psi is not None else prop.f_psi_z(z)
        l_psi_bar = 1 - delta * ((1 - psi) / rho_p + alpha * (psi - eta))

        p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

        return p_bar

    def peak_bar(tag: Point, z: float, state: State, p_bar: float) -> float:
        if tag == Point.PEAK_AVG:
            return p_bar
        _, l_bar, v_bar, eta, tau = state
        ps, p0, pb, _, _ = _to_ps_p0_pb_vb_stag(
            params, l_bar * params.l_0, v_bar * params.v_j, p_bar * p_scale, tau, eta
        )
        if tag == Point.PEAK_BREECH:
            return pb / p_scale
        elif tag == Point.PEAK_STAG:
            return p0 / p_scale
        elif tag == Point.PEAK_SHOT:
            return ps / p_scale
        raise ValueError("tag not handled.")

    def build_row(row: BarRow, trace_steps: int) -> tuple[TableEntry, PressureTraceEntry]:
        dtag, t_bar, l_bar, z, v_bar, p_bar, eta, tau = row

        t = t_bar * (params.l_0 / params.v_j)
        l = l_bar * params.l_0
        psi = prop.f_psi_z(z)
        v = v_bar * params.v_j
        p = p_bar * p_scale
        temp = tau * prop.temp_v if prop.temp_v else None

        ps, p0, pb, vb, stag = _to_ps_p0_pb_vb_stag(params, l, v, p, tau, eta)

        p_line = []
        for i in range(trace_steps):
            x = i / trace_steps * (l + params.l_c)
            px = _to_px(params, l, v, vb, ps, eta, x)
            p_line.append(PressureProbePoint(x, px))
        p_line.append(PressureProbePoint(l + params.l_c, ps))

        return TableEntry(
            tag=dtag,
            time=t,
            travel=l,
            burnup=psi,
            velocity=v,
            breech_pressure=pb,
            avg_pressure=p,
            shot_pressure=ps,
            temperature=temp,
            outflow_velocity=vb,
            stag_pressure=p0,
            outflow_fraction=eta,
            rel_stag_point=stag,
        ), PressureTraceEntry(dtag, temp, p_line)

    return _Model(
        f_p_bar=f_p_bar,
        rhs=make_recoilless_rhs(
            theta=theta,
            n=n,
            b=b,
            c_a=c_a,
            s_j_bar=s_j_bar,
            f_psi_z=prop.f_psi_z,
            f_sigma_z=prop.f_sigma_z,
            f_p_bar=f_p_bar,
        ),
        peak_tags=(Point.PEAK_AVG, Point.PEAK_SHOT, Point.PEAK_BREECH, Point.PEAK_STAG),
        peak_bar=peak_bar,
        build_row=build_row,
    )


def _to_ps_pb(params: Params, l: float, p: float) -> tuple[float, float]:
    """
    Convert average chamber pressure at certain travel to
    shot base pressure, and breech face pressure

    l: travel of the projectile
    p: average pressure

    Ps: pressure at shot
    Pb: pressure at breech
    """
    labda_g = l / params.l_0
    labda_1_prime = params.labda_1 * (1 / params.chi_k + labda_g) / (1 + labda_g)
    labda_2_prime = params.labda_2 * (1 / params.chi_k + labda_g) / (1 + labda_g)

    factor_s = 1 + labda_2_prime * (params.w / (params.phi_1 * params.m))
    factor_b = (params.phi_1 * params.m + labda_2_prime * params.w) / (
        params.phi_1 * params.m + labda_1_prime * params.w
    )

    return p / factor_s, p / factor_b


def _to_px_u(params: Params, l: float, p_s: float, p_b: float, v: float, x: float) -> tuple[float, float]:
    """
    Convert the average chamber to pressure and gas flow speed
    at arbitrary point x for projectile travel of l and average pressure
    of p, **assuming the Lagrangian distribution**.

    Note that with the current state of research, only characteristic point
    values are available for other distributions, use _to_ps_pb() instead for that.

    l: projectile travel
    p_s: pressure of shot
    p_b: pressure of breech
    x: probe point, start from the breech face.
    """
    l_c, l_0, chi_k = params.l_c, params.l_0, params.chi_k
    r = chi_k * x if x < l_c else (x - l_c) + l_0
    k = (r / (l_0 + l)) ** 2
    p_x = p_s * k + p_b * (1 - k)

    if x < l_c:
        u = x * v / (l_0 + l)
    else:
        u = (x + (chi_k - 1) * l_c) * v / (l_0 + l)

    return p_x, u


def _get_temperature(params: Params, psi: float, l: float, p: float) -> float | None:
    """
    given pressure and travel, return temperature
    using the Nobel-Abel EOS
    """
    prop = params.propellant
    if not prop.temp_v:
        return None

    if psi:
        r = prop.f / prop.temp_v
        l_psi = params.l_0 * (1 - params.delta / prop.rho_p - params.delta * (prop.alpha - 1 / prop.rho_p) * psi)
        return params.s * p * (l + l_psi) / (params.w * psi * r)
    else:
        return prop.temp_v


def _to_ps_p0_pb_vb_stag(
    params: Params, l: float, v: float, p: float, tau: float, eta: float
) -> tuple[float, float, float, float, float]:
    """
    Diagrammatic explanation of the calculated values:
        ----\\___     __________________
                |---|                  |_________________________________
                   ==                                       |        \
       Nozzle Breech|    Chamber                  Barrel    |  Shot  |>
                   ==                   ____________________|________/____
             ___|---|__________________|
        ----//
    Cross-section of the breech face:
       _
     *- -*
    | 0 0 |   0+0: Nozzle throat area, S_j
     *-_-*

    returns:
        Ps: Shot base pressure
        Pb: Breech pressure, chamber pressure at rearward of chamber
        P0: Stagnation point pressure
        vb: Rearward flow velocity at the rear of chamber
        stag: relative location of the stagnation point compared to the volume behind projectile.
    """
    prop = params.propellant
    w, m, phi_1 = params.w, params.m, params.phi_1
    y = w * eta
    m_dot = params.c_a * params.v_j * params.s_j * p / (prop.f * tau**0.5)
    sb = params.s * params.chi_k
    vb = m_dot * (params.vol_0 + params.s * l) / (sb * (w - y))

    h_lim = 2 * phi_1 * m / (w * (1 - eta)) + 1
    h = h_lim if v == 0 else min(h_lim, vb / v)

    ps = p / (1 + (w - y) / (3 * phi_1 * m) * (1 - 0.5 * h))
    p0 = ps * (1 + (w - y) / (2 * phi_1 * m) * (1 + h) ** -1)
    pb = 0 if h == h_lim else ps * (1 + (w - y) / (2 * phi_1 * m) * (1 - h))
    stag = h / (1 + h)
    return ps, p0, pb, vb, stag


def _to_px(params: Params, l: float, v: float, vb: float, ps: float, eta: float, x: float) -> float:
    """
    convert x, the physical displacement from breech bottom, to
    effective length in the equivalent gun.
    """
    w, m, phi_1 = params.w, params.m, params.phi_1
    l_k = params.l_0 / params.chi_k
    s_k = params.s * params.chi_k
    if x < l_k:
        z = x * s_k / (l_k * s_k + l * params.s)
    else:
        z = (l_k * s_k + (x - l_k) * params.s) / (l_k * s_k + l * params.s)

    h_lim = 2 * phi_1 * m / (w * (1 - eta)) + 1
    h = h_lim if v == 0 else min(h_lim, vb / v)

    px = ps * (1 + w * (1 - eta) / (2 * phi_1 * m) * (1 + h) * (1 - z**2) - w * (1 - eta) / (phi_1 * m) * h * (1 - z))
    return px


def _integrate(params: Params, model: _Model, step: int, dom: Domain, logger: logging.Logger) -> SimResult:
    prop = params.propellant
    tol = params.tol

    ode_t = in_domain(model.rhs, T)
    ode_z = in_domain(model.rhs, Z)
    ode_l = in_domain(model.rhs, L)

    bar_data: list[BarRow] = []
    l_g_bar = params.l_g / params.l_0
    z_0, z_b = params.z_0, prop.z_b

    def append_bar_data(tag: Point, t_bar: float, l_bar: float, z: float, v_bar: float, eta: float, tau: float):
        bar_data.append((tag, t_bar, l_bar, z, v_bar, model.f_p_bar(z, l_bar, v_bar, eta, tau), eta, tau))

    append_bar_data(Point.START, 0, 0, z_0, 0, 0.0, 1.0)

    p_bar_max = 1e9 / (prop.f * params.delta)  # 1 GPa

    def abort_condition(_z: float, state: State, _) -> bool:
        _t_bar, _l_bar, _v_bar, _eta, _tau = state
        return _l_bar > l_g_bar or model.f_p_bar(_z, _l_bar, _v_bar, _eta, _tau) > p_bar_max

    z_record: list[tuple[float, State]] = [(z_0, (0.0, 0.0, 0.0, 0.0, 1.0))]

    # Phase 1: integrate to burnout or exit
    z_end, (t_bar_end, l_bar_end, v_bar_end, eta_end, tau_end), aborted = rkf(
        ode_z,
        (0.0, 0.0, 0.0, 0.0, 1.0),
        z_0,
        z_b,
        rel_tol=tol,
        abort_func=abort_condition,
        record=z_record,
    )

    # check if integration exited due to excess pressure
    if model.f_p_bar(z_end, l_bar_end, v_bar_end, eta_end, tau_end) > p_bar_max:
        raise ValueError(
            "excessive pressure encountered during integration (>1GPa mean). results cannot be expected to "
            "be accurate due to gross violation of the applicable domain of Nobel-Abel equation-of-state. "
        )

    # Determine if exit happened before burnout
    is_burn_out_contained = not (aborted and l_bar_end >= l_g_bar)

    if is_burn_out_contained:
        append_bar_data(Point.BURNOUT, t_bar_end, l_bar_end, z_b, v_bar_end, eta_end, tau_end)
        logger.info("integrated to burnout point")
    else:
        logger.warning("shot exited barrel before burnout")

    # integrate to actual exit point for both cases
    _, (t_bar_exit, z_exit, v_bar_exit, eta_exit, tau_exit), _ = rkf(
        ode_l,
        (t_bar_end, z_end, v_bar_end, eta_end, tau_end),
        l_bar_end,
        l_g_bar,
        rel_tol=tol,
    )
    append_bar_data(Point.EXIT, t_bar_exit, l_g_bar, z_exit, v_bar_exit, eta_exit, tau_exit)

    # Populate fracture point entry at Z = 1
    if z_b > 1.0 and z_exit >= 1.0:
        t_bar_f, l_bar_f, v_bar_f, eta_f, tau_f = rkf(ode_z, (0.0, 0.0, 0.0, 0.0, 1.0), z_0, 1, rel_tol=tol)[1]
        append_bar_data(Point.FRACTURE, t_bar_f, l_bar_f, 1, v_bar_f, eta_f, tau_f)

    # Peak finding
    def func_p_z(_z: float, tag: Point) -> tuple[float, tuple[float, float, float, float, float, float]]:
        _i = find_last_le(z_record, _z)
        _x = z_record[_i][0]
        ys = z_record[_i][1]

        r = []
        _t_bar, _l_bar, _v_bar, _eta, _tau = rkf(ode_z, ys, _x, _z, rel_tol=tol, record=r)[1]
        merge_sorted_records(z_record, r)

        _p_bar = model.f_p_bar(_z, _l_bar, _v_bar, _eta, _tau)
        return model.peak_bar(tag, _z, (_t_bar, _l_bar, _v_bar, _eta, _tau), _p_bar), (
            _z,
            _t_bar,
            _l_bar,
            _v_bar,
            _eta,
            _tau,
        )

    def find_peak(tag: Point) -> None:
        z_p = gss(lambda _z: func_p_z(_z, tag)[0], z_0, z_end, x_tol=tol, find_min=False)
        _, (z_p, t_bar_p, l_bar_p, v_bar_p, eta_p, tau_p) = func_p_z(z_p, tag)
        append_bar_data(tag, t_bar_p, l_bar_p, z_p, v_bar_p, eta_p, tau_p)

    for peak in model.peak_tags:
        find_peak(peak)

    # Sampling
    if dom == Domain.TIME:
        z_j, l_bar_j, v_bar_j, t_bar_j, eta_j, tau_j = z_0, 0.0, 0.0, 0.0, 0.0, 1.0
    else:  # length domain ODE requires starting from some point
        t_bar_j, (z_j, l_bar_j, v_bar_j, eta_j, tau_j), _ = rkf(
            ode_t, (z_0, 0.0, 0.0, 0.0, 1.0), 0, 0.5 * t_bar_exit, rel_tol=tol
        )

    for j in range(step):
        if dom == Domain.TIME:
            t_bar_k = t_bar_exit / (step + 1) * (j + 1)
            z_j, l_bar_j, v_bar_j, eta_j, tau_j = rkf(
                ode_t, (z_j, l_bar_j, v_bar_j, eta_j, tau_j), t_bar_j, t_bar_k, rel_tol=tol
            )[1]
            t_bar_j = t_bar_k
        else:
            l_bar_k = l_g_bar / (step + 1) * (j + 1)
            t_bar_j, z_j, v_bar_j, eta_j, tau_j = rkf(
                ode_l, (t_bar_j, z_j, v_bar_j, eta_j, tau_j), l_bar_j, l_bar_k, rel_tol=tol
            )[1]
            l_bar_j = l_bar_k

        append_bar_data(Point.SAMPLE, t_bar_j, l_bar_j, z_j, v_bar_j, eta_j, tau_j)

    logger.info(f"sampled {step} points")

    # Data processing
    trace_steps = max(step, 1)
    data, p_trace = zip(*(model.build_row(row, trace_steps) for row in bar_data))
    data, p_trace = zip(*sorted(zip(data, p_trace), key=lambda e: e[0].time))

    return SimResult(params=params, table_data=list(data), pressure_trace=list(p_trace))


def barrel_monoblock(
    xs: list[float],
    ps: list[float],
    ss: list[float],
    yield_strength: float,
) -> tuple[float, list[float], list[float]]:
    v, ks, ms = 0.0, [], []

    for p in ps:
        if p > yield_strength * 0.5:
            raise ValueError(
                f"Limit to conventional construction ({yield_strength * 0.5 * 1e-6:.3f} MPa)" + " exceeded in section."
            )
        k = (1 - 2 * p / yield_strength) ** -0.5
        ks.append(k)
        ms.append(1)

    for i in range(len(xs) - 1):
        x_0, x_1 = xs[i], xs[i + 1]
        s_0, s_1 = ss[i], ss[i + 1]
        k_0, k_1 = ks[i], ks[i + 1]
        dv = 0.5 * ((k_0**2 - 1) * s_0 + (k_1**2 - 1) * s_1) * (x_1 - x_0)
        v += dv

    return v, ks, ms


def barrel_autofrettage(
    xs: list[float],
    ps: list[float],
    ss: list[float],
    yield_strength: float,
) -> tuple[float, list[float], list[float]]:
    v, ks, ms = 0.0, [], []

    for _p in ps:
        m_opt = exp(_p / yield_strength)
        k = m_opt
        ks.append(k)
        ms.append(m_opt)

    for i in range(len(xs) - 1):
        x_0, x_1 = xs[i], xs[i + 1]
        s_0, s_1 = ss[i], ss[i + 1]
        k_0, k_1 = ks[i], ks[i + 1]
        dv = 0.5 * ((k_0**2 - 1) * s_0 + (k_1**2 - 1) * s_1) * (x_1 - x_0)
        v += dv
    return v, ks, ms


def gun_structure(
    params: Params,
    material: Material,
    safety_factor: float,
    autofrettage: bool,
    table_data: list[TableEntry],
    step: int,
    logger: logging.Logger,
) -> tuple[list[OutlineEntry], float]:
    tol = params.tol
    step = max(step, 1)

    r_b = 0.5 * params.caliber
    r_c = r_b * params.chi_k**0.5
    x_probes = (
        [i / step * params.l_c for i in range(step)]
        + [params.l_c * (1 - tol)]
        + [i / step * params.l_g + params.l_c for i in range(step)]
        + [params.l_g + params.l_c]
    )
    p_probes = [0.0 for _ in range(len(x_probes))]

    for entry in table_data:
        l = entry.travel
        v = entry.velocity
        p_s = entry.shot_pressure
        p_b = entry.breech_pressure
        for i, x in enumerate(x_probes):
            if (x - params.l_c) <= l:
                p_x, _ = _to_px_u(params, l, p_s, p_b, v, x)
                p_probes[i] = max(p_probes[i], p_x)
            else:
                break

    for i in range(len(p_probes)):
        p_probes[i] *= safety_factor

    barrel_fn = barrel_autofrettage if autofrettage else barrel_monoblock
    i = step + 1
    x_c, p_c = x_probes[:i], p_probes[:i]
    x_b, p_b = x_probes[i:], p_probes[i:]

    v_c, k_c, m_c = barrel_fn(x_c, p_c, [params.s * params.chi_k for _ in x_c], material.yield_strength)
    v_b, k_b, m_b = barrel_fn(x_b, p_b, [params.s for _ in x_b], material.yield_strength)

    v = v_c + v_b
    k_probes = k_c + k_b
    m_probes = m_c + m_b

    hull = []
    for x, k, m in zip(x_probes, k_probes, m_probes):
        if x < params.l_c:
            hull.append(OutlineEntry(x, r_c, k * r_c, m * r_c))
        else:
            hull.append(OutlineEntry(x, r_b, k * r_b, m * r_b))

    logger.info("structural calculation done")
    return hull, v * material.density


def recoilless_structure(
    params: Params,
    material: Material,
    safety_factor: float,
    autofrettage: bool,
    table_data: list[TableEntry],
    step: int,
    logger: logging.Logger,
) -> tuple[list[OutlineEntry], float]:
    tol = params.tol
    l_c = params.l_c
    l_g = params.l_g

    sigma = material.yield_strength
    s = params.s
    gamma = params.propellant.theta + 1
    r_b = 0.5 * params.caliber

    r_c = r_b * params.chi_k**0.5
    r_t = r_c * params.s_j_bar**0.5
    r_n = r_t * params.a_bar**0.5

    x_probes = (
        [i / step * l_c for i in range(step)]
        + [l_c * (1 - tol)]
        + [i / step * l_g + l_c for i in range(step)]
        + [l_g + l_c]
    )

    p_probes = [0.0 for _ in range(len(x_probes))]
    for entry in table_data:
        l = entry.travel
        v = entry.velocity
        vb = entry.outflow_velocity
        ps = entry.shot_pressure
        eta = entry.outflow_fraction

        for i, x in enumerate(x_probes):
            if (x - l_c) <= l:
                p_x = _to_px(params, l, v, vb, ps, eta, x)
                p_probes[i] = max(p_probes[i], p_x)
            else:
                break

    for i in range(len(p_probes)):
        p_probes[i] *= safety_factor

    """
    subscript k denote the end of the chamber
    subscript b denote the throat of nozzle
    subscript a denote the nozzle base.

    beta is the half angle of the constriction section
    alpha is the half angle of the expansion section
    """
    beta = 30
    l_b = (r_c - r_t) / tan(beta * pi / 180)
    alpha = 15
    l_a = (r_n - r_t) / tan(alpha * pi / 180)

    p_entry = p_probes[0]
    neg_x_probes = (
        [-l_b - l_a + l_a * i / step for i in range(step)]
        + [-l_b - tol * l_a]
        + [-l_b + l_b * i / step for i in range(step)]
        + [0]
    )

    r_probes, neg_ps, neg_ss = [], [], []

    for x in neg_x_probes:
        if x < -l_b:
            k = (x + l_b + l_a) / l_a
            r = k * r_t + (1 - k) * r_n
            p = get_pr(gamma, (r / r_t) ** 2, tol)[1] * p_entry
        else:
            k = (x + l_b) / l_b
            r = k * r_c + (1 - k) * r_t
            p = get_pr(gamma, (r / r_t) ** 2, tol)[0] * p_entry

        r_probes.append(r)
        neg_ps.append(p)
        neg_ss.append(pi * r**2)

    barrel_fn = barrel_autofrettage if autofrettage else barrel_monoblock
    v_n, k_n, m_n = barrel_fn(neg_x_probes, neg_ps, neg_ss, sigma)

    hull = []
    for x, r, k, m in zip(neg_x_probes, r_probes, k_n, m_n):
        hull.append(OutlineEntry(x, r, r * k, r * m))

    nozzle_mass = v_n * material.density

    i = step + 1
    x_c, p_c = x_probes[:i], p_probes[:i]
    x_b, p_b = x_probes[i:], p_probes[i:]

    v_c, k_c, m_c = barrel_fn(x_c, p_c, [s * params.chi_k for _ in x_c], sigma)
    v_b, k_b, m_b = barrel_fn(x_b, p_b, [s for _ in x_b], sigma)

    v = v_c + v_b
    k_probes = k_c + k_b
    m_probes = m_c + m_b

    tube_mass = v * material.density
    for x, k, m in zip(x_probes, k_probes, m_probes):
        if x < l_c:
            hull.append(OutlineEntry(x, r_c, k * r_c, m * r_c))
        else:
            hull.append(OutlineEntry(x, r_b, k * r_b, m * r_b))

    logger.info("structural calculation done")
    return hull, tube_mass + nozzle_mass


def simulate(cfg: SimConfig, step: int | None = None, logger: logging.Logger | None = None) -> SimResult:
    """Run a full interior-ballistics simulation: derive parameters, integrate
    to burnout/exit, find peak pressures, sample the trajectory, and (if a
    structural material is configured) size the barrel."""
    logger = logger or cfg.logger or logging.getLogger(__name__)
    params = derive_params(cfg)
    model = _recoilless_model(params) if params.is_recoilless else _gun_model(params)
    result = _integrate(params, model, cfg.step if step is None else step, cfg.domain, logger)

    if cfg.structural_material is not None:
        structure_fn = recoilless_structure if params.is_recoilless else gun_structure
        result.outline, result.tube_mass = structure_fn(
            params,
            cfg.structural_material,
            cfg.structural_safety_factor,
            cfg.autofrettage,
            result.table_data,
            cfg.step if step is None else step,
            logger,
        )

    return result
