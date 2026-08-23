from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass
from math import inf, pi, tan
from typing import Callable

from . import (
    DOMAIN_TIME,
    POINT_BURNOUT,
    POINT_EXIT,
    POINT_FRACTURE,
    POINT_PEAK_AVG,
    POINT_PEAK_BREECH,
    POINT_PEAK_SHOT,
    POINT_PEAK_STAG,
    POINT_START,
    SAMPLE,
    Domains,
    Points,
    RecoillessPeakPoints,
)
from .base_gun import BaseGun
from .config import GunGeometry, Nozzle, PropellantLoad, Solver, Structural
from .gun import GenericEntry, GenericResult, OutlineEntry, PressureProbePoint, PressureTraceEntry
from .material import Material
from .num import dekker, gss, rkf


@dataclass
class RecoillessTableEntry(GenericEntry):
    outflow_velocity: float
    stag_pressure: float
    outflow_fraction: float
    rel_stag_point: float


@dataclass
class RecoillessResult(GenericResult):
    gun: Recoilless
    table_data: list[RecoillessTableEntry]


class Recoilless(BaseGun):
    def __init__(
        self,
        geometry: GunGeometry,
        load: PropellantLoad,
        nozzle: Nozzle,
        solver: Solver | None = None,
        logger=None,
    ):
        super().__init__(
            geometry=geometry,
            load=load,
            solver=solver,
            logger=logger,
        )

        self.nozzle = nozzle
        self.chi_0 = nozzle.efficiency
        self.a_bar = nozzle.expansion_ratio

        gamma = self.theta + 1
        """
        Enforce the "recoilless condition" by setting the size of the
        throat.
        """
        phi_2 = 1
        self.C_a = (
            (0.5 * self.theta * self.phi * self.m / self.w) ** 0.5
            * gamma**0.5
            * (2 / (gamma + 1)) ** (0.5 * (gamma + 1) / self.theta)
            * phi_2
        )
        self.c_f = self.get_cf(gamma, self.a_bar, self.tol)
        self.s_j_bar = 1 / (self.c_f * self.chi_0)
        if self.s_j_bar > self.chi_k:
            raise ValueError(
                "Achieving recoilless condition necessitates a larger throat area than could be fit into breech face."
            )
        self.s_j = self.s_j_bar * self.s

    def to_json(self) -> str:
        return json.dumps(
            {
                **json.loads(super().to_json()),
                "nozzle": asdict(self.nozzle),
            },
            ensure_ascii=False,
        )

    def f_p_bar(self, z: float, l_bar: float, eta: float, tau: float, psi: None | float = None) -> float:
        psi = psi if psi else self.f_psi_z(z)
        l_psi_bar = 1 - self.delta * ((1 - psi) / self.rho_p + self.alpha * (psi - eta))

        p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

        return p_bar

    def ode_t(
        self, t: float, z_t_l_v_eta_tau: tuple[float, float, float, float, float]
    ) -> tuple[float, float, float, float, float]:
        z, l_bar, v_bar, eta, tau = z_t_l_v_eta_tau
        psi = self.f_psi_z(z)
        d_psi = self.f_sigma_z(z)
        p_bar = self.f_p_bar(z, l_bar, eta, tau, psi)

        dz = (0.5 * self.theta / self.b) ** 0.5 * p_bar**self.n

        dl_bar = v_bar
        dv_bar = self.theta * 0.5 * p_bar

        d_eta = self.C_a * self.s_j_bar * p_bar * tau**-0.5
        d_tau = ((1 - tau) * (d_psi * dz) - 2 * v_bar * dv_bar - self.theta * tau * d_eta) / (psi - eta)

        return dz, dl_bar, dv_bar, d_eta, d_tau

    def ode_l(
        self, l_bar: float, t_z_v_eta_tau: tuple[float, float, float, float, float]
    ) -> tuple[float, float, float, float, float]:
        """length domain ode of internal ballistics
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition.
        in general, d/dl_bar = d/dt_bar * dt_bar/dl_bar
        """

        t, z, v_bar, eta, tau = t_z_v_eta_tau
        psi = self.f_psi_z(z)
        d_psi__d_z = self.f_sigma_z(z)
        p_bar = self.f_p_bar(z, l_bar, eta, tau, psi)
        dz = (0.5 * self.theta / self.b) ** 0.5 * p_bar**self.n / v_bar
        dv_bar = self.theta * 0.5 * p_bar / v_bar
        dt_bar = 1 / v_bar

        d_eta = self.C_a * self.s_j_bar * p_bar * tau**-0.5 * dt_bar
        d_tau = ((1 - tau) * (d_psi__d_z * dz) - 2 * v_bar * dv_bar - self.theta * tau * d_eta) / (psi - eta)

        return dt_bar, dz, dv_bar, d_eta, d_tau

    def ode_z(
        self, z: float, t_l_v_eta_tau: tuple[float, float, float, float, float]
    ) -> tuple[float, float, float, float, float]:
        t, l_bar, v_bar, eta, tau = t_l_v_eta_tau
        psi = self.f_psi_z(z)
        d_psi__d_z = self.f_sigma_z(z)
        p_bar = self.f_p_bar(z, l_bar, eta, tau, psi)

        dt_bar = (2 * self.b / self.theta) ** 0.5 * p_bar**-self.n
        dl_bar = v_bar * dt_bar
        dv_bar = 0.5 * self.theta * p_bar * dt_bar
        d_eta = self.C_a * self.s_j_bar * p_bar * tau**-0.5 * dt_bar
        d_tau = ((1 - tau) * d_psi__d_z - 2 * v_bar * dv_bar - self.theta * tau * d_eta) / (psi - eta)

        return dt_bar, dl_bar, dv_bar, d_eta, d_tau

    def integrate(
        self,
        step: int = 33,
        tol: float | None = None,
        dom: Domains = DOMAIN_TIME,
    ) -> RecoillessResult:
        tol = tol if tol is not None else self.tol

        bar_data = []
        t_scale = self.l_0 / self.v_j
        p_scale = self.f * self.delta
        l_g_bar = self.l_g / self.l_0
        z_0, z_b = self.z_0, self.z_b

        self.append_bar_data(bar_data, tag=POINT_START, t_bar=0, l_bar=0, z=z_0, v_bar=0, eta=0, tau=1)

        # Phase 1: Integrate to burnout or exit
        p_bar_max = 1e9 / p_scale  # 1 GPa

        def abort_condition(z, tlv_eta_tau, _):
            _, l_bar, v_bar, eta, tau = tlv_eta_tau
            p_bar = self.f_p_bar(z, l_bar, eta, tau)
            return l_bar > l_g_bar or p_bar > p_bar_max  # or p_bar < 0

        z_end, (t_bar_end, l_bar_end, v_bar_end, eta_end, tau_end), aborted = rkf(
            self.ode_z,
            (0, 0, 0, 0, 1),
            z_0,
            z_b,
            rel_tol=tol,
            abort_func=abort_condition,
        )

        # Check for excessive pressure
        p_bar_end = self.f_p_bar(z_end, l_bar_end, eta_end, tau_end)
        if p_bar_end > p_bar_max:
            raise ValueError(
                "excessive pressure encountered during integration (>1GPa mean). results cannot be expected to "
                + "be accurate due to gross violation of the applicable domain of Nobel-Abel equation-of-state."
            )

        # Determine if exit happened before burnout
        is_burn_out_contained = not (aborted and l_bar_end >= l_g_bar)

        if is_burn_out_contained:
            self.append_bar_data(
                bar_data,
                tag=POINT_BURNOUT,
                t_bar=t_bar_end,
                l_bar=l_bar_end,
                z=z_b,
                v_bar=v_bar_end,
                eta=eta_end,
                tau=tau_end,
            )
            self.logger.info("integrated to burnout point.")
        else:
            self.logger.warning("shot exited barrel before burnout.")

        # Integrate to actual exit point
        l_bar_exit, (t_bar_exit, z_exit, v_bar_exit, eta_exit, tau_exit), _ = rkf(
            self.ode_l, (t_bar_end, z_end, v_bar_end, eta_end, tau_end), l_bar_end, l_g_bar, rel_tol=tol
        )
        self.append_bar_data(
            bar_data,
            tag=POINT_EXIT,
            t_bar=t_bar_exit,
            l_bar=l_g_bar,
            z=z_exit,
            v_bar=v_bar_exit,
            eta=eta_exit,
            tau=tau_exit,
        )

        # Fracture point
        if z_b > 1.0 and z_exit >= 1.0:
            t_bar_f, l_bar_f, v_bar_f, eta_f, tau_f = rkf(self.ode_z, (0, 0, 0, 0, 1), z_0, 1, rel_tol=tol)[1]
            self.append_bar_data(
                bar_data,
                tag=POINT_FRACTURE,
                t_bar=t_bar_f,
                l_bar=l_bar_f,
                z=1.0,
                v_bar=v_bar_f,
                eta=eta_f,
                tau=tau_f,
            )

        # Peak finding
        def find_peak(g: Callable[[float], float], tag: RecoillessPeakPoints) -> None:
            t_bar_p = 0.5 * sum(gss(g, 0, t_bar_exit, x_tol=t_bar_exit * tol, find_min=False))
            z_p, l_bar_p, v_bar_p, eta_p, tau_p = self.g(t_bar_p, tag, tol)[1]
            self.append_bar_data(
                bar_data, tag=tag, t_bar=t_bar_p, l_bar=l_bar_p, z=z_p, v_bar=v_bar_p, eta=eta_p, tau=tau_p
            )

        for peak in [POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH, POINT_PEAK_STAG]:
            find_peak(lambda _t_bar: self.g(_t_bar, peak, tol)[0], peak)

        # Sampling
        if dom == DOMAIN_TIME:
            z_j, l_bar_j, v_bar_j, t_bar_j, eta_j, tau_j = z_0, 0.0, 0.0, 0.0, 0.0, 1.0
        else:  # length domain ODE requires starting from some point
            t_bar_j, (z_j, l_bar_j, v_bar_j, eta_j, tau_j), _ = rkf(
                self.ode_t, (z_0, 0.0, 0.0, 0.0, 1.0), 0, 0.5 * t_bar_exit, rel_tol=tol
            )

        for j in range(step):
            if dom == DOMAIN_TIME:
                t_bar_k = t_bar_exit / (step + 1) * (j + 1)
                z_j, l_bar_j, v_bar_j, eta_j, tau_j = rkf(
                    self.ode_t, (z_j, l_bar_j, v_bar_j, eta_j, tau_j), t_bar_j, t_bar_k, rel_tol=tol
                )[1]
                t_bar_j = t_bar_k
            else:
                l_bar_k = l_g_bar / (step + 1) * (j + 1)
                t_bar_j, z_j, v_bar_j, eta_j, tau_j = rkf(
                    self.ode_l, (t_bar_j, z_j, v_bar_j, eta_j, tau_j), l_bar_j, l_bar_k, rel_tol=tol
                )[1]
                l_bar_j = l_bar_k

            self.append_bar_data(
                bar_data, tag=SAMPLE, t_bar=t_bar_j, l_bar=l_bar_j, z=z_j, v_bar=v_bar_j, eta=eta_j, tau=tau_j
            )

        self.logger.info(f"sampled for {step} points.")

        # Data processing
        data, p_trace = [], []
        l_c = self.l_c
        trace_step = max(step, 1)

        for bar_data_line in bar_data:
            dtag, t_bar, l_bar, z, v_bar, p_bar, eta, tau = bar_data_line

            t = t_bar * t_scale
            l = l_bar * self.l_0
            psi = self.f_psi_z(z)
            v = v_bar * self.v_j
            p = p_bar * p_scale
            temp = tau * self.temp_v if self.temp_v else None

            ps, p0, pb, vb, stag = self.to_ps_p0_pb_vb_stag(l, v, p, tau, eta)

            p_line = []
            for i in range(trace_step):
                x = i / trace_step * (l + l_c)
                px = self.to_px(l, v, vb, ps, eta, x)
                p_line.append(PressureProbePoint(x, px))

            p_line.append(PressureProbePoint(l + l_c, ps))
            p_trace.append(PressureTraceEntry(dtag, temp, p_line))

            data.append(
                RecoillessTableEntry(
                    tag=dtag,
                    time=t,
                    travel=l,
                    burnup=psi,
                    velocity=v,
                    outflow_velocity=vb,
                    breech_pressure=pb,
                    stag_pressure=p0,
                    avg_pressure=p,
                    shot_pressure=ps,
                    temperature=temp,
                    outflow_fraction=eta,
                    rel_stag_point=stag,
                )
            )

        data, p_trace = zip(*sorted(zip(data, p_trace), key=lambda e: e[0].time))
        return RecoillessResult(self, data, p_trace)

    def g(
        self, t: float, tag: RecoillessPeakPoints, tol: float
    ) -> tuple[float, tuple[float, float, float, float, float]]:
        z, l_bar, v_bar, eta, tau = rkf(self.ode_t, (self.z_0, 0, 0, 0, 1), 0, t, rel_tol=tol)[1]
        p_bar = self.f_p_bar(z, l_bar, eta, tau)

        if tag == POINT_PEAK_AVG:
            return p_bar, (z, l_bar, v_bar, eta, tau)

        ps, p0, pb, _, _ = self.to_ps_p0_pb_vb_stag(
            l_bar * self.l_0, v_bar * self.v_j, p_bar * self.f * self.delta, tau, eta
        )

        ps_bar = ps / (self.f * self.delta)
        p0_bar = p0 / (self.f * self.delta)
        pb_bar = pb / (self.f * self.delta)
        if tag == POINT_PEAK_BREECH:
            return pb_bar, (z, l_bar, v_bar, eta, tau)
        elif tag == POINT_PEAK_STAG:
            return p0_bar, (z, l_bar, v_bar, eta, tau)
        elif tag == POINT_PEAK_SHOT:
            return ps_bar, (z, l_bar, v_bar, eta, tau)
        else:
            raise ValueError("tag not handled.")

    def abort_z(
        self,
        x: float,
        ys: tuple[float, float, float, float, float],
        _: tuple[float, tuple[float, float, float, float, float]],
        p_bar_max: float,
        l_g_bar: float,
    ) -> bool:

        z, (_, l_bar, v_bar, eta, tau) = x, ys
        p_bar = self.f_p_bar(z, l_bar, eta, tau)

        return any((l_bar > l_g_bar, p_bar > p_bar_max, v_bar <= 0, p_bar < 0))

    def append_bar_data(
        self,
        bar_data: list[tuple[str, float, float, float, float, float, float, float]],
        tag: Points,
        t_bar: float,
        l_bar: float,
        z: float,
        v_bar: float,
        eta: float,
        tau: float,
    ) -> None:
        bar_data.append((tag, t_bar, l_bar, z, v_bar, self.f_p_bar(z, l_bar, eta, tau), eta, tau))

    def to_ps_p0_pb_vb_stag(
        self, l: float, v: float, p: float, tau: float, eta: float
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
        y = self.w * eta
        m_dot = self.C_a * self.v_j * self.s_j * p / (self.f * tau**0.5)
        sb = self.s * self.chi_k
        vb = m_dot * (self.vol_0 + self.s * l) / (sb * (self.w - y))

        h_lim = 2 * self.phi_1 * self.m / (self.w * (1 - eta)) + 1
        h = h_lim if v == 0 else min(h_lim, vb / v)

        ps = p / (1 + (self.w - y) / (3 * self.phi_1 * self.m) * (1 - 0.5 * h))
        p0 = ps * (1 + (self.w - y) / (2 * self.phi_1 * self.m) * (1 + h) ** -1)
        pb = 0 if h == h_lim else ps * (1 + (self.w - y) / (2 * self.phi_1 * self.m) * (1 - h))
        stag = h / (1 + h)
        return ps, p0, pb, vb, stag

    def to_px(self, l: float, v: float, vb: float, ps: float, eta: float, x: float) -> float:
        """
        convert x, the physical displacement from breech bottom, to
        effective length in the equivalent gun.
        """
        y = self.w * eta
        l_k = self.l_0 / self.chi_k
        s_k = self.s * self.chi_k
        if x < l_k:
            z = x * s_k / (l_k * s_k + l * self.s)
        else:
            z = (l_k * s_k + (x - l_k) * self.s) / (l_k * s_k + l * self.s)

        h_lim = 2 * self.phi_1 * self.m / (self.w * (1 - eta)) + 1
        h = h_lim if v == 0 else min(h_lim, vb / v)

        px = ps * (
            1
            + self.w * (1 - eta) / (2 * self.phi_1 * self.m) * (1 + h) * (1 - z**2)
            - self.w * (1 - eta) / (self.phi_1 * self.m) * h * (1 - z)
        )
        return px

    @staticmethod
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

    def structure(
        self,
        recoilless_result: RecoillessResult,
        structural: Structural,
        step: int = 33,
        tol: float | None = None,
    ) -> None:
        tol = tol if tol is not None else self.tol
        step = max(step, 1)
        l_c = self.l_c
        l_g = self.l_g

        structural_material = structural.material
        if structural_material is None:
            raise ValueError("No structural material provided")

        sigma = structural_material.yield_strength
        s = self.s
        gamma = self.theta + 1
        r_b = 0.5 * self.caliber

        r_c = r_b * self.chi_k**0.5
        r_t = r_c * self.s_j_bar**0.5
        r_n = r_t * self.a_bar**0.5

        x_probes = (
            [i / step * l_c for i in range(step)]
            + [l_c * (1 - tol)]
            + [i / step * l_g + l_c for i in range(step)]
            + [l_g + l_c]
        )

        p_probes = [0.0 for _ in range(len(x_probes))]
        for recoillessTableEntry in recoilless_result.table_data:
            l = recoillessTableEntry.travel
            v = recoillessTableEntry.velocity
            vb = recoillessTableEntry.outflow_velocity
            ps = recoillessTableEntry.shot_pressure
            eta = recoillessTableEntry.outflow_fraction

            for i, x in enumerate(x_probes):
                if (x - l_c) <= l:
                    p_x = self.to_px(l, v, vb, ps, eta, x)
                    p_probes[i] = max(p_probes[i], p_x)
                else:
                    break

        for i in range(len(p_probes)):
            p_probes[i] *= structural.safety_factor

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
                p = self.get_pr(gamma, (r / r_t) ** 2, tol)[1] * p_entry
            else:
                k = (x + l_b) / l_b
                r = k * r_c + (1 - k) * r_t
                p = self.get_pr(gamma, (r / r_t) ** 2, tol)[0] * p_entry

            r_probes.append(r)
            neg_ps.append(p)
            neg_ss.append(pi * r**2)

        if structural.autofrettage:
            v_n, k_n, m_n = self.barrel_autofrettage(neg_x_probes, neg_ps, neg_ss, sigma)
        else:
            v_n, k_n, m_n = self.barrel_monoblock(neg_x_probes, neg_ps, neg_ss, sigma)

        hull = []
        for x, r, k, m in zip(neg_x_probes, r_probes, k_n, m_n):
            hull.append(OutlineEntry(x, r, r * k, r * m))

        nozzle_mass = v_n * structural_material.density

        i = step + 1
        x_c, p_c = x_probes[:i], p_probes[:i]
        x_b, p_b = x_probes[i:], p_probes[i:]

        if structural.autofrettage:
            v_c, k_c, m_c = self.barrel_autofrettage(x_c, p_c, [s * self.chi_k for _ in x_c], sigma)
            v_b, k_b, m_b = self.barrel_autofrettage(x_b, p_b, [s for _ in x_b], sigma)
        else:
            v_c, k_c, m_c = self.barrel_monoblock(x_c, p_c, [s * self.chi_k for _ in x_c], sigma)
            v_b, k_b, m_b = self.barrel_monoblock(x_b, p_b, [s for _ in x_b], sigma)

        v = v_c + v_b
        k_probes = k_c + k_b
        m_probes = m_c + m_b

        tube_mass = v * structural_material.density
        for x, k, m in zip(x_probes, k_probes, m_probes):
            if x < l_c:
                hull.append(OutlineEntry(x, r_c, k * r_c, m * r_c))
            else:
                hull.append(OutlineEntry(x, r_b, k * r_b, m * r_b))

        recoilless_result.outline = hull
        recoilless_result.tube_mass = tube_mass + nozzle_mass

        self.logger.info("conducted structural calculation.")

    @staticmethod
    def get_ar(gamma: float, pr: float) -> float:
        if pr <= 0 or pr >= 1:
            return inf
        return (
            (0.5 * (gamma + 1)) ** (1 / (gamma - 1))
            * pr ** (1 / gamma)
            * ((gamma + 1) / (gamma - 1) * (1 - pr ** ((gamma - 1) / gamma))) ** 0.5
        ) ** -1

    @staticmethod
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
            pr_sup, _ = dekker(lambda pr: Recoilless.get_ar(gamma, pr), 0, pr_c, y=ar, y_rel_tol=tol)
            pr_sub, _ = dekker(lambda pr: Recoilless.get_ar(gamma, pr), pr_c, 1, y=ar, y_rel_tol=tol)
            return pr_sub, pr_sup
