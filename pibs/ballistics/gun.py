from __future__ import annotations

import json
import logging
import math
from dataclasses import asdict, dataclass
from math import log
from typing import Callable, Literal

from . import (
    DOMAIN_TIME,
    POINT_BURNOUT,
    POINT_EXIT,
    POINT_FRACTURE,
    POINT_PEAK_AVG,
    POINT_PEAK_BREECH,
    POINT_PEAK_SHOT,
    POINT_START,
    SAMPLE,
    SOL_LAGRANGE,
    SOL_MAMONTOV,
    SOL_PIDDUCK,
    Domains,
    Points,
    Solutions,
)
from .base_gun import BaseGun
from .config import Environment, GunGeometry, PropellantLoad, Solver, Structural
from .generics import (
    GenericEntry,
    GenericResult,
    OutlineEntry,
    PressureProbePoint,
    PressureTraceEntry,
)
from .material import Material
from .num import dekker, gss, integrate, rkf


@dataclass
class GunResult(GenericResult):
    gun: Gun
    table_data: list[GunTableEntry]


@dataclass
class GunTableEntry(GenericEntry):
    pass


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

        i, _ = integrate(lambda x: f(om, x), 0, 1, tol)

        if k == 1:
            return i - 0.5 * wpm * math.exp(-om) / om
        else:
            return i - 0.5 * ((k - 1) / k) * wpm * ((1 - om) ** (k / (k - 1)) / om)

    omega = 0.5 * sum(dekker(f_omega, 0, 1, x_tol=tol))

    if k == 1:
        labda_1 = (math.exp(omega) - 1) / wpm
    else:
        labda_1 = ((1 - omega) ** (k / (1 - k)) - 1) / wpm

    i_u, _ = integrate(lambda x: g(omega, x), 0, 1, tol)
    i_l, _ = integrate(lambda x: f(omega, x), 0, 1, tol)
    labda_2 = i_u / i_l

    return labda_1, labda_2


class Gun(BaseGun):
    def __init__(
        self,
        geometry: GunGeometry,
        load: PropellantLoad,
        solver: Solver | None = None,
        environment: Environment | None = None,
        logger: logging.Logger | None = None,
    ):
        super().__init__(
            geometry=geometry,
            load=load,
            environment=environment,
            solver=solver,
            logger=logger,
        )

        self.sol = self.solver.solution_method

        if self.sol == SOL_LAGRANGE:
            self.labda_1, self.labda_2 = 1 / 2, 1 / 3
        elif self.sol == SOL_PIDDUCK:
            self.labda_1, self.labda_2 = pidduck(self.w / (self.phi_1 * self.m), self.theta + 1, self.tol)
        elif self.sol == SOL_MAMONTOV:
            self.labda_1, self.labda_2 = pidduck(self.w / (self.phi_1 * self.m), 1, self.tol)
        else:
            raise ValueError("Unknown Solution")

        labda = self.l_g / self.l_0
        cc = 1 - (1 - 1 / self.chi_k) * log(labda + 1) / labda

        self.phi = self.phi_1 + self.labda_2 * cc * self.w / self.m

        """
        见《枪炮内弹道学》（金，2014）p.70 式
        """
        self.b = (
            self.s**2
            * self.e_1**2
            / (self.f * self.phi * self.w * self.m * self.u_1**2)
            * (self.f * self.delta) ** (2 * (1 - self.n))
        )

        self.v_j = (2 * self.f * self.w / (self.theta * self.phi * self.m)) ** 0.5

    def to_json(self) -> str:
        return json.dumps(
            {**json.loads(super().to_json()), "sol": self.sol},
            ensure_ascii=False,
        )

    def f_p_bar(self, z: float, l_bar: float, v_bar: float) -> float:
        psi = self.f_psi_z(z)
        l_psi_bar = 1 - self.delta * ((1 - psi) / self.rho_p + (self.alpha * psi))
        p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)
        return max(p_bar, self.p_a_bar)

    def ode_t(self, _: float, zlv: tuple[float, float, float]) -> tuple[float, float, float]:
        z, l_bar, v_bar = zlv
        p_bar = self.f_p_bar(z, l_bar, v_bar)
        dz = (0.5 * self.theta / self.b) ** 0.5 * p_bar**self.n
        dl_bar = v_bar
        dv_bar = self.theta * 0.5 * (p_bar - self.func_p_ad_bar(v_bar))

        return dz, dl_bar, dv_bar

    def ode_l(self, l_bar: float, tzv: tuple[float, float, float]) -> tuple[float, float, float]:
        """length domain ode of internal ballistics
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition."""
        t_bar, z, v_bar = tzv

        p_bar = self.f_p_bar(z, l_bar, v_bar)

        dz = (0.5 * self.theta / self.b) ** 0.5 * p_bar**self.n / v_bar

        dv_bar = self.theta * 0.5 * (p_bar - self.func_p_ad_bar(v_bar)) / v_bar
        dt_bar = 1 / v_bar

        return dt_bar, dz, dv_bar

    def ode_z(self, z: float, tlv: tuple[float, float, float]) -> tuple[float, float, float]:
        t_bar, l_bar, v_bar = tlv
        p_bar = self.f_p_bar(z, l_bar, v_bar)

        dt_bar = (2 * self.b / self.theta) ** 0.5 * p_bar**-self.n
        dl_bar = v_bar * dt_bar
        dv_bar = 0.5 * self.theta * (p_bar - self.func_p_ad_bar(v_bar)) * dt_bar

        return dt_bar, dl_bar, dv_bar

    def get_temperature(self, psi: float, l: float, p: float) -> float | None:
        """
        given pressure and travel, return temperature
        using the Nobel-Abel EOS
        """
        if not self.temp_v:
            return None

        if psi:
            r = self.f / self.temp_v
            l_psi = self.l_0 * (1 - self.delta / self.rho_p - self.delta * (self.alpha * -1 / self.rho_p) * psi)
            return self.s * p * (l + l_psi) / (self.w * psi * r)
        else:
            return self.temp_v

    def dp_dz(self, z: float, l_bar: float, v_bar: float) -> float:
        psi = self.f_psi_z(z)
        p_bar = self.f_p_bar(z, l_bar, v_bar)

        dz = (0.5 * self.theta / self.b) ** 0.5 * p_bar**self.n

        l_psi_bar = 1 - self.delta * ((1 - psi) / self.rho_p + (self.alpha * psi))
        dp_bar = (
            (
                (1 + p_bar * self.delta * (self.alpha - 1 / self.rho_p)) * self.f_sigma_z(z) * dz
                - p_bar * v_bar * (1 + self.theta)
            )
            / (l_bar + l_psi_bar)
            / dz
        )

        return dp_bar

    def integrate(
        self,
        step: int = 33,
        tol: float | None = None,
        dom: Domains = DOMAIN_TIME,
    ) -> GunResult:
        tol = tol if tol is not None else self.tol

        bar_data = []
        t_scale, p_scale = self.l_0 / self.v_j, self.f * self.delta
        l_g_bar = self.l_g / self.l_0
        z_0, z_b = self.z_0, self.z_b

        self.append_bar_data(bar_data, tag=POINT_START, t_bar=0, l_bar=0, z=z_0, v_bar=0)

        p_bar_max = 1e9 / p_scale  # 1 GPa

        def abort_condition(z, tlv, _):
            t_bar, l_bar, v_bar = tlv
            p_bar = self.f_p_bar(z, l_bar, v_bar)
            return l_bar > l_g_bar or p_bar > p_bar_max

        z_end, (t_bar_end, l_bar_end, v_bar_end), aborted = rkf(
            self.ode_z,
            (0, 0, 0),
            z_0,
            z_b,
            rel_tol=tol,
            abort_func=abort_condition,
        )

        # check if integration exited due to excess pressure
        p_bar_end = self.f_p_bar(z_end, l_bar_end, v_bar_end)
        if p_bar_end > p_bar_max:
            raise ValueError(
                "excessive pressure encountered during integration (>1GPa mean). results cannot be expected to \
be accurate due to gross violation of the applicable domain of Nobel-Abel equation-of-state. "
            )

        # Determine if exit happened before burnout
        is_burn_out_contained = not (aborted and l_bar_end >= l_g_bar)

        if is_burn_out_contained:
            self.append_bar_data(bar_data, tag=POINT_BURNOUT, t_bar=t_bar_end, l_bar=l_bar_end, z=z_b, v_bar=v_bar_end)
            self.logger.info("integrated to burnout point.")
        else:
            self.logger.warning("shot exited barrel before burnout.")

        # integrate to actual exit point for both cases
        l_bar_exit, (t_bar_exit, z_exit, v_bar_exit), _ = rkf(
            self.ode_l,
            (t_bar_end, z_end, v_bar_end),
            l_bar_end,
            l_g_bar,
            rel_tol=tol,
        )
        # Populate exit point
        self.append_bar_data(bar_data, tag=POINT_EXIT, t_bar=t_bar_exit, l_bar=l_g_bar, z=z_exit, v_bar=v_bar_exit)

        # Populate fracture point entry at Z = 1
        if z_b > 1.0 and z_exit >= 1.0:
            t_bar_f, l_bar_f, v_bar_f = rkf(self.ode_z, (0, 0, 0), z_0, 1, rel_tol=tol)[1]
            self.append_bar_data(bar_data, tag=POINT_FRACTURE, t_bar=t_bar_f, l_bar=l_bar_f, z=1, v_bar=v_bar_f)

        def find_peak(g: Callable[[float], float], tag: Points) -> None:
            t_bar_p = 0.5 * sum(gss(g, 0, t_bar_exit, x_tol=t_bar_exit * tol, find_min=False))
            z_p, l_bar_p, v_bar_p = self.g(t_bar_p, tag, tol)[1]
            self.append_bar_data(bar_data, tag=tag, t_bar=t_bar_p, l_bar=l_bar_p, z=z_p, v_bar=v_bar_p)

        # Find peak pressure in time domain
        for peak in [POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH]:
            find_peak(lambda _t_bar: self.g(_t_bar, peak, z_0)[0], peak)

        # Sampling
        if dom == DOMAIN_TIME:
            z_j, l_bar_j, v_bar_j, t_bar_j = z_0, 0, 0, 0
        else:  # length domain ODE requires starting from some point
            t_bar_j, (z_j, l_bar_j, v_bar_j), _ = rkf(self.ode_t, (z_0, 0, 0), 0, 0.5 * t_bar_exit, rel_tol=tol)

        for j in range(step):
            if dom == DOMAIN_TIME:
                t_bar_k = t_bar_exit / (step + 1) * (j + 1)
                z_j, l_bar_j, v_bar_j = rkf(self.ode_t, (z_j, l_bar_j, v_bar_j), t_bar_j, t_bar_k, rel_tol=tol)[1]
                t_bar_j = t_bar_k
            else:
                l_bar_k = l_g_bar / (step + 1) * (j + 1)
                t_bar_j, z_j, v_bar_j = rkf(self.ode_l, (t_bar_j, z_j, v_bar_j), l_bar_j, l_bar_k, rel_tol=tol)[1]
                l_bar_j = l_bar_k

            self.append_bar_data(bar_data, tag=SAMPLE, t_bar=t_bar_j, l_bar=l_bar_j, z=z_j, v_bar=v_bar_j)

        self.logger.info(f"sampled for {step} points.")

        # Data processing
        data = []
        p_trace = []
        trace_steps = max(step, 1)

        for bar_data_line in bar_data:
            dtag, t_bar, l_bar, z, v_bar, p_bar = bar_data_line

            t = t_bar * t_scale
            l = l_bar * self.l_0
            psi = self.f_psi_z(z)
            v = v_bar * self.v_j
            p = p_bar * p_scale
            ps, pb = self.to_ps_pb(l, p)
            temp = self.get_temperature(psi, l, p)

            p_line = []
            for i in range(trace_steps):
                x = i / trace_steps * (l + self.l_c)
                p_x, _ = self.to_px_u(l, ps, pb, v, x)
                p_line.append(PressureProbePoint(x, p_x))

            p_line.append(PressureProbePoint(l + self.l_c, ps))
            p_trace.append(PressureTraceEntry(dtag, temp, p_line))

            data.append(
                GunTableEntry(
                    tag=dtag,
                    time=t,
                    travel=l,
                    burnup=psi,
                    velocity=v,
                    breech_pressure=pb,
                    avg_pressure=p,
                    shot_pressure=ps,
                    temperature=temp,
                )
            )

        data, p_trace = zip(*sorted(zip(data, p_trace), key=lambda e: e[0].time))
        return GunResult(self, data, p_trace)

    def append_bar_data(
        self,
        bar_data: list[tuple[str, float, float, float, float, float]],
        tag: str,
        t_bar: float,
        l_bar: float,
        z: float,
        v_bar: float,
    ) -> None:
        bar_data.append((tag, t_bar, l_bar, z, v_bar, self.f_p_bar(z, l_bar, v_bar)))

    def g(self, t_bar: float, tag: str, tol: float) -> tuple[float, tuple[float, float, float]]:
        z, l_bar, v_bar = rkf(self.ode_t, (self.z_0, 0, 0), 0, t_bar, rel_tol=tol)[1]
        p_bar = self.f_p_bar(z, l_bar, v_bar)
        if tag == POINT_PEAK_AVG:
            return p_bar, (z, l_bar, v_bar)
        else:
            ps_bar, pb_bar = self.to_ps_pb(l_bar * self.l_0, p_bar)
            if tag == POINT_PEAK_SHOT:
                return ps_bar, (z, l_bar, v_bar)
            elif tag == POINT_PEAK_BREECH:
                return pb_bar, (z, l_bar, v_bar)

        raise ValueError(f"tag {tag} not handled.")

    def abort_z(self, z: float, tlv: tuple[float, float, float], _, p_bar_max: float, l_g_bar: float) -> bool:
        t_bar, l_bar, v_bar = tlv
        p_bar = self.f_p_bar(z, l_bar, v_bar)

        return l_bar > l_g_bar or p_bar > p_bar_max or v_bar < 0

    def to_ps_pb(self, l: float, p: float) -> tuple[float, float]:
        """
        Convert average chamber pressure at certain travel to
        shot base pressure, and breech face pressure

        l: travel of the projectile
        p: average pressure

        Ps: pressure at shot
        Pb: pressure at breech
        """
        labda_g = l / self.l_0
        labda_1_prime = self.labda_1 * (1 / self.chi_k + labda_g) / (1 + labda_g)
        labda_2_prime = self.labda_2 * (1 / self.chi_k + labda_g) / (1 + labda_g)

        factor_s = 1 + labda_2_prime * (self.w / (self.phi_1 * self.m))
        factor_b = (self.phi_1 * self.m + labda_2_prime * self.w) / (self.phi_1 * self.m + labda_1_prime * self.w)

        return p / factor_s, p / factor_b

    def to_px_u(self, l: float, p_s: float, p_b: float, v: float, x: float) -> tuple[float, float]:
        """
        Convert the average chamber to pressure and gas flow speed
        at arbitrary point x for projectile travel of l and average pressure
        of p, **assuming the Lagrangian distribution**.

        Note that with the current state of research, only characteristic point
        values are available for other distributions, use to_ps_pb() instead for that.

        l: projectile travel
        p_s: pressure of shot
        p_b: pressure of breech
        x: probe point, start from the breech face.
        """
        r = self.chi_k * x if x < self.l_c else (x - self.l_c) + self.l_0
        k = (r / (self.l_0 + l)) ** 2
        p_x = p_s * k + p_b * (1 - k)

        if x < self.l_c:
            u = x * v / (self.l_0 + l)
        else:
            u = (x + (self.chi_k - 1) * self.l_c) * v / (self.l_0 + l)

        return p_x, u

    def structure(
        self,
        gun_result: GunResult,
        structural: Structural,
        step: int = 33,
        tol: float | None = None,
    ) -> None:
        tol = tol if tol is not None else self.tol
        step = max(step, 1)

        structural_material = structural.material
        if structural_material is None:
            raise ValueError("No structural material provided")

        structural_safety_factor = structural.safety_factor

        r_b = 0.5 * self.caliber
        r_c = r_b * self.chi_k**0.5
        x_probes = (
            [i / step * self.l_c for i in range(step)]
            + [self.l_c * (1 - tol)]
            + [i / step * self.l_g + self.l_c for i in range(step)]
            + [self.l_g + self.l_c]
        )
        p_probes = [0.0 for _ in range(len(x_probes))]

        for gun_table_entry in gun_result.table_data:
            l = gun_table_entry.travel
            v = gun_table_entry.velocity
            p_s = gun_table_entry.shot_pressure
            p_b = gun_table_entry.breech_pressure
            for i, x in enumerate(x_probes):
                if (x - self.l_c) <= l:
                    p_x, _ = self.to_px_u(l, p_s, p_b, v, x)
                    p_probes[i] = max(p_probes[i], p_x)
                else:
                    break

        for i in range(len(p_probes)):
            p_probes[i] *= structural_safety_factor

        i = step + 1
        x_c, p_c = x_probes[:i], p_probes[:i]
        x_b, p_b = x_probes[i:], p_probes[i:]

        if structural.autofrettage:
            v_c, k_c, m_c = self.barrel_autofrettage(
                x_c, p_c, [self.s * self.chi_k for _ in x_c], structural_material.yield_strength
            )
            v_b, k_b, m_b = self.barrel_autofrettage(
                x_b, p_b, [self.s for _ in x_b], structural_material.yield_strength
            )
        else:
            v_c, k_c, m_c = self.barrel_monoblock(
                x_c, p_c, [self.s * self.chi_k for _ in x_c], structural_material.yield_strength
            )
            v_b, k_b, m_b = self.barrel_monoblock(x_b, p_b, [self.s for _ in x_b], structural_material.yield_strength)

        v = v_c + v_b
        k_probes = k_c + k_b
        m_probes = m_c + m_b

        hull = []
        for x, k, m in zip(x_probes, k_probes, m_probes):
            if x < self.l_c:
                hull.append(OutlineEntry(x, r_c, k * r_c, m * r_c))
            else:
                hull.append(OutlineEntry(x, r_b, k * r_b, m * r_b))

        gun_result.outline = hull
        gun_result.tube_mass = v * structural_material.density

        self.logger.info("conducted structural calculation.")
