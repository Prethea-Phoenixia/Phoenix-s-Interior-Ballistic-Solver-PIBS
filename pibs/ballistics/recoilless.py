from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from math import inf, pi, tan
from typing import Callable


from . import (
    COMPUTE,
    DOMAIN_LEN,
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
)
from .gun import Gun
from .num import dekker, gss, rkf
from .prop import Propellant
from .material import Material
from . import JSONable

logger = logging.getLogger(__name__)


from .generics import (
    DelegatesPropellant,
    GenericEntry,
    GenericResult,
    OutlineEntry,
    PressureProbePoint,
    PressureTraceEntry,
)


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


class Recoilless(DelegatesPropellant, JSONable):
    def __init__(
        self,
        caliber: float,
        shot_mass: float,
        propellant: Propellant,
        web: float,  # 2e_1
        charge_mass: float,
        chamber_volume: float,
        start_pressure: float,
        length_gun: float,
        chambrage: float,
        nozzle_expansion: float,
        tol: float,
        drag_coefficient: float = 0.0,
        nozzle_efficiency: float = 0.92,
        ambient_density: float = 1.204,
        ambient_pressure: float = 101.325e3,
        ambient_adb_index: float = 1.4,
        **_,
    ):
        super().__init__(propellant=propellant)

        self.propellant = propellant
        self.caliber = caliber
        self.s = (caliber / 2) ** 2 * pi
        self.m = shot_mass
        self.w = charge_mass
        self.vol_0 = chamber_volume
        self.p_0 = start_pressure
        self.l_g = length_gun
        self.chi_0 = nozzle_efficiency
        self.a_bar = nozzle_expansion
        self.e_1 = web / 2
        self.tol = tol
        self.chi_k = chambrage  # ration of l_0 / l_chamber
        self.l_0 = self.vol_0 / self.s
        self.l_c = self.l_0 / self.chi_k

        self.ambient_density = ambient_density
        self.ambient_pressure = ambient_pressure
        self.ambient_adb_index = ambient_adb_index

        self.Delta = self.w / self.vol_0
        self.phi_1 = 1 / (1 - drag_coefficient)  # drag work coefficient
        self.phi = self.phi_1 + self.w / (3 * self.m)

        self.v_j = (2 * self.f * self.w / (self.theta * self.phi * self.m)) ** 0.5
        if self.p_0 == 0:
            raise NotImplementedError(
                "Current implementation use exponential burn rate and does not"
                + " allow for solving the case with 0 shot start pressure."
            )
        else:
            self.psi_0 = (1 / self.Delta - 1 / self.rho_p) / (self.f / self.p_0 + self.alpha - 1 / self.rho_p)
            if self.psi_0 <= 0:
                raise ValueError(
                    "Initial burnup fraction is solved to be negative."
                    + " This indicate an excessively high load density for the"
                    + " start-pressure target."
                )
            elif self.psi_0 >= 1:
                raise ValueError(
                    "Initial burnup fraction is solved to be greater than unity."
                    + " This indicate an excessively low loading density for the"
                    + " start-pressure target."
                )

        self.B = (
            self.s**2
            * self.e_1**2
            / (self.f * self.phi * self.w * self.m * self.u_1**2)
            * (self.f * self.Delta) ** (2 * (1 - self.n))
        )

        # additional calculation for recoilless weapons:
        gamma = self.theta + 1
        """
        Enforce the "recoilless condition" by setting the size of the
        throat.
        """

        phi_2 = 1
        self.c_a = (
            (0.5 * self.theta * self.phi * self.m / self.w) ** 0.5
            * gamma**0.5
            * (2 / (gamma + 1)) ** (0.5 * (gamma + 1) / self.theta)
            * phi_2
        )

        ambient_pressure = max(ambient_pressure, 1)
        ambient_adb_index = max(ambient_adb_index, 1)

        self.psi_0 = (1 / self.Delta - 1 / self.rho_p) / (self.f / self.p_0 + self.alpha - 1 / self.rho_p)
        if self.psi_0 <= 0:
            raise ValueError(
                "Initial burnup fraction is solved to be negative."
                + " This indicate an excessively high load density for the"
                + " start-pressure target."
            )
        elif self.psi_0 >= 1:
            raise ValueError(
                "Initial burnup fraction is solved to be greater than unity."
                + " This indicate an excessively low loading density for the"
                + " start-pressure target."
            )
        self.z_0, _ = dekker(
            lambda _z: self.propellant.f_psi_z(_z), 0, 1, y=self.psi_0, x_tol=tol, y_rel_tol=tol, y_abs_tol=tol**2
        )

        if ambient_density < 0:
            raise ValueError("Invalid ambient condition")

        gamma = self.theta + 1
        """
        self.S_j_bar = self.get_cf(gamma, 1, tol) / (
            self.get_cf(gamma, self.A_bar, tol) * self.chi_0
        )  # = S_j/S
        """
        self.c_f = self.get_cf(gamma, self.a_bar, tol)
        self.s_j_bar = 1 / (self.c_f * self.chi_0)
        if self.s_j_bar > self.chi_k:
            raise ValueError(
                "Achieving recoilless condition necessitates a larger throat area than could be fit into breech face."
            )
        self.s_j = self.s_j_bar * self.s

        self.p_a_bar = ambient_pressure / (self.f * self.Delta)
        if ambient_density != 0:
            self.c_a_bar = (ambient_adb_index * ambient_pressure / ambient_density) ** 0.5 / self.v_j
        else:
            self.c_a_bar = 0

        self.k_1 = ambient_adb_index

    def to_json(self) -> str:
        return json.dumps(
            {
                "caliber": self.caliber,
                "shot_mass": self.m,
                "propellant": json.loads(self.propellant.to_json()),
                "web": self.e_1 * 2,
                "charge_mass": self.w,
                "chamber_volume": self.vol_0,
                "start_pressure": self.p_0,
                "length_gun": self.l_g,
                "chambrage": self.chi_k,
                "tol": self.tol,
                "drag_coefficient": 1 - 1 / self.phi_1,
                "nozzle_expansion": self.a_bar,
                "nozzle_efficiency": self.chi_0,
                "ambient_density": self.ambient_density,
                "ambient_pressure": self.ambient_pressure,
                "ambient_adb_index": self.ambient_adb_index,
            },
            ensure_ascii=False,
        )

    @classmethod
    def from_json(cls, json_dict: dict) -> Recoilless:
        deserialized_dict = {}
        for key, value in json_dict.items():
            if key == "propellant":
                value = Propellant.from_json(value)
            deserialized_dict[key] = value

        return cls(**deserialized_dict)

    def f_p_bar(self, z: float, l_bar: float, eta: float, tau: float, psi: None | float = None) -> float:
        psi = psi if psi else self.f_psi_z(z)
        l_psi_bar = 1 - self.Delta * ((1 - psi) / self.rho_p + self.alpha * (psi - eta))

        p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

        return max(p_bar, self.p_a_bar)

    def f_p_d_bar(self, v_bar: float) -> float:
        if self.c_a_bar and v_bar > 0:
            v_r = v_bar / self.c_a_bar
            return (
                0.25 * self.k_1 * (self.k_1 + 1) * v_r**2
                + self.k_1 * v_r * (1 + (0.25 * (self.k_1 + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a_bar
        else:
            return 0.0

    def ode_t(
        self, _: float, zlvetatau: tuple[float, float, float, float, float], __: float
    ) -> tuple[float, float, float, float, float]:
        z, l_bar, v_bar, eta, tau = zlvetatau
        psi = self.f_psi_z(z)
        dpsi = self.f_sigma_z(z)  # dpsi/dZ
        p_bar = self.f_p_bar(z, l_bar, eta, tau, psi)

        p_d_bar = self.f_p_d_bar(v_bar)
        dz = (0.5 * self.theta / self.B) ** 0.5 * p_bar**self.n

        dl_bar = v_bar
        dv_bar = self.theta * 0.5 * (p_bar - p_d_bar)

        deta = self.c_a * self.s_j_bar * p_bar * tau**-0.5  # deta / dt_bar
        dtau = ((1 - tau) * (dpsi * dz) - 2 * v_bar * dv_bar - self.theta * tau * deta) / (psi - eta)  # dtau/dt_bar

        return dz, dl_bar, dv_bar, deta, dtau

    def ode_l(
        self, l_bar: float, tzvetatau: tuple[float, float, float, float, float], _: float
    ) -> tuple[float, float, float, float, float]:
        """length domain ode of internal ballistics
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition.

        in general, d/dl_bar = d/dt_bar * dt_bar/dl_bar

        """
        t, z, v_bar, eta, tau = tzvetatau

        psi = self.f_psi_z(z)
        dpsi = self.f_sigma_z(z)  # dpsi/dZ
        p_bar = self.f_p_bar(z, l_bar, eta, tau, psi)
        p_d_bar = self.f_p_d_bar(v_bar)

        dz = (0.5 * self.theta / self.B) ** 0.5 * p_bar**self.n / v_bar
        dv_bar = self.theta * 0.5 * (p_bar - p_d_bar) / v_bar
        dt_bar = 1 / v_bar  # dt_bar / dl_bar

        deta = self.c_a * self.s_j_bar * p_bar * tau**-0.5 * dt_bar  # deta / dl_bar

        dtau = ((1 - tau) * (dpsi * dz) - 2 * v_bar * dv_bar - self.theta * tau * deta) / (psi - eta)  # dtau/dl_bar

        return dt_bar, dz, dv_bar, deta, dtau

    def ode_z(
        self, z: float, tlvetatau: tuple[float, float, float, float, float], _: float
    ) -> tuple[float, float, float, float, float]:
        t, l_bar, v_bar, eta, tau = tlvetatau
        psi = self.f_psi_z(z)
        dpsi = self.f_sigma_z(z)  # dpsi/dZ
        p_bar = self.f_p_bar(z, l_bar, eta, tau, psi)
        p_d_bar = self.f_p_d_bar(v_bar)

        dt_bar = (2 * self.B / self.theta) ** 0.5 * p_bar**-self.n
        dl_bar = v_bar * dt_bar
        dv_bar = 0.5 * self.theta * (p_bar - p_d_bar) * dt_bar
        deta = self.c_a * self.s_j_bar * p_bar * tau**-0.5 * dt_bar
        dtau = ((1 - tau) * dpsi - 2 * v_bar * dv_bar - self.theta * tau * deta) / (psi - eta)

        return dt_bar, dl_bar, dv_bar, deta, dtau

    def integrate(
        self,
        step: int = 10,
        tol: float = 1e-5,
        dom: Domains = DOMAIN_TIME,
        **_,
    ) -> RecoillessResult:
        """
        Runs a full numerical solution for the gun in the specified domain sampled
        evenly at specified number of steps, using a scaled numerical tolerance as
        specified.

        tolerance is meant to be interpreted as the maximum relative deviation each
        component is allowed to have, at each step of integration point.

        Through significant trials and errors, it was determined that for this particular
        system, the error due to compounding does not appear to be significant,
        usually on the order of 1e-16 - 1e-14 as compared to much larger for component
        errors.
        """

        record, bar_data = [], []

        t_scale = self.l_0 / self.v_j
        p_scale = self.f * self.Delta

        l_g_bar = self.l_g / self.l_0
        p_bar_0 = self.p_0 / p_scale

        self.append_bar_data(bar_data, tag=POINT_START, t_bar=0, l_bar=0, z=self.z_0, v_bar=0, eta=0, tau=1)

        record.append((0, (0, self.psi_0, 0, p_bar_0, 0, 1)))

        z_i = self.z_0
        z_j = self.z_b
        n = 1
        delta_z = self.z_b - self.z_0

        t_bar_i, l_bar_i, v_bar_i, p_bar_i, eta_i, tau_i = 0, 0, 0, p_bar_0, 0, 1

        is_burn_out_contained = True

        ztlvet_record = [(self.z_0, (0, 0, 0, 0, 1))]
        p_max = 1e9
        p_bar_max = p_max / p_scale

        while z_i < self.z_b:  # terminates if burnout is achieved
            ztlvet_record_i = []
            if z_j == z_i:
                raise ValueError("Numerical accuracy exhausted in search of exit/burnout point.")
            try:
                if z_j > self.z_b:
                    z_j = self.z_b

                z, (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j) = rkf(
                    self.ode_z,
                    (t_bar_i, l_bar_i, v_bar_i, eta_i, tau_i),
                    z_i,
                    z_j,
                    rel_tol=tol,
                    abs_tol=tol**2,
                    abort_func=lambda _x, _ys, _records: self.abort_z(
                        _x, _ys, _records, p_bar_max=p_bar_max, l_g_bar=l_g_bar
                    ),
                    record=ztlvet_record_i,
                )

                p_bar_j = self.f_p_bar(z_j, l_bar_j, eta_j, tau_j)

            except ValueError as e:
                ztlvet_record.extend(ztlvet_record_i)
                z, (t_bar, l_bar, v_bar, eta, tau) = ztlvet_record[-1]
                dt_bar, dl_bar, dv_bar, deta, dtau = self.ode_z(z, (t_bar, l_bar, v_bar, eta, tau), 0)

                if all((dt_bar > 0, dl_bar > 0, dv_bar < 0)):
                    raise ValueError(
                        "Extremely low propulsive effort exerted on shot,"
                        + " impossible to integrate down to numerical precision.\n"
                        + "Shot last calculated at {:.0f} mm with velocity {:.0f} mm/s after {:.0f} ms\n".format(
                            l_bar * self.l_0 * 1e3,
                            v_bar * self.v_j * 1e3,
                            t_bar * t_scale * 1e3,
                        )
                    )
                else:
                    raise e  # unknown issues

            if l_bar_j >= l_g_bar:
                if abs(l_bar_i - l_g_bar) / l_g_bar > tol or l_bar_i == 0:
                    n *= 2
                    z_j = z_i + delta_z / n
                else:
                    is_burn_out_contained = False
                    break  # l_bar_i is solved to within a tol of l_bar_g

            else:
                ztlvet_record.extend(ztlvet_record_i)
                if p_bar_j > p_bar_max:
                    raise ValueError(
                        "Nobel-Abel EoS is generally accurate enough below 600MPa. However,"
                        + " Unreasonably high pressure (>{:.0f} MPa) was encountered.".format(p_max / 1e6)
                    )

                if v_bar_j <= 0:
                    z, (t_bar, l_bar, v_bar, eta, tau) = ztlvet_record[-1]

                    raise ValueError(
                        "Squib load condition detected: Shot stopped in bore.\n"
                        + "Shot is last calculated at {:.0f} mm at {:.0f} mm/s after {:.0f} ms".format(
                            l_bar * self.l_0 * 1e3, v_bar * self.v_j * 1e3, t_bar * t_scale * 1e3
                        )
                    )

                if any(v < 0 for v in (t_bar_j, l_bar_j, p_bar_j)):
                    raise ValueError(
                        "Numerical Integration diverged: negative"
                        + " values encountered in results.\n"
                        + "{:.0f} ms, {:.0f} mm, {:.0f} m/s, {:.0f} MPa".format(
                            t_bar_j * t_scale * 1e3,
                            l_bar_j * self.l_0 * 1e3,
                            v_bar_j * self.v_j,
                            p_bar_j * p_scale * 1e-6,
                        )
                    )  # this will catch any case where t, l, p are negative

                (t_bar_i, l_bar_i, v_bar_i, eta_i, tau_i) = (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j)
                z_i = z_j
                """
                this way the group of values denoted by _i is always updated
                as a group.
                """
                z_j += delta_z / n

        if t_bar_i == 0:
            raise ValueError("burnout point found to be at the origin.")

        record.extend(
            (t_bar, (l_bar, self.f_psi_z(z), v_bar, self.f_p_bar(z, l_bar, eta, tau), eta, tau))
            for (z, (t_bar, l_bar, v_bar, eta, tau)) in ztlvet_record
        )

        if is_burn_out_contained:
            logger.info("integrated to burnout point.")
        else:
            logger.warning("shot exited barrel before burnout.")

        ltzvet_record = []
        t_bar_e, z_e, v_bar_e, eta_e, tau_e = rkf(
            self.ode_l,
            (t_bar_i, z_i, v_bar_i, eta_i, tau_i),
            l_bar_i,
            l_g_bar,
            rel_tol=tol,
            abs_tol=tol**2,
            record=ltzvet_record,
        )[1]

        record.extend(
            (t_bar, (l_bar, self.f_psi_z(z), v_bar, self.f_p_bar(z, l_bar, eta, tau), eta, tau))
            for (l_bar, (t_bar, z, v_bar, eta, tau)) in ltzvet_record
        )

        self.append_bar_data(
            bar_data, tag=POINT_EXIT, t_bar=t_bar_e, l_bar=l_g_bar, z=z_e, v_bar=v_bar_e, eta=eta_e, tau=tau_e
        )

        t_bar_f = None
        if self.z_b > 1.0 and z_e >= 1.0:  # fracture point exist and is contained
            """
            Subscript f indicate fracture condition
            ODE w.r.t Z is integrated from Z_0 to 1, from onset of projectile
            movement to charge fracture
            """
            t_bar_f, l_bar_f, v_bar_f, eta_f, tau_f = rkf(
                self.ode_z, (0, 0, 0, 0, 1), self.z_0, 1, rel_tol=tol, abs_tol=tol**2
            )[1]

            self.append_bar_data(
                bar_data,
                tag=POINT_FRACTURE,
                t_bar=t_bar_f,
                l_bar=l_bar_f,
                z=1,
                eta=eta_f,
                tau=tau_f,
                v_bar=v_bar_f,
            )

        t_bar_b = None
        if is_burn_out_contained:
            """
            Subscript b indicated burnout condition
            ODE w.r.t Z is integrated from Z_0 to z_b, from onset of projectile
            movement to charge burnout.
            """

            t_bar_b, l_bar_b, v_bar_b, eta_b, tau_b = rkf(
                self.ode_z, (0, 0, 0, 0, 1), self.z_0, self.z_b, rel_tol=tol, abs_tol=tol**2
            )[1]

            self.append_bar_data(
                bar_data,
                tag=POINT_BURNOUT,
                t_bar=t_bar_b,
                l_bar=l_bar_b,
                z=self.z_b,
                v_bar=v_bar_b,
                eta=eta_b,
                tau=tau_b,
            )

        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the
        peak pressure point. As well as, not having starting issues.

        we hereby simply solve p golden section searching it from origin
        to point e, i.e. inside the barrel.
        """

        def find_peak(g: Callable[[float], float], tag: Points) -> None:
            """
            tolerance is specified a bit differently for gold section search
            GSS tol is the length between the upper bound and lower bound
            of the maxima/minima, thus by our definition of tolerance (one-sided),
            we take the median value.
            """
            t_bar_tol = tol * min(_t_bar for _t_bar in (t_bar_e, t_bar_b, t_bar_f) if _t_bar)
            t_bar_p = 0.5 * sum(gss(g, 0, t_bar_e if t_bar_b is None else t_bar_b, x_tol=t_bar_tol, find_min=False))

            z_p, l_bar_p, v_bar_p, eta_p, tau_p = self.g(t_bar_p, tag, tol)[1]
            self.append_bar_data(
                bar_data, tag=tag, t_bar=t_bar_p, l_bar=l_bar_p, z=z_p, v_bar=v_bar_p, eta=eta_p, tau=tau_p
            )

        for i, peak in enumerate([POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH, POINT_PEAK_STAG]):
            find_peak(lambda _t_bar: self.g(_t_bar, peak, tol)[0], peak)

        """
        populate data for output purposes
        """
        if dom == DOMAIN_TIME:
            z_j, l_bar_j, v_bar_j, t_bar_j, eta_j, tau_j = self.z_0, 0, 0, 0, 0, 1
            for j in range(step):
                t_bar_k = t_bar_e / (step + 1) * (j + 1)
                z_j, l_bar_j, v_bar_j, eta_j, tau_j = rkf(
                    self.ode_t, (z_j, l_bar_j, v_bar_j, eta_j, tau_j), t_bar_j, t_bar_k, rel_tol=tol, abs_tol=tol**2
                )[1]
                t_bar_j = t_bar_k

                self.append_bar_data(
                    bar_data, tag=SAMPLE, t_bar=t_bar_j, l_bar=l_bar_j, z=z_j, v_bar=v_bar_j, eta=eta_j, tau=tau_j
                )

        elif dom == DOMAIN_LEN:
            """
            Due to two issues, i.e. 1.the length domain ODE
            cannot be integrated from the origin point, and 2.the
            correct behaviour can only be expected when starting from
            a point with active burning else dZ flat lines.
            we do another Z domain integration to seed the initial values
            to a point where ongoing burning is guaranteed.
            (in the case of gun barrel length >= burn length, the group
             of value by subscript i will not guarantee burning is still
             ongoing).
            """
            t_bar_j = 0.5 * t_bar_i
            z_j, l_bar_j, v_bar_j, eta_j, tau_j = rkf(
                self.ode_t, (self.z_0, 0, 0, 0, 1), 0, t_bar_j, rel_tol=tol, abs_tol=tol**2
            )[1]

            for j in range(step):
                l_bar_k = l_g_bar / (step + 1) * (j + 1)

                t_bar_j, z_j, v_bar_j, eta_j, tau_j = rkf(
                    self.ode_l, (t_bar_j, z_j, v_bar_j, eta_j, tau_j), l_bar_j, l_bar_k, rel_tol=tol, abs_tol=tol**2
                )[1]

                l_bar_j = l_bar_k

                self.append_bar_data(
                    bar_data, tag=SAMPLE, t_bar=t_bar_j, l_bar=l_bar_j, z=z_j, v_bar=v_bar_j, eta=eta_j, tau=tau_j
                )

        logger.info(f"sampled for {step} points.")

        """
        sort the data points
        """

        data, p_trace = [], []
        l_c = self.l_c

        trace_step = max(step, 1)

        for bar_dataLine in bar_data:
            dtag, t_bar, l_bar, z, v_bar, p_bar, eta, tau = bar_dataLine

            t = t_bar * t_scale
            l = l_bar * self.l_0
            psi = self.f_psi_z(z)
            v = v_bar * self.v_j
            p = p_bar * p_scale
            temp = tau * self.temp_v if self.temp_v else None

            ps, p0, pb, vb, stag = self.to_ps_p0_pb_vb_stag(l, v, p, tau, eta)

            p_line = []
            for i in range(trace_step):  # 0....step -1
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

        for t_bar, (l_bar, psi, v_bar, p_bar, eta, tau) in record:
            t = t_bar * t_scale
            if t in [tableEntry.time for tableEntry in data]:
                continue
            l, v, p, temp = (
                l_bar * self.l_0,
                v_bar * self.v_j,
                p_bar * p_scale,
                tau * self.temp_v if self.temp_v else None,
            )

            ps, p0, pb, vb, stag = self.to_ps_p0_pb_vb_stag(l, v, p, tau, eta)

            p_line = []
            for i in range(trace_step):  # 0....step -1
                x = i / trace_step * (l + l_c)
                px = self.to_px(l, v, vb, ps, eta, x)
                pp = PressureProbePoint(x, px)
                p_line.append(pp)

            p_line.append(PressureProbePoint(l + l_c, ps))
            p_trace.append(PressureTraceEntry(COMPUTE, temp, p_line))

            data.append(
                RecoillessTableEntry(
                    tag=COMPUTE,
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

        data, p_trace = zip(*sorted(zip(data, p_trace), key=lambda entries: entries[0].time))

        recoilless_result = RecoillessResult(self, data, p_trace)

        return recoilless_result

    def g(self, t: float, tag: Points, tol: float) -> tuple[float, tuple[float, float, float, float, float]]:
        z, l_bar, v_bar, eta, tau = rkf(self.ode_t, (self.z_0, 0, 0, 0, 1), 0, t, rel_tol=tol, abs_tol=tol**2)[1]
        p_bar = self.f_p_bar(z, l_bar, eta, tau)

        if tag == POINT_PEAK_AVG:
            return p_bar, (z, l_bar, v_bar, eta, tau)

        ps_bar, p0_bar, pb_bar, _, _ = self.to_ps_p0_pb_vb_stag(l_bar * self.l_0, v_bar * self.v_j, p_bar, tau, eta)

        if tag == POINT_PEAK_BREECH:
            return ps_bar, (z, l_bar, v_bar, eta, tau)
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
        m_dot = self.c_a * self.v_j * self.s_j * p / (self.f * tau**0.5)
        # mass flow rate, rearward
        sb = self.s * self.chi_k
        vb = m_dot * (self.vol_0 + self.s * l) / (sb * (self.w - y))
        # velocity impinging upon the rear of the breech before nozzle constriction

        h = min(inf if v == 0 else (vb / v), 2 * self.phi_1 * self.m / (self.w - y) + 1)

        ps = p / (1 + (self.w - y) / (3 * self.phi_1 * self.m) * (1 - 0.5 * h))  # shot base pressure
        p0 = ps * (1 + (self.w - y) / (2 * self.phi_1 * self.m) * (1 + h) ** -1)  # stagnation point pressure
        pb = ps * (1 + (self.w - y) / (2 * self.phi_1 * self.m) * (1 - h))  # breech pressure
        stag = h / (1 + h)  # relative location of the stagnation point
        return ps, p0, pb, vb, stag

    def to_px(self, l: float, v: float, vb: float, ps: float, eta: float, x: float) -> float:
        """
        convert x, the physical displacement from breech bottom, to
        effective length in the equivalent gun.
        """
        l_k = self.l_0 / self.chi_k  # physical length of the chamber.
        s_k = self.s * self.chi_k
        if x < l_k:
            z = x * s_k / (l_k * s_k + l * self.s)
        else:
            z = (l_k * s_k + (x - l_k) * self.s) / (l_k * s_k + l * self.s)

        h = min(inf if v == 0 else (vb / v), 2 * self.phi_1 * self.m / (self.w * (1 - eta)) + 1)
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
            vr = (
                (gamma + 1) / (gamma - 1) * (1 - 2 / (gamma + 1) * (vr * sr) ** (1 - gamma))
            ) ** 0.5  # Eq. A.23 with neglect of covolume term

        cf = (2 / (gamma + 1)) ** (gamma / (gamma - 1)) * (
            gamma * vr + sr ** (1 - gamma) * vr ** (-gamma)
        )  # Eq. A.28 in source with neglect of covolume

        return cf

    def structure(
        self,
        recoilless_result: RecoillessResult,
        step: int,
        tol: float,
        structural_material: Material,
        structural_safety_factor: float = 1.1,
        autofrettage: bool = True,
        **_,
    ) -> None:
        step = max(step, 1)
        l_c = self.l_c
        l_g = self.l_g

        sigma = structural_material.yield_strength
        s = self.s
        gamma = self.theta + 1

        # a_bar = self.a_bar
        r_b = 0.5 * self.caliber  # radius of the shot.

        r_c = r_b * self.chi_k**0.5  # chamber/breech radius
        r_t = r_c * self.s_j_bar**0.5  # throat radius
        r_n = r_t * self.a_bar**0.5  # nozzle exit radius

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

        # strength requirement given structural safety factor.
        for i in range(len(p_probes)):
            p_probes[i] *= structural_safety_factor

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
                p = Recoilless.get_pr(gamma, (r / r_t) ** 2, tol)[1] * p_entry
            else:
                k = (x + l_b) / l_b
                r = k * r_c + (1 - k) * r_t
                p = Recoilless.get_pr(gamma, (r / r_t) ** 2, tol)[0] * p_entry

            r_probes.append(r)
            neg_ps.append(p)
            neg_ss.append(pi * r**2)

        if autofrettage:
            v_n, k_n, m_n = Gun.barrel_autofrettage(neg_x_probes, neg_ps, neg_ss, sigma)
        else:
            v_n, k_n, m_n = Gun.barrel_monoblock(neg_x_probes, neg_ps, neg_ss, sigma)

        hull = []
        for x, r, k, m in zip(neg_x_probes, r_probes, k_n, m_n):
            hull.append(OutlineEntry(x, r, r * k, r * m))

        nozzle_mass = v_n * structural_material.density

        i = step + 1
        x_c, p_c = x_probes[:i], p_probes[:i]  # c for chamber
        x_b, p_b = x_probes[i:], p_probes[i:]  # b for barrel

        if autofrettage:
            v_c, k_c, m_c = Gun.barrel_autofrettage(x_c, p_c, [s * self.chi_k for _ in x_c], sigma)
            v_b, k_b, m_b = Gun.barrel_autofrettage(x_b, p_b, [s for _ in x_b], sigma)
        else:
            v_c, k_c, m_c = Gun.barrel_monoblock(x_c, p_c, [s * self.chi_k for _ in x_c], sigma)
            v_b, k_b, m_b = Gun.barrel_monoblock(x_b, p_b, [s for _ in x_b], sigma)

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
        recoilless_result.tubeMass = tube_mass + nozzle_mass

        logger.info("conducted structural calculation.")

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
        # pressure ratio of nozzle throat over the upstream chamber pressure
        pr_c = (2 / (gamma + 1)) ** (gamma / (gamma - 1))

        if ar == 1:
            return pr_c, pr_c
        else:
            pr_sup, _ = dekker(lambda pr: Recoilless.get_ar(gamma, pr), 0, pr_c, y=ar, y_rel_tol=tol)
            pr_sub, _ = dekker(lambda pr: Recoilless.get_ar(gamma, pr), pr_c, 1, y=ar, y_rel_tol=tol)
            return pr_sub, pr_sup


if __name__ == "__main__":
    pass
