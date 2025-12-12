from __future__ import annotations

import json
import logging
from math import inf, pi
from random import uniform
from .optimize_gun import probe_func
from . import JSONable
from . import (
    MAX_GUESSES,
    POINT_PEAK_AVG,
    POINT_PEAK_BREECH,
    POINT_PEAK_SHOT,
    POINT_PEAK_STAG,
    MIN_BORE_VOLUME,
    MIN_PROJ_TRAVEL,
)
from . import Points, Optimization_Targets
from .generics import DelegatesPropellant
from .num import dekker, gss, rkf
from .prop import Propellant
from .recoilless import Recoilless

logger = logging.getLogger(__name__)


class ConstrainedRecoilless(DelegatesPropellant, JSONable):
    def __init__(
        self,
        caliber: float,
        shot_mass: float,
        propellant: Propellant,
        start_pressure: float,
        drag_coefficient: float,
        nozzle_expansion: float,
        nozzle_efficiency: float,
        chambrage: float,
        tol: float,
        design_pressure: float,
        design_velocity: float,
        min_web: float = 1e-6,
        max_length: float = 1e3,
        ambient_density: float = 1.204,
        ambient_pressure: float = 101.325e3,
        ambient_adb_index: float = 1.4,
        control: Points = POINT_PEAK_AVG,
        **_,
    ):
        super().__init__(propellant=propellant)

        if any(
            (
                caliber <= 0,
                shot_mass <= 0,
                start_pressure <= 0,
                drag_coefficient < 0,
                drag_coefficient >= 1,
                nozzle_expansion < 1,
                nozzle_efficiency > 1,
                nozzle_efficiency <= 0,
                chambrage < 1,
            )
        ):
            raise ValueError("Invalid parameters for constrained design")

        if any((design_pressure <= 0, design_velocity <= 0)):
            raise ValueError("Invalid design constraint")

        ambient_pressure = max(ambient_pressure, 1)
        ambient_adb_index = max(ambient_adb_index, 1)

        self.caliber = caliber
        self.s = (0.5 * caliber) ** 2 * pi
        self.m = shot_mass
        self.propellant = propellant
        self.p_0 = start_pressure
        self.phi_1 = 1 / (1 - drag_coefficient)
        self.p_d = design_pressure
        self.v_d = design_velocity

        self.chi_0 = nozzle_efficiency
        self.a_bar = nozzle_expansion

        self.chi_k = chambrage

        self.tol = tol
        self.min_web = min_web
        self.max_length = max_length
        self.ambient_density = ambient_density
        self.ambient_pressure = ambient_pressure
        self.ambient_adb_index = ambient_adb_index
        self.control = control

    def to_json(self) -> str:
        return json.dumps(
            {
                "caliber": self.caliber,
                "shot_mass": self.m,
                "propellant": json.loads(self.propellant.to_json()),
                "start_pressure": self.p_0,
                "drag_coefficient": 1 - 1 / self.phi_1,
                "nozzle_expansion": self.a_bar,
                "nozzle_efficiency": self.chi_0,
                "chambrage": self.chi_k,
                "tol": self.tol,
                "design_pressure": self.p_d,
                "design_velocity": self.v_d,
                "min_web": self.min_web,
                "max_length": self.max_length,
                "ambient_density": self.ambient_density,
                "ambient_pressure": self.ambient_pressure,
                "ambient_adb_index": self.ambient_adb_index,
                "control": self.control,
            },
            ensure_ascii=False,
        )

    @classmethod
    def from_json(cls, json_dict: dict) -> ConstrainedRecoilless:
        deserialized_dict = {}
        for key, value in json_dict.items():
            if key == "propellant":
                value = Propellant.from_json(value)

            deserialized_dict[key] = value

        return cls(**deserialized_dict)

    def _f_p_bar_ad(self, v_bar: float, c_a_bar: float, p_a_bar: float) -> float:
        if c_a_bar and v_bar > 0.0:
            v_r = v_bar / c_a_bar
            return (
                +0.25 * self.ambient_adb_index * (self.ambient_adb_index + 1) * v_r**2
                + self.ambient_adb_index * v_r * (1 + (0.25 * (self.ambient_adb_index + 1)) ** 2 * v_r**2) ** 0.5
            ) * p_a_bar

        else:
            return 0

    def solve(
        self,
        load_fraction: float,
        charge_mass_ratio: float,
        length_gun: float | None = None,
        known_bore: bool = False,
        suppress: bool = False,
        **_,
    ) -> tuple[float, float]:
        if any((charge_mass_ratio <= 0, load_fraction <= 0, load_fraction > 1)):
            raise ValueError("Invalid parameters to solve constrained design problem")
        """
        minWeb  : represents minimum possible grain size
        """
        w = self.m * charge_mass_ratio
        vol_0 = w / (self.rho_p * load_fraction)
        delta = w / vol_0
        l_0 = vol_0 / self.s
        gamma = self.theta + 1
        phi = self.phi_1 + w / (3 * self.m)

        s_j_bar = 1 / (Recoilless.get_cf(gamma, self.a_bar, self.tol) * self.chi_0)
        if s_j_bar > self.chi_k:
            raise ValueError(
                "Achieving recoilless condition necessitates a larger throat area than could be fit into breech face."
            )
        s_j = s_j_bar * self.s

        k_0 = (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1))) * gamma**0.5

        phi_2 = 1
        c_a = (0.5 * self.theta * phi * self.m / w) ** 0.5 * k_0 * phi_2  # flow rate value

        """
        it is impossible to account for the chambrage effect given unspecified
        barrel length, in our formulation
        """
        v_j = (2 * self.f * w / (self.theta * phi * self.m)) ** 0.5

        if self.ambient_density != 0:
            c_a_bar = (self.ambient_adb_index * self.ambient_pressure / self.ambient_density) ** 0.5 / v_j
            p_a_bar = self.ambient_pressure / (self.f * delta)
        else:
            c_a_bar = 0
            p_a_bar = 0

        if v_j < self.v_d and not known_bore:
            raise ValueError(
                "Propellant load too low to achieve design velocity. "
                + " The 2nd ballistic limit for this loading conditions is"
                + " {:.4g} m/s,".format(v_j)
                + " and recoilless guns only achieve a part of that as well."
            )

        psi_0 = (1 / delta - 1 / self.rho_p) / (self.f / self.p_0 + self.alpha - 1 / self.rho_p)
        z_0, _ = dekker(
            lambda _z: self.propellant.f_psi_z(_z), 0, 1, y=psi_0, y_rel_tol=self.tol, y_abs_tol=self.tol**2
        )

        def _f_p_bar(z: float, l_bar: float, v_bar: float, eta: float, tau: float) -> float:
            psi = self.f_psi_z(z)
            l_psi_bar = 1 - delta * ((1 - psi) / self.rho_p + self.alpha * (psi - eta))
            p_bar = max(tau / (l_bar + l_psi_bar) * (psi - eta), p_a_bar)

            if self.control == POINT_PEAK_AVG:
                return p_bar
            else:
                m_dot = c_a * v_j * s_j * p_bar * delta / (tau**0.5)
                vb = m_dot * (vol_0 + self.s * l_bar * l_0) / (self.s * self.chi_k * w * (1 - eta))

                h = min(inf if v_bar == 0 else (vb / (v_j * v_bar)), 2 * self.phi_1 * self.m / (w * (1 - eta)) + 1)
                p_s_bar = p_bar / (1 + w * (1 - eta) / (3 * self.phi_1 * self.m) * (1 - 0.5 * h))
                if self.control == POINT_PEAK_SHOT:
                    return p_s_bar
                elif self.control == POINT_PEAK_STAG:
                    return p_s_bar * (1 + w * (1 - eta) / (2 * self.phi_1 * self.m) * (1 + h) ** -1)
                elif self.control == POINT_PEAK_BREECH:
                    return p_s_bar * (1 + w * (1 - eta) / (2 * self.phi_1 * self.m) * (1 - h))
                else:
                    raise ValueError(f"unknown control {self.control}")

        p_bar_d = self.p_d / (self.f * delta)  # convert to unitless
        l_bar_d = self.max_length / l_0

        """
        step 1, find grain size that satisfies design pressure
        """

        def _abort_z(
            x: float,
            ys: tuple[float, float, float, float, float],
            _: list[tuple[float, tuple[float, float, float, float, float]]],
        ) -> bool:
            z, (_, l_bar, v_bar, eta, tau) = x, ys
            p_bar = _f_p_bar(z, l_bar, v_bar, eta, tau)
            return (p_bar > 2 * p_bar_d) or l_bar > l_bar_d

        def _f_p_e_1(e_1: float) -> tuple[float, float, float, float, float, float, float]:
            b_e_1 = (
                (self.s**2 * e_1**2)
                / (self.f * phi * w * self.m * self.u_1**2)
                * (self.f * delta) ** (2 * (1 - self.n))
            )

            def _ode_z(
                z: float, tlvetatau: tuple[float, float, float, float, float], _: float
            ) -> tuple[float, float, float, float, float]:
                t_bar, l_bar, v_bar, eta, tau = tlvetatau
                psi = self.f_psi_z(z)
                dpsi = self.f_sigma_z(z)

                l_psi_bar = 1 - delta * ((1 - psi) / self.rho_p + self.alpha * (psi - eta))
                p_bar = max(tau / (l_bar + l_psi_bar) * (psi - eta), p_a_bar)

                dt_bar = (2 * b_e_1 / self.theta) ** 0.5 * p_bar**-self.n
                dl_bar = v_bar * dt_bar
                dv_bar = 0.5 * self.theta * (p_bar - self._f_p_bar_ad(v_bar, c_a_bar, p_a_bar)) * dt_bar

                deta = c_a * s_j_bar * p_bar / tau**0.5 * dt_bar
                dtau = ((1 - tau) * dpsi - 2 * v_bar * dv_bar - self.theta * tau * deta) / (psi - eta)

                return dt_bar, dl_bar, dv_bar, deta, dtau

            record = [(z_0, (0.0, 0.0, 0.0, 0.0, 1.0))]
            try:
                z_k, (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j) = rkf(
                    d_func=_ode_z,
                    ini_val=(0.0, 0.0, 0.0, 0.0, 1.0),
                    x_0=z_0,
                    x_1=self.z_b,
                    rel_tol=self.tol,
                    abs_tol=self.tol**2,
                    abort_func=_abort_z,
                    record=record,
                )

                if z_k not in [line[0] for line in record]:
                    record.append((z_k, (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j)))
            except ValueError:
                z_k, (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j) = record[-1]

            p_bar_j = _f_p_bar(z_k, l_bar_j, v_bar_j, eta_j, tau_j)

            peak = None

            # find the last peak
            if len(record) > 2:
                p_bars = [
                    _f_p_bar(_z, _l_bar, _v_bar, _eta, _tau) for (_z, (_t_bar, _l_bar, _v_bar, _eta, _tau)) in record
                ]

                for _i, (_l, _c, _r) in enumerate(zip(p_bars[:-2], p_bars[1:-1], p_bars[2:])):
                    if _l < _c and _c > _r:
                        peak = _i + 1

            if peak is None:
                # no peak, so it suffice to compare the end points.
                if p_bar_j == 0:
                    p_bar_j = inf
                z_j, (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j) = record[-1]
                return p_bar_j - p_bar_d, z_j, t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j
            else:  # peak exist, must compare the peak and the two end points.
                z_j = record[peak - 1][0]
                z_k = record[peak + 1][0]

                def _f_p_z(z: float) -> tuple[float, float, float, float, float, float, float]:
                    i = record.index([v for v in record if v[0] <= z][-1])
                    x = record[i][0]
                    ys = record[i][1]

                    r = []
                    t_bar, l_bar, v_bar, eta, tau = rkf(
                        d_func=_ode_z, ini_val=ys, x_0=x, x_1=z, rel_tol=self.tol, abs_tol=self.tol**2, record=r
                    )[1]

                    xs = [v[0] for v in record]
                    record.extend(v for v in r if v[0] not in xs)
                    record.sort()
                    return _f_p_bar(z, l_bar, v_bar, eta, tau) - p_bar_d, z, t_bar, l_bar, v_bar, eta, tau

                z_p = 0.5 * sum(gss(lambda _z: _f_p_z(_z)[0], z_j, z_k, y_rel_tol=self.tol, find_min=False))

                return _f_p_z(z_p)

        probe_web = 0.5 * self.min_web
        dp_bar_probe, z_probe, *_ = _f_p_e_1(probe_web)

        if dp_bar_probe < 0:
            raise ValueError(
                "Design pressure cannot be achieved by varying web down to minimum. "
                + "Peak pressure found at phi = {:.4g} at {:.4g} MPa".format(
                    self.f_psi_z(z_probe), (dp_bar_probe + p_bar_d) * 1e-6 * self.f * delta
                )
            )

        while dp_bar_probe > 0:
            probe_web *= 2
            dp_bar_probe = _f_p_e_1(probe_web)[0]

        e_1_solved, _ = dekker(lambda _e_1: _f_p_e_1(_e_1)[0], probe_web, 0.5 * probe_web, y_rel_tol=self.tol)

        dp_bar_i, *vals_1 = _f_p_e_1(e_1_solved)
        z_i, t_bar_i, l_bar_i, v_bar_i, eta_i, tau_i = vals_1

        if known_bore:
            return e_1_solved, length_gun

        """
        step 2, find the requisite muzzle length to achieve design velocity
        """
        v_bar_d = self.v_d / v_j

        if v_bar_i > v_bar_d:
            if suppress:
                logger.warning("velocity target point occurred before peak pressure point.")
                logger.warning("this is currently being suppressed due to program control.")
            else:
                raise ValueError(f"Design velocity exceeded before peak pressure point (V = {v_bar_i * v_j:.4g} m/s).")

        b = (
            self.s**2
            * e_1_solved**2
            / (self.f * phi * w * self.m * self.u_1**2)
            * (self.f * delta) ** (2 * (1 - self.n))
        )

        def _ode_v(
            v_bar: float, tzletatau: tuple[float, float, float, float, float], __: float
        ) -> tuple[float, float, float, float, float]:
            _, z, l_bar, eta, tau = tzletatau
            psi = self.f_psi_z(z)
            dpsi = self.f_sigma_z(z)
            l_psi_bar = 1 - delta * ((1 - psi) / self.rho_p + self.alpha * (psi - eta))
            p_bar = max(tau / (l_bar + l_psi_bar) * (psi - eta), p_a_bar)

            dt_bar = 2 / (self.theta * (p_bar - self._f_p_bar_ad(v_bar, c_a_bar, p_a_bar)))

            dz = dt_bar * (0.5 * self.theta / b) ** 0.5 * p_bar**self.n

            dl_bar = v_bar * dt_bar

            deta = c_a * s_j_bar * p_bar / tau**0.5 * dt_bar  # deta / dv_bar
            dtau = (
                (1 - tau) * (dpsi * dz)
                - 2 * v_bar
                - self.theta * tau * deta  # dZ/dt_bar  # dv_bar/dt_bar  # deta/dt_bar
            ) / (psi - eta)

            return dt_bar, dz, dl_bar, deta, dtau

        def _abort_v(
            _: float,
            ys: tuple[float, float, float, float, float],
            record: list[tuple[float, tuple[float, float, float, float, float]]],
        ) -> bool:
            t_bar, _, l_bar, _, _ = ys
            ot_bar, *_ = record[-1][-1]
            return l_bar > l_bar_d or t_bar < ot_bar

        vtzlet_record = [(v_bar_i, (t_bar_i, z_i, l_bar_i, eta_i, tau_i))]
        try:
            v_bar_g, (t_bar_g, z_g, l_bar_g, eta_g, tau_g) = rkf(
                d_func=_ode_v,
                ini_val=(t_bar_i, z_i, l_bar_i, eta_i, tau_i),
                x_0=v_bar_i,
                x_1=v_bar_d,
                rel_tol=self.tol,
                abs_tol=self.tol**2,
                abort_func=_abort_v,
                record=vtzlet_record,
            )

        except ValueError:
            v_bar_m, (t_bar_m, z_m, l_bar_m, eta_m, tau_m) = vtzlet_record[-1]
            pmax = _f_p_bar(z_m, l_bar_m, v_bar_m, eta_m, tau_m) * self.f * delta
            vmax = v_bar_m * v_j
            lmax = l_bar_m * l_0
            raise ValueError(
                "Integration appears to be approaching asymptote, "
                + "last calculated to v = {:.4g} m/s, ".format(vmax)
                + "x = {:.4g} m, p = {:.4g} MPa. ".format(lmax, pmax * 1e-6)
                + "This indicates an excessive velocity target relative to pressure developed."
            )

        p_bar_g = _f_p_bar(z_g, l_bar_g, v_bar_g, eta_g, tau_g)

        v_g = v_bar_g * v_j
        l_g = l_bar_g * l_0
        p_g = p_bar_g * self.f * delta

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
            f"ω/m = {charge_mass_ratio:.2f}, Δ/ρ = {load_fraction:.2f} -> e_1 = {e_1_solved * 1e3:.2f} mm, l_g = {l_g * 1e3:.0f} mm"
        )
        return e_1_solved, l_bar_g * l_0

    def find_min_v(
        self,
        charge_mass_ratio: float,
        max_guess: int = MAX_GUESSES,
        target: Optimization_Targets = MIN_BORE_VOLUME,
        **_,
    ) -> tuple[float, float, float]:
        """
        find the minimum volume solution.
        """

        def _f(load_fraction: float) -> tuple[float, float, float]:
            e_1_delta, l_g_delta = self.solve(
                load_fraction=load_fraction, charge_mass_ratio=charge_mass_ratio, known_bore=False, suppress=True
            )
            return (
                e_1_delta,
                l_g_delta,
                l_g_delta + (self.m * charge_mass_ratio / (self.rho_p * load_fraction) / self.s),
            )

        records = []
        for i in range(max_guess):
            start_probe = uniform(self.tol, 1 - self.tol)
            try:
                records.append((start_probe, _f(start_probe)))
                break
            except ValueError:
                pass
        else:
            raise ValueError(f"Unable to find any valid load fraction with {max_guess} random samples.")

        logger.info(f"valid Δ/ρ = {start_probe:.3%}.")

        low, low_record = probe_func(_f, start_probe, self.tol, self.tol)
        records.extend(low_record)
        logger.info(f"min Δ/ρ = {low:.3%}.")

        high, high_record = probe_func(_f, start_probe, 1 - self.tol, self.tol)
        records.extend(high_record)
        logger.info(f"max Δ/ρ = {high:.3%}.")

        if target == MIN_PROJ_TRAVEL:
            _f_index = 1
        elif target == MIN_BORE_VOLUME:
            _f_index = 2
        else:
            raise ValueError(f"unknown target {target}")

        if len(records) > 2:
            records.sort(key=lambda x: x[0])  # sort by e_1
            for l, m, h in zip(records[:-2], records[1:-1], records[2:]):
                if l[1][_f_index] > m[1][_f_index] and h[1][_f_index] > m[1][_f_index]:
                    low, high = l[0], h[0]

        logger.info(f"solution constrained Δ/ρ : {low:.3%} - {high:.3%}")
        lf_low, lf_high = gss(lambda _lf: _f(_lf)[1], low, high, x_tol=self.tol, find_min=True)
        lf = 0.5 * (lf_high + lf_low)
        e_1, l_g, _ = _f(lf)
        logger.info(f"Optimal Δ/ρ = {lf:.2f}")
        return lf, e_1, l_g


if __name__ == "__main__":
    pass
