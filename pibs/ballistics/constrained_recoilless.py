from __future__ import annotations

import json
import logging

from .constrained import Constrained
from . import POINT_PEAK_AVG, POINT_PEAK_BREECH, POINT_PEAK_SHOT, POINT_PEAK_STAG
from . import Points
from .num import dekker, gss, rkf
from .prop import Propellant
from .recoilless import Recoilless

logger = logging.getLogger(__name__)


class ConstrainedRecoilless(Constrained):
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
        super().__init__(
            caliber=caliber,
            shot_mass=shot_mass,
            propellant=propellant,
            start_pressure=start_pressure,
            drag_coefficient=drag_coefficient,
            chambrage=chambrage,
            tol=tol,
            design_pressure=design_pressure,
            design_velocity=design_velocity,
            min_web=min_web,
            max_length=max_length,
            ambient_density=ambient_density,
            ambient_pressure=ambient_pressure,
            ambient_adb_index=ambient_adb_index,
            control=control,
        )

        if any((nozzle_expansion < 1, nozzle_efficiency > 1, nozzle_efficiency <= 0)):
            raise ValueError("Invalid parameters for constrained design")

        self.chi_0 = nozzle_efficiency
        self.a_bar = nozzle_expansion

    def to_json(self) -> str:
        return json.dumps(
            {
                **super().to_json(),
                "nozzle_expansion": self.a_bar,
                "nozzle_efficiency": self.chi_0,
            },
            ensure_ascii=False,
        )

    @Constrained.validate_solve_inputs
    def solve(
        self,
        load_fraction: float,
        charge_mass_ratio: float,
        length_gun: float | None = None,
        known_bore: bool = False,
        **_,
    ) -> tuple[float, float]:

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
        v_j = (2 * self.f * w / (self.theta * phi * self.m)) ** 0.5

        t_scale = l_0 / v_j

        if self.ambient_density != 0:
            c_a_bar = (self.ambient_adb_index * self.ambient_pressure / self.ambient_density) ** 0.5 / v_j
            p_a_bar = self.ambient_pressure / (self.f * delta)
        else:
            c_a_bar = 0
            p_a_bar = 0

        if v_j < self.v_d and not known_bore:
            raise ValueError(
                "propellant load too low to achieve design velocity, "
                + f"the 2nd ballistic limit for this loading conditions is {v_j:.4g} m/s,"
                + " and recoilless guns only achieve a part of that as well.",
            )

        psi_0 = (1 / delta - 1 / self.rho_p) / (self.f / self.p_0 + self.alpha - 1 / self.rho_p)
        z_0, _ = dekker(lambda _z: self.propellant.f_psi_z(_z), 0, 1, y=psi_0, y_rel_tol=self.tol)

        def func_p_avg_bar(psi: float, l_bar: float, eta: float, tau: float) -> float:
            l_psi_bar = 1 - delta * ((1 - psi) / self.rho_p + self.alpha * (psi - eta))
            return max(tau / (l_bar + l_psi_bar) * (psi - eta), p_a_bar)

        def f_p_bar_control(z: float, l_bar: float, v_bar: float, eta: float, tau: float) -> float:
            psi = self.f_psi_z(z)
            p_bar = func_p_avg_bar(psi=psi, l_bar=l_bar, eta=eta, tau=tau)

            m_dot = c_a * v_j * s_j * p_bar * delta * tau**-0.5
            vb = m_dot * (vol_0 + self.s * l_bar * l_0) / (self.s * self.chi_k * w * (1 - eta))

            h_lim = 2 * self.phi_1 * self.m / (w * (1 - eta)) + 1
            h = h_lim if v_bar == 0 else min(h_lim, vb / (v_j * v_bar))

            if self.control == POINT_PEAK_AVG:
                return p_bar
            else:
                p_s_bar = p_bar / (1 + w * (1 - eta) / (3 * self.phi_1 * self.m) * (1 - 0.5 * h))
                if self.control == POINT_PEAK_SHOT:
                    return p_s_bar
                elif self.control == POINT_PEAK_STAG:
                    return p_s_bar * (1 + w * (1 - eta) / (2 * self.phi_1 * self.m) * (1 + h) ** -1)
                elif self.control == POINT_PEAK_BREECH:
                    if h == h_lim:
                        return 0.0
                    else:
                        return p_s_bar * (1 + w * (1 - eta) / (2 * self.phi_1 * self.m) * (1 - h))
            raise ValueError(f"unknown control {self.control}")

        p_bar_d = self.p_d / (self.f * delta)  # convert to unitless
        l_bar_d = self.max_length / l_0

        def abort_t(
            x: float,
            ys: tuple[float, float, float, float, float],
            record: list[tuple[float, tuple[float, float, float, float, float]]],
        ) -> bool:
            t_bar, (z, l_bar, v_bar, eta, tau) = x, ys
            p_bar_control = f_p_bar_control(z=z, l_bar=l_bar, v_bar=v_bar, eta=eta, tau=tau)
            p_avg_bar = func_p_avg_bar(psi=self.f_psi_z(z), l_bar=l_bar, eta=eta, tau=tau)

            o_t_bar, (o_z, o_l_bar, o_v_bar, o_eta, o_tau) = record[-1]
            op_bar_control = f_p_bar_control(z=o_z, l_bar=o_l_bar, v_bar=o_v_bar, eta=o_eta, tau=o_tau)
            op_avg_bar = func_p_avg_bar(psi=self.f_psi_z(o_z), l_bar=o_l_bar, eta=o_eta, tau=o_tau)

            return (p_bar_control > p_bar_d * (1 + 2 * self.tol)) or (
                (p_bar_control < op_bar_control) and (p_avg_bar < op_avg_bar)
            )

        def func_p_e1(e_1: float) -> tuple[float, float, tuple[float, float, float, float, float]]:
            b_e_1 = (
                (self.s**2 * e_1**2)
                / (self.f * phi * w * self.m * self.u_1**2)
                * (self.f * delta) ** (2 * (1 - self.n))
            )

            def ode_t(
                t: float, z_l_v_eta_tau: tuple[float, float, float, float, float], dt: float
            ) -> tuple[float, float, float, float, float]:
                z, l_bar, v_bar, eta, tau = z_l_v_eta_tau
                psi = self.f_psi_z(z)
                d_psi__d_z = self.f_sigma_z(z)
                p_bar = func_p_avg_bar(psi=psi, l_bar=l_bar, eta=eta, tau=tau)
                dz = (0.5 * self.theta / b_e_1) ** 0.5 * p_bar**self.n

                dl_bar = v_bar
                dv_bar = self.theta * 0.5 * (p_bar - self.func_p_ad_bar(v_bar, c_a_bar, p_a_bar))

                d_eta = c_a * s_j_bar * p_bar * tau**-0.5
                d_tau = ((1 - tau) * (d_psi__d_z * dz) - 2 * v_bar * dv_bar - self.theta * tau * d_eta) / (psi - eta)

                return dz, dl_bar, dv_bar, d_eta, d_tau

            record = [(0.0, (z_0, 0.0, 0.0, 0.0, 1.0))]

            aborted = False
            while not aborted:
                t_bar_next = record[-1][0] + 10 * t_scale
                t_bar_k, (z_k, l_bar_k, v_bar_k, eta_k, tau_k), aborted = rkf(
                    d_func=ode_t,
                    ini_val=record[-1][1],
                    x_0=record[-1][0],
                    x_1=t_bar_next,
                    rel_tol=self.tol,
                    abort_func=abort_t,
                    record=record,
                )
                record.append((t_bar_k, (z_k, l_bar_k, v_bar_k, eta_k, tau_k)))

            t_bar_k, (z_k, l_bar_k, v_bar_k, eta_k, tau_k) = record[-1]
            p_bar_k = f_p_bar_control(z_k, l_bar_k, v_bar_k, eta_k, tau_k)

            if p_bar_k > (1 + 2 * self.tol) * p_bar_d:  # case for abort due to excessive pressure
                return p_bar_k - p_bar_d, t_bar_k, (z_k, l_bar_k, v_bar_k, eta_k, tau_k)

            t_bar_j = 0.0

            def func_p_t_bar(t_bar: float) -> tuple[float, float, tuple[float, float, float, float, float]]:
                i = record.index([v for v in record if v[0] <= t_bar][-1])
                x = record[i][0]
                ys = record[i][1]

                z, l_bar, v_bar, eta, tau = rkf(
                    d_func=ode_t, ini_val=ys, x_0=x, x_1=t_bar, rel_tol=self.tol, record=record
                )[1]

                record.sort()
                return f_p_bar_control(z, l_bar, v_bar, eta, tau) - p_bar_d, t_bar, (z, l_bar, v_bar, eta, tau)

            t_p = 0.5 * sum(gss(lambda t_bar: func_p_t_bar(t_bar)[0], t_bar_j, t_bar_k, x_tol=self.tol, find_min=False))
            return func_p_t_bar(t_p)

        probe_web = next_probe_web = 0.5 * self.min_web
        dp_bar_probe = next_dp_bar_probe = func_p_e1(probe_web)[0]

        while dp_bar_probe * next_dp_bar_probe > 0:  # is same sign
            probe_web = next_probe_web
            dp_bar_probe = next_dp_bar_probe
            next_probe_web = (probe_web * 2) if dp_bar_probe > 0 else (probe_web * 0.5)
            next_dp_bar_probe = func_p_e1(next_probe_web)[0]

        e_1_solved, _ = dekker(
            lambda _e_1: func_p_e1(_e_1)[0],
            probe_web,
            next_probe_web,
            y_abs_tol=self.tol * p_bar_d,  # debug=True
        )

        if known_bore:
            return e_1_solved, length_gun

        v_bar_d = self.v_d / v_j
        dp_bar_i, *vals_1 = func_p_e1(e_1_solved)
        t_bar_i, (z_i, l_bar_i, v_bar_i, eta_i, tau_i) = vals_1

        if v_bar_i > v_bar_d:
            raise ValueError(f"Design velocity exceeded before peak pressure point (V = {v_bar_i * v_j:.4g} m/s).")

        b = (
            (self.s**2 * e_1_solved**2)
            / (self.f * phi * w * self.m * self.u_1**2)
            * (self.f * delta) ** (2 * (1 - self.n))
        )

        def ode_v(
            v_bar: float, t_z_l_eta_tau: tuple[float, float, float, float, float], d_v_bar: float
        ) -> tuple[float, float, float, float, float]:
            t_bar, z, l_bar, eta, tau = t_z_l_eta_tau
            psi = self.f_psi_z(z)
            d_psi = self.f_sigma_z(z)

            p_bar = func_p_avg_bar(psi=psi, l_bar=l_bar, eta=eta, tau=tau)
            dt_bar = 2 / (self.theta * (p_bar - self.func_p_ad_bar(v_bar, c_a_bar, p_a_bar)))
            dz = dt_bar * (0.5 * self.theta / b) ** 0.5 * p_bar**self.n
            dl_bar = v_bar * dt_bar
            d_eta = c_a * s_j_bar * p_bar / tau**0.5 * dt_bar
            d_tau = ((1 - tau) * (d_psi * dz) - 2 * v_bar - self.theta * tau * d_eta) / (psi - eta)
            return dt_bar, dz, dl_bar, d_eta, d_tau

        def abort_v(
            x: float,
            ys: tuple[float, float, float, float, float],
            record: list[tuple[float, tuple[float, float, float, float, float]]],
        ) -> bool:
            t_bar, z, l_bar, eta, tau = ys
            ot_bar, *_ = record[-1][-1]
            return l_bar > l_bar_d or t_bar < ot_bar

        v_t_z_l_eta_tau_record = [(v_bar_i, (t_bar_i, z_i, l_bar_i, eta_i, tau_i))]
        try:
            v_bar_g, (t_bar_g, z_g, l_bar_g, eta_g, tau_g), _ = rkf(
                d_func=ode_v,
                ini_val=(t_bar_i, z_i, l_bar_i, eta_i, tau_i),
                x_0=v_bar_i,
                x_1=v_bar_d,
                rel_tol=self.tol,
                abort_func=abort_v,
                record=v_t_z_l_eta_tau_record,
            )

        except ValueError:
            v_bar_m, (t_bar_m, z_m, l_bar_m, eta_m, tau_m) = v_t_z_l_eta_tau_record[-1]
            p_max = f_p_bar_control(z_m, l_bar_m, v_bar_m, eta_m, tau_m) * self.f * delta
            v_max = v_bar_m * v_j
            lmax = l_bar_m * l_0
            raise ValueError(
                "Integration appears to be approaching asymptote, "
                + "last calculated to v = {:.4g} m/s, ".format(v_max)
                + "x = {:.4g} m, p = {:.4g} MPa. ".format(lmax, p_max * 1e-6)
                + "This indicates an excessive velocity target relative to pressure developed."
            )

        p_bar_g = f_p_bar_control(z_g, l_bar_g, v_bar_g, eta_g, tau_g)

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


if __name__ == "__main__":
    pass
