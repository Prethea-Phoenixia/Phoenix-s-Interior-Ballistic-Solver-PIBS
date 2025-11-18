from __future__ import annotations

import json
import logging
from math import log


from . import MAX_ITER, POINT_PEAK_AVG, POINT_PEAK_BREECH, POINT_PEAK_SHOT, SOL_LAGRANGE, SOL_MAMONTOV, SOL_PIDDUCK
from . import Points, Solutions
from .gun import pidduck
from .num import dekker, gss, rkf
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .prop import Propellant

from .constrained import Constrained


logger = logging.getLogger(__name__)


class ConstrainedGun(Constrained):
    def __init__(
        self,
        caliber: float,
        shot_mass: float,
        propellant: Propellant,
        start_pressure: float,
        drag_coefficient: float,
        design_pressure: float,
        design_velocity: float,
        chambrage: float,
        tol: float,
        min_web: float = 1e-6,
        max_length: float = 1e3,
        sol: Solutions = SOL_LAGRANGE,
        ambient_density: float = 1.204,
        ambient_pressure: float = 101.325e3,
        ambient_adb_index: float = 1.4,
        control: Points = POINT_PEAK_AVG,
        max_iteration: int = MAX_ITER,
        **_,
    ):
        super().__init__(
            caliber=caliber,
            shot_mass=shot_mass,
            propellant=propellant,
            start_pressure=start_pressure,
            drag_coefficient=drag_coefficient,
            design_pressure=design_pressure,
            design_velocity=design_velocity,
            tol=tol,
            chambrage=chambrage,
            min_web=min_web,
            max_length=max_length,
            ambient_density=ambient_density,
            ambient_pressure=ambient_pressure,
            ambient_adb_index=ambient_adb_index,
            control=control,
        )

        self.sol = sol
        self.max_iteration = max_iteration

    def to_json(self) -> str:
        return json.dumps(
            {**json.loads(super().to_json()), "sol": self.sol},
            ensure_ascii=False,
        )

    @Constrained.validate_solve_inputs
    def solve(
        self,
        *_,
        load_fraction: float,
        charge_mass_ratio: float,
        length_gun: float | None = None,
        known_bore: bool = False,
        labda_1: float | None = None,
        labda_2: float | None = None,
        cc: float | None = None,
        it: int = 0,
        **__,
    ) -> tuple[float, float]:
        w = self.m * charge_mass_ratio
        v_0 = w / (self.rho_p * load_fraction)
        delta = w / v_0

        l_0 = v_0 / self.s

        if length_gun is None:
            l_bar_g_0 = self.max_length / l_0
        else:
            l_bar_g_0 = length_gun / l_0

        if cc is None:
            cc = 1 - (1 - 1 / self.chi_k) * log(l_bar_g_0 + 1) / l_bar_g_0

        if any((labda_1 is None, labda_2 is None)):
            if self.sol == SOL_LAGRANGE:
                labda_1, labda_2 = 1 / 2, 1 / 3
            elif self.sol == SOL_PIDDUCK:
                labda_1, labda_2 = pidduck(w / (self.phi_1 * self.m), self.theta + 1, self.tol)
            elif self.sol == SOL_MAMONTOV:
                labda_1, labda_2 = pidduck(w / (self.phi_1 * self.m), 1, self.tol)
            else:
                raise ValueError("Unknown Solution")

        phi = self.phi_1 + labda_2 * w / self.m * cc
        v_j = (2 * self.f * w / (self.theta * phi * self.m)) ** 0.5
        v_bar_d = self.v_d / v_j

        if self.ambient_density:
            c_a_bar = (self.ambient_adb_index * self.ambient_pressure / self.ambient_density) ** 0.5 / v_j
            p_a_bar = self.ambient_pressure / (self.f * delta)
        else:
            c_a_bar, p_a_bar = 0, 0

        if v_j < self.v_d and not known_bore:
            raise ValueError(
                f"Propellant load too low to achieve design velocity. The 2nd ballistic limit for this loading \
conditions is {v_j:.4g} m/s."
            )

        psi_0 = (1 / delta - 1 / self.rho_p) / (self.f / self.p_0 + self.alpha - 1 / self.rho_p)
        z_0, _ = dekker(self.propellant.f_psi_z, 0, 1, y=psi_0, y_rel_tol=self.tol)

        def func_p_control_bar(z: float, l_bar: float, v_bar: float) -> float:
            psi = self.f_psi_z(z)
            l_psi_bar = 1 - delta / self.rho_p - delta * (self.alpha - 1 / self.rho_p) * psi
            p_bar = max((psi - v_bar**2) / (l_bar + l_psi_bar), p_a_bar)
            if self.control == POINT_PEAK_AVG:
                return p_bar
            else:
                cc_prime = (1 / self.chi_k + l_bar) / (1 + l_bar)
                labda_1_prime = labda_1 * cc_prime
                labda_2_prime = labda_2 * cc_prime

                factor_s = 1 + labda_2_prime * (w / (self.phi_1 * self.m))
                factor_b = (self.phi_1 * self.m + labda_2_prime * w) / (self.phi_1 * self.m + labda_1_prime * w)

                if self.control == POINT_PEAK_SHOT:
                    return p_bar / factor_s
                elif self.control == POINT_PEAK_BREECH:
                    return p_bar / factor_b
                else:
                    raise ValueError(f"Unknown control {self.control}")

        p_bar_s = func_p_control_bar(z_0, 0, 0)
        p_bar_d = self.p_d / (self.f * delta)

        if p_bar_d < p_bar_s:
            raise ValueError("Design pressure cannot be lower than shot start pressure.")
        l_bar_d = self.max_length / l_0

        def abort_z(
            x: float, ys: tuple[float, float, float], record: list[tuple[float, tuple[float, float, float]]]
        ) -> bool:
            z, (t_bar, l_bar, v_bar) = x, ys
            p_bar = func_p_control_bar(z, l_bar, v_bar)
            o_x, o_ys = record[-1]
            o_z = o_x
            ot_bar, ol_bar, ov_bar = o_ys
            op_bar = func_p_control_bar(o_z, ol_bar, ov_bar)

            return (p_bar < op_bar) or (p_bar > 2 * p_bar_d)

        def func_p_e_1(e_1: float) -> tuple[float, float, float, float, float]:
            b_e_1 = (
                self.s**2 * e_1**2 / (self.f * phi * w * self.m * self.u_1**2) * (self.f * delta) ** (2 * (1 - self.n))
            )

            def ode_z(z: float, t_l_v: tuple[float, float, float], __: float) -> tuple[float, float, float]:
                """burnup domain ode of internal ballistics"""
                t_bar, l_bar, v_bar = t_l_v
                psi = self.f_psi_z(z)
                l_psi_bar = 1 - delta / self.rho_p - delta * (self.alpha - 1 / self.rho_p) * psi

                p_bar = max((psi - v_bar**2) / (l_bar + l_psi_bar), p_a_bar)

                dt_bar = (2 * b_e_1 / self.theta) ** 0.5 * p_bar**-self.n  # dt_bar/dz
                dl_bar = v_bar * dt_bar
                dv_bar = 0.5 * self.theta * (p_bar - self.func_p_ad_bar(v_bar, c_a_bar, p_a_bar)) * dt_bar

                return dt_bar, dl_bar, dv_bar

            record = [(z_0, (0, 0, 0))]

            z_k, (t_bar_k, l_bar_k, v_bar_k), _ = rkf(
                d_func=ode_z,
                ini_val=(0, 0, 0),
                x_0=z_0,
                x_1=self.z_b,
                rel_tol=self.tol,
                abort_func=abort_z,
                record=record,
            )

            p_bar_j = func_p_control_bar(z_k, l_bar_k, v_bar_k)

            if p_bar_j >= 2 * p_bar_d:  # case for abort due to excessive pressure
                return p_bar_j - p_bar_d, z_k, t_bar_k, l_bar_k, v_bar_k

            if len(record) > 1:
                z_j = record[-2][0]
            else:
                z_j = z_0

            def func_p_z(z: float) -> tuple[float, float, float, float, float]:
                i = record.index([v for v in record if v[0] <= z][-1])
                x = record[i][0]
                ys = record[i][1]

                r = []
                t_bar, l_bar, v_bar = rkf(d_func=ode_z, ini_val=ys, x_0=x, x_1=z, rel_tol=self.tol, record=r)[1]
                xs = [v[0] for v in record]
                record.extend(v for v in r if v[0] not in xs)
                record.sort()
                return func_p_control_bar(z, l_bar, v_bar) - p_bar_d, z, t_bar, l_bar, v_bar

            z_p = 0.5 * sum(gss(lambda z: func_p_z(z)[0], z_j, z_k, x_tol=self.tol, find_min=False))
            return func_p_z(z_p)

        probe_web = 0.5 * self.min_web

        dp_bar_probe = next_dp_bar_probe = func_p_e_1(probe_web)[0]
        next_probe_web = probe_web

        while dp_bar_probe * next_dp_bar_probe > 0:  # is same sign
            probe_web = next_probe_web
            dp_bar_probe = next_dp_bar_probe
            next_probe_web = probe_web * 2 if dp_bar_probe > 0 else probe_web * 0.5
            next_dp_bar_probe = func_p_e_1(next_probe_web)[0]

        e_1_solved, _ = dekker(
            lambda _e_1: func_p_e_1(_e_1)[0], probe_web, next_probe_web, y_rel_tol=self.tol  # >0  # ?0
        )
        p_bar_dev, z_i, t_bar_i, l_bar_i, v_bar_i = func_p_e_1(e_1_solved)

        if known_bore:
            return e_1_solved, length_gun

        if v_j * v_bar_i > self.v_d:
            raise ValueError(f"Design velocity exceeded before peak pressure point (V = {v_bar_i * v_j:.4g} m/s).")

        b = (
            self.s**2
            * e_1_solved**2
            / (self.f * phi * w * self.m * self.u_1**2)
            * (self.f * delta) ** (2 * (1 - self.n))
        )

        def ode_v(v_bar: float, t_z_l: tuple[float, float, float], __: float) -> tuple[float, float, float]:
            t_bar, z, l_bar = t_z_l
            psi = self.f_psi_z(z)

            l_psi_bar = 1 - delta / self.rho_p - delta * (self.alpha - 1 / self.rho_p) * psi
            p_bar = max((psi - v_bar**2) / (l_bar + l_psi_bar), p_a_bar)
            dt_bar = 2 / (self.theta * (p_bar - self.func_p_ad_bar(v_bar, c_a_bar, p_a_bar)))
            dz = dt_bar * (0.5 * self.theta / b) ** 0.5 * p_bar**self.n
            dl_bar = v_bar * dt_bar

            return dt_bar, dz, dl_bar

        def abort_v(
            _: float, ys: tuple[float, float, float], record: list[tuple[float, tuple[float, float, float]]]
        ) -> bool:
            t_bar, _, l_bar = ys
            _, (ot_bar, _, _) = record[-1]
            return l_bar > l_bar_d or t_bar < ot_bar

        v_t_z_l_record = [(v_bar_i, (t_bar_i, z_i, l_bar_i))]
        try:
            v_bar_g, (t_bar_g, z_g, l_bar_g), _ = rkf(
                d_func=ode_v,
                ini_val=(t_bar_i, z_i, l_bar_i),
                x_0=v_bar_i,
                x_1=v_bar_d,
                rel_tol=self.tol,
                abort_func=abort_v,
                record=v_t_z_l_record,
            )
        except ValueError:
            v_bar_m, (t_bar_m, z_m, l_bar_m) = v_t_z_l_record[-1]
            p_max = func_p_control_bar(z_m, l_bar_m, v_bar_m) * self.f * delta
            v_max = v_bar_m * v_j
            lmax = l_bar_m * l_0
            raise ValueError(
                "Integration appears to be approaching asymptote, "
                + f"last calculated to v = {v_max:.4g} m/s, x = {lmax:.4g} m, p = {p_max * 1e-6:.4g} MPa. "
                + "This indicates an excessive velocity target relative to pressure developed."
            )

        p_bar_g = func_p_control_bar(z_g, l_bar_g, v_bar_g)

        v_g = v_bar_g * v_j
        l_g = l_bar_g * l_0
        p_g = p_bar_g * self.f * delta
        if l_bar_g > l_bar_d:
            raise ValueError(
                f"Solution requires excessive tube length, last calculated to v = {v_g:.4g} m/s, x = {l_g:.4g} m, p = {p_g * 1e-6:.4g} MPa."
            )
        cc_n = 1 - (1 - 1 / self.chi_k) * log(l_bar_g + 1) / l_bar_g
        l_g = l_bar_g * l_0
        if abs((l_bar_g - l_bar_g_0) / min(l_bar_g, l_bar_g_0)) > self.tol and it < self.max_iteration:
            return self.solve(
                load_fraction=load_fraction,
                charge_mass_ratio=charge_mass_ratio,
                labda_1=labda_1,
                labda_2=labda_2,
                cc=cc_n,
                it=it + 1,
                length_gun=l_g,
                known_bore=known_bore,
            )
        else:
            logger.info(
                f"ω/m = {charge_mass_ratio:.2f}, Δ/ρ = {load_fraction:.2f} -> e_1 = {e_1_solved * 1e3:.2f} mm, l_g = {l_g * 1e3:.0f} mm ({it+1} iterations)"
            )
            return e_1_solved, l_g


if __name__ == "__main__":
    pass
