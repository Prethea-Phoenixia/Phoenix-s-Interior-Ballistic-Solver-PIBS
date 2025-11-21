from __future__ import annotations

import logging
from math import floor, log, pi
from random import uniform
from typing import Callable, TypeVar

from . import (
    MAX_GUESSES,
    MAX_ITER,
    POINT_PEAK_AVG,
    POINT_PEAK_BREECH,
    POINT_PEAK_SHOT,
    SOL_LAGRANGE,
    SOL_MAMONTOV,
    SOL_PIDDUCK,
    MIN_BORE_VOLUME,
    MIN_PROJ_TRAVEL,
)
from . import Points, Solutions, Optimization_Targets
from .generics import DelegatesPropellant
from .gun import pidduck
from .num import dekker, gss, rkf
from .prop import Propellant

"""
Machine-accuracy factor, determines that, if a numerical method
is used within another, then how much more accurate should the
inner method be called as compared to the outer. A small value
is necessary to prevent numerical instability form causing sporadic
appearance of outlier, at the cost of increased computation times.
"""


logger = logging.getLogger(__name__)


T = TypeVar("T")


def probe_func(
    func: Callable[[float], T], start: float, stop: float, tol: float, exceptions: tuple[Exception] = (ValueError,)
) -> tuple[float, list[tuple[float, T]]]:

    delta = stop - start
    k, n = 0, floor(log(abs(delta) / tol, 2)) + 1
    records = []

    probe = start
    new_probe = probe + delta
    while abs(2 * delta) > tol:
        try:
            records.append((new_probe, func(new_probe)))
            probe = new_probe
        except exceptions:
            delta *= 0.5
            k += 1
        finally:
            new_probe = probe + delta

    return probe, records


class Constrained(DelegatesPropellant):
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
        ambient_rho: float = 1.204,
        ambient_p: float = 101.325e3,
        ambient_gamma: float = 1.4,
        control: Points = POINT_PEAK_AVG,
        **_,
    ):
        super().__init__(propellant=propellant)
        if any((caliber <= 0, shot_mass <= 0, start_pressure <= 0, drag_coefficient < 0, drag_coefficient >= 1)):
            raise ValueError("Invalid parameters for constrained design")

        if any((design_pressure <= 0, design_velocity <= 0)):
            raise ValueError("Invalid design constraint")

        ambient_p = max(ambient_p, 1)
        ambient_gamma = max(ambient_gamma, 1)

        self.s = (caliber / 2) ** 2 * pi
        self.m = shot_mass
        self.p_0 = start_pressure
        self.phi_1 = 1 / (1 - drag_coefficient)

        # design limits
        self.p_d = design_pressure
        self.v_d = design_velocity
        self.chi_k = chambrage

        self.sol = sol
        self.ambient_rho = ambient_rho
        self.ambient_p = ambient_p
        self.ambient_gamma = ambient_gamma
        self.control = control

        self.min_web = min_web
        self.max_length = max_length

        self.tol = tol

    def _f_p_bar_ad(self, v_bar: float, c_a_bar: float, p_a_bar: float) -> float:
        if c_a_bar and v_bar > 0.0:
            v_r = v_bar / c_a_bar
            return (
                +0.25 * self.ambient_gamma * (self.ambient_gamma + 1) * v_r**2
                + self.ambient_gamma * v_r * (1 + (0.25 * (self.ambient_gamma + 1)) ** 2 * v_r**2) ** 0.5
            ) * p_a_bar

        else:
            return 0

    def solve(
        self,
        load_fraction: float,
        charge_mass_ratio: float,
        labda_1: float | None = None,
        labda_2: float | None = None,
        cc: float | None = None,
        it: int = 0,
        length_gun: float | None = None,
        known_bore: bool = False,
        suppress: bool = False,  # suppress design velocity exceeded before peak pressure check
        max_iteration: int = MAX_ITER,
        **_,
    ) -> tuple[float, float]:
        if any((charge_mass_ratio <= 0, load_fraction <= 0, load_fraction > 1)):
            raise ValueError("Invalid parameters to solve constrained design problem")

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

        if self.ambient_rho:
            c_a_bar = (self.ambient_gamma * self.ambient_p / self.ambient_rho) ** 0.5 / v_j
            p_a_bar = self.ambient_p / (self.f * delta)
        else:
            c_a_bar, p_a_bar = 0, 0

        if v_j < self.v_d and not known_bore:
            raise ValueError(
                f"Propellant load too low to achieve design velocity. The 2nd ballistic limit for this loading \
conditions is {v_j:.4g} m/s."
            )

        psi_0 = (1 / delta - 1 / self.rho_p) / (self.f / self.p_0 + self.alpha - 1 / self.rho_p)
        z_0, _ = dekker(lambda _z: self.propellant.f_psi_z(_z) - psi_0, 0, 1, y_rel_tol=self.tol, y_abs_tol=self.tol**2)

        def _f_p_bar(z: float, l_bar: float, v_bar: float) -> float:
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
                    raise ValueError(f"unknown control {self.control}")

        """
        step 1, find grain size that satisfies design pressure
        """
        p_bar_s = _f_p_bar(z_0, 0, 0)
        p_bar_d = self.p_d / (self.f * delta)

        if p_bar_d < p_bar_s:
            raise ValueError(
                "Interior ballistics of conventional gun precludes" + " a pressure lower than that at shot start."
            )
        l_bar_d = self.max_length / l_0

        def _abort_z(
            x: float, ys: tuple[float, float, float], record: list[tuple[float, tuple[float, float, float]]]
        ) -> bool:
            z, (t_bar, l_bar, v_bar) = x, ys
            p_bar = _f_p_bar(z, l_bar, v_bar)
            o_x, o_ys = record[-1]
            o_z = o_x
            ot_bar, ol_bar, ov_bar = o_ys
            op_bar = _f_p_bar(o_z, ol_bar, ov_bar)

            return (p_bar < op_bar) or (p_bar > 2 * p_bar_d)

        def _f_p_e_1(e_1: float) -> tuple[float, float, float, float, float]:
            """
            calculate either the peak pressure, given the arc thickness,
            or until the system develops 2x design pressure.
            """

            b_e_1 = (
                self.s**2 * e_1**2 / (self.f * phi * w * self.m * self.u_1**2) * (self.f * delta) ** (2 * (1 - self.n))
            )

            def _ode_z(z: float, tlv: tuple[float, float, float], __: float) -> tuple[float, float, float]:
                """burnup domain ode of internal ballistics"""
                t_bar, l_bar, v_bar = tlv
                psi = self.f_psi_z(z)
                l_psi_bar = 1 - delta / self.rho_p - delta * (self.alpha - 1 / self.rho_p) * psi

                p_bar = max((psi - v_bar**2) / (l_bar + l_psi_bar), p_a_bar)

                dt_bar = (2 * b_e_1 / self.theta) ** 0.5 * p_bar**-self.n  # dt_bar/dz
                dl_bar = v_bar * dt_bar
                dv_bar = 0.5 * self.theta * (p_bar - self._f_p_bar_ad(v_bar, c_a_bar, p_a_bar)) * dt_bar

                return dt_bar, dl_bar, dv_bar

            record = [(z_0, (0, 0, 0))]

            z_k, (t_bar_j, l_bar_j, v_bar_j) = rkf(
                d_func=_ode_z,
                ini_val=(0, 0, 0),
                x_0=z_0,
                x_1=self.z_b,
                rel_tol=self.tol,
                abs_tol=self.tol**2,
                abort_func=_abort_z,
                record=record,
            )

            p_bar_j = _f_p_bar(z_k, l_bar_j, v_bar_j)

            if p_bar_j >= 2 * p_bar_d:  # case for abort due to excessive pressure
                return p_bar_j - p_bar_d, z_k, t_bar_j, l_bar_j, v_bar_j

            def _f_p_z(z: float) -> tuple[float, float, float, float, float]:
                j = record.index([v for v in record if v[0] <= z][-1])
                x = record[j][0]
                ys = record[j][1]

                r = []
                t_bar, l_bar, v_bar = rkf(
                    d_func=_ode_z, ini_val=ys, x_0=x, x_1=z, rel_tol=self.tol, abs_tol=self.tol**2, record=r
                )[1]
                xs = [v[0] for v in record]
                record.extend(v for v in r if v[0] not in xs)
                record.sort()
                return _f_p_bar(z, l_bar, v_bar) - p_bar_d, z, t_bar, l_bar, v_bar

            if len(record) > 1:
                z_j = record[-2][0]
            else:
                z_j = z_0

            z_p = 0.5 * sum(
                gss(lambda z: _f_p_z(z)[0], z_j, z_k, y_rel_tol=self.tol, y_abs_tol=self.tol**2, find_min=False)
            )

            return _f_p_z(z_p)

        probe_web = 0.5 * self.min_web
        dp_bar_probe = _f_p_e_1(probe_web)[0]

        if dp_bar_probe < 0:
            raise ValueError("Design pressure cannot be achieved by varying web down to minimum")

        while dp_bar_probe > 0:
            probe_web *= 2
            dp_bar_probe = _f_p_e_1(probe_web)[0]

        e_1_solved, _ = dekker(
            lambda _e_1: _f_p_e_1(_e_1)[0],
            probe_web,  # >0
            0.5 * probe_web,  # ?0
            y_rel_tol=self.tol,
            y_abs_tol=p_bar_d * self.tol,
        )  # this is the e_1 that satisfies the pressure specification.

        p_bar_dev, z_i, t_bar_i, l_bar_i, v_bar_i = _f_p_e_1(e_1_solved)

        if known_bore:
            return e_1_solved, length_gun

        if v_j * v_bar_i > self.v_d:
            if suppress:
                logger.warning("velocity target point occurred before peak pressure point.")
                logger.warning("this is currently being suppressed due to program control.")
            else:
                raise ValueError(f"Design velocity exceeded before peak pressure point (V = {v_bar_i * v_j:.4g} m/s).")

        """
        step 2, find the requisite muzzle length to achieve design velocity
        """

        b = (
            self.s**2
            * e_1_solved**2
            / (self.f * phi * w * self.m * self.u_1**2)
            * (self.f * delta) ** (2 * (1 - self.n))
        )

        def _ode_v(v_bar: float, tzl: tuple[float, float, float], __: float) -> tuple[float, float, float]:
            t_bar, z, l_bar = tzl
            psi = self.f_psi_z(z)

            l_psi_bar = 1 - delta / self.rho_p - delta * (self.alpha - 1 / self.rho_p) * psi
            p_bar = max((psi - v_bar**2) / (l_bar + l_psi_bar), p_a_bar)
            dt_bar = 2 / (self.theta * (p_bar - self._f_p_bar_ad(v_bar, c_a_bar, p_a_bar)))
            dz = dt_bar * (0.5 * self.theta / b) ** 0.5 * p_bar**self.n
            dl_bar = v_bar * dt_bar

            return dt_bar, dz, dl_bar

        def _abort_v(
            _: float, ys: tuple[float, float, float], record: list[tuple[float, tuple[float, float, float]]]
        ) -> bool:
            t_bar, _, l_bar = ys
            _, (ot_bar, _, _) = record[-1]
            return l_bar > l_bar_d or t_bar < ot_bar

        vtzl_record = [(v_bar_i, (t_bar_i, z_i, l_bar_i))]
        try:
            v_bar_g, (t_bar_g, z_g, l_bar_g) = rkf(
                d_func=_ode_v,
                ini_val=(t_bar_i, z_i, l_bar_i),
                x_0=v_bar_i,
                x_1=v_bar_d,
                rel_tol=self.tol,
                abs_tol=self.tol**2,
                abort_func=_abort_v,
                record=vtzl_record,
            )

        except ValueError:
            v_bar_m, (t_bar_m, z_m, l_bar_m) = vtzl_record[-1]
            pmax = _f_p_bar(z_m, l_bar_m, v_bar_m) * self.f * delta
            vmax = v_bar_m * v_j
            lmax = l_bar_m * l_0
            raise ValueError(
                "Integration appears to be approaching asymptote, "
                + f"last calculated to v = {vmax:.4g} m/s, x = {lmax:.4g} m, p = {pmax * 1e-6:.4g} MPa. "
                + "This indicates an excessive velocity target relative to pressure developed."
            )

        p_bar_g = _f_p_bar(z_g, l_bar_g, v_bar_g)

        v_g = v_bar_g * v_j
        l_g = l_bar_g * l_0
        p_g = p_bar_g * self.f * delta
        if l_bar_g > l_bar_d:
            raise ValueError(
                "Solution requires excessive tube length, last calculated to "
                + f"v = {v_g:.4g} m/s, x = {l_g:.4g} m, p = {p_g * 1e-6:.4g} MPa."
            )

        cc_n = 1 - (1 - 1 / self.chi_k) * log(l_bar_g + 1) / l_bar_g
        l_g = l_bar_g * l_0
        if abs((l_bar_g - l_bar_g_0) / min(l_bar_g, l_bar_g_0)) > self.tol and it < max_iteration:
            return self.solve(
                load_fraction=load_fraction,
                charge_mass_ratio=charge_mass_ratio,
                labda_1=labda_1,
                labda_2=labda_2,
                cc=cc_n,
                it=it + 1,
                length_gun=l_g,
                known_bore=known_bore,
                suppress=suppress,
            )
        else:
            logger.info(
                f"ω/m = {charge_mass_ratio:.2f}, Δ/ρ = {load_fraction:.2f} -> e_1 = {e_1_solved * 1e3:.2f} mm, l_g = {l_g * 1e3:.0f} mm ({it} iterations)"
            )
            return e_1_solved, l_g

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
        logger.info("solving minimum chamber volume for constraint.")

        if self.sol == SOL_LAGRANGE:
            labda_1, labda_2 = 1 / 2, 1 / 3
        elif self.sol == SOL_PIDDUCK:
            labda_1, labda_2 = pidduck(charge_mass_ratio / self.phi_1, self.theta + 1, self.tol)
        elif self.sol == SOL_MAMONTOV:
            labda_1, labda_2 = pidduck(charge_mass_ratio / self.phi_1, 1, self.tol)

        def _f(load_fraction: float) -> tuple[float, float, float]:

            e_1_delta, l_g_delta = self.solve(
                load_fraction=load_fraction,
                charge_mass_ratio=charge_mass_ratio,
                labda_1=labda_1,
                labda_2=labda_2,
                known_bore=False,
                suppress=True,
            )

            return (
                e_1_delta,
                l_g_delta,
                l_g_delta + (self.m * charge_mass_ratio / (self.rho_p * load_fraction)) / self.s,
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

        logger.info(f"solution constrained to Δ/ρ : {low:.3%} - {high:.3%}")
        lf_low, lf_high = gss(
            lambda load_fraction: _f(load_fraction)[_f_index], low, high, x_tol=self.tol, find_min=True
        )
        lf = 0.5 * (lf_high + lf_low)
        e_1, l_g, _ = _f(lf)
        logger.info(f"Optimal Δ/ρ = {lf:.2f}")
        return lf, e_1, l_g


if __name__ == "__main__":
    pass
