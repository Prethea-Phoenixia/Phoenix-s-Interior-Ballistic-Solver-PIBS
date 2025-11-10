from __future__ import annotations

import logging
from math import floor, inf, log, pi
from random import uniform

from . import (
    MAX_GUESSES,
    POINT_PEAK_AVG,
    POINT_PEAK_BREECH,
    POINT_PEAK_SHOT,
    POINT_PEAK_STAG,
    Points,
)
from .generics import DelegatesPropellant
from .num import dekker, gss, rkf
from .prop import Propellant
from .recoilless import Recoilless

logger = logging.getLogger(__name__)


class ConstrainedRecoilless(DelegatesPropellant):
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
        ambient_rho: float = 1.204,
        ambient_p: float = 101.325e3,
        ambient_gamma: float = 1.4,
        control: Points = POINT_PEAK_AVG,
        traveling_charge: bool = False,
        **_,
    ):
        # constants for constrained designs

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

        ambient_p = max(ambient_p, 1)
        ambient_gamma = max(ambient_gamma, 1)

        self.s = (0.5 * caliber) ** 2 * pi
        self.m = shot_mass
        self.propellant = propellant
        self.p_0 = start_pressure
        self.phi_1 = 1 / (1 - drag_coefficient)

        # design limits
        self.p_d = design_pressure
        self.v_d = design_velocity

        self.chi_0 = nozzle_efficiency
        self.a_bar = nozzle_expansion

        self.chi_k = chambrage

        self.tol = tol
        self.min_web = min_web
        self.max_length = max_length
        self.ambient_rho = ambient_rho
        self.ambient_p = ambient_p
        self.ambient_gamma = ambient_gamma
        self.control = control

        self.traveling_charge = traveling_charge

    def solve(
        self,
        load_fraction: float,
        charge_mass_ratio: float,
        length_gun: float | None = None,
        known_bore: bool = False,
        suppress: bool = False,
        **_,
    ):
        if any((charge_mass_ratio <= 0, load_fraction <= 0, load_fraction > 1)):
            raise ValueError("Invalid parameters to solve constrained design problem")
        """
        minWeb  : represents minimum possible grain size
        """
        m = self.m
        rho_p = self.rho_p
        theta = self.theta
        f = self.f
        chi_k = self.chi_k
        s = self.s
        phi_1 = self.phi_1
        p_0 = self.p_0
        v_d = self.v_d
        p_d = self.p_d
        u_1 = self.u_1
        n = self.n
        alpha = self.alpha
        z_b = self.z_b

        f_psi_z = self.f_psi_z
        f_sigma_z = self.f_sigma_z

        chi_0 = self.chi_0
        a_bar = self.a_bar

        sb = s * chi_k
        tol = self.tol

        is_tc = self.traveling_charge

        w = m * charge_mass_ratio
        vol_0 = w / (rho_p * load_fraction)
        delta = w / vol_0
        l_0 = vol_0 / s

        gamma = theta + 1

        phi = phi_1 + w / (3 * m)

        s_j_bar = 1 / (Recoilless.get_cf(gamma, a_bar, tol) * chi_0)
        if s_j_bar > chi_k:
            raise ValueError(
                "Achieving recoilless condition necessitates a larger throat area than could be fit into breech face."
            )
        s_j = s_j_bar * s

        k_0 = (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1))) * gamma**0.5

        phi_2 = 1
        c_a = (0.5 * theta * phi * m / w) ** 0.5 * k_0 * phi_2  # flow rate value

        """
        it is impossible to account for the chambrage effect given unspecified
        barrel length, in our formulation
        """
        v_j = (2 * f * w / (theta * phi * m)) ** 0.5

        if self.ambient_rho != 0:
            c_a_bar = (self.ambient_gamma * self.ambient_p / self.ambient_rho) ** 0.5 / v_j
            p_a_bar = self.ambient_p / (f * delta)
        else:
            c_a_bar = 0
            p_a_bar = 0

        gamma_1 = self.ambient_gamma

        if v_j < v_d and not known_bore:
            raise ValueError(
                "Propellant load too low to achieve design velocity. "
                + " The 2nd ballistic limit for this loading conditions is"
                + " {:.4g} m/s,".format(v_j)
                + " and recoilless guns only achieve a part of that as well."
            )

        psi_0 = (1 / delta - 1 / rho_p) / (f / p_0 + alpha - 1 / rho_p)
        z_0, _ = dekker(lambda _z: self.propellant.f_psi_z(_z), 0, 1, y=psi_0, y_rel_tol=tol, y_abs_tol=tol**2)

        def _f_p_bar(z, l_bar, v_bar, eta, tau):
            psi = f_psi_z(z)
            l_psi_bar = 1 - delta * ((1 - psi) / rho_p + alpha * (psi - eta))
            p_bar = max(tau / (l_bar + l_psi_bar) * (psi - eta), p_a_bar)

            if self.control == POINT_PEAK_AVG:
                return p_bar

            else:
                y = w * eta
                m_dot = c_a * v_j * s_j * p_bar * delta / (tau**0.5)
                vb = m_dot * (vol_0 + s * l_bar * l_0) / (sb * (w - y))

                h_1, h_2 = vb / (v_j * v_bar) if v_bar != 0 else inf, 2 * phi_1 * m / (w - y) + 1

                h = min(h_1, h_2)

                p_s_bar = p_bar / (1 + (w - y) / (3 * phi_1 * m) * (1 - 0.5 * h))
                if self.control == POINT_PEAK_SHOT:
                    return p_s_bar
                elif self.control == POINT_PEAK_STAG:
                    return p_s_bar * (1 + (w - y) / (2 * phi_1 * m) * (1 + h) ** -1)
                elif self.control == POINT_PEAK_BREECH:
                    return p_s_bar * (1 + (w - y) / (2 * phi_1 * m) * (1 - h)) if h == h_1 else 0
                else:
                    raise ValueError("tag unhandled.")

        p_bar_d = p_d / (f * delta)  # convert to unitless
        l_bar_d = self.max_length / l_0

        """
        step 1, find grain size that satisfies design pressure
        """

        def _abort_z(x, ys, _):
            z, (_, l_bar, v_bar, eta, tau) = x, ys
            p_bar = _f_p_bar(z, l_bar, v_bar, eta, tau)
            return (p_bar > 2 * p_bar_d) or l_bar > l_bar_d

        def _f_p_e_1(e_1):
            """
            Find pressure maximum bracketed by the range of
            Z_0 < Z < Z_d
            l_bar < l_bar_d,
            p_bar < 2 * p_bar_d.
            """
            b_e_1 = (s**2 * e_1**2) / (f * phi * w * m * u_1**2) * (f * delta) ** (2 * (1 - n))

            # integrate this to end of burn

            def _ode_z(
                z: float, tlvetatau: tuple[float, float, float, float, float], _: float
            ) -> tuple[float, float, float, float, float]:
                """burnout domain ode of internal ballistics"""
                t_bar, l_bar, v_bar, eta, tau = tlvetatau
                psi = f_psi_z(z)
                dpsi = f_sigma_z(z)  # dpsi/dZ

                l_psi_bar = 1 - delta * ((1 - psi) / rho_p + alpha * (psi - eta))
                p_bar = max(tau / (l_bar + l_psi_bar) * (psi - eta), p_a_bar)

                if c_a_bar != 0 and v_bar > 0:
                    v_r = v_bar / c_a_bar
                    p_d_bar = (
                        +0.25 * gamma_1 * (gamma_1 + 1) * v_r**2
                        + gamma_1 * v_r * (1 + (0.25 * (gamma_1 + 1)) ** 2 * v_r**2) ** 0.5
                    ) * p_a_bar
                else:
                    p_d_bar = 0

                if z <= z_b:
                    dt_bar = (2 * b_e_1 / theta) ** 0.5 * p_bar**-n
                    dl_bar = v_bar * dt_bar
                    dv_bar = 0.5 * theta * (p_bar - p_d_bar) * dt_bar
                    dv_bar /= (1 + w / m * (1 - psi)) if is_tc else 1

                else:
                    # technically speaking it is undefined in this area
                    dt_bar = 0  # dt_bar/dZ
                    dl_bar = 0  # dl_bar/dZ
                    dv_bar = 0  # dv_bar/dZ

                deta = c_a * s_j_bar * p_bar / tau**0.5 * dt_bar  # deta / dZ
                dtau = ((1 - tau) * dpsi - 2 * v_bar * dv_bar - theta * tau * deta) / (psi - eta)

                return dt_bar, dl_bar, dv_bar, deta, dtau

            # stepVanished = False
            record = [(z_0, (0.0, 0.0, 0.0, 0.0, 1.0))]
            try:
                z_k, (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j) = rkf(
                    d_func=_ode_z,
                    ini_val=(0.0, 0.0, 0.0, 0.0, 1.0),
                    x_0=z_0,
                    x_1=z_b,
                    rel_tol=tol,
                    abs_tol=tol**2,
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
                return p_bar_j - p_bar_d, record[-1][0], *record[-1][-1]
            else:  # peak exist, must compare the peak and the two end points.
                z_j = record[peak - 1][0]
                z_k = record[peak + 1][0]

                def _f_p_z(z):
                    i = record.index([v for v in record if v[0] <= z][-1])
                    x = record[i][0]
                    ys = record[i][1]

                    r = []
                    t_bar, l_bar, v_bar, eta, tau = rkf(
                        d_func=_ode_z, ini_val=ys, x_0=x, x_1=z, rel_tol=tol, abs_tol=tol**2, record=r
                    )[1]

                    xs = [v[0] for v in record]
                    record.extend(v for v in r if v[0] not in xs)
                    record.sort()
                    return _f_p_bar(z, l_bar, v_bar, eta, tau), z, t_bar, l_bar, v_bar, eta, tau

                z_1, z_2 = gss(lambda _z: _f_p_z(_z)[0], z_j, z_k, y_rel_tol=tol, find_min=False)
                z_p = 0.5 * (z_1 + z_2)
                p_bar_p, *vals = _f_p_z(z_p)
                return p_bar_p - p_bar_d, *vals

        probe_web = 0.5 * self.min_web
        dp_bar_probe, z_probe, *_ = _f_p_e_1(probe_web)

        if dp_bar_probe < 0:
            raise ValueError(
                "Design pressure cannot be achieved by varying web down to minimum. "
                + "Peak pressure found at phi = {:.4g} at {:.4g} MPa".format(
                    f_psi_z(z_probe), (dp_bar_probe + p_bar_d) * 1e-6 * f * delta
                )
            )

        while dp_bar_probe > 0:
            probe_web *= 2
            dp_bar_probe = _f_p_e_1(probe_web)[0]

        e_1_solved, _ = dekker(lambda _e_1: _f_p_e_1(_e_1)[0], probe_web, 0.5 * probe_web, y_rel_tol=tol)

        dp_bar_i, *vals_1 = _f_p_e_1(e_1_solved)
        z_i, t_bar_i, l_bar_i, v_bar_i, eta_i, tau_i = vals_1

        if known_bore:
            return e_1_solved, length_gun

        """
        step 2, find the requisite muzzle length to achieve design velocity
        """
        v_bar_d = v_d / v_j

        if v_bar_i > v_bar_d:
            if suppress:
                logger.warning("velocity target point occurred before peak pressure point.")
                logger.warning("this is currently being suppressed due to program control.")
            else:
                raise ValueError(f"Design velocity exceeded before peak pressure point (V = {v_bar_i * v_j:.4g} m/s).")

        b = s**2 * e_1_solved**2 / (f * phi * w * m * u_1**2) * (f * delta) ** (2 * (1 - n))

        def _ode_v(
            v_bar: float, tzletatau: tuple[float, float, float, float, float], __: float
        ) -> tuple[float, float, float, float, float]:
            _, z, l_bar, eta, tau = tzletatau
            psi = f_psi_z(z)

            dpsi = f_sigma_z(z)  # dpsi/dZ

            l_psi_bar = 1 - delta * ((1 - psi) / rho_p + alpha * (psi - eta))
            p_bar = max(tau / (l_bar + l_psi_bar) * (psi - eta), p_a_bar)

            if c_a_bar != 0 and v_bar > 0:
                v_r = v_bar / c_a_bar
                p_d_bar = (
                    +0.25 * gamma_1 * (gamma_1 + 1) * v_r**2
                    + gamma_1 * v_r * (1 + (0.25 * (gamma_1 + 1)) ** 2 * v_r**2) ** 0.5
                ) * p_a_bar
            else:
                p_d_bar = 0

            dt_bar = 2 / (theta * (p_bar - p_d_bar))
            dt_bar *= (1 + w / m * (1 - psi)) if is_tc else 1

            dz = dt_bar * (0.5 * theta / b) ** 0.5 * p_bar**n

            dl_bar = v_bar * dt_bar

            deta = c_a * s_j_bar * p_bar / tau**0.5 * dt_bar  # deta / dv_bar
            dtau = (
                (1 - tau) * (dpsi * dz) - 2 * v_bar - theta * tau * deta  # dZ/dt_bar  # dv_bar/dt_bar  # deta/dt_bar
            ) / (psi - eta)

            return dt_bar, dz, dl_bar, deta, dtau

        def _abort_v(_, ys, record):
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
                rel_tol=tol,
                abs_tol=tol**2,
                abort_func=_abort_v,
                record=vtzlet_record,
            )

        except ValueError:
            v_bar_m, (t_bar_m, z_m, l_bar_m, eta_m, tau_m) = vtzlet_record[-1]
            pmax = _f_p_bar(z_m, l_bar_m, v_bar_m, eta_m, tau_m) * f * delta
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
            f"ω/m = {charge_mass_ratio:.2f}, Δ/ρ = {load_fraction:.2f} -> e_1 = {e_1_solved * 1e3:.2f} mm, l_g = {l_g * 1e3:.0f} mm"
        )
        return e_1_solved, l_bar_g * l_0

    def find_min_v(self, charge_mass_ratio: float, max_guess: int = MAX_GUESSES, **_):
        """
        find the minimum volume solution.
        """

        w = self.m * charge_mass_ratio
        rho_p = self.rho_p
        s = self.s
        solve = self.solve
        tol = self.tol

        def f(load_fraction):
            vol_0 = w / (rho_p * load_fraction)
            l_0 = vol_0 / s

            e_1_delta, l_g_delta = solve(
                load_fraction=load_fraction, charge_mass_ratio=charge_mass_ratio, known_bore=False, suppress=True
            )
            return e_1_delta, (l_g_delta + l_0), l_g_delta

        records = []
        for i in range(max_guess):
            start_probe = uniform(tol, 1 - tol)
            try:
                _, lt_i, lg_i = f(start_probe)
                records.append((start_probe, lt_i))
                break
            except ValueError:
                pass
        else:
            raise ValueError(f"Unable to find any valid load fraction with {max_guess:d} random samples.")

        logger.info(f"valid Δ/ρ = {start_probe:.3%}.")

        low = tol
        probe = start_probe
        delta_low = low - probe
        new_low = probe + delta_low

        k, n = 0, floor(log(abs(delta_low) / tol, 2)) + 1
        while abs(2 * delta_low) > tol:
            try:
                _, lt_i, lg_i = f(new_low)
                records.append((new_low, lt_i))
                probe = new_low
            except ValueError:
                delta_low *= 0.5
                k += 1
            finally:
                new_low = probe + delta_low

        low = probe

        logger.info(f"min Δ/ρ = {low:.3%}.")

        high = 1 - tol
        probe = start_probe
        delta_high = high - probe
        new_high = probe + delta_high

        k, n = 0, floor(log(abs(delta_high) / tol, 2)) + 1
        while abs(2 * delta_high) > tol:
            try:
                _, lt_i, lg_i = f(new_high)
                records.append((new_high, lt_i))
                probe = new_high
            except ValueError:
                delta_high *= 0.5
                k += 1
            finally:
                new_high = probe + delta_high

        high = probe

        logger.info(f"max Δ/ρ = {high:.3%}.")

        if abs(high - low) < tol:
            raise ValueError("No range of values satisfying constraint.")

        if len(records) > 2:
            records.sort(key=lambda x: x[0])
            for l, m, h in zip(records[:-2], records[1:-1], records[2:]):
                if l[1] > m[1] and h[1] > m[1]:
                    low = l[0]
                    high = h[0]

        """
        Step 2, gss to min.
        """
        logger.info(f"solution constrained Δ/ρ : {low:.3%} - {high:.3%}")
        lf_low, lf_high = gss(lambda _lf: f(_lf)[1], low, high, y_rel_tol=tol, find_min=True)
        lf = 0.5 * (lf_high + lf_low)
        e_1, l_t, l_g = f(lf)
        logger.info(f"Optimal Δ/ρ = {lf:.2f}")
        return lf, e_1, l_g


if __name__ == "__main__":
    pass
