import logging
from math import floor, inf, log, pi
from random import uniform

from .generics import DelegatesPropellant
from . import POINT_PEAK_AVG, POINT_PEAK_BREECH, POINT_PEAK_SHOT, POINT_PEAK_STAG, Points, MAX_GUESSES
from .num import RKF78, dekker, gss
from .prop import Propellant
from .recoilless import Recoilless

logger = logging.getLogger(__name__)


class ConstrainedRecoilless(DelegatesPropellant):
    def __init__(
        self,
        caliber,
        shot_mass,
        propellant,
        start_pressure,
        drag_coefficient,
        nozzle_expansion,
        nozzle_efficiency,
        chambrage,
        tol,
        design_pressure,
        design_velocity,
        min_web=1e-6,
        max_length=1e3,
        ambient_rho=1.204,
        ambient_p=101.325e3,
        ambient_gamma=1.4,
        control: Points = POINT_PEAK_AVG,
        **_,
    ):
        # constants for constrained designs

        super().__init__(propellant=propellant)
        logger.info("initializing constrained recoiless gun object.")

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

        logger.info("constraints successfully initialized.")

    def solve(
        self,
        load_fraction,
        charge_mass_ratio,
        length_gun=None,
        known_bore=False,
        suppress=False,  # suppress design velocity exceeded before peak pressure check
        progress_queue=None,
        **_,
    ):
        logger.info("solving constraint at specified load fraction")
        if any((charge_mass_ratio <= 0, load_fraction <= 0, load_fraction > 1)):
            raise ValueError("Invalid parameters to solve constrained design problem")
        """
        minWeb  : represents minimum possible grain size
        """
        if progress_queue is not None:
            progress_queue.put(1)

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

        omega = m * charge_mass_ratio
        v_0 = omega / (rho_p * load_fraction)
        delta = omega / v_0
        l_0 = v_0 / s

        gamma = theta + 1

        phi = phi_1 + omega / (3 * m)

        s_j_bar = 1 / (Recoilless.get_cf(gamma, a_bar, tol) * chi_0)
        if s_j_bar > chi_k:
            raise ValueError(
                "Achieving recoiless condition necessitates"
                + " a larger throat area than could be fit into"
                + " breech face."
            )
        s_j = s_j_bar * s

        k_0 = (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1))) * gamma**0.5

        phi_2 = 1
        c_a = (0.5 * theta * phi * m / omega) ** 0.5 * k_0 * phi_2  # flow rate value

        """
        it is impossible to account for the chambrage effect given unspecified
        barrel length, in our formulation
        """
        v_j = (2 * f * omega / (theta * phi * m)) ** 0.5

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
                + " and recoiless guns only achieve a part of that as well."
            )

        psi_0 = (1 / delta - 1 / rho_p) / (f / p_0 + alpha - 1 / rho_p)
        z_0, _ = dekker(lambda z: self.propellant.f_psi_z(z) - psi_0, 0, 1, y_rel_tol=tol, y_abs_tol=tol**2)
        logger.info("solved starting burnup.")

        # p_bar_0 = p_0 / (f * Delta)

        def _f_p_bar(z, l_bar, v_bar, eta, tau):
            psi = f_psi_z(z)
            l_psi_bar = 1 - delta * ((1 - psi) / rho_p + alpha * (psi - eta))
            p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

            if self.control == POINT_PEAK_AVG:
                return p_bar

            else:
                y = omega * eta
                m_dot = c_a * v_j * s_j * p_bar * delta / (tau**0.5)
                vb = m_dot * (v_0 + s * l_bar * l_0) / (sb * (omega - y))

                h_1, h_2 = vb / (v_j * v_bar) if v_bar != 0 else inf, 2 * phi_1 * m / (omega - y) + 1

                h = min(h_1, h_2)

                p_s_bar = p_bar / (1 + (omega - y) / (3 * phi_1 * m) * (1 - 0.5 * h))
                if self.control == POINT_PEAK_SHOT:
                    return p_s_bar
                elif self.control == POINT_PEAK_STAG:
                    return p_s_bar * (1 + (omega - y) / (2 * phi_1 * m) * (1 + h) ** -1)
                elif self.control == POINT_PEAK_BREECH:
                    return p_s_bar * (1 + (omega - y) / (2 * phi_1 * m) * (1 - h)) if h == h_1 else 0
                else:
                    raise ValueError("tag unhandled.")

        p_bar_d = p_d / (f * delta)  # convert to unitless
        l_bar_d = self.max_length / l_0

        """
        step 1, find grain size that satisifies design pressure
        """

        # noinspection PyUnusedLocal
        def abort_z(x, ys, record):
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
            b = (s**2 * e_1**2) / (f * phi * omega * m * u_1**2) * (f * delta) ** (2 * (1 - n))

            # integrate this to end of burn

            def _ode_z(z, t_bar, l_bar, v_bar, eta, tau, _):
                """burnout domain ode of internal ballistics"""
                psi = f_psi_z(z)
                dpsi = f_sigma_z(z)  # dpsi/dZ

                l_psi_bar = 1 - delta * ((1 - psi) / rho_p + alpha * (psi - eta))
                p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

                if c_a_bar != 0 and v_bar > 0:
                    v_r = v_bar / c_a_bar
                    p_d_bar = (
                        +0.25 * gamma_1 * (gamma_1 + 1) * v_r**2
                        + gamma_1 * v_r * (1 + (0.25 * (gamma_1 + 1)) ** 2 * v_r**2) ** 0.5
                    ) * p_a_bar
                else:
                    p_d_bar = 0

                if z <= z_b:
                    dt_bar = (2 * b / theta) ** 0.5 * p_bar**-n
                    dl_bar = v_bar * dt_bar
                    dv_bar = 0.5 * theta * (p_bar - p_d_bar) * dt_bar

                else:
                    # technically speaking it is undefined in this area
                    dt_bar = 0  # dt_bar/dZ
                    dl_bar = 0  # dl_bar/dZ
                    dv_bar = 0  # dv_bar/dZ

                deta = c_a * s_j_bar * p_bar / tau**0.5 * dt_bar  # deta / dZ
                dtau = ((1 - tau) * dpsi - 2 * v_bar * dv_bar - theta * tau * deta) / (psi - eta)

                return dt_bar, dl_bar, dv_bar, deta, dtau

            # stepVanished = False
            record = [[z_0, [0, 0, 0, 0, 1]]]
            try:
                (z_j, (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j), _) = RKF78(
                    dFunc=_ode_z,
                    iniVal=(0, 0, 0, 0, 1),
                    x_0=z_0,
                    x_1=z_b,
                    relTol=tol,
                    absTol=tol**2,
                    abortFunc=abort_z,
                    record=record,
                )

                if z_j not in [line[0] for line in record]:
                    record.append([z_j, [t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j]])
            except ValueError:
                z_j, (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j) = record[-1]

            p_bar_j = _f_p_bar(z_j, l_bar_j, v_bar_j, eta_j, tau_j)

            peak = None

            # find the last peak
            if len(record) > 2:
                p_bars = [_f_p_bar(z, l_bar, v_bar, eta, tau) for (z, (t_bar, l_bar, v_bar, eta, tau)) in record]

                for i, (l, c, r) in enumerate(zip(p_bars[:-2], p_bars[1:-1], p_bars[2:])):
                    if l < c and c > r:
                        peak = i + 1

            if peak is None:
                # no peak, so it suffice to compare the end points.
                if p_bar_j == 0:
                    p_bar_j = inf
                return p_bar_j - p_bar_d, record[-1][0], *record[-1][-1]
            else:  # peak exist, must compare the peak and the two end points.
                z_i = record[peak - 1][0]
                z_j = record[peak + 1][0]

                def _f_p_z(z):
                    i = record.index([v for v in record if v[0] <= z][-1])
                    x = record[i][0]
                    ys = record[i][1]

                    r = []
                    _, (t_bar, l_bar, v_bar, eta, tau), _ = RKF78(
                        dFunc=_ode_z,
                        iniVal=ys,
                        x_0=x,
                        x_1=z,
                        relTol=tol,
                        absTol=tol**2,
                        record=r,
                    )

                    xs = [v[0] for v in record]
                    record.extend(v for v in r if v[0] not in xs)
                    record.sort()
                    return _f_p_bar(z, l_bar, v_bar, eta, tau)

                z_1, z_2 = gss(_f_p_z, z_i, z_j, y_rel_tol=tol, findMin=False)
                z_p = 0.5 * (z_1 + z_2)

                p_bar_p = _f_p_z(z_p)
                i = [line[0] for line in record].index(z_p)

                return p_bar_p - p_bar_d, record[i][0], *record[i][-1]

        dp_bar_probe, z, *_ = _f_p_e_1(self.min_web)
        probe_web = self.min_web

        if dp_bar_probe < 0:
            raise ValueError(
                "Design pressure cannot be achieved by varying web down to minimum. "
                + "Peak pressure found at phi = {:.4g} at {:.4g} MPa".format(
                    f_psi_z(z), (dp_bar_probe + p_bar_d) * 1e-6 * f * delta
                )
            )

        while dp_bar_probe > 0:
            probe_web *= 2
            dp_bar_probe = _f_p_e_1(probe_web)[0]

        def fr(x):
            progress_queue.put(round(x * 100))

        e_1, e_1_2 = dekker(
            lambda web: _f_p_e_1(web)[0],
            probe_web,
            0.5 * probe_web,
            y_abs_tol=p_bar_d * tol,
            f_report=fr if progress_queue is not None else None,
        )  # this is the e_1 that satisifies the pressure specification.

        """
        e_1 and e_2 brackets the true solution
        """

        dp_bar_i, *vals_1 = _f_p_e_1(e_1)

        z_i, t_bar_i, l_bar_i, v_bar_i, eta_i, tau_i = vals_1

        logger.info("solved web satisfying pressure constraint.")

        if known_bore:
            if progress_queue is not None:
                progress_queue.put(100)
            return e_1, length_gun

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

        b = s**2 * e_1**2 / (f * phi * omega * m * u_1**2) * (f * delta) ** (2 * (1 - n))

        def _ode_v(v_bar, _, z, l_bar, eta, tau, __):
            psi = f_psi_z(z)
            dpsi = f_sigma_z(z)  # dpsi/dZ

            l_psi_bar = 1 - delta * ((1 - psi) / rho_p + alpha * (psi - eta))
            p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

            if c_a_bar != 0 and v_bar > 0:
                v_r = v_bar / c_a_bar
                p_d_bar = (
                    +0.25 * gamma_1 * (gamma_1 + 1) * v_r**2
                    + gamma_1 * v_r * (1 + (0.25 * (gamma_1 + 1)) ** 2 * v_r**2) ** 0.5
                ) * p_a_bar
            else:
                p_d_bar = 0

            dt_bar = 2 / (theta * (p_bar - p_d_bar))

            dz = dt_bar * (0.5 * theta / b) ** 0.5 * p_bar**n

            dl_bar = v_bar * dt_bar

            deta = c_a * s_j_bar * p_bar / tau**0.5 * dt_bar  # deta / dv_bar
            dtau = (
                (1 - tau) * (dpsi * dz) - 2 * v_bar - theta * tau * deta  # dZ/dt_bar  # dv_bar/dt_bar  # deta/dt_bar
            ) / (
                psi - eta
            )  # dtau/dt_bar

            return dt_bar, dz, dl_bar, deta, dtau

        def abort_v(x, ys, record):
            # v_bar = x
            t_bar, _, l_bar, _, _ = ys

            ot_bar, *_ = record[-1][-1]
            return l_bar > l_bar_d or t_bar < ot_bar

        vtzlet_record = [[v_bar_i, (t_bar_i, z_i, l_bar_i, eta_i, tau_i)]]
        try:
            (
                v_bar_g,
                (t_bar_g, z_g, l_bar_g, eta_g, tau_g),
                _,
            ) = RKF78(
                dFunc=_ode_v,
                iniVal=(t_bar_i, z_i, l_bar_i, eta_i, tau_i),
                x_0=v_bar_i,
                x_1=v_bar_d,
                relTol=tol,
                absTol=tol**2,
                abortFunc=abort_v,
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

        if progress_queue is not None:
            progress_queue.put(100)
        logger.info("solved tube length to satisfy design velocity.")

        return e_1, l_bar_g * l_0

    def find_min_v(self, charge_mass_ratio, progress_queue=None, **_):
        """
        find the minimum volume solution.
        """
        if progress_queue is not None:
            progress_queue.put(1)

        logger.info("solving minimum chamber volume for constraint.")

        """
        Step 1, find a valid range of values for load fraction,
        using psuedo-bisection.

        high lf -> high web
        low lf -> low web
        """
        omega = self.m * charge_mass_ratio
        rho_p = self.rho_p
        s = self.s
        solve = self.solve
        tol = self.tol

        def f(lf):
            logger.info(f"dispatching calculation at load factor = {lf:.3%}")
            v_0 = omega / (rho_p * lf)
            l_0 = v_0 / s

            e_1, l_g = solve(
                load_fraction=lf,
                charge_mass_ratio=charge_mass_ratio,
                known_bore=False,
                suppress=True,
            )
            return e_1, (l_g + l_0), l_g

        logger.info(f"attempting to find valid load fraction with {MAX_GUESSES} guesses.")

        records = []
        for i in range(MAX_GUESSES):
            start_probe = uniform(tol, 1 - tol)
            try:
                _, lt_i, lg_i = f(start_probe)
                records.append((start_probe, lt_i))
                break
            except ValueError:
                if progress_queue is not None:
                    progress_queue.put(round(i / MAX_GUESSES * 33))
        else:
            raise ValueError(
                "Unable to find any valid load" + " fraction with {:d} random samples.".format(MAX_GUESSES)
            )

        if progress_queue is not None:
            progress_queue.put(33)

        logger.info(f"found valid load fraction at {start_probe:.3%}.")
        logger.info("attempting to find minimum valid load fraction.")

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
                if progress_queue is not None:
                    progress_queue.put(round(k / n * 17) + 33)
                k += 1
            finally:
                new_low = probe + delta_low

        low = probe

        logger.info(f"found minimum valid load fraction at {low:.3%}.")
        logger.info("attempting to find maximum valid load fraction.")

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
                if progress_queue is not None:
                    progress_queue.put(round(k / n * 16) + 50)
                k += 1
            finally:
                new_high = probe + delta_high

        high = probe

        logger.info(f"found maximum valid load fraction at {high:.3%}.")

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
        logger.info(f"found maximum valid load fraction at {high:.3%}.")

        def fr(x):
            progress_queue.put(round(x * 33) + 66)

        lf_low, lf_high = gss(
            lambda lf: f(lf)[1],
            low,
            high,
            y_rel_tol=tol,
            findMin=True,
            f_report=fr if progress_queue is not None else None,
        )

        lf = 0.5 * (lf_high + lf_low)
        logger.info(f"validating found optimal at {lf:.3%}")
        e_1, l_t, l_g = f(lf)

        if progress_queue is not None:
            progress_queue.put(100)
        logger.info("minimum chamber solution found.")
        return lf, e_1, l_g


if __name__ == "__main__":
    from prop import Composition, MultPerfGeometry, SimpleGeometry

    compositions = Composition.read_file("data/propellants.csv")
    S22 = compositions["ATK PRD(S)22"]
    M8 = compositions["M8"]
    M1 = compositions["M1"]

    S22S = Propellant(S22, SimpleGeometry.TUBE, 1, 2.5)
    M8S = Propellant(M8, SimpleGeometry.TUBE, 1, 2.5)
    M17P = Propellant(M1, MultPerfGeometry.SEVEN_PERF_CYLINDER, 1.0, 2.50)
    test = ConstrainedRecoilless(
        caliber=93e-3,
        shot_mass=4,
        propellant=M17P,
        start_pressure=10e6,
        drag_coefficient=5e-2,
        design_pressure=50e6,
        design_velocity=120,
        nozzle_expansion=4,
        nozzle_efficiency=0.92,
        chambrage=1.5,
        tol=1e-3,
        max_length=100,
        min_web=1e-6,
    )

    datas = []
    for i in range(10):
        datas.append(test.find_min_v(charge_mass_ratio=0.309 / 4))

    from tabulate import tabulate

    means = [sum(x) / len(datas) for x in zip(*datas)]

    delta = []
    for line in datas:
        delta.append((v - m) / m for v, m in zip(line, means))

    print(tabulate(datas, headers=("load fract.", "web", "length")))
    print(*means)
    print(tabulate(delta, headers=("load fract.", "web", "length"), floatfmt=".3e"))
