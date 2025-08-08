import logging
from math import floor, log, pi
from random import uniform

from . import (
    MAX_GUESSES,
    MAX_ITER,
    POINT_PEAK_AVG,
    POINT_PEAK_BREECH,
    POINT_PEAK_SHOT,
    SOL_LAGRANGE,
    SOL_MAMONTOV,
    SOL_PIDDUCK,
    Points,
    Solutions,
)
from .generics import DelegatesPropellant
from .gun import pidduck
from .num import RKF78, dekker, gss
from .prop import Propellant

"""
Machine-accuracy factor, determines that, if a numerical method
is used within another, then how much more accurate should the
inner method be called as compared to the outter. A small value
is necessary to prevent numerical instability form causing sporadic
appearance of outlier, at the cost of increased computation times.
"""


logger = logging.getLogger(__name__)


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
        # constants for constrained designs
        super().__init__(propellant=propellant)
        logger.info("initializing constrained gun object.")

        if any((caliber <= 0, shot_mass <= 0, start_pressure <= 0, drag_coefficient < 0, drag_coefficient >= 1)):
            raise ValueError("Invalid parameters for constrained design")

        if any((design_pressure <= 0, design_velocity <= 0)):
            raise ValueError("Invalid design constraint")

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

        logger.info("constraints successfully initialized.")

    def solve(
        self,
        load_fraction,
        charge_mass_ratio,
        labda_1=None,
        labda_2=None,
        cc=None,
        it=0,
        length_gun=None,
        known_bore=False,
        suppress=False,
        progress_queue=None,
        **_,
    ):  # suppress design velocity exceeded before peak pressure check

        logger.info(f"solving constraint at specified load fraction, cycle {it}/{MAX_ITER}.")
        if any((charge_mass_ratio <= 0, load_fraction <= 0, load_fraction > 1)):
            raise ValueError("Invalid parameters to solve constrained design problem")

        if progress_queue is not None:
            progress_queue.put(1)

        """
        minWeb  : represents minimum possible grain size
        """
        m = self.m
        rho_p = self.rho_p
        theta = self.theta
        f = self.f
        s = self.s
        phi_1 = self.phi_1
        p_0 = self.p_0
        v_d = self.v_d
        p_d = self.p_d
        u_1 = self.u_1
        n = self.n
        alpha = self.alpha
        z_b = self.z_b
        chi_k = self.chi_k
        f_psi_z = self.f_psi_z
        tol = self.tol

        omega = m * charge_mass_ratio
        v_0 = omega / (rho_p * load_fraction)
        delta = omega / v_0
        l_0 = v_0 / s
        gamma = theta + 1

        if length_gun is None:
            l_bar_g_0 = self.max_length / l_0
        else:
            l_bar_g_0 = length_gun / l_0

        if cc is None:
            cc = 1 - (1 - 1 / chi_k) * log(l_bar_g_0 + 1) / l_bar_g_0

        if any((labda_1 is None, labda_2 is None)):
            if self.sol == SOL_LAGRANGE:
                labda_1, labda_2 = 1 / 2, 1 / 3
            elif self.sol == SOL_PIDDUCK:
                labda_1, labda_2 = pidduck(omega / (phi_1 * m), gamma, tol)
            elif self.sol == SOL_MAMONTOV:
                labda_1, labda_2 = pidduck(omega / (phi_1 * m), 1, tol)
            else:
                raise ValueError("Unknown Solution")

        phi = phi_1 + labda_2 * omega / m * cc
        v_j = (2 * f * omega / (theta * phi * m)) ** 0.5
        v_bar_d = v_d / v_j

        if self.ambient_rho != 0:
            c_a_bar = (self.ambient_gamma * self.ambient_p / self.ambient_rho) ** 0.5 / v_j
            p_a_bar = self.ambient_p / (f * delta)
        else:
            c_a_bar, p_a_bar = 0, 0

        gamma_1 = self.ambient_gamma

        if v_j < v_d and not known_bore:
            raise ValueError(
                "Propellant load too low to achieve design velocity. "
                + " The 2nd ballistic limit for this loading conditions is"
                + " {:.4g} m/s.".format(v_j)
            )

        psi_0 = (1 / delta - 1 / rho_p) / (f / p_0 + alpha - 1 / rho_p)
        z_0, _ = dekker(lambda z: self.propellant.f_psi_z(z) - psi_0, 0, 1, y_rel_tol=tol, y_abs_tol=tol**2)
        logger.info("solved starting burnup.")

        def _f_p_bar(z, l_bar, v_bar):
            psi = f_psi_z(z)
            l_psi_bar = 1 - delta / rho_p - delta * (alpha - 1 / rho_p) * psi
            p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

            if self.control == POINT_PEAK_AVG:
                return p_bar
            else:
                prime = (1 / chi_k + l_bar) / (1 + l_bar)
                labda_1_prime = labda_1 * prime
                labda_2_prime = labda_2 * prime

                factor_s = 1 + labda_2_prime * (omega / (phi_1 * m))
                factor_b = (phi_1 * m + labda_2_prime * omega) / (phi_1 * m + labda_1_prime * omega)

                if self.control == POINT_PEAK_SHOT:
                    return p_bar / factor_s
                elif self.control == POINT_PEAK_BREECH:
                    return p_bar / factor_b
                else:
                    raise ValueError("tag unhandled.")

        """
        step 1, find grain size that satisfies design pressure
        """
        p_bar_s = _f_p_bar(z_0, 0, 0)
        p_bar_d = p_d / (f * delta)  # convert to unitless

        if p_bar_d < p_bar_s:
            raise ValueError(
                "Interior ballistics of conventional gun precludes" + " a pressure lower than that at shot start."
            )
        l_bar_d = self.max_length / l_0

        def abort_z(x, ys, record):
            z = x
            t_bar, l_bar, v_bar = ys
            p_bar = _f_p_bar(z, l_bar, v_bar)

            o_x, o_ys = record[-1]

            o_z = o_x
            ot_bar, ol_bar, ov_bar = o_ys
            op_bar = _f_p_bar(o_z, ol_bar, ov_bar)

            return (p_bar < op_bar) or (p_bar > 2 * p_bar_d)

        def _f_p_e_1(e_1: float):
            """
            calculate either the peak pressure, given the arc thickness,
            or until the system develops 2x design pressure.
            """

            b = s**2 * e_1**2 / (f * phi * omega * m * u_1**2) * (f * delta) ** (2 * (1 - n))

            def _ode_z(z, _, l_bar, v_bar, __):
                """burnup domain ode of internal ballistics"""
                psi = f_psi_z(z)
                l_psi_bar = 1 - delta / rho_p - delta * (alpha - 1 / rho_p) * psi

                p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)
                if c_a_bar and v_bar > 0:
                    v_r = v_bar / c_a_bar
                    p_d_bar = (
                        +0.25 * gamma_1 * (gamma_1 + 1) * v_r**2
                        + gamma_1 * v_r * (1 + (0.25 * (gamma_1 + 1)) ** 2 * v_r**2) ** 0.5
                    ) * p_a_bar

                else:
                    p_d_bar = 0

                if z <= z_b:
                    dt_bar = (2 * b / theta) ** 0.5 * p_bar**-n  # dt_bar/dz
                    dl_bar = v_bar * dt_bar
                    dv_bar = 0.5 * theta * (p_bar - p_d_bar) * dt_bar

                else:
                    dt_bar, dl_bar, dv_bar = 0, 0, 0

                return [dt_bar, dl_bar, dv_bar]

            record = [[z_0, [0, 0, 0]]]

            z_j, (t_bar_j, l_bar_j, v_bar_j), e = RKF78(
                dFunc=_ode_z,
                iniVal=(0, 0, 0),
                x_0=z_0,
                x_1=z_b,
                relTol=tol,
                absTol=tol**2,
                abortFunc=abort_z,
                record=record,
            )

            p_bar_j = _f_p_bar(z_j, l_bar_j, v_bar_j)

            if p_bar_j >= 2 * p_bar_d:  # case for abort due to excessive pressure
                return p_bar_j - p_bar_d, z_j, t_bar_j, l_bar_j, v_bar_j

            # case for abort due to decreasing pressure

            def _f_p_z(z):
                i = record.index([v for v in record if v[0] <= z][-1])
                x = record[i][0]
                ys = record[i][1]

                r = []
                _, (t_bar, l_bar, v_bar), _ = RKF78(
                    dFunc=_ode_z, iniVal=ys, x_0=x, x_1=z, relTol=tol, absTol=tol**2, record=r
                )
                xs = [v[0] for v in record]
                record.extend(v for v in r if v[0] not in xs)
                record.sort()
                return _f_p_bar(z, l_bar, v_bar), z, t_bar, l_bar, v_bar

            """
            find the peak pressure point.
            """

            if len(record) > 1:
                z_i = record[-2][0]
            else:
                z_i = z_0

            z_1, z_2 = gss(lambda z: _f_p_z(z)[0], z_i, z_j, y_rel_tol=tol, y_abs_tol=tol**2, findMin=False)
            z_p = 0.5 * (z_1 + z_2)

            p_bar_p, *vals = _f_p_z(z_p)

            return p_bar_p - p_bar_d, *vals

        dp_bar_probe = _f_p_e_1(self.min_web)[0]
        probe_web = self.min_web

        if dp_bar_probe < 0:
            raise ValueError("Design pressure cannot be achieved by varying web down to minimum")

        while dp_bar_probe > 0:
            probe_web *= 2
            dp_bar_probe = _f_p_e_1(probe_web)[0]

        def fr(x):
            progress_queue.put(round(x * 100))

        e_1, _ = dekker(
            lambda web: _f_p_e_1(web)[0],
            probe_web,  # >0
            0.5 * probe_web,  # ?0
            y_rel_tol=tol,
            y_abs_tol=p_bar_d * tol,
            f_report=fr if progress_queue is not None else None,
        )  # this is the e_1 that satisfies the pressure specification.

        p_bar_dev, z_i, t_bar_i, l_bar_i, v_bar_i = _f_p_e_1(e_1)

        logger.info("solved web satisfying pressure constraint.")

        if known_bore:
            if progress_queue:
                progress_queue.put(100)
            return e_1, length_gun

        if v_j * v_bar_i > v_d:
            if suppress:
                logger.warning("velocity target point occured before peak pressure point.")
                logger.warning("this is currently being suppressed due to program control.")
            else:
                raise ValueError(f"Design velocity exceeded before peak pressure point (V = {v_bar_i * v_j:.4g} m/s).")

        """
        step 2, find the requisite muzzle length to achieve design velocity
        """

        b = s**2 * e_1**2 / (f * phi * omega * m * u_1**2) * (f * delta) ** (2 * (1 - n))

        def _ode_v(v_bar, _, z, l_bar, __):
            psi = f_psi_z(z)

            l_psi_bar = 1 - delta / rho_p - delta * (alpha - 1 / rho_p) * psi
            p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

            if c_a_bar != 0 and v_bar > 0:
                v_r = v_bar / c_a_bar
                p_d_bar = (
                    0.25 * gamma_1 * (gamma_1 + 1) * v_r**2
                    + gamma_1 * v_r * (1 + (0.25 * (gamma_1 + 1)) ** 2 * v_r**2) ** 0.5
                ) * p_a_bar

            else:
                p_d_bar = 0

            dt_bar = 2 / (theta * (p_bar - p_d_bar))
            dz = dt_bar * (0.5 * theta / b) ** 0.5 * p_bar**n
            dl_bar = v_bar * dt_bar

            return [dt_bar, dz, dl_bar]

        def abort_v(x, ys, record):
            _, _, l_bar = ys
            return l_bar > l_bar_d

        vtzl_record = [[v_bar_i, (t_bar_i, z_i, l_bar_i)]]
        try:
            (v_bar_g, (t_bar_g, z_g, l_bar_g), _) = RKF78(
                dFunc=_ode_v,
                iniVal=(t_bar_i, z_i, l_bar_i),
                x_0=v_bar_i,
                x_1=v_bar_d,
                relTol=tol,
                absTol=tol**2,
                abortFunc=abort_v,
                record=vtzl_record,
            )

        except ValueError:
            v_bar_m, (t_bar_m, z_m, l_bar_m) = vtzl_record[-1]
            pmax = _f_p_bar(z_m, l_bar_m, v_bar_m) * f * delta
            vmax = v_bar_m * v_j
            lmax = l_bar_m * l_0
            raise ValueError(
                "Integration appears to be approaching asymptote, "
                + "last calculated to v = {:.4g} m/s, ".format(vmax)
                + "x = {:.4g} m, p = {:.4g} MPa. ".format(lmax, pmax * 1e-6)
                + "This indicates an excessive velocity target relative to pressure developed."
            )

        p_bar_g = _f_p_bar(z_g, l_bar_g, v_bar_g)

        v_g = v_bar_g * v_j
        l_g = l_bar_g * l_0
        p_g = p_bar_g * f * delta
        if l_bar_g > l_bar_d:
            raise ValueError(
                "Solution requires excessive tube length, last calculated to "
                + "v = {:.4g} m/s, x = {:.4g} m, ".format(v_g, l_g)
                + "p = {:.4g} MPa.".format(p_g * 1e-6)
            )

        logger.info("solved tube length to satisfy design velocity.")

        # calculate the averaged chambrage correction factor
        # implied by this solution
        cc_n = 1 - (1 - 1 / chi_k) * log(l_bar_g + 1) / l_bar_g

        if progress_queue is not None:
            progress_queue.put(100)

        if abs((l_bar_g - l_bar_g_0) / min(l_bar_g, l_bar_g_0)) > tol and it < MAX_ITER:
            # successive better approximations will eventually
            # result in value within tolerance.

            logger.info("recursively approaching stable solution under chambrage specification.")
            return self.solve(
                load_fraction=load_fraction,
                charge_mass_ratio=charge_mass_ratio,
                labda_1=labda_1,
                labda_2=labda_2,
                cc=cc_n,
                it=it + 1,
                length_gun=l_bar_g * l_0,
                known_bore=known_bore,
                suppress=suppress,
                progress_queue=progress_queue,
            )
            # TODO: Maximum recursion depth exceeded in comparison is occasionally thrown here. Investigate why.
        else:
            logger.info(f"solution satisfying chambrage found after {it} iterations.")
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
        m = self.m
        omega = m * charge_mass_ratio
        rho_p = self.rho_p
        s = self.s
        phi_1 = self.phi_1
        theta = self.theta
        gamma = theta + 1
        tol = self.tol

        if self.sol == SOL_LAGRANGE:
            labda_1, labda_2 = 1 / 2, 1 / 3
        elif self.sol == SOL_PIDDUCK:
            labda_1, labda_2 = pidduck(omega / (phi_1 * m), gamma, tol)
        elif self.sol == SOL_MAMONTOV:
            labda_1, labda_2 = pidduck(omega / (phi_1 * m), 1, tol)
        else:
            raise ValueError("Unknown Solution")

        logger.info(f"solved parameters for gas distribution {self.sol}")

        def f(lf: float) -> tuple[float, float, float]:
            logger.info(f"dispatching calculation at load factor = {lf:.3%}")
            v_0 = omega / (rho_p * lf)
            l_0 = v_0 / s

            e_1, l_g = self.solve(
                load_fraction=lf,
                charge_mass_ratio=charge_mass_ratio,
                labda_1=labda_1,
                labda_2=labda_2,
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
                "Unable to find any valid load fraction" + " with {:d} random samples.".format(MAX_GUESSES)
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

        It was found that at this step, setting the accuracy metric
        on the x-value (or the load fraction) gives more consistent
        result than requriing a relative tolerance on the function
        values.
        """

        logger.info(f"finding minimum volume solution between {low:.3%} and {high:.3%}")

        def fr(x):
            progress_queue.put(round(x * 33) + 66)

        lf_low, lf_high = gss(
            lambda lf: f(lf)[1],
            low,
            high,
            x_tol=tol,
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
    pass
