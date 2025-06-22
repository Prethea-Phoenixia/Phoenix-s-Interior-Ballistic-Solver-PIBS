from math import pi, log, floor
from random import uniform
from .num import cubic, RKF78, dekker, gss
from . import POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BLEED
from .optGun import MAX_GUESSES
from time import time

import logging
import sys
import traceback


class LPCPMaxBelowStartError(ValueError):
    pass


class ConstrainedHighlow:
    def __init__(
        self,
        caliber,
        propellant,
        shotMass,
        burstPressure,
        startPressure,
        dragCoefficient,
        chambrage,
        tol,
        designHighPressure,
        designLowPressure,
        designVelocity,
        minWeb=1e-6,
        maxLength=1e3,
        maxEV=1,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        control=POINT_PEAK_AVG,
        **_,
    ):

        logger = logging.getLogger("opt.highlow")
        logger.info("initializing constrained highlow object.")

        # cache the constants:
        self.S = (0.5 * caliber) ** 2 * pi
        self.propellant = propellant
        self.m = shotMass
        self.p_0_e = burstPressure
        self.p_0_s = startPressure

        self.p_d_h = designHighPressure
        self.p_d_l = designLowPressure
        self.v_d = designVelocity

        self.chi_k = chambrage

        self.phi_1 = 1 / (1 - dragCoefficient)  # drag work coefficient

        self.tol = tol

        self.maxLength = maxLength
        self.maxEV = maxEV

        self.minWeb = minWeb
        self.ambientRho = ambientRho
        self.ambientP = ambientP
        self.ambientGamma = ambientGamma
        self.control = control

        theta = self.theta
        gamma = theta + 1
        self.phi_2 = 0.15  # for small ports 1.5mm to 2mm in size

        self.cfpr = (2 / (gamma + 1)) ** (gamma / theta)
        self.K_0 = gamma**0.5 * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * theta))

    def __getattr__(self, attrName):
        if "propellant" in vars(self) and not (
            attrName.startswith("__") and attrName.endswith("__")
        ):
            try:
                return getattr(self.propellant, attrName)
            except AttributeError:
                raise AttributeError(
                    "%r object has no attribute %r"
                    % (self.__class__.__name__, attrName)
                )
        else:
            raise AttributeError

    def solve(
        self,
        loadFraction,
        chargeMassRatio,
        portArea,
        lengthGun=None,
        knownBore=False,
        suppress=False,
        progressQueue=None,
        **_,
    ):
        logger = logging.getLogger("opt.highlow.solve")
        logger.info("solving constraint at specified load fraction.")
        if any((chargeMassRatio <= 0, loadFraction <= 0, loadFraction > 1)):
            raise ValueError("Invalid parameters to solve constrained design problem")

        if progressQueue is not None:
            progressQueue.put(1)

        p_max = 1e9

        m = self.m
        rho_p = self.rho_p
        theta = self.theta
        f = self.f
        chi = self.chi
        # chi_k = self.chi_k
        labda = self.labda
        mu = self.mu
        S = self.S
        maxLF = self.maxLF
        phi_1 = self.phi_1
        p_d_l = self.p_d_l
        p_d_h = self.p_d_h
        u_1, n = self.u_1, self.n
        alpha = self.alpha
        Z_b = self.Z_b
        f_psi_Z = self.f_psi_Z
        f_sigma_Z = self.f_sigma_Z
        p_0_e = self.p_0_e
        p_0_s = self.p_0_s
        cfpr = self.cfpr
        phi_2 = self.phi_2
        K_0 = self.K_0
        Z_b = self.Z_b
        tol = self.tol
        p_a = self.ambientP
        k_1 = self.ambientGamma
        control = self.control
        l_d = self.maxLength

        if self.ambientRho != 0:
            c_a = (self.ambientGamma * self.ambientP / self.ambientRho) ** 0.5
        else:
            c_a = 0

        if loadFraction > maxLF:
            raise ValueError("Specified Load Fraction Violates Geometrical Constraint")

        omega = m * chargeMassRatio
        V_0 = omega / (rho_p * loadFraction)

        labda_1, labda_2 = 1 / 2, 1 / 3

        psi_0 = (
            p_0_e
            * (V_0 - omega / rho_p)
            / (f * omega - p_0_e * (omega / rho_p - alpha * omega))
        )

        if psi_0 > 1 or psi_0 < 0:
            raise ValueError(
                "Initial burnup solved to be greater than unity, or negative"
            )

        Zs = cubic(chi * mu, chi * labda, chi, -psi_0)
        # pick a valid solution between 0 and 1
        Zs = sorted(
            Z for Z in Zs if not isinstance(Z, complex) and (Z > 0.0 and Z < 1.0)
        )  # evaluated from left to right, guards against complex >/< float
        if len(Zs) < 1:
            raise ValueError(
                "Propellant either could not develop enough pressure to overcome"
                + " port open pressure, or has burnt to post fracture."
            )

        logger.info("solved shot start burnup for current constraint.")

        def _f_p_1(Z, eta, tau_1, psi=None):
            psi = psi if psi else f_psi_Z(Z)
            V_psi = V_0 - omega / rho_p * (1 - psi) - alpha * omega * (psi - eta)
            return f * omega * tau_1 / V_psi * (psi - eta)

        S_j = portArea
        phi = phi_1 + labda_2 * omega / m

        logger_e_1 = logging.getLogger("opt.highlow.solve.highchamber")

        def f_e_1(e_1):
            logger_e_1.info(
                f"independent high chamber condition calculated at web of {e_1*1e3:.3f} mm"
            )
            Z_0, t_0, eta_0, tau_1_0 = Zs[0], 0, 0, 1

            def _ode_Zi(Z, t, eta, tau_1, _):
                psi = f_psi_Z(Z)
                p_1 = _f_p_1(Z, eta, tau_1, psi)
                dt = 1 / (u_1 / e_1 * p_1**n)  # dt / dZ
                deta = (phi_2 * K_0 * p_1 * S_j / ((f * tau_1) ** 0.5 * omega)) * dt
                dpsi = f_sigma_Z(Z)  # dpsi / dZ
                dtau_1 = ((1 - tau_1) * dpsi - theta * tau_1 * deta) / (psi - eta)

                return dt, deta, dtau_1

            record_0 = [[Z_0, [t_0, eta_0, tau_1_0]]]

            def abort_Zi(x, ys, record):
                Z, _, eta, tau_1 = x, *ys
                p_1 = _f_p_1(Z, eta, tau_1)
                o_x, o_ys = record[-1]
                oZ, _, oeta, otau_1 = o_x, *o_ys
                op_1 = _f_p_1(oZ, oeta, otau_1)

                return p_1 < op_1 or p_1 > p_max

            def g(Z):
                i = record_0.index([v for v in record_0 if v[0] <= Z][-1])
                x = record_0[i][0]
                ys = record_0[i][1]
                r = [(x, ys)]
                Z, (t, eta, tau_1), _ = RKF78(
                    _ode_Zi,
                    ys,
                    x,
                    Z,
                    relTol=tol,
                    absTol=tol**2,
                    abortFunc=abort_Zi,
                    record=r,
                )
                p_1 = _f_p_1(Z, eta, tau_1)
                if p_1 > p_max:
                    raise ValueError
                if Z not in [line[0] for line in r]:
                    r.append([Z, (t, eta, tau_1)])

                xs = [v[0] for v in record_0]
                record_0.extend(v for v in r if v[0] not in xs)
                record_0.sort()
                return p_1

            try:
                Z_p_i = sum(gss(g, Z_0, Z_b, y_rel_tol=0.5 * tol, findMin=False)) * 0.5
            except ValueError:
                return p_max

            p = g(Z_p_i)
            logger_e_1.info(
                f"peak pressure developed is greater or equal to {p*1e-6:.3f} MPa."
            )
            return p

        logger.info(
            f"probing web to satisfy peak high chamber pressure from {self.minWeb * 1e3:.3f} mm"
        )
        p_probe = f_e_1(self.minWeb)
        probeWeb = self.minWeb
        if p_probe < p_d_h:
            raise ValueError(
                "design pressure cannot be achieved by varying web down to minimum"
            )
        else:
            logger.info(
                "high chamber pressure can be achieved with web greater than minimum, "
                + "proceeding to probe."
            )

        while p_probe > p_d_h:
            probeWeb *= 2
            p_probe = f_e_1(probeWeb)

        def fr(x):
            progressQueue.put(round(x * 33))

        logger.info("probing loop found valid range.")
        e_1, _ = dekker(
            f_e_1,
            probeWeb,  # >0
            0.5 * probeWeb,  # ?0
            y=p_d_h,
            y_abs_tol=tol * p_d_h,
            f_report=fr if progressQueue is not None else None,
        )  # this is the e_1 that satisifies the pressure specification.

        logger.info("solved web satisfying pressure constraint.")
        ################################################################################
        ################################################################################
        f_ev_records = dict()
        logger_ev = logging.getLogger("opt.highlow.solve.lowchamber")

        def f_ev(expansionVolume):
            if expansionVolume in f_ev_records:
                return f_ev_records[expansionVolume]
            logger_ev.info(
                "coupled low chamber condition calculated at expansion volume of"
                + f" {expansionVolume*1e3:.6f} L."
            )

            V_1 = expansionVolume

            def _f_p_2(l, eta, tau_2, point=control):
                l_star = (V_1 - alpha * omega * eta) / S
                p_avg = f * omega * tau_2 * eta / (S * (l_star + l))

                if point == POINT_PEAK_AVG:
                    return p_avg
                elif point == POINT_PEAK_SHOT:
                    factor_s = 1 + labda_2 * (omega * eta / (phi_1 * m))
                    return p_avg / factor_s
                elif point == POINT_PEAK_BLEED:
                    factor_b = (phi_1 * m + labda_2 * omega * eta) / (
                        phi_1 * m + labda_1 * omega * eta
                    )
                    return p_avg / factor_b
                else:
                    raise ValueError("Unknown control point.")

            Z_0, t_0, eta_0, tau_1_0, tau_2_0 = Zs[0], 0, 0, 1, 1 + theta

            def _ode_Zs(Z, t, eta, tau_1, tau_2, delta):
                psi = f_psi_Z(Z)
                p_1 = _f_p_1(Z, eta, tau_1, psi)
                p_2 = _f_p_2(0, eta, tau_2, POINT_PEAK_AVG)

                dt = 1 / (u_1 / e_1 * p_1**n)  # dt / dZ

                pr = p_2 / p_1
                if pr <= cfpr:
                    deta = (phi_2 * K_0 * p_1 * S_j / ((f * tau_1) ** 0.5 * omega)) * dt
                else:
                    gamma = theta + 1
                    deta = (
                        (phi_2 * p_1 * S_j)
                        / ((f * tau_1) ** 0.5 * omega)
                        * (
                            (2 * gamma / theta)
                            * (pr ** (2 / gamma) - pr ** ((gamma + 1) / (gamma)))
                        )
                        ** 0.5
                    ) * dt

                dpsi = f_sigma_Z(Z)  # dpsi / dZ
                dtau_1 = ((1 - tau_1) * dpsi - theta * tau_1 * deta) / (psi - eta)

                # internal half-stepping
                Z += delta
                t += dt * delta
                eta += deta * delta
                tau_1 += dtau_1 * delta

                dtau_2 = (((1 + theta) * tau_1 - tau_2) * deta) / eta

                return dt, deta, dtau_1, dtau_2

            record_Zs = [[Z_0, [t_0, eta_0, tau_1_0, tau_2_0]]]

            def abort_Zs(x, ys, record):
                _, _, eta, _, tau_2 = x, *ys
                p_2 = _f_p_2(0, eta, tau_2, POINT_PEAK_AVG)

                o_x, o_ys = record[-1]
                _, _, oeta, _, otau_2 = o_x, *o_ys
                op_2 = _f_p_2(0, oeta, otau_2, POINT_PEAK_AVG)

                return p_2 > p_0_s or (p_2 - op_2) < p_2 * tol

            Z_k, (t, eta, tau_1, tau_2), _ = RKF78(
                _ode_Zs,
                record_Zs[0][1],
                record_Zs[0][0],
                Z_b,
                relTol=tol,
                absTol=tol**2,
                abortFunc=abort_Zs,
                record=record_Zs,
            )

            def g(Z):
                i = record_Zs.index([v for v in record_Zs if v[0] <= Z][-1])
                x, ys = record_Zs[i][0], record_Zs[i][1]
                r = []
                Z, (t, eta, tau_1, tau_2), _ = RKF78(
                    _ode_Zs, ys, x, Z, relTol=tol, absTol=tol**2, record=r
                )
                xs = [v[0] for v in record_Zs]
                record_Zs.extend(v for v in r if v[0] not in xs)
                record_Zs.sort()
                p_2 = _f_p_2(0, eta, tau_2, POINT_PEAK_AVG)
                return p_2, (t, eta, tau_1, tau_2)

            Z_j = record_Zs[-1][0]
            Z_p = (
                sum(
                    gss(
                        lambda x: g(x)[0],
                        Z_j,
                        Z_k,
                        y_abs_tol=p_0_s * tol,
                        findMin=False,
                    )
                )
                * 0.5
            )
            p_p_avg, (t, eta, tau_1, tau_2) = g(Z_p)
            # p_p = _f_p_2(0, eta, tau_2)
            if p_p_avg < p_0_s:
                logger_ev.warning(
                    f"maximum low chamber pressure not enough to start shot, "
                    + f"{p_p_avg * 1e-6:.3f} MPa",
                )
                raise LPCPMaxBelowStartError(
                    f"Maximum pressure developed in low-chamber ({p_p_avg * 1e-6:.3f} MPa) "
                    + "not enough to start the shot."
                )

            Z_1, _ = dekker(lambda x: g(x)[0] - p_0_s, Z_j, Z_p, y_abs_tol=p_0_s * tol)

            _, (t_1, eta_1, tau_1_1, tau_2_1) = g(Z_1)

            def abort_t(x, ys, record):
                _, _, l, _, eta, _, tau_2 = x, *ys
                p_l = _f_p_2(l, eta, tau_2)

                o_x, o_ys = record[-1]
                _, _, ol, _, oeta, _, otau_2 = o_x, *o_ys
                op_l = _f_p_2(ol, oeta, otau_2)

                return p_l < op_l

            def _ode_t(t, Z, l, v, eta, tau_1, tau_2, _):
                psi = f_psi_Z(Z)
                p_1 = _f_p_1(Z, eta, tau_1, psi)
                p_2 = _f_p_2(l, eta, tau_2, POINT_PEAK_AVG)

                dZ = u_1 / e_1 * p_1**n if Z <= Z_b else 0

                pr = p_2 / p_1
                if pr <= cfpr:
                    deta = (phi_2 * K_0 * p_1 * S_j) / ((f * tau_1) ** 0.5 * omega)
                else:
                    gamma = theta + 1
                    deta = (
                        (phi_2 * p_1 * S_j)
                        / ((f * tau_1) ** 0.5 * omega)
                        * (
                            (2 * gamma / theta)
                            * (pr ** (2 / gamma) - pr ** ((gamma + 1) / (gamma)))
                        )
                        ** 0.5
                    )

                dpsi = f_sigma_Z(Z) * dZ
                dtau_1 = ((1 - tau_1) * dpsi - theta * tau_1 * deta) / (psi - eta)

                if c_a != 0 and v > 0:
                    k = k_1  # gamma
                    v_r = v / c_a
                    p_d = (
                        0.25 * k * (k + 1) * v_r**2
                        + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
                    ) * p_a
                else:
                    p_d = 0

                dv = S / (phi * m) * (p_2 - p_d)
                dl = v

                dtau_2 = (
                    ((1 + theta) * tau_1 - tau_2) * deta
                    - (theta * phi * m) / (f * omega) * v * dv
                ) / eta

                return dZ, dl, dv, deta, dtau_1, dtau_2

            record_t = [[t_1, (Z_1, 0, 0, eta_1, tau_1_1, tau_2_1)]]

            t_probe, t_j = t_1, t_1
            while t_probe == t_j:  # has not aborted
                logger_ev.info("extended time of integration")
                t_probe = t_j * 2
                t_j, (Z_j, l_j, v_j, eta_j, tau_1_j, tau_2_j), _ = RKF78(
                    _ode_t,
                    record_t[-1][1],
                    record_t[-1][0],
                    t_probe,
                    relTol=tol,
                    absTol=tol**2,
                    abortFunc=abort_t,
                    record=record_t,
                )

            def _f_p_2_t(t):
                i = record_t.index([v for v in record_t if v[0] <= t][-1])
                x, ys = record_t[i]
                r = []
                t, (Z, l, v, eta, tau_1, tau_2), _ = RKF78(
                    dFunc=_ode_t,
                    iniVal=ys,
                    x_0=x,
                    x_1=t,
                    relTol=tol,
                    absTol=tol**2,
                    record=r,
                )
                xs = [v[0] for v in record_t]
                record_t.extend(v for v in r if v[0] not in xs)
                record_t.sort()

                return _f_p_2(l, eta, tau_2), (t, Z, l, v, eta, tau_1, tau_2)

            pressures = []
            for line in record_t:
                t, (Z, l, v, eta, tau_1, tau_2) = line
                p_2 = _f_p_2(l, eta, tau_2)
                pressures.append(p_2)

            t_i, t_k = None, None

            if len(pressures) > 3:  # barely 30ms improvement, *sigh*.
                for i in range(len(pressures) - 2):
                    p_i, p_j, p_k = pressures[i], pressures[i + 1], pressures[i + 2]
                    if p_i <= p_j and p_j >= p_k:
                        t_i, t_k = record_t[i][0], record_t[i + 2][0]

            if t_i is None and t_k is None:
                t_i, t_k = t_1, t_j

            t_p_2 = (
                sum(
                    gss(
                        lambda t: _f_p_2_t(t)[0],
                        t_i,
                        t_k,
                        y_rel_tol=0.5 * tol,
                        findMin=False,
                    )
                )
                * 0.5
            )

            p_p_l_t, vals = _f_p_2_t(t_p_2)

            f_ev_records.update({expansionVolume: (p_p_l_t, vals)})
            logger_ev.info(
                f"peak low chamber pressure calculated is {p_p_l_t*1e-6:.3f}MPa."
            )
            return p_p_l_t, vals

        def findBound(func, x_probe, x_bound, tol, record=None, exception=ValueError):
            if record is None:
                record = []
            try:
                record.append((x_bound, func(x_bound)))
                return x_bound
            except exception:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                logger.warning("exception raised in bound finding outer bound")
                logger.warning(
                    "".join(
                        traceback.format_exception(exc_type, exc_value, exc_traceback)
                    )
                )
                x_valid = x_probe

                delta = x_bound - x_probe
                up = x_bound > x_probe

                while abs(2 * delta) > tol and (
                    x_probe <= x_bound + tol if up else x_probe >= x_bound - tol
                ):
                    try:
                        record.append((x_probe, func(x_probe)))
                        x_valid = x_probe
                    except exception:
                        exc_type, exc_value, exc_traceback = sys.exc_info()
                        logger.warning("exception raised in bound finding loop")
                        logger.warning(
                            "".join(
                                traceback.format_exception(
                                    exc_type, exc_value, exc_traceback
                                )
                            )
                        )
                        delta *= 0.5
                    finally:
                        x_probe = x_valid + delta

                return x_valid

        probeEV = 2 * tol * V_0
        lastEV = tol * V_0

        logger.info("probing for expansion volume lower bound.")
        while True:
            try:
                p = f_ev(probeEV)[0]
                break

            except LPCPMaxBelowStartError:  # overshoot
                logger.warning(
                    "exception raised: Low Pressure Chamber Pressure Max Below Start Error",
                )
                delta = probeEV - lastEV
                if delta < tol * V_0:
                    raise ValueError(
                        "No valid expansion volume can be found that "
                        + "starts the shot without the pressures coalescing"
                        + " between high and low chamber."
                    )
                probeEV = lastEV + 0.5 * delta
            except ValueError:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                logger.warning("exception raised: ValueError")
                logger.warning(
                    "".join(
                        traceback.format_exception(exc_type, exc_value, exc_traceback)
                    )
                )
                lastEV = probeEV
                probeEV *= 2

        logger.info(
            f"finding lower bound for f_ev between {probeEV * 1e3:.6f} L and {lastEV * 1e3:.6f} L",
        )
        evLowerBound = findBound(f_ev, probeEV, lastEV, tol * V_0)
        logger.info(f"expansion volume lower bound found as {evLowerBound * 1e3:.6f} L")

        if progressQueue is not None:
            progressQueue.put(50)

        probeEV = evLowerBound * 2
        lastEV = evLowerBound

        logger.info("probing for expansion volume upper bound.")
        while True:
            try:
                f_ev(probeEV)[0]
                lastEV = probeEV
                probeEV *= 2
            except ValueError:
                break

        logger.info(
            f"finding upper bound for f_ev between {probeEV * 1e3:.3f} L and {lastEV * 1e3:.3f} L",
        )
        evUpperBound = findBound(f_ev, lastEV, probeEV, tol * V_0)
        logger.info(f"expansion volume upper bound found as {evUpperBound * 1e3:.3f} L")

        if progressQueue is not None:
            progressQueue.put(66)

        p_l_max = f_ev(evLowerBound)[0]
        p_l_min = f_ev(evUpperBound)[0]

        if p_l_min > p_d_l or p_d_l > p_l_max:
            raise ValueError(
                "Range of valid solution does not allow low chamber design pressure "
                + "to be met, "
                + f"between V > {evLowerBound * 1e3:.3f} L, P < {p_l_max * 1e-6:.3f} MPa, "
                + f"and V < {evUpperBound * 1e3:.3f} L, P > {p_l_min * 1e-6:.3f} MPa."
            )

        def fr(x):
            progressQueue.put(round(x * 34) + 66)

        logger.info(
            f"dekker method called on f_ev between {evLowerBound * 1e3:.3f} L and {evUpperBound*1e3:.3f} L",
        )
        # ascending volume, descending pressure
        f_ev_evs, f_ev_results = zip(*sorted(f_ev_records.items(), key=lambda v: v[0]))

        index_maxVol = [p > p_d_l for p, _ in f_ev_results].index(False)
        index_minVol = index_maxVol - 1

        V_1, _ = dekker(
            lambda ev: f_ev(ev)[0],
            f_ev_evs[index_minVol],
            f_ev_evs[index_maxVol],
            y=p_d_l,
            y_rel_tol=tol,
            f_report=fr if progressQueue is not None else None,
        )
        logger.info(
            f"expansion volume to match low chamber pressure solved at {V_1 * 1e3:.3f} L"
        )
        p_l_act, (t_i, Z_i, l_i, v_i, eta_i, tau_1_i, tau_2_i) = f_ev(V_1)
        logger.info(f"check: actual low chamber pressure {p_l_act * 1e-6:.3f} MPa")

        if knownBore:
            if progressQueue is not None:
                progressQueue.put(100)
            return e_1, V_1, lengthGun

        def _f_p_2(l, eta, tau_2, point=control):
            l_star = (V_1 - alpha * omega * eta) / S
            p_avg = f * omega * tau_2 * eta / (S * (l_star + l))

            if point == POINT_PEAK_AVG:
                return p_avg

            elif point == POINT_PEAK_SHOT:
                factor_s = 1 + labda_2 * (omega * eta / (phi_1 * m))
                return p_avg / factor_s

            elif point == POINT_PEAK_BLEED:
                factor_b = (phi_1 * m + labda_2 * omega * eta) / (
                    phi_1 * m + labda_1 * omega * eta
                )
                return p_avg / factor_b
            else:
                raise ValueError("Unknown control point.")

        def _ode_v(v, t, Z, l, eta, tau_1, tau_2, _):
            psi = f_psi_Z(Z)
            p_1 = _f_p_1(Z, eta, tau_1, psi)
            p_2 = _f_p_2(l, eta, tau_2, POINT_PEAK_AVG)

            if c_a != 0 and v > 0:
                k = k_1  # gamma
                v_r = v / c_a
                p_d = (
                    0.25 * k * (k + 1) * v_r**2
                    + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
                ) * p_a
            else:
                p_d = 0

            dt = (phi * m) / (S * (p_2 - p_d))  # dt / dv
            dl = v * dt  # dl / dv

            dZ = u_1 / e_1 * p_1**n * dt if Z <= Z_b else 0

            pr = p_2 / p_1
            if pr <= cfpr:
                deta = (phi_2 * K_0 * p_1 * S_j) / ((f * tau_1) ** 0.5 * omega) * dt
            else:
                gamma = theta + 1
                deta = (
                    (phi_2 * p_1 * S_j)
                    / ((f * tau_1) ** 0.5 * omega)
                    * (
                        (2 * gamma / theta)
                        * (pr ** (2 / gamma) - pr ** ((gamma + 1) / (gamma)))
                    )
                    ** 0.5
                ) * dt

            dpsi = f_sigma_Z(Z) * dZ  # dpsi / dv
            dtau_1 = ((1 - tau_1) * dpsi - theta * tau_1 * deta) / (
                psi - eta
            )  # dtau_1 / dv

            dtau_2 = (
                ((1 + theta) * tau_1 - tau_2) * deta
                - (theta * phi * m) / (f * omega) * v
            ) / eta  # dtau_2 / dv

            return dt, dZ, dl, deta, dtau_1, dtau_2

        v_d = self.v_d
        if v_i > v_d and not suppress:
            if suppress:
                logger.warning(
                    "velocity target point occured before peak pressure point."
                )
                logger.warning(
                    "this is currently being suppressed due to program control."
                )
            else:
                raise ValueError(
                    f"Design velocity exceeded before peak pressure point (V = {v_i:.4g} m/s)."
                )

        def abort_v(x, ys, record):
            _, _, l, _, _, _ = ys
            return l > l_d

        vtzlett_record = [v_i, (t_i, Z_i, l_i, eta_i, tau_1_i, tau_2_i)]
        try:
            v_g, (t_g, Z_g, l_g, eta_g, tau_1_g, tau_2_g), _ = RKF78(
                _ode_v,
                (t_i, Z_i, l_i, eta_i, tau_1_i, tau_2_i),
                v_i,
                v_d,
                relTol=tol,
                absTol=tol**2,
                record=vtzlett_record,
                abortFunc=abort_v,
            )

            p_1_g, p_2_g = _f_p_1(Z_g, eta_g, tau_1_g), _f_p_2(l_g, eta_g, tau_2_g)

        except ValueError:
            v_m, (t_m, Z_m, l_m, eta_m, tau_1_m, tau_2_m) = vtzlett_record[-1]
            p_1_m, p_2_m = _f_p_1(Z_m, eta_m, tau_1_m), _f_p_2(l_m, eta_m, tau_2_m)

            raise ValueError(
                "Integration appears to be approaching asymptote, "
                + f"last calculated to v = {v_m:.4g} m/s, "
                + f"x = {l_m:.4g} m, p = {p_1_m * 1e-6:.4g} MPa (high), {p_2_m * 1e-6:.4g} MPa (low) . "
                + "This indicates an excessive velocity target relative to pressure developed."
            )

        if l_g > l_d:
            raise ValueError(
                "Solution requires excessive tube length, last calculated to "
                + f"v = {v_g:.4g} m/s, x = {l_g:.4g} m, "
                + f"p = {p_1_g * 1e-6:.4g} MPa (high), {p_2_g * 1e-6:.4g} (low)."
            )

        if abs(v_g - v_d) > (tol * v_d):
            raise ValueError(
                "Velocity target is not met, last calculated to "
                + f"v = {v_g:.4g} m/s ({(v_g - v_d) / v_d:+.3g} %), x = {l_g:.4g} m, "
                + f"p = {p_1_g * 1e-6:.4g} MPa (high), {p_2_g * 1e-6:.4g} (low)."
            )

        if progressQueue is not None:
            progressQueue.put(100)

        logger.info("solution satisfying constraint found.")

        return e_1, V_1, l_g

    def findMinV(self, chargeMassRatio, portArea, progressQueue=None, **_):
        """
        find the minimum volume solution.
        """
        if progressQueue is not None:
            progressQueue.put(1)

        """
        Step 1, find a valid range of values for load fraction,
        using psuedo-bisection.

        high lf -> high web
        low lf -> low web
        """
        m = self.m
        omega = m * chargeMassRatio
        rho_p = self.rho_p
        S = self.S
        solve = self.solve
        tol = self.tol

        def f(lf):
            V_0 = omega / (rho_p * lf)
            l_0 = V_0 / S

            e_1, V_1, l_g = solve(
                loadFraction=lf,
                chargeMassRatio=chargeMassRatio,
                portArea=portArea,
                knownBore=False,
                suppress=True,
            )
            l_1 = V_1 / S

            return l_g + l_0 + l_1, e_1, V_1, l_g

        records = []
        for i in range(MAX_GUESSES):
            startProbe = uniform(tol, 1 - tol)

            try:
                lt_i, *_ = f(startProbe)
                records.append((startProbe, lt_i))
                break
            except ValueError:
                if progressQueue is not None:
                    progressQueue.put(round(i / MAX_GUESSES * 33))

        else:
            raise ValueError(
                "Unable to find any valid load fraction"
                + " with {:d} random samples.".format(MAX_GUESSES)
            )

        if progressQueue is not None:
            progressQueue.put(33)

        # find the lower limit

        low = tol
        cap_low = low
        probe = startProbe
        delta_low = low - probe
        new_low = probe + delta_low

        k, n = -1, floor(log(abs(delta_low) / tol, 2)) + 1
        while abs(2 * delta_low) > tol:
            try:
                lt_i, *_ = f(new_low)
                records.append((new_low, lt_i))
                probe = new_low
            except ValueError:
                cap_low = new_low
                delta_low *= 0.5
                k += 1
            finally:
                if probe + delta_low < cap_low + tol:
                    delta_low *= 0.5
                    k += 1
                new_low = probe + delta_low

                if progressQueue is not None:
                    progressQueue.put(round(k / n * 17) + 33)

        low = probe

        # find the upper limit

        high = 1 - tol
        cap_high = high
        probe = startProbe
        delta_high = high - probe
        new_high = probe + delta_high

        k, n = -1, floor(log(abs(delta_high) / tol, 2)) + 1
        while abs(2 * delta_high) > tol:
            try:
                lt_i, *_ = f(new_high)
                records.append((new_high, lt_i))
                probe = new_high
            except ValueError:
                cap_high = new_high
                delta_high *= 0.5
                k += 1
            finally:
                if probe + delta_high > cap_high - tol:
                    delta_high *= 0.5
                    k += 1

                new_high = probe + delta_high
                if progressQueue is not None:
                    progressQueue.put(round(k / n * 17) + 50)

        high = probe

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

        def fr(x):
            progressQueue.put(round(x * 33) + 67)

        lf_low, lf_high = gss(
            lambda lf: f(lf)[0],
            low,
            high,
            x_tol=tol,
            findMin=True,
            f_report=fr if progressQueue is not None else None,
        )

        lf = 0.5 * (lf_high + lf_low)
        _, e_1, V_1, l_g = f(lf)

        if progressQueue is not None:
            progressQueue.put(100)

        return lf, e_1, V_1, l_g
