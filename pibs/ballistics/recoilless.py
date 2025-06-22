from __future__ import annotations

import logging
import sys
import traceback
from dataclasses import dataclass
from math import inf, pi, tan
from typing import List

from pibs.ballistics import (
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
    Domains,
)
from .gun import (
    GenericEntry,
    GenericResult,
    Gun,
    OutlineEntry,
    PressureProbePoint,
    PressureTraceEntry,
)
from .num import RKF78, cubic, dekker, gss
from .prop import GrainComp, Propellant

logger = logging.getLogger(__name__)


@dataclass
class RecoillessTableEntry(GenericEntry):
    tag: str
    time: float
    travel: float
    burnup: float
    velocity: float
    outflowVelocity: float
    breechPressure: float
    stagPressure: float
    avgPressure: float
    shotPressure: float
    temperature: float
    outflowFraction: float


@dataclass
class RecoillessErrorEntry(GenericEntry):
    tag: str
    time: float | None = None
    travel: float | None = None
    burnup: float | None = None
    velocity: float | None = None
    outflowVelocity: float | None = None
    breechPressure: float | None = None
    stagPressure: float | None = None
    avgPressure: float | None = None
    shotPressure: float | None = None
    temperature: float | None = None
    outflowFraction: float | None = None


@dataclass
class RecoillessResult(GenericResult):
    gun: "Recoilless"
    tableData: List[RecoillessTableEntry]
    errorData: List[RecoillessErrorEntry]


class Recoilless:
    def __init__(
        self,
        caliber,
        shotMass,
        propellant,
        grainSize,
        chargeMass,
        chamberVolume,
        startPressure,
        lengthGun,
        chambrage,
        nozzleExpansion,
        dragCoefficient=0,
        nozzleEfficiency=0.92,
        structuralMaterial=None,
        structuralSafetyFactor=1.1,
        autofrettage=True,
        **_,
    ):

        logger.info("initializing recoilless object.")
        if any(
            (
                caliber <= 0,
                shotMass <= 0,
                chargeMass <= 0,
                grainSize <= 0,
                chamberVolume <= 0,
                lengthGun <= 0,
                nozzleExpansion < 1,
                nozzleEfficiency > 1,
                nozzleEfficiency <= 0,
                dragCoefficient < 0,
                dragCoefficient >= 1,
                startPressure < 0,
                structuralSafetyFactor <= 1,
            )
        ):
            raise ValueError("Invalid gun parameters")

        if chargeMass > (propellant.maxLF * propellant.rho_p * chamberVolume):
            raise ValueError("Specified Load Fraction Violates Geometrical Constraint")

        self.propellant = propellant
        self.caliber = caliber

        e_1 = 0.5 * grainSize
        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.omega = chargeMass
        self.V_0 = chamberVolume
        self.p_0 = startPressure
        self.l_g = lengthGun
        self.chi_0 = nozzleEfficiency
        self.A_bar = nozzleExpansion

        self.chi_k = chambrage  # ration of l_0 / l_chamber
        self.l_0 = self.V_0 / self.S
        self.l_c = self.l_0 / self.chi_k

        self.Delta = self.omega / self.V_0
        self.phi_1 = 1 / (1 - dragCoefficient)  # drag work coefficient
        self.phi = self.phi_1 + self.omega / (3 * self.m)

        self.v_j = (2 * self.f * self.omega / (self.theta * self.phi * self.m)) ** 0.5

        self.material = structuralMaterial
        self.ssf = structuralSafetyFactor
        self.is_af = autofrettage

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
            self.S**2
            * e_1**2
            / (self.f * self.phi * self.omega * self.m * self.u_1**2)
            * (self.f * self.Delta) ** (2 * (1 - self.n))
        )

        Zs = cubic(self.chi * self.mu, self.chi * self.labda, self.chi, -self.psi_0)
        # pick a valid solution between 0 and 1

        Zs = sorted(
            Z for Z in Zs if not isinstance(Z, complex) and (0.0 < Z < 1.0)
        )  # evaluated from left to right, guards against complex >/< float
        if len(Zs) < 1:
            raise ValueError(
                "Propellant either could not develop enough pressure to overcome"
                + " start pressure, or has burnt to post fracture."
            )

        self.Z_0 = Zs[0]  # pick the smallest solution

        # additional calculation for recoilless weapons:
        gamma = self.theta + 1
        """
        Enforce the "recoilless condition" by setting the size of the
        throat.
        """

        phi_2 = 1
        self.C_A = (
            (0.5 * self.theta * self.phi * self.m / self.omega) ** 0.5
            * gamma**0.5
            * (2 / (gamma + 1)) ** (0.5 * (gamma + 1) / self.theta)
            * phi_2
        )  # flow rate value
        self.k_1 = 0.0
        self.c_a_bar = 0.0
        self.p_a_bar = 0.0
        self.S_j = 0.0
        self.S_j_bar = 0.0
        self.C_f = 0.0
        logger.info("recoilless successfully initialized.")

    def __getattr__(self, attrName):
        if "propellant" in vars(self) and not (attrName.startswith("__") and attrName.endswith("__")):
            try:
                return getattr(self.propellant, attrName)
            except AttributeError:
                raise AttributeError("%r object has no attribute %r" % (self.__class__.__name__, attrName))
        else:
            raise AttributeError

    def _f_p_bar(self, Z, l_bar, eta, tau, psi=None):
        psi = psi if psi else self.f_psi_Z(Z)
        l_psi_bar = 1 - self.Delta * ((1 - psi) / self.rho_p + self.alpha * (psi - eta))

        p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

        return p_bar

    def _ode_t(self, _, Z, l_bar, v_bar, eta, tau, __):
        psi = self.f_psi_Z(Z)
        dpsi = self.f_sigma_Z(Z)  # dpsi/dZ
        p_bar = self._f_p_bar(Z, l_bar, eta, tau, psi)

        if self.c_a_bar != 0 and v_bar > 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (
                +0.25 * k * (k + 1) * v_r**2 + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a_bar
        else:
            p_d_bar = 0

        if Z <= self.Z_b:
            dZ = (0.5 * self.theta / self.B) ** 0.5 * p_bar**self.n
        else:
            dZ = 0  # dZ/dt_bar

        dl_bar = v_bar
        dv_bar = self.theta * 0.5 * (p_bar - p_d_bar)

        deta = self.C_A * self.S_j_bar * p_bar * tau**-0.5  # deta / dt_bar
        dtau = ((1 - tau) * (dpsi * dZ) - 2 * v_bar * dv_bar - self.theta * tau * deta) / (psi - eta)  # dtau/dt_bar

        return dZ, dl_bar, dv_bar, deta, dtau

    def _ode_l(self, l_bar, _, Z, v_bar, eta, tau, __):
        """length domain ode of internal ballistics
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition.

        in general, d/dl_bar = d/dt_bar * dt_bar/dl_bar

        """

        psi = self.f_psi_Z(Z)
        dpsi = self.f_sigma_Z(Z)  # dpsi/dZ
        p_bar = self._f_p_bar(Z, l_bar, eta, tau, psi)

        if self.c_a_bar != 0 and v_bar > 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (+0.25 * k * (k + 1) * v_r**2 + k * v_r * (1 + (0.25 * (k + 1) * v_r) ** 2) ** 0.5) * self.p_a_bar
        else:
            p_d_bar = 0

        if Z <= self.Z_b:
            dZ = (0.5 * self.theta / self.B) ** 0.5 * p_bar**self.n / v_bar
        else:
            dZ = 0  # dZ /dl_bar
        dv_bar = self.theta * 0.5 * (p_bar - p_d_bar) / v_bar  # dv_bar/dl_bar
        dt_bar = 1 / v_bar  # dt_bar / dl_bar

        deta = self.C_A * self.S_j_bar * p_bar * tau**-0.5 * dt_bar  # deta / dl_bar

        dtau = ((1 - tau) * (dpsi * dZ) - 2 * v_bar * dv_bar - self.theta * tau * deta) / (psi - eta)  # dtau/dl_bar

        return dt_bar, dZ, dv_bar, deta, dtau

    def _ode_Z(self, Z, _, l_bar, v_bar, eta, tau, __):
        psi = self.f_psi_Z(Z)
        dpsi = self.f_sigma_Z(Z)  # dpsi/dZ
        p_bar = self._f_p_bar(Z, l_bar, eta, tau, psi)

        if self.c_a_bar != 0 and v_bar > 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (
                +0.25 * k * (k + 1) * v_r**2 + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a_bar
        else:
            p_d_bar = 0

        dt_bar = (2 * self.B / self.theta) ** 0.5 * p_bar**-self.n
        dl_bar = v_bar * dt_bar
        dv_bar = 0.5 * self.theta * (p_bar - p_d_bar) * dt_bar

        deta = self.C_A * self.S_j_bar * p_bar * tau**-0.5 * dt_bar  # deta / dZ

        dtau = ((1 - tau) * dpsi - 2 * v_bar * dv_bar - self.theta * tau * deta) / (psi - eta)

        return dt_bar, dl_bar, dv_bar, deta, dtau

    def integrate(
        self,
        step=10,
        tol=1e-5,
        dom: Domains = DOMAIN_TIME,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        progressQueue=None,
        **_,
    ):
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

        if progressQueue is not None:
            progressQueue.put(1)
        logger.info("commencing integration for characteristic points.")

        record = []
        if any((step < 0, tol < 0)):
            raise ValueError("Invalid integration specification")

        if any((ambientP < 0, ambientRho < 0, ambientGamma < 1)):
            raise ValueError("Invalid ambient condition")

        gamma = self.theta + 1
        """
        self.S_j_bar = self._getCf(gamma, 1, tol) / (
            self._getCf(gamma, self.A_bar, tol) * self.chi_0
        )  # = S_j/S
        """
        self.C_f = self._getCf(gamma, self.A_bar, tol)
        self.S_j_bar = 1 / (self.C_f * self.chi_0)
        if self.S_j_bar > self.chi_k:
            raise ValueError(
                "Achieving recoiless condition necessitates"
                + " a larger throat area than could be fit into"
                + " breech face."
            )
        self.S_j = self.S_j_bar * self.S

        tScale = self.l_0 / self.v_j
        pScale = self.f * self.Delta

        # ambient conditions
        self.p_a_bar = ambientP / pScale
        if ambientRho != 0:
            self.c_a_bar = (ambientGamma * ambientP / ambientRho) ** 0.5 / self.v_j
        else:
            self.c_a_bar = 0

        self.k_1 = ambientGamma

        l_g_bar = self.l_g / self.l_0
        p_bar_0 = self.p_0 / pScale
        Z_b = self.Z_b
        Z_0 = self.Z_0

        bar_data = []
        bar_err = []

        def updBarData(tag, t_bar, l_bar, Z, v_bar, eta, tau, t_bar_err, l_bar_err, Z_err, v_bar_err, eta_err, tau_err):
            p_bar = self._f_p_bar(Z, l_bar, eta, tau)
            p_bar_err = "N/A"
            bar_data.append((tag, t_bar, l_bar, Z, v_bar, p_bar, eta, tau))
            """
            Worst case scenario for pressure deviation is when:
            Z higher than actual, l lower than actual, v lower than actual
            more burnt propellant, less volume, lower speed of projectile.
            """
            bar_err.append(("L", t_bar_err, l_bar_err, Z_err, v_bar_err, p_bar_err, eta_err, tau_err))

        updBarData(
            tag=POINT_START,
            t_bar=0,
            l_bar=0,
            Z=Z_0,
            v_bar=0,
            eta=0,
            tau=1,
            t_bar_err=0,
            l_bar_err=0,
            Z_err=0,
            v_bar_err=0,
            eta_err=0,
            tau_err=0,
        )

        record.append((0, (0, self.psi_0, 0, p_bar_0, 0, 1)))

        Z_i = Z_0
        Z_j = Z_b
        N = 1
        Delta_Z = Z_b - Z_0

        t_bar_i, l_bar_i, v_bar_i, p_bar_i, eta_i, tau_i = 0, 0, 0, p_bar_0, 0, 1

        isBurnOutContained = True

        """
        Instead of letting the integrator handle the heavy lifting, we
        partition Z and integrate upwards until either barrel exit or
        burnout has been achieved. This seek-and-validate inspired
        from bisection puts the group of value denoted by subscript
        i within or on the muzzle, with the propellant either still
        burning or right on the burnout point..
        """
        ztlvet_record = [(Z_0, (0, 0, 0, 0, 1))]
        p_max = 1e9  # 1GPa
        p_bar_max = p_max / pScale

        # noinspection PyUnusedLocal
        def abort(x, ys, record):
            Z, (_, l_bar, v_bar, eta, tau) = x, ys
            p_bar = self._f_p_bar(Z, l_bar, eta, tau)

            return any((l_bar > l_g_bar, p_bar > p_bar_max, v_bar <= 0, p_bar < 0))

        while Z_i < Z_b:  # terminates if burnout is achieved
            ztlvet_record_i = []
            if Z_j == Z_i:
                raise ValueError("Numerical accuracy exhausted in search of exit/burnout point.")
            try:
                if Z_j > Z_b:
                    Z_j = Z_b

                Z, (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j), _ = RKF78(
                    self._ode_Z,
                    (t_bar_i, l_bar_i, v_bar_i, eta_i, tau_i),
                    Z_i,
                    Z_j,
                    relTol=tol,
                    absTol=tol**2,
                    abortFunc=abort,
                    record=ztlvet_record_i,
                )

                p_bar_j = self._f_p_bar(Z_j, l_bar_j, eta_j, tau_j)

            except ValueError as e:
                ztlvet_record.extend(ztlvet_record_i)
                Z, (t_bar, l_bar, v_bar, eta, tau) = ztlvet_record[-1]
                dt_bar, dl_bar, dv_bar, deta, dtau = self._ode_Z(Z, t_bar, l_bar, v_bar, eta, tau, _)
                # p_bar = self._f_p_bar(Z, l_bar, eta, tau)

                if all((dt_bar > 0, dl_bar > 0, dv_bar < 0)):
                    raise ValueError(
                        "Extremely low propulsive effort exerted on shot,"
                        + " impossible to integrate down to numerical precision.\n"
                        + "Shot last calculated at {:.0f} mm with velocity {:.0f} mm/s after {:.0f} ms\n".format(
                            l_bar * self.l_0 * 1e3,
                            v_bar * self.v_j * 1e3,
                            t_bar * tScale * 1e3,
                        )
                    )
                else:
                    raise e  # unknown issues

            if l_bar_j >= l_g_bar:
                if abs(l_bar_i - l_g_bar) / (l_g_bar) > tol or l_bar_i == 0:
                    N *= 2
                    Z_j = Z_i + Delta_Z / N
                else:
                    isBurnOutContained = False
                    break  # l_bar_i is solved to within a tol of l_bar_g

            else:
                ztlvet_record.extend(ztlvet_record_i)
                if p_bar_j > p_bar_max:
                    raise ValueError(
                        "Nobel-Abel EoS is generally accurate enough below 600MPa. However,"
                        + " Unreasonably high pressure (>{:.0f} MPa) was encountered.".format(
                            p_max / 1e6
                        )  # in practice most of the pressure-realted spikes are captured here.
                    )

                if v_bar_j <= 0:
                    Z, (t_bar, l_bar, v_bar, eta, tau) = ztlvet_record[-1]

                    raise ValueError(
                        "Squib load condition detected: Shot stopped in bore.\n"
                        + "Shot is last calculated at {:.0f} mm at {:.0f} mm/s after {:.0f} ms".format(
                            l_bar * self.l_0 * 1e3,
                            v_bar * self.v_j * 1e3,
                            t_bar * tScale * 1e3,
                        )
                    )

                if any(v < 0 for v in (t_bar_j, l_bar_j, p_bar_j)):
                    raise ValueError(
                        "Numerical Integration diverged: negative"
                        + " values encountered in results.\n"
                        + "{:.0f} ms, {:.0f} mm, {:.0f} m/s, {:.0f} MPa".format(
                            t_bar_j * tScale * 1e3,
                            l_bar_j * self.l_0 * 1e3,
                            v_bar_j * self.v_j,
                            p_bar_j * pScale / 1e6,
                        )
                    )  # this will catch any case where t, l, p are negative

                (t_bar_i, l_bar_i, v_bar_i, eta_i, tau_i) = (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j)
                Z_i = Z_j
                """
                this way the group of values denoted by _i is always updated
                as a group.
                """
                Z_j += Delta_Z / N

        if t_bar_i == 0:
            raise ValueError("burnout point found to be at the origin.")

        """
        Cludge code to force the SoE past the discontinuity at Z = Z_b, since
        we wrote the SoE to be be piecewise continous from (0, Z_b] and (Z_b, +inf)
        it is necessary to do this to prevent the RKF integrator coming up with
        irreducible error estimates and driving the step size to 0 around Z = Z_b
        """
        if isBurnOutContained:
            Z_i = Z_b + tol

        record.extend(
            (t_bar, (l_bar, self.f_psi_Z(Z), v_bar, self._f_p_bar(Z, l_bar, eta, tau), eta, tau))
            for (Z, (t_bar, l_bar, v_bar, eta, tau)) in ztlvet_record
        )

        if progressQueue is not None:
            progressQueue.put(10)

        if isBurnOutContained:
            logger.info("integrated to burnout point.")
        else:
            logger.warning("shot exited barrel before burnout.")

        """
        Subscript e indicate exit condition.
        At this point, since its guaranteed that point i will be further
        towards the chamber of the firearm than point e, we do not have
        to worry about the dependency of the correct behaviour of this ODE
        on the positive direction of integration.
        """

        ltzvet_record = []
        (_, (t_bar_e, Z_e, v_bar_e, eta_e, tau_e), (t_bar_err, Z_err, v_bar_err, eta_err, tau_err)) = RKF78(
            self._ode_l,
            (t_bar_i, Z_i, v_bar_i, eta_i, tau_i),
            l_bar_i,
            l_g_bar,
            relTol=tol,
            absTol=tol**2,
            record=ltzvet_record,
        )

        record.extend(
            (t_bar, (l_bar, self.f_psi_Z(Z), v_bar, self._f_p_bar(Z, l_bar, eta, tau), eta, tau))
            for (l_bar, (t_bar, Z, v_bar, eta, tau)) in ltzvet_record
        )

        updBarData(
            tag=POINT_EXIT,
            t_bar=t_bar_e,
            l_bar=l_g_bar,
            Z=Z_e,
            v_bar=v_bar_e,
            eta=eta_e,
            tau=tau_e,
            t_bar_err=t_bar_err,
            l_bar_err=0,
            Z_err=Z_err,
            v_bar_err=v_bar_err,
            eta_err=eta_err,
            tau_err=tau_err,
        )

        if progressQueue is not None:
            progressQueue.put(20)

        logger.info("integrated and determined conditions at exit point.")

        t_bar_f = None
        if Z_b > 1.0 and Z_e >= 1.0:  # fracture point exist and is contained
            """
            Subscript f indicate fracture condition
            ODE w.r.t Z is integrated from Z_0 to 1, from onset of projectile
            movement to charge fracture
            """
            (
                _,
                (t_bar_f, l_bar_f, v_bar_f, eta_f, tau_f),
                (t_bar_err_f, l_bar_err_f, v_bar_err_f, eta_err_f, tau_err_f),
            ) = RKF78(self._ode_Z, (0, 0, 0, 0, 1), Z_0, 1, relTol=tol, absTol=tol**2)

            updBarData(
                tag=POINT_FRACTURE,
                t_bar=t_bar_f,
                l_bar=l_bar_f,
                Z=1,
                eta=eta_f,
                tau=tau_f,
                v_bar=v_bar_f,
                t_bar_err=t_bar_err_f,
                l_bar_err=l_bar_err_f,
                Z_err=0,
                v_bar_err=v_bar_err_f,
                eta_err=eta_err_f,
                tau_err=tau_err_f,
            )

            logger.info("integrated and determined conditions at propellant frature point.")

        t_bar_b = None
        if isBurnOutContained:
            """
            Subscript b indicated burnout condition
            ODE w.r.t Z is integrated from Z_0 to Z_b, from onset of projectile
            movement to charge burnout.
            """

            (
                _,
                (t_bar_b, l_bar_b, v_bar_b, eta_b, tau_b),
                (t_bar_err_b, l_bar_err_b, v_bar_err_b, eta_err_b, tau_err_b),
            ) = RKF78(self._ode_Z, (0, 0, 0, 0, 1), Z_0, Z_b, relTol=tol, absTol=tol**2)

            updBarData(
                tag=POINT_BURNOUT,
                t_bar=t_bar_b,
                l_bar=l_bar_b,
                Z=Z_b,
                v_bar=v_bar_b,
                eta=eta_b,
                tau=tau_b,
                t_bar_err=t_bar_err_b,
                l_bar_err=l_bar_err_b,
                Z_err=0,
                v_bar_err=v_bar_err_b,
                eta_err=eta_err_b,
                tau_err=tau_err_b,
            )
            logger.info("integrated and determined conditions at propellant burnout point.")

        if progressQueue is not None:
            progressQueue.put(30)

        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the
        peak pressure point. As well as, not having starting issues.

        we hereby simply solve p golden section searching it from origin
        to point e, i.e. inside the barrel.
        """

        def g(t, tag):
            Z, l_bar, v_bar, eta, tau = RKF78(self._ode_t, (Z_0, 0, 0, 0, 1), 0, t, relTol=tol, absTol=tol**2)[1]
            p_bar = self._f_p_bar(Z, l_bar, eta, tau)

            if tag == POINT_PEAK_AVG:
                return p_bar

            Ps, P0, Pb, Vb = self._toPsP0PbVb(l_bar * self.l_0, v_bar * self.v_j, p_bar * pScale, tau * self.T_v, eta)

            if tag == POINT_PEAK_BREECH:
                return Pb / pScale
            elif tag == POINT_PEAK_STAG:
                return P0 / pScale
            elif tag == POINT_PEAK_SHOT:
                return Ps / pScale
            else:
                raise ValueError("tag not handled.")

        def findPeak(g, tag):
            """
            tolerance is specified a bit differently for gold section search
            GSS tol is the length between the upper bound and lower bound
            of the maxima/minima, thus by our definition of tolerance (one-sided),
            we take the median value.
            """

            t_bar_tol = tol * min(t for t in (t_bar_e, t_bar_b, t_bar_f) if t is not None)

            t_bar_1, t_bar_2 = gss(
                g,
                0,
                t_bar_e if t_bar_b is None else t_bar_b,
                x_tol=t_bar_tol,
                findMin=False,
            )

            t_bar = 0.5 * (t_bar_1 + t_bar_2)

            (_, (Z, l_bar, v_bar, eta, tau), (Z_err, l_bar_err, v_bar_err, eta_err, tau_err)) = RKF78(
                self._ode_t, (Z_0, 0, 0, 0, 1), 0, t_bar, relTol=tol, absTol=tol**2
            )
            t_bar_err = 0.5 * t_bar_tol

            updBarData(
                tag=tag,
                t_bar=t_bar,
                l_bar=l_bar,
                Z=Z,
                v_bar=v_bar,
                eta=eta,
                tau=tau,
                t_bar_err=t_bar_err,
                l_bar_err=l_bar_err,
                Z_err=Z_err,
                v_bar_err=v_bar_err,
                eta_err=eta_err,
                tau_err=tau_err,
            )

        for i, peak in enumerate([POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH, POINT_PEAK_STAG]):
            findPeak(lambda x: g(x, peak), peak)
            logger.info(f"found peak conditions {peak}.")
            if progressQueue is not None:
                progressQueue.put(30 + (i + 1) * 12.5)  # 42.5, 55, 67.5, 80

        """
        populate data for output purposes
        """
        logger.info(f"sampling for {step} points.")
        if dom == DOMAIN_TIME:
            # fmt:off
            (Z_j, l_bar_j, v_bar_j, t_bar_j, eta_j, tau_j) = (
                Z_0, 0, 0, 0, 0, 1
            )
            # fmt: on
            for j in range(step):
                t_bar_k = t_bar_e / (step + 1) * (j + 1)
                (_, (Z_j, l_bar_j, v_bar_j, eta_j, tau_j), (Z_err, l_bar_err, v_bar_err, eta_err, tau_err)) = RKF78(
                    self._ode_t, (Z_j, l_bar_j, v_bar_j, eta_j, tau_j), t_bar_j, t_bar_k, relTol=tol, absTol=tol**2
                )
                t_bar_j = t_bar_k

                updBarData(
                    tag="",
                    t_bar=t_bar_j,
                    l_bar=l_bar_j,
                    Z=Z_j,
                    v_bar=v_bar_j,
                    eta=eta_j,
                    tau=tau_j,
                    t_bar_err=0,
                    l_bar_err=l_bar_err,
                    Z_err=Z_err,
                    v_bar_err=v_bar_err,
                    eta_err=eta_err,
                    tau_err=tau_err,
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
             of value by subscipt i will not guarantee burning is still
             ongoing).
            """
            t_bar_j = 0.5 * t_bar_i
            Z_j, l_bar_j, v_bar_j, eta_j, tau_j = RKF78(
                self._ode_t,
                (Z_0, 0, 0, 0, 1),
                0,
                t_bar_j,
                relTol=tol,
                absTol=tol**2,
            )[1]

            for j in range(step):
                l_bar_k = l_g_bar / (step + 1) * (j + 1)

                (_, (t_bar_j, Z_j, v_bar_j, eta_j, tau_j), (t_bar_err, Z_err, v_bar_err, eta_err, tau_err)) = RKF78(
                    self._ode_l, (t_bar_j, Z_j, v_bar_j, eta_j, tau_j), l_bar_j, l_bar_k, relTol=tol, absTol=tol**2
                )

                l_bar_j = l_bar_k

                updBarData(
                    tag="",
                    t_bar=t_bar_j,
                    l_bar=l_bar_j,
                    Z=Z_j,
                    v_bar=v_bar_j,
                    eta=eta_j,
                    tau=tau_j,
                    t_bar_err=t_bar_err,
                    l_bar_err=0,
                    Z_err=Z_err,
                    v_bar_err=v_bar_err,
                    eta_err=eta_err,
                    tau_err=tau_err,
                )

        logger.info("sampling completed.")

        if progressQueue is not None:
            progressQueue.put(100)

        """
        sort the data points
        """

        data = []
        error = []
        p_trace = []
        l_c = self.l_c
        # l_g = self.l_g

        for bar_dataLine, bar_errorLine in zip(bar_data, bar_err):
            dtag, t_bar, l_bar, Z, v_bar, p_bar, eta, tau = bar_dataLine
            (etag, t_bar_err, l_bar_err, Z_err, v_bar_err, p_bar_err, eta_err, tau_err) = bar_errorLine

            t, t_err = (t_bar * tScale, t_bar_err * tScale)
            l, l_err = (l_bar * self.l_0, l_bar_err * self.l_0)
            psi, psi_err = (self.f_psi_Z(Z), abs(self.f_sigma_Z(Z) * Z_err))
            v, v_err = v_bar * self.v_j, v_bar_err * self.v_j

            p, p_err = p_bar * pScale, (p_bar_err if isinstance(p_bar_err, str) else p_bar_err * pScale)
            T = tau * self.T_v
            T_err = tau_err * self.T_v

            ps, p0, pb, vb = self._toPsP0PbVb(l, v, p, T, eta)

            p_line = []
            for i in range(step):  # 0....step -1
                x = i / step * (l + l_c)
                px = self._toPx(l, v, vb, ps, T, eta, x)

                p_line.append(PressureProbePoint(x, px))

            p_line.append(PressureProbePoint(l + l_c, ps))
            p_trace.append(PressureTraceEntry(dtag, T, p_line))

            data.append(RecoillessTableEntry(dtag, t, l, psi, v, vb, pb, p0, p, ps, T, eta))
            error.append(
                RecoillessErrorEntry(etag, t_err, l_err, psi_err, v_err, None, None, None, p_err, None, T_err, eta_err)
            )

        for t_bar, (l_bar, psi, v_bar, p_bar, eta, tau) in record:
            t = t_bar * tScale
            if t in [tableEntry.time for tableEntry in data]:
                continue
            l, v, p, T = (l_bar * self.l_0, v_bar * self.v_j, p_bar * pScale, tau * self.T_v)

            ps, p0, pb, vb = self._toPsP0PbVb(l, v, p, T, eta)

            p_line = []
            for i in range(step):  # 0....step -1
                x = i / step * (l + l_c)
                px = self._toPx(l, v, vb, ps, T, eta, x)
                pp = PressureProbePoint(x, px)
                p_line.append(pp)

            p_line.append(PressureProbePoint(l + l_c, ps))
            p_trace.append(PressureTraceEntry("*", T, p_line))

            data.append(RecoillessTableEntry("*", t, l, psi, v, vb, pb, p0, p, ps, T, eta))
            error.append(RecoillessErrorEntry("L"))

        data, error, p_trace = zip(
            *((a, b, c) for a, b, c in sorted(zip(data, error, p_trace), key=lambda entries: entries[0].time))
        )

        logger.info("scaled intermediate results to SI.")

        recoilessResult = RecoillessResult(self, data, error, p_trace)

        if self.material is None:
            logger.warning("material is not specified, skipping structural calculation.")

        else:
            try:
                self._getStructural(recoilessResult, step, tol)
            except Exception:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                logger.error("exception occured during structural calculation:")
                logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
                logger.info("structural calculation skipped.")

        logger.info("integration concluded.")

        return recoilessResult

    def _toPsP0PbVb(self, l, v, p, T, eta):
        """
        Diagramatic explanation of the calculated values:
            ----\___     __________________
                    |---|                  |_________________________________
                       ==                                       |        \
           Nozzle Breech|    Chamber                  Barrel    |  Shot  |>
                       ==                   ____________________|________/____
                 ___|---|__________________|
            ----/
        Cross-section of the breech face:
           _
         *- -*
        | 0 0 |   0+0: Nozzle throat area, S_j
         *-_-*

        Ps  : Shot base pressure
        Pb  : Breech pressure, chamber pressure at rearward of chamber
        P0  : Stagnation point pressure
        vb  : Rearward flow velocity at the rear of chamber
        """
        tau = T / self.T_v
        y = self.omega * eta
        m_dot = self.C_A * self.v_j * self.S_j * p / (self.f * tau**0.5)
        # mass flow rate, rearward
        Sb = self.S * self.chi_k
        vb = m_dot * (self.V_0 + self.S * l) / (Sb * (self.omega - y))
        # velocity impinging upon the rear of the breech before nozzle constriction

        H_1 = vb / v if v != 0 else inf
        H_2 = 2 * self.phi_1 * self.m / (self.omega - y) + 1
        H = min(H_1, H_2)

        ps = p / (1 + (self.omega - y) / (3 * self.phi_1 * self.m) * (1 - 0.5 * H))  # shot base pressure
        p0 = ps * (1 + (self.omega - y) / (2 * self.phi_1 * self.m) * (1 + H) ** -1)  # stagnation point pressure
        pb = ps * (1 + (self.omega - y) / (2 * self.phi_1 * self.m) * (1 - H)) if H == H_1 else 0  # breech pressure
        # l0 = H / (1 + H) * l  # location of the stagnation point
        return ps, p0, pb, vb

    def _toPx(self, l, v, vb, ps, T, eta, x):
        m = self.m
        omega = self.omega
        phi_1 = self.phi_1
        y = self.omega * eta

        """
        convert x, the physical displacement from breech bottom, to
        effective length in the equivalent gun.
        """

        L_1 = l
        L_0 = self.l_0 / self.chi_k  # physical length of the chamber.

        A_1 = self.S
        A_0 = A_1 * self.chi_k

        if x < L_0:
            z = x * A_0 / (L_0 * A_0 + L_1 * A_1)
        else:
            z = (L_0 * A_0 + (x - L_0) * A_1) / (L_0 * A_0 + L_1 * A_1)

        H_1 = vb / v if v != 0 else inf
        H_2 = 2 * phi_1 * self.m / (self.omega - y) + 1
        H = min(H_1, H_2)

        px = ps * (1 + (omega - y) / (2 * phi_1 * m) * (1 + H) * (1 - z**2) - (omega - y) / (phi_1 * m) * H * (1 - z))
        return px

    @staticmethod
    def _getCf(gamma, Sr, tol=1e-5):
        """
        takes the adiabatic index and area ration between constriction throat
        and the exit, calculate the thrust factor Cf
        See Hunt (1953) appendix I.A01-A03
        Sr = Se/St
        Vr = V/Vt
        """
        Vr_old = 0
        Vr = 1
        while abs(Vr - Vr_old) / Vr > tol:
            Vr_old = Vr
            Vr = ((gamma + 1) / (gamma - 1) * (1 - 2 / (gamma + 1) * (Vr * Sr) ** (1 - gamma))) ** 0.5

        Cf = (2 / (gamma + 1)) ** (gamma / (gamma - 1)) * (gamma * Vr + Sr ** (1 - gamma) * Vr ** (-gamma))

        return Cf

    def _getStructural(self, recoillessResult: RecoillessResult, step: int, tol: float):
        logger = logging.getLogger("recoiless.structure")
        logger.info("commencing structural calculation")
        l_c = self.l_c
        l_g = self.l_g
        chi_k = self.chi_k
        sigma = self.material.Y
        S = self.S
        gamma = self.theta + 1

        S_j_bar = self.S_j_bar
        A_bar = self.A_bar
        r_s = 0.5 * self.caliber  # radius of the shot.
        r_b = r_s * chi_k**0.5  # chamber/breech radius
        r_t = r_b * S_j_bar**0.5  # throat radius
        r_n = r_t * A_bar**0.5  # nozzle exit radius

        x_probes = (
            [i / step * l_c for i in range(step)]
            + [l_c * (1 - tol)]
            + [i / step * l_g + l_c for i in range(step)]
            + [l_g + l_c]
        )

        p_probes = [0] * len(x_probes)
        for recoillessTableEntry in recoillessResult.tableData:
            l = recoillessTableEntry.travel
            v = recoillessTableEntry.velocity
            vb = recoillessTableEntry.outflowVelocity
            ps = recoillessTableEntry.shotPressure
            T = recoillessTableEntry.temperature
            eta = recoillessTableEntry.outflowFraction

            for i, x in enumerate(x_probes):
                if (x - l_c) <= l:
                    p_x = self._toPx(l, v, vb, ps, T, eta, x)
                    p_probes[i] = max(p_probes[i], p_x)
                else:
                    break

        # strength requirement given structural safety factor.
        for i in range(len(p_probes)):
            p_probes[i] *= self.ssf

        """
        subscript k denote the end of the chamber
        subscript b denote the throat of nozzle
        subscript a denote the nozzle base.

        beta is the half angle of the constriction section
        alpha is the half angle of the expansion section
        """
        beta = 45
        l_b = (r_b - r_t) / tan(beta * pi / 180)
        alpha = 15
        l_a = (r_n - r_t) / tan(alpha * pi / 180)

        p_entry = p_probes[0]

        neg_x_probes = (
            [-l_b - l_a + l_a * i / step for i in range(step)]
            + [-l_b - tol * l_a]
            + [-l_b + l_b * i / step for i in range(step)]
            + [0]
        )
        neg_p_probes = []
        neg_S_probes = []

        r_probes = []

        for x in neg_x_probes:
            if x < -l_b:
                k = (x + l_b + l_a) / l_a
                r = k * r_t + (1 - k) * r_n
                p = Recoilless._getPr(gamma, (r / r_t) ** 2, tol)[1] * p_entry

            else:
                k = (x + l_b) / l_b
                r = k * r_b + (1 - k) * r_t
                p = Recoilless._getPr(gamma, (r / r_t) ** 2, tol)[0] * p_entry

            if p > 3**-0.5 * sigma:
                raise ValueError(
                    f"Limit to conventional construction ({sigma * 1e-6:.3f} MPa)" + " exceeded in nozzle."
                )

            r_probes.append(r)
            neg_p_probes.append(p)
            neg_S_probes.append(pi * r**2)

        if self.is_af:
            V_n, rho_n = Gun._Vrho_k(neg_x_probes, neg_p_probes, neg_S_probes, sigma, tol)

        else:
            V_n, rho_n = 0, []

            for p in neg_p_probes:
                y = p / sigma
                rho = ((1 + y * (4 - 3 * y**2) ** 0.5) / (1 - 3 * y**2)) ** 0.5
                rho_n.append(rho)

            for i in range(len(neg_x_probes) - 1):
                x_0, x_1 = neg_x_probes[i], neg_x_probes[i + 1]
                h = x_1 - x_0
                S_0, S_1 = neg_S_probes[i], neg_S_probes[i + 1]
                rho_0, rho_1 = rho_n[i], rho_n[i + 1]
                dV = 0.5 * ((rho_0**2 - 1) * S_0 + (rho_1**2 - 1) * S_1) * h
                V_n += dV

        hull = []

        for x, r, rho in zip(neg_x_probes, r_probes, rho_n):
            hull.append(OutlineEntry(x, r, r * rho))

        nozzle_mass = V_n * self.material.rho
        logger.info("structural calculation of nozzle section complete.")

        rho_probes = []
        V = 0
        if self.is_af:
            """
            m : r_n / r_i
            k : r_o / r_i
            n : p_vM_max / sigma

            1 < m < k

            The point of optimum autofrettage, or the minimum autofrettage
            necessary to contain the working pressure, is
            """
            i = step + 1
            x_c, p_c = x_probes[:i], p_probes[:i]  # c for chamber
            x_b, p_b = x_probes[i:], p_probes[i:]  # b for barrel
            V_c, rho_c = Gun._Vrho_k(x_c, p_c, [S * chi_k for _ in x_c], sigma, tol, k_min=rho_n[-1])  # c for chamber
            V_b, rho_b = Gun._Vrho_k(
                x_b,
                p_b,
                [S for _ in x_b],
                sigma,
                tol,
                k_max=rho_c[-1] * chi_k**0.5,
                p_ref=max(p_c),
            )  # b for bore
            V = V_c + V_b
            rho_probes = rho_c + rho_b

        else:
            """
            The yield criterion chosen here is the fourth strength
            theory (von Mises) as it is generally accepted to be the most
            applicable for this application.

            The limiting stress points circumferentially along the circum-
            ference of the barrel.

            P_4 = sigma_e * (rho^2-1)/ (3*rho**4 + 1) ** 0.5
            lim (x->inf) (x^2-1)/sqrt(3*x**4+1) = 1/sqrt(3)

            the inverse of (x^2-1)/sqrt(3*x**4+1) is:
            sqrt(
                [-sqrt(-x**2*(3*x**2-4)) - 1]/(3 * x**2 - 1)
            )
            (x < -1 or x > 1)
            and
            sqrt(
                [sqrt(-x**2*(3*x**2-4)) - 1]/(3 * x**2 - 1)
            )
            (x from -1 to 1)
            """
            for p in p_probes:
                y = p / sigma
                if y > 3**-0.5:
                    raise ValueError(
                        f"Limit to conventional construction ({sigma * 3*1e-6:.3f} MPa)" + " exceeded in tube section."
                    )
                rho = ((1 + y * (4 - 3 * y**2) ** 0.5) / (1 - 3 * y**2)) ** 0.5
                rho_probes.append(rho)

            for i in range(len(x_probes) - 1):
                x_0 = x_probes[i]
                x_1 = x_probes[i + 1]
                rho_0 = rho_probes[i]
                rho_1 = rho_probes[i + 1]
                dV = (rho_1**2 + rho_0**2 - 2) * 0.5 * S * (x_1 - x_0)
                if x_1 < l_c:
                    V += dV * chi_k
                else:
                    V += dV

        tube_mass = V * self.material.rho
        logger.info("structural calculation of tube section complete.")

        for x, rho in zip(x_probes, rho_probes):
            if x < l_c:
                hull.append(OutlineEntry(x, r_b, rho * r_b))
            else:
                hull.append(OutlineEntry(x, r_s, rho * r_s))

        recoillessResult.outline = hull
        recoillessResult.tubeMass = tube_mass
        recoillessResult.breechMass = nozzle_mass

        logger.info("structural calculation results attached.")

    @staticmethod
    def _getAr(gamma, Pr):
        if Pr == 0 or Pr == 1:
            return inf
        return (
            (0.5 * (gamma + 1)) ** (1 / (gamma - 1))
            * Pr ** (1 / gamma)
            * ((gamma + 1) / (gamma - 1) * (1 - Pr ** ((gamma - 1) / gamma))) ** 0.5
        ) ** -1

    @staticmethod
    def _getPr(gamma, Ar, tol):
        """
        Given an area ratio Ar, calculate the pressure ratio Pr for a
        converging-diverging nozzle. One area ratio corresponds to two
        pressure ratio, for subsonic and supersonic, respectively.

        Ar: Nozzle cs area at probe point x over throat area
            Ar = Ax / At
        Pr: pressure ratio at probe point x over upstream pressure
            Pr = Px / Pc

        returns:
            Pr_sub: subsonic solution
            Pr_sup: supersonic solution

        """

        # pressure ratio of nozzle throat over the upstream pressure
        Pr_c = (2 / (gamma + 1)) ** (gamma / (gamma - 1))

        if Ar == 1:
            return Pr_c, Pr_c
        else:
            Pr_sup, _ = dekker(lambda Pr: Recoilless._getAr(gamma, Pr), 0, Pr_c, y=Ar, y_rel_tol=tol)

            Pr_sub, _ = dekker(lambda Pr: Recoilless._getAr(gamma, Pr), Pr_c, 1, y=Ar, y_rel_tol=tol)

            return Pr_sub, Pr_sup


if __name__ == "__main__":
    """standard 7 hole cylinder has d_0=e_1, hole dia = 0.5 * arc width
    d_0 = 4*2e_1+3*d_0 = 11 * e_1
    """

    from tabulate import tabulate

    compositions = GrainComp.readFile("data/propellants.csv")

    M17 = compositions["M17"]
    M1 = compositions["M1"]
    from prop import SimpleGeometry

    M17C = Propellant(M17, SimpleGeometry.CYLINDER, None, 2.5)
    M1C = Propellant(M1, SimpleGeometry.CYLINDER, None, 10)
    lf = 0.6
    print("DELTA/rho:", lf)
    test = Recoilless(
        caliber=0.082,
        shotMass=2,
        propellant=M1C,
        grainSize=4e-4,
        chargeMass=0.3,
        chamberVolume=0.3 / M1C.rho_p / lf,
        startPressure=10e6,
        lengthGun=3.5,
        nozzleExpansion=2.0,
        chambrage=1.0,
    )
    record = []
    # fmt: off
    print("\nnumerical: time")
    print(
        tabulate(
            test.integrate(10, 1e-3, dom=DOMAIN_TIME).getRawTableData(),
            headers=("tag", "t", "l", "psi", "v", "vb", "pb", "p0", "p", "ps",
                     "T", "eta")
        )
    )
    print("\nnumerical: length")
    print(
        tabulate(
            test.integrate(10, 1e-3, dom=DOMAIN_LEN).getRawTableData(),
            headers=("tag", "t", "l", "psi", "v", "vb", "pb", "p0", "p", "ps",
                     "T", "eta")
        )
    )
    # fmt:on
