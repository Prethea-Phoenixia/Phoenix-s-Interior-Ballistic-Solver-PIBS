from __future__ import annotations

import logging
import sys
import traceback
from dataclasses import dataclass
from math import inf, pi

from pibs.ballistics.material import Material

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
)
from .gun import Gun
from .num import dekker, gss, rkf78
from .prop import Composition, Propellant

logger = logging.getLogger(__name__)


from .generics import (
    DelegatesPropellant,
    GenericEntry,
    GenericErrorEntry,
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


@dataclass
class RecoillessErrorEntry(GenericErrorEntry):
    outflow_velocity: float | None = None
    stag_pressure: float | None = None
    outflow_fraction: float | None = None


@dataclass
class RecoillessResult(GenericResult):
    gun: Recoilless
    table_data: list[RecoillessTableEntry]
    error_data: list[RecoillessErrorEntry]


class Recoilless(DelegatesPropellant):
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
        drag_coefficient: float = 0.0,
        nozzle_efficiency: float = 0.92,
        structural_material: Material = None,
        structural_safety_factor: float = 1.1,
        autofrettage: bool = True,
        # traveling_charge: bool = True,
        **_,
    ):
        super().__init__(propellant=propellant)

        if any(
            (
                caliber <= 0,
                shot_mass <= 0,
                charge_mass <= 0,
                web <= 0,
                chamber_volume <= 0,
                length_gun <= 0,
                nozzle_expansion < 1,
                nozzle_efficiency > 1,
                nozzle_efficiency <= 0,
                drag_coefficient < 0,
                drag_coefficient >= 1,
                start_pressure < 0,
                structural_safety_factor <= 1,
            )
        ):
            raise ValueError("Invalid gun parameters")

        self.propellant = propellant
        self.caliber = caliber
        # self.is_tc = traveling_charge

        e_1 = 0.5 * web
        self.s = (caliber / 2) ** 2 * pi
        self.m = shot_mass
        self.w = charge_mass
        self.V_0 = chamber_volume
        self.p_0 = start_pressure
        self.l_g = length_gun
        self.chi_0 = nozzle_efficiency
        self.a_bar = nozzle_expansion

        self.chi_k = chambrage  # ration of l_0 / l_chamber
        self.l_0 = self.V_0 / self.s
        self.l_c = self.l_0 / self.chi_k

        self.Delta = self.w / self.V_0
        self.phi_1 = 1 / (1 - drag_coefficient)  # drag work coefficient
        self.phi = self.phi_1 + self.w / (3 * self.m)

        self.v_j = (2 * self.f * self.w / (self.theta * self.phi * self.m)) ** 0.5

        self.material = structural_material
        self.ssf = structural_safety_factor
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
            self.s**2
            * e_1**2
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
        )  # flow rate value

        self.k_1, self.c_a_bar, self.p_a_bar, self.S_j, self.s_j_bar, self.c_f = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        self.z_0, self.psi_0 = 0, 0

    def f_p_bar(self, z, l_bar, eta, tau, psi=None):
        psi = psi if psi else self.f_psi_z(z)
        l_psi_bar = 1 - self.Delta * ((1 - psi) / self.rho_p + self.alpha * (psi - eta))

        p_bar = tau / (l_bar + l_psi_bar) * (psi - eta)

        return max(p_bar, self.p_a_bar)

    def ode_t(self, _, z, l_bar, v_bar, eta, tau, __):
        psi = self.f_psi_z(z)
        dpsi = self.f_sigma_z(z)  # dpsi/dZ
        p_bar = self.f_p_bar(z, l_bar, eta, tau, psi)

        if self.c_a_bar != 0 and v_bar > 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (
                0.25 * k * (k + 1) * v_r**2 + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a_bar
        else:
            p_d_bar = 0

        dz = (0.5 * self.theta / self.B) ** 0.5 * p_bar**self.n

        dl_bar = v_bar
        dv_bar = self.theta * 0.5 * (p_bar - p_d_bar)

        # dv_bar /= (1 + self.w / self.m * (1 - psi)) if self.is_tc else 1

        deta = self.c_a * self.s_j_bar * p_bar * tau**-0.5  # deta / dt_bar
        dtau = ((1 - tau) * (dpsi * dz) - 2 * v_bar * dv_bar - self.theta * tau * deta) / (psi - eta)  # dtau/dt_bar

        return dz, dl_bar, dv_bar, deta, dtau

    def ode_l(self, l_bar, _, z, v_bar, eta, tau, __):
        """length domain ode of internal ballistics
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition.

        in general, d/dl_bar = d/dt_bar * dt_bar/dl_bar

        """

        psi = self.f_psi_z(z)
        dpsi = self.f_sigma_z(z)  # dpsi/dZ
        p_bar = self.f_p_bar(z, l_bar, eta, tau, psi)

        if self.c_a_bar != 0 and v_bar > 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (+0.25 * k * (k + 1) * v_r**2 + k * v_r * (1 + (0.25 * (k + 1) * v_r) ** 2) ** 0.5) * self.p_a_bar
        else:
            p_d_bar = 0

        dz = (0.5 * self.theta / self.B) ** 0.5 * p_bar**self.n / v_bar
        dv_bar = self.theta * 0.5 * (p_bar - p_d_bar) / v_bar
        # dv_bar /= (1 + self.w / self.m * (1 - psi)) if self.is_tc else 1
        # dv_bar/dl_bar
        dt_bar = 1 / v_bar  # dt_bar / dl_bar

        deta = self.c_a * self.s_j_bar * p_bar * tau**-0.5 * dt_bar  # deta / dl_bar

        dtau = ((1 - tau) * (dpsi * dz) - 2 * v_bar * dv_bar - self.theta * tau * deta) / (psi - eta)  # dtau/dl_bar

        return dt_bar, dz, dv_bar, deta, dtau

    def ode_z(self, z, _, l_bar, v_bar, eta, tau, __):
        psi = self.f_psi_z(z)
        dpsi = self.f_sigma_z(z)  # dpsi/dZ
        p_bar = self.f_p_bar(z, l_bar, eta, tau, psi)

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
        # dv_bar /= (1 + self.w / self.m * (1 - psi)) if self.is_tc else 1
        deta = self.c_a * self.s_j_bar * p_bar * tau**-0.5 * dt_bar  # deta / dZ
        dtau = ((1 - tau) * dpsi - 2 * v_bar * dv_bar - self.theta * tau * deta) / (psi - eta)

        return dt_bar, dl_bar, dv_bar, deta, dtau

    def integrate(
        self,
        step: int = 10,
        tol: float = 1e-5,
        dom: Domains = DOMAIN_TIME,
        ambient_rho: float = 1.204,
        ambient_p: float = 101.325e3,
        ambient_gamma: float = 1.4,
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

        ambient_p = max(ambient_p, 1)
        ambient_gamma = max(ambient_gamma, 1)

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
            lambda z: self.propellant.f_psi_z(z) - self.psi_0, 0, 1, x_tol=tol, y_rel_tol=tol, y_abs_tol=tol**2
        )

        record = []
        if any((step < 0, tol < 0)):
            raise ValueError("Invalid integration specification")

        if ambient_rho < 0:
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
                "Achieving recoilless condition necessitates"
                + " a larger throat area than could be fit into"
                + " breech face."
            )
        self.S_j = self.s_j_bar * self.s

        t_scale = self.l_0 / self.v_j
        p_scale = self.f * self.Delta

        # ambient conditions
        self.p_a_bar = ambient_p / p_scale
        if ambient_rho != 0:
            self.c_a_bar = (ambient_gamma * ambient_p / ambient_rho) ** 0.5 / self.v_j
        else:
            self.c_a_bar = 0

        self.k_1 = ambient_gamma

        l_g_bar = self.l_g / self.l_0
        p_bar_0 = self.p_0 / p_scale
        z_b = self.z_b
        z_0 = self.z_0

        bar_data = []
        bar_err = []

        def upd_bar_data(
            tag, t_bar, l_bar, z, v_bar, eta, tau, t_bar_err, l_bar_err, z_err, v_bar_err, eta_err, tau_err
        ):
            p_bar = self.f_p_bar(z, l_bar, eta, tau)
            p_bar_err = "N/A"
            bar_data.append((tag, t_bar, l_bar, z, v_bar, p_bar, eta, tau))
            """
            Worst case scenario for pressure deviation is when:
            Z higher than actual, l lower than actual, v lower than actual
            more burnt propellant, less volume, lower speed of projectile.
            """
            bar_err.append(("L", t_bar_err, l_bar_err, z_err, v_bar_err, p_bar_err, eta_err, tau_err))

        upd_bar_data(
            tag=POINT_START,
            t_bar=0,
            l_bar=0,
            z=z_0,
            v_bar=0,
            eta=0,
            tau=1,
            t_bar_err=0,
            l_bar_err=0,
            z_err=0,
            v_bar_err=0,
            eta_err=0,
            tau_err=0,
        )

        record.append((0, (0, self.psi_0, 0, p_bar_0, 0, 1)))

        z_i = z_0
        z_j = z_b
        n = 1
        delta_z = z_b - z_0

        t_bar_i, l_bar_i, v_bar_i, p_bar_i, eta_i, tau_i = 0, 0, 0, p_bar_0, 0, 1

        is_burn_out_contained = True

        """
        Instead of letting the integrator handle the heavy lifting, we
        partition Z and integrate upwards until either barrel exit or
        burnout has been achieved. This seek-and-validate inspired
        from bisection puts the group of value denoted by subscript
        i within or on the muzzle, with the propellant either still
        burning or right on the burnout point..
        """
        ztlvet_record = [(z_0, (0, 0, 0, 0, 1))]
        p_max = 1e9  # 1GPa
        p_bar_max = p_max / p_scale

        # noinspection PyUnusedLocal
        def abort(x, ys, record):

            z, (_, l_bar, v_bar, eta, tau) = x, ys
            p_bar = self.f_p_bar(z, l_bar, eta, tau)

            return any((l_bar > l_g_bar, p_bar > p_bar_max, v_bar <= 0, p_bar < 0))

        while z_i < z_b:  # terminates if burnout is achieved
            ztlvet_record_i = []
            if z_j == z_i:
                raise ValueError("Numerical accuracy exhausted in search of exit/burnout point.")
            try:
                if z_j > z_b:
                    z_j = z_b

                z, (t_bar_j, l_bar_j, v_bar_j, eta_j, tau_j), _ = rkf78(
                    self.ode_z,
                    (t_bar_i, l_bar_i, v_bar_i, eta_i, tau_i),
                    z_i,
                    z_j,
                    rel_tol=tol,
                    abs_tol=tol**2,
                    abort_func=abort,
                    record=ztlvet_record_i,
                    debug=True,
                )

                p_bar_j = self.f_p_bar(z_j, l_bar_j, eta_j, tau_j)

            except ValueError as e:
                ztlvet_record.extend(ztlvet_record_i)
                z, (t_bar, l_bar, v_bar, eta, tau) = ztlvet_record[-1]
                dt_bar, dl_bar, dv_bar, deta, dtau = self.ode_z(z, t_bar, l_bar, v_bar, eta, tau, _)

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
                        + " Unreasonably high pressure (>{:.0f} MPa) was encountered.".format(
                            p_max / 1e6
                        )  # in practice most of the pressure-realted spikes are captured here.
                    )

                if v_bar_j <= 0:
                    z, (t_bar, l_bar, v_bar, eta, tau) = ztlvet_record[-1]

                    raise ValueError(
                        "Squib load condition detected: Shot stopped in bore.\n"
                        + "Shot is last calculated at {:.0f} mm at {:.0f} mm/s after {:.0f} ms".format(
                            l_bar * self.l_0 * 1e3,
                            v_bar * self.v_j * 1e3,
                            t_bar * t_scale * 1e3,
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
                            p_bar_j * p_scale / 1e6,
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

        """
        Cludge code to force the SoE past the discontinuity at Z = z_b, since
        we wrote the SoE to be be piecewise continous from (0, z_b] and (z_b, +inf)
        it is necessary to do this to prevent the RKF integrator coming up with
        irreducible error estimates and driving the step size to 0 around Z = z_b
        """
        if is_burn_out_contained:
            z_i = z_b + tol

        record.extend(
            (t_bar, (l_bar, self.f_psi_z(z), v_bar, self.f_p_bar(z, l_bar, eta, tau), eta, tau))
            for (z, (t_bar, l_bar, v_bar, eta, tau)) in ztlvet_record
        )

        if is_burn_out_contained:
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
        (_, (t_bar_e, z_e, v_bar_e, eta_e, tau_e), (t_bar_err, z_err, v_bar_err, eta_err, tau_err)) = rkf78(
            self.ode_l,
            (t_bar_i, z_i, v_bar_i, eta_i, tau_i),
            l_bar_i,
            l_g_bar,
            rel_tol=tol,
            abs_tol=tol**2,
            record=ltzvet_record,
        )

        record.extend(
            (t_bar, (l_bar, self.f_psi_z(z), v_bar, self.f_p_bar(z, l_bar, eta, tau), eta, tau))
            for (l_bar, (t_bar, z, v_bar, eta, tau)) in ltzvet_record
        )

        upd_bar_data(
            tag=POINT_EXIT,
            t_bar=t_bar_e,
            l_bar=l_g_bar,
            z=z_e,
            v_bar=v_bar_e,
            eta=eta_e,
            tau=tau_e,
            t_bar_err=t_bar_err,
            l_bar_err=0,
            z_err=z_err,
            v_bar_err=v_bar_err,
            eta_err=eta_err,
            tau_err=tau_err,
        )

        t_bar_f = None
        if z_b > 1.0 and z_e >= 1.0:  # fracture point exist and is contained
            """
            Subscript f indicate fracture condition
            ODE w.r.t Z is integrated from Z_0 to 1, from onset of projectile
            movement to charge fracture
            """
            (
                _,
                (t_bar_f, l_bar_f, v_bar_f, eta_f, tau_f),
                (t_bar_err_f, l_bar_err_f, v_bar_err_f, eta_err_f, tau_err_f),
            ) = rkf78(self.ode_z, (0, 0, 0, 0, 1), z_0, 1, rel_tol=tol, abs_tol=tol**2)

            upd_bar_data(
                tag=POINT_FRACTURE,
                t_bar=t_bar_f,
                l_bar=l_bar_f,
                z=1,
                eta=eta_f,
                tau=tau_f,
                v_bar=v_bar_f,
                t_bar_err=t_bar_err_f,
                l_bar_err=l_bar_err_f,
                z_err=0,
                v_bar_err=v_bar_err_f,
                eta_err=eta_err_f,
                tau_err=tau_err_f,
            )

        t_bar_b = None
        if is_burn_out_contained:
            """
            Subscript b indicated burnout condition
            ODE w.r.t Z is integrated from Z_0 to z_b, from onset of projectile
            movement to charge burnout.
            """

            (
                _,
                (t_bar_b, l_bar_b, v_bar_b, eta_b, tau_b),
                (t_bar_err_b, l_bar_err_b, v_bar_err_b, eta_err_b, tau_err_b),
            ) = rkf78(self.ode_z, (0, 0, 0, 0, 1), z_0, z_b, rel_tol=tol, abs_tol=tol**2)

            upd_bar_data(
                tag=POINT_BURNOUT,
                t_bar=t_bar_b,
                l_bar=l_bar_b,
                z=z_b,
                v_bar=v_bar_b,
                eta=eta_b,
                tau=tau_b,
                t_bar_err=t_bar_err_b,
                l_bar_err=l_bar_err_b,
                z_err=0,
                v_bar_err=v_bar_err_b,
                eta_err=eta_err_b,
                tau_err=tau_err_b,
            )

        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the
        peak pressure point. As well as, not having starting issues.

        we hereby simply solve p golden section searching it from origin
        to point e, i.e. inside the barrel.
        """

        def g(t, tag):
            z, l_bar, v_bar, eta, tau = rkf78(self.ode_t, (z_0, 0, 0, 0, 1), 0, t, rel_tol=tol, abs_tol=tol**2)[1]
            p_bar = self.f_p_bar(z, l_bar, eta, tau)

            if tag == POINT_PEAK_AVG:
                return p_bar

            ps, p0, pb, vb = self.to_ps_p0_pb_vb(l_bar * self.l_0, v_bar * self.v_j, p_bar * p_scale, tau, eta)

            if tag == POINT_PEAK_BREECH:
                return pb / p_scale
            elif tag == POINT_PEAK_STAG:
                return p0 / p_scale
            elif tag == POINT_PEAK_SHOT:
                return ps / p_scale
            else:
                raise ValueError("tag not handled.")

        def find_peak(g, tag):
            """
            tolerance is specified a bit differently for gold section search
            GSS tol is the length between the upper bound and lower bound
            of the maxima/minima, thus by our definition of tolerance (one-sided),
            we take the median value.
            """

            t_bar_tol = tol * min(t for t in (t_bar_e, t_bar_b, t_bar_f) if t)

            t_bar_1, t_bar_2 = gss(g, 0, t_bar_e if t_bar_b is None else t_bar_b, x_tol=t_bar_tol, find_min=False)

            t_bar = 0.5 * (t_bar_1 + t_bar_2)

            (_, (z, l_bar, v_bar, eta, tau), (z_err, l_bar_err, v_bar_err, eta_err, tau_err)) = rkf78(
                self.ode_t, (z_0, 0, 0, 0, 1), 0, t_bar, rel_tol=tol, abs_tol=tol**2
            )
            t_bar_err = 0.5 * t_bar_tol

            upd_bar_data(
                tag=tag,
                t_bar=t_bar,
                l_bar=l_bar,
                z=z,
                v_bar=v_bar,
                eta=eta,
                tau=tau,
                t_bar_err=t_bar_err,
                l_bar_err=l_bar_err,
                z_err=z_err,
                v_bar_err=v_bar_err,
                eta_err=eta_err,
                tau_err=tau_err,
            )

        for i, peak in enumerate([POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH, POINT_PEAK_STAG]):
            find_peak(lambda x: g(x, peak), peak)

        """
        populate data for output purposes
        """
        if dom == DOMAIN_TIME:
            # fmt:off
            (z_j, l_bar_j, v_bar_j, t_bar_j, eta_j, tau_j) = (
                z_0, 0, 0, 0, 0, 1
            )
            # fmt: on
            for j in range(step):
                t_bar_k = t_bar_e / (step + 1) * (j + 1)
                (_, (z_j, l_bar_j, v_bar_j, eta_j, tau_j), (z_err, l_bar_err, v_bar_err, eta_err, tau_err)) = rkf78(
                    self.ode_t, (z_j, l_bar_j, v_bar_j, eta_j, tau_j), t_bar_j, t_bar_k, rel_tol=tol, abs_tol=tol**2
                )
                t_bar_j = t_bar_k

                upd_bar_data(
                    tag=COMPUTE,
                    t_bar=t_bar_j,
                    l_bar=l_bar_j,
                    z=z_j,
                    v_bar=v_bar_j,
                    eta=eta_j,
                    tau=tau_j,
                    t_bar_err=0,
                    l_bar_err=l_bar_err,
                    z_err=z_err,
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
            z_j, l_bar_j, v_bar_j, eta_j, tau_j = rkf78(
                self.ode_t, (z_0, 0, 0, 0, 1), 0, t_bar_j, rel_tol=tol, abs_tol=tol**2
            )[1]

            for j in range(step):
                l_bar_k = l_g_bar / (step + 1) * (j + 1)

                (_, (t_bar_j, z_j, v_bar_j, eta_j, tau_j), (t_bar_err, z_err, v_bar_err, eta_err, tau_err)) = rkf78(
                    self.ode_l, (t_bar_j, z_j, v_bar_j, eta_j, tau_j), l_bar_j, l_bar_k, rel_tol=tol, abs_tol=tol**2
                )

                l_bar_j = l_bar_k

                upd_bar_data(
                    tag=COMPUTE,
                    t_bar=t_bar_j,
                    l_bar=l_bar_j,
                    z=z_j,
                    v_bar=v_bar_j,
                    eta=eta_j,
                    tau=tau_j,
                    t_bar_err=t_bar_err,
                    l_bar_err=0,
                    z_err=z_err,
                    v_bar_err=v_bar_err,
                    eta_err=eta_err,
                    tau_err=tau_err,
                )

        logger.info(f"sampled for {step} points.")

        """
        sort the data points
        """

        data, error, p_trace = [], [], []
        l_c = self.l_c
        # l_g = self.l_g

        trace_step = max(step, 1)

        for bar_dataLine, bar_errorLine in zip(bar_data, bar_err):
            dtag, t_bar, l_bar, z, v_bar, p_bar, eta, tau = bar_dataLine
            (etag, t_bar_err, l_bar_err, z_err, v_bar_err, p_bar_err, eta_err, tau_err) = bar_errorLine

            t, t_err = t_bar * t_scale, t_bar_err * t_scale
            l, l_err = l_bar * self.l_0, l_bar_err * self.l_0
            psi, psi_err = self.f_psi_z(z), abs(self.f_sigma_z(z) * z_err)
            v, v_err = v_bar * self.v_j, v_bar_err * self.v_j

            p, p_err = p_bar * p_scale, (p_bar_err if isinstance(p_bar_err, str) else p_bar_err * p_scale)
            temp = tau * self.T_v if self.T_v else None
            temp_err = tau_err * self.T_v if self.T_v else None

            ps, p0, pb, vb = self.to_ps_p0_pb_vb(l, v, p, tau, eta)

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
                )
            )
            error.append(
                RecoillessErrorEntry(
                    tag=etag,
                    time=t_err,
                    travel=l_err,
                    burnup=psi_err,
                    velocity=v_err,
                    avg_pressure=p_err,
                    temperature=temp_err,
                    outflow_fraction=eta_err,
                )
            )

        for t_bar, (l_bar, psi, v_bar, p_bar, eta, tau) in record:
            t = t_bar * t_scale
            if t in [tableEntry.time for tableEntry in data]:
                continue
            l, v, p, temp = (l_bar * self.l_0, v_bar * self.v_j, p_bar * p_scale, tau * self.T_v if self.T_v else None)

            ps, p0, pb, vb = self.to_ps_p0_pb_vb(l, v, p, tau, eta)

            p_line = []
            for i in range(trace_step):  # 0....step -1
                x = i / trace_step * (l + l_c)
                px = self.to_px(l, v, vb, ps, eta, x)
                pp = PressureProbePoint(x, px)
                p_line.append(pp)

            p_line.append(PressureProbePoint(l + l_c, ps))
            p_trace.append(PressureTraceEntry(SAMPLE, temp, p_line))

            data.append(
                RecoillessTableEntry(
                    tag=SAMPLE,
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
                )
            )
            error.append(RecoillessErrorEntry("L"))

        data, error, p_trace = zip(
            *((a, b, c) for a, b, c in sorted(zip(data, error, p_trace), key=lambda entries: entries[0].time))
        )

        recoilless_result = RecoillessResult(self, data, error, p_trace)

        if self.material is None:
            logger.warning("material is not specified, skipping structural calculation.")
        else:

            try:
                self.get_structural(recoilless_result, step, tol)

            except Exception:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                logger.error("exception occured during structural calculation:")
                logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
                logger.info("structural calculation skipped.")

        return recoilless_result

    def to_ps_p0_pb_vb(self, l, v, p, tau, eta):
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

        Ps  : Shot base pressure
        Pb  : Breech pressure, chamber pressure at rearward of chamber
        P0  : Stagnation point pressure
        vb  : Rearward flow velocity at the rear of chamber
        """
        y = self.w * eta
        m_dot = self.c_a * self.v_j * self.S_j * p / (self.f * tau**0.5)
        # mass flow rate, rearward
        sb = self.s * self.chi_k
        vb = m_dot * (self.V_0 + self.s * l) / (sb * (self.w - y))
        # velocity impinging upon the rear of the breech before nozzle constriction

        h_1 = vb / v if v != 0 else inf
        h_2 = 2 * self.phi_1 * self.m / (self.w - y) + 1
        h = min(h_1, h_2)

        ps = p / (1 + (self.w - y) / (3 * self.phi_1 * self.m) * (1 - 0.5 * h))  # shot base pressure
        p0 = ps * (1 + (self.w - y) / (2 * self.phi_1 * self.m) * (1 + h) ** -1)  # stagnation point pressure
        pb = ps * (1 + (self.w - y) / (2 * self.phi_1 * self.m) * (1 - h)) if h == h_1 else 0  # breech pressure
        # l0 = H / (1 + H) * l  # location of the stagnation point
        return ps, p0, pb, vb

    def to_px(self, l, v, vb, ps, eta, x):
        m = self.m
        w = self.w
        phi_1 = self.phi_1
        y = self.w * eta

        """
        convert x, the physical displacement from breech bottom, to
        effective length in the equivalent gun.
        """

        l_1 = l
        l_0 = self.l_0 / self.chi_k  # physical length of the chamber.

        a_1 = self.s
        a_0 = a_1 * self.chi_k

        if x < l_0:
            z = x * a_0 / (l_0 * a_0 + l_1 * a_1)
        else:
            z = (l_0 * a_0 + (x - l_0) * a_1) / (l_0 * a_0 + l_1 * a_1)

        h_1 = vb / v if v != 0 else inf
        h_2 = 2 * phi_1 * self.m / (self.w - y) + 1
        h = min(h_1, h_2)

        px = ps * (1 + (w - y) / (2 * phi_1 * m) * (1 + h) * (1 - z**2) - (w - y) / (phi_1 * m) * h * (1 - z))
        return px

    @staticmethod
    def get_cf(gamma, sr, tol=1e-5):
        """
        takes the adiabatic index and area ration between constriction throat
        and the exit, calculate the thrust factor Cf
        See Hunt (1953) appendix I.A01-A03
        Sr = Se/St
        Vr = V/Vt
        """
        vr_old = 0
        vr = 1
        while abs(vr - vr_old) / vr > tol:
            vr_old = vr
            vr = ((gamma + 1) / (gamma - 1) * (1 - 2 / (gamma + 1) * (vr * sr) ** (1 - gamma))) ** 0.5

        cf = (2 / (gamma + 1)) ** (gamma / (gamma - 1)) * (gamma * vr + sr ** (1 - gamma) * vr ** (-gamma))

        return cf

    def get_structural(self, recoilless_result: RecoillessResult, step: int, tol: float):
        if not self.material:
            raise ValueError("Material must be supplied for structural calculation.")

        step = max(step, 1)
        l_c = self.l_c
        l_g = self.l_g
        chi_k = self.chi_k
        sigma = self.material.yield_strength
        s = self.s
        gamma = self.theta + 1

        s_j_bar = self.s_j_bar
        a_bar = self.a_bar
        r_s = 0.5 * self.caliber  # radius of the shot.
        r_b = r_s * chi_k**0.5  # chamber/breech radius
        r_t = r_b * s_j_bar**0.5  # throat radius
        r_n = r_t * a_bar**0.5  # nozzle exit radius

        x_probes = (
            [i / step * l_c for i in range(step)]
            + [l_c * (1 - tol)]
            + [i / step * l_g + l_c for i in range(step)]
            + [l_g + l_c]
        )

        p_probes = [0] * len(x_probes)
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
            p_probes[i] *= self.ssf

        # """
        # subscript k denote the end of the chamber
        # subscript b denote the throat of nozzle
        # subscript a denote the nozzle base.
        #
        # beta is the half angle of the constriction section
        # alpha is the half angle of the expansion section
        # """
        # beta = 45
        # l_b = (r_b - r_t) / tan(beta * pi / 180)
        # alpha = 15
        # l_a = (r_n - r_t) / tan(alpha * pi / 180)
        #
        # p_entry = p_probes[0]
        #
        # neg_x_probes = (
        #     [-l_b - l_a + l_a * i / step for i in range(step)]
        #     + [-l_b - tol * l_a]
        #     + [-l_b + l_b * i / step for i in range(step)]
        #     + [0]
        # )
        # neg_p_probes = []
        # neg_s_probes = []
        #
        # r_probes = []
        #
        # for x in neg_x_probes:
        #     if x < -l_b:
        #         k = (x + l_b + l_a) / l_a
        #         r = k * r_t + (1 - k) * r_n
        #         p = Recoilless.get_pr(gamma, (r / r_t) ** 2, tol)[1] * p_entry
        #
        #     else:
        #         k = (x + l_b) / l_b
        #         r = k * r_b + (1 - k) * r_t
        #         p = Recoilless.get_pr(gamma, (r / r_t) ** 2, tol)[0] * p_entry
        #
        #     if p > 3**-0.5 * sigma:
        #         raise ValueError(
        #             f"Limit to conventional construction ({sigma * 1e-6:.3f} MPa)" + " exceeded in nozzle."
        #         )
        #
        #     r_probes.append(r)
        #     neg_p_probes.append(p)
        #     neg_s_probes.append(pi * r**2)
        #
        # if self.is_af:
        #     v_n, rho_n = Gun.vrho_k(neg_x_probes, neg_p_probes, neg_s_probes, sigma, tol)
        #
        # else:
        #     v_n, rho_n = 0, []
        #
        #     for p in neg_p_probes:
        #         y = p / sigma
        #         rho = ((1 + y * (4 - 3 * y**2) ** 0.5) / (1 - 3 * y**2)) ** 0.5
        #         rho_n.append(rho)
        #
        #     for i in range(len(neg_x_probes) - 1):
        #         x_0, x_1 = neg_x_probes[i], neg_x_probes[i + 1]
        #         h = x_1 - x_0
        #         s_0, s_1 = neg_s_probes[i], neg_s_probes[i + 1]
        #         rho_0, rho_1 = rho_n[i], rho_n[i + 1]
        #         dv = 0.5 * ((rho_0**2 - 1) * s_0 + (rho_1**2 - 1) * s_1) * h
        #         v_n += dv

        hull = []
        #
        # for x, r, rho in zip(neg_x_probes, r_probes, rho_n):
        #     hull.append(OutlineEntry(x, r, r * rho))

        # nozzle_mass = v_n * self.material.rho

        rho_probes = []
        v = 0
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
            # v_c, rho_c = Gun.vrho_k(x_c, p_c, [s * chi_k for _ in x_c], sigma, tol, k_min=rho_n[-1])  # c for chamber

            v_c, rho_c = Gun.vrho_k(x_c, p_c, [s * chi_k for _ in x_c], sigma, tol)
            v_b, rho_b = Gun.vrho_k(
                x_b, p_b, [s for _ in x_b], sigma, tol, k_max=rho_c[-1] * chi_k**0.5, p_ref=max(p_c)
            )  # b for bore
            v = v_c + v_b
            rho_probes = rho_c + rho_b

        else:
            """
            The yield criterion chosen here is the fourth strength
            theory (von Mises) as it is generally accepted to be the most
            applicable for this application.

            The limiting stress points circumferentially along the circumference of the barrel.

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
                dv = (rho_1**2 + rho_0**2 - 2) * 0.5 * s * (x_1 - x_0)
                if x_1 < l_c:
                    v += dv * chi_k
                else:
                    v += dv

        tube_mass = v * self.material.rho

        for x, rho in zip(x_probes, rho_probes):
            if x < l_c:
                hull.append(OutlineEntry(x, r_b, rho * r_b))
            else:
                hull.append(OutlineEntry(x, r_s, rho * r_s))

        recoilless_result.outline = hull
        recoilless_result.tubeMass = tube_mass  # + nozzle_mass

        logger.info("conducted structural calculation.")

    @staticmethod
    def get_ar(gamma, pr):
        if pr == 0 or pr == 1:
            return inf
        return (
            (0.5 * (gamma + 1)) ** (1 / (gamma - 1))
            * pr ** (1 / gamma)
            * ((gamma + 1) / (gamma - 1) * (1 - pr ** ((gamma - 1) / gamma))) ** 0.5
        ) ** -1

    @staticmethod
    def get_pr(gamma, ar, tol):
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
        pr_c = (2 / (gamma + 1)) ** (gamma / (gamma - 1))

        if ar == 1:
            return pr_c, pr_c
        else:
            pr_sup, _ = dekker(lambda pr: Recoilless.get_ar(gamma, pr), 0, pr_c, y=ar, y_rel_tol=tol)
            pr_sub, _ = dekker(lambda pr: Recoilless.get_ar(gamma, pr), pr_c, 1, y=ar, y_rel_tol=tol)

            return pr_sub, pr_sup


if __name__ == "__main__":
    """standard 7 hole cylinder has d_0=e_1, hole dia = 0.5 * arc width
    d_0 = 4*2e_1+3*d_0 = 11 * e_1
    """

    compositions = Composition.read_file("data/propellants.csv")

    M17 = compositions["M17"]
    M1 = compositions["M1"]
    from prop import SimpleGeometry

    M17C = Propellant(M17, SimpleGeometry.CYLINDER, 0.0, 2.5)
    M1C = Propellant(M1, SimpleGeometry.CYLINDER, 0.0, 10)
    lf = 0.6
    print("DELTA/rho:", lf)
    test = Recoilless(
        caliber=0.082,
        shot_mass=2,
        propellant=M1C,
        web=4e-4,
        charge_mass=0.3,
        chamber_volume=0.3 / M1C.rho_p / lf,
        start_pressure=10e6,
        length_gun=3.5,
        nozzle_expansion=2.0,
        chambrage=1.0,
    )
    record = []
