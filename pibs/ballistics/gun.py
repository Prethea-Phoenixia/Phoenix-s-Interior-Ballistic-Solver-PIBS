from __future__ import annotations

import logging
import sys
import traceback
from dataclasses import dataclass
from math import exp, inf, log, pi

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
    POINT_START,
    SAMPLE,
    SOL_LAGRANGE,
    SOL_MAMONTOV,
    SOL_PIDDUCK,
    Domains,
    Solutions,
)
from .generics import (
    DelegatesPropellant,
    GenericEntry,
    GenericErrorEntry,
    GenericResult,
    OutlineEntry,
    PressureProbePoint,
    PressureTraceEntry,
)
from .material import Material
from .num import dekker, gss, intg, rkf78, secant
from .prop import Composition, Propellant

logger = logging.getLogger(__name__)


@dataclass
class GunResult(GenericResult):
    gun: Gun
    table_data: list[GunTableEntry]
    error_data: list[GunErrorEntry]


@dataclass
class GunTableEntry(GenericEntry):
    pass


@dataclass
class GunErrorEntry(GenericErrorEntry):
    pass


def pidduck(wpm, k, tol):
    """
    Pidduck's limiting solution to the Lagrange problem.
    wpm : w/(phi_1 * m), charge mass to equivalent corrected (fictitious) shot
          weight
    k   : adiabatic index of the gas, in practice this is not a great influence
    tol : numerical tolerance

    Pidduck's solution is reduced to that of M.A.Mamontov's solution at k -> 1,
    however numerical difficulty necessitates taking the limit.
    """
    if k < 1:
        raise ValueError("Invalid adiabatic index passed", k)

    def f(om, x):
        if k == 1:
            return exp(-om * x**2)
        else:
            return (1 - om * x**2) ** (1 / (k - 1))

    def g(om, x):
        return f(om, x) * x**2

    def f_omega(om):
        if om == 0:
            return -inf

        i, _ = intg(lambda x: f(om, x), 0, 1, tol)

        if k == 1:
            return i - 0.5 * wpm * exp(-om) / om
        else:
            return i - 0.5 * ((k - 1) / k) * wpm * ((1 - om) ** (k / (k - 1)) / om)

    a, b = dekker(f_omega, 0, 1, tol, y_abs_tol=tol)
    omega = 0.5 * (a + b)

    if k == 1:
        labda_1 = (exp(omega) - 1) / wpm
    else:
        labda_1 = ((1 - omega) ** (k / (1 - k)) - 1) / wpm

    # Pidduck's solution
    i_u, _ = intg(lambda x: g(omega, x), 0, 1, tol)
    i_l, _ = intg(lambda x: f(omega, x), 0, 1, tol)
    labda_2 = i_u / i_l

    return labda_1, labda_2


class Gun(DelegatesPropellant):
    def __init__(
        self,
        caliber,
        shot_mass: float,
        propellant: Propellant,
        web: float,  # 2e_1
        charge_mass: float,
        chamber_volume: float,
        start_pressure: float,
        length_gun: float,
        chambrage: float,
        structural_material: Material | None = None,
        structural_safety_factor: float = 1.1,
        drag_coefficient: float = 0.0,
        autofrettage: bool = True,
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
                chambrage < 1,
                drag_coefficient < 0,
                drag_coefficient >= 1,
                structural_safety_factor <= 1,
            )
        ):
            raise ValueError("Invalid gun parameters")

        # if chargeMass > (propellant.maxLF * propellant.rho_p * chamberVolume):
        #     raise ValueError("Specified Load Fraction Violates Geometrical Constraint")

        self.caliber = caliber

        self.e_1 = 0.5 * web
        self.s = (0.5 * caliber) ** 2 * pi
        self.m = shot_mass
        self.omega = charge_mass
        self.vol_0 = chamber_volume
        self.p_0 = start_pressure
        self.l_g = length_gun
        self.chi_k = chambrage  # ration of l_0 / l_chamber
        self.l_0 = self.vol_0 / self.s
        self.l_c = self.l_0 / self.chi_k
        self.delta = self.omega / self.vol_0

        self.phi_1 = 1 / (1 - drag_coefficient)  # drag work coefficient
        self.material = structural_material
        self.ssf = structural_safety_factor
        self.is_af = autofrettage

        if self.p_0 == 0:
            raise NotImplementedError(
                "Current implementation use exponential burn rate and does not"
                + " allow for solving the case with 0 shot start pressure."
            )

        self.p_a_bar, self.c_a_bar, self.labda_1, self.labda_2 = 0.0, 0.0, 0.0, 0.0
        self.v_j, self.b, self.phi, self.k_1 = 0.0, 0.0, 0.0, 0.0
        self.z_0, self.psi_0 = 0, 0

    def f_p_bar(self, z, l_bar, v_bar) -> float:
        psi = self.f_psi_z(z)

        l_psi_bar = 1 - self.delta * ((1 - psi) / self.rho_p + (self.alpha * psi))

        p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

        return p_bar

    def ode_t(self, _, z, l_bar, v_bar, __) -> tuple[float, float, float]:
        psi = self.f_psi_z(z)

        l_psi_bar = 1 - self.delta / self.rho_p - self.delta * (self.alpha - 1 / self.rho_p) * psi
        p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

        if self.c_a_bar != 0 and v_bar > 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (
                0.25 * k * (k + 1) * v_r**2 + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a_bar
        else:
            p_d_bar = 0

        dz = (0.5 * self.theta / self.b) ** 0.5 * p_bar**self.n

        dl_bar = v_bar
        dv_bar = self.theta * 0.5 * (p_bar - p_d_bar)

        return dz, dl_bar, dv_bar

    def ode_l(self, l_bar, _, z, v_bar, __) -> tuple[float, float, float]:
        """length domain ode of internal ballistics
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition."""

        psi = self.f_psi_z(z)

        l_psi_bar = 1 - self.delta / self.rho_p - self.delta * (self.alpha - 1 / self.rho_p) * psi
        p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

        if self.c_a_bar != 0 and v_bar > 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (0.25 * k * (k + 1) * v_r**2 + k * v_r * (1 + (0.25 * (k + 1) * v_r) ** 2) ** 0.5) * self.p_a_bar
        else:
            p_d_bar = 0

        dz = (0.5 * self.theta / self.b) ** 0.5 * p_bar**self.n / v_bar

        dv_bar = self.theta * 0.5 * (p_bar - p_d_bar) / v_bar  # dv_bar/dl_bar
        dt_bar = 1 / v_bar  # dt_bar / dl_bar

        return dt_bar, dz, dv_bar

    def ode_z(self, z, t_bar, l_bar, v_bar, _):
        psi = self.f_psi_z(z)

        l_psi_bar = 1 - self.delta / self.rho_p - self.delta * (self.alpha - 1 / self.rho_p) * psi
        p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

        if self.c_a_bar != 0 and v_bar > 0:
            k = self.k_1  # gamma
            v_r = v_bar / self.c_a_bar
            p_d_bar = (
                0.25 * k * (k + 1) * v_r**2 + k * v_r * (1 + (0.25 * (k + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a_bar

        else:
            p_d_bar = 0

        if z <= self.z_b:
            dt_bar = (2 * self.b / self.theta) ** 0.5 * p_bar**-self.n  # dt_bar/dZ
            dl_bar = v_bar * dt_bar  # dv_bar/dZ
            dv_bar = 0.5 * self.theta * (p_bar - p_d_bar) * dt_bar
        else:
            # technically speaking it is undefined in this area
            dt_bar = 0  # dt_bar/dZ
            dl_bar = 0  # dl_bar/dZ
            dv_bar = 0  # dv_bar/dZ

        return dt_bar, dl_bar, dv_bar

    def get_temperature(self, psi, l, p):
        """
        given pressure and travel, return temperature
        using the Nobel-Abel EOS
        """
        if self.T_v:
            r = self.f / self.T_v
            l_psi = self.l_0 * (1 - self.delta / self.rho_p - self.delta * (self.alpha * -1 / self.rho_p) * psi)
            return self.s * p * (l + l_psi) / (self.omega * psi * r)
        else:
            return None

    def dp_dz(self, z, l_bar, v_bar):
        psi = self.f_psi_z(z)
        p_bar = self.f_p_bar(z, l_bar, v_bar)

        dz = (0.5 * self.theta / self.b) ** 0.5 * p_bar**self.n

        l_psi_bar = 1 - self.delta * ((1 - psi) / self.rho_p + (self.alpha * psi))
        dp_bar = (
            (
                (1 + p_bar * self.delta * (self.alpha - 1 / self.rho_p))
                * self.f_sigma_z(z)  # dpsi/dz
                * dz  # dz/dt_bar
                - p_bar * v_bar * (1 + self.theta)
            )
            / (l_bar + l_psi_bar)  # dp_bar/dt_bar
            / dz
        )  # dp_bar/dz

        return dp_bar

    def integrate(
        self,
        step: int = 10,
        tol: float = 1e-5,
        dom: Domains = DOMAIN_TIME,
        sol: Solutions = SOL_PIDDUCK,
        ambient_rho: float = 1.204,
        ambient_p: float = 101.325e3,
        ambient_gamma: float = 1.4,
        **_,
    ) -> GunResult:
        """
        Runs a full numerical solution for the gun in the specified domain
        sampled evenly at specified number of step, using a scaled numerical
        tolerance as specified.

        tolerance is meant to be interpreted as the maximum relative deviation
        each component is allowed to have, at each step of integration point.

        Through significant trials and errors, it was determined that for this
        particular system, the error due to compounding does not appear to be
        significant, usually on the order of 1e-16 - 1e-14 as compared to much
        larger for component errors.

        Partial analytical approximation to the Lagrangian problem:
        d rho / dx = 0: Lagrange
          d S / dx = 0: Pidduck
          d T / dx = 0: M.A.Mamontov
        All solutions assumes gas velocity increasing linearlly from 0
        at breech face and shot velocity at shot base.
        """

        self.psi_0 = (1 / self.delta - 1 / self.rho_p) / (self.f / self.p_0 + self.alpha - 1 / self.rho_p)
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

        if any((ambient_p < 0, ambient_rho < 0, ambient_gamma < 1)):
            raise ValueError("Invalid ambient condition")

        if sol == SOL_LAGRANGE:
            labda_1, labda_2 = 0.5, 1 / 3
        elif sol == SOL_PIDDUCK:
            labda_1, labda_2 = pidduck(self.omega / (self.phi_1 * self.m), self.theta + 1, tol)
        elif sol == SOL_MAMONTOV:
            labda_1, labda_2 = pidduck(self.omega / (self.phi_1 * self.m), 1, tol)
        else:
            raise ValueError("Unknown Solution")

        self.labda_1, self.labda_2 = labda_1, labda_2

        labda = self.l_g / self.l_0
        cc = 1 - (1 - 1 / self.chi_k) * log(labda + 1) / labda  # chambrage correction factor

        self.phi = self.phi_1 + labda_2 * self.omega / self.m * cc
        """
        见《枪炮内弹道学》（金，2014）p.70 式
        """

        self.b = (
            self.s**2
            * self.e_1**2
            / (self.f * self.phi * self.omega * self.m * self.u_1**2)
            * (self.f * self.delta) ** (2 * (1 - self.n))
        )

        self.v_j = (2 * self.f * self.omega / (self.theta * self.phi * self.m)) ** 0.5

        t_scale = self.l_0 / self.v_j
        p_scale = self.f * self.delta

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

        def upd_bar_data(tag, t_bar, l_bar, z, v_bar, t_bar_err, l_bar_err, z_err, v_bar_err):
            p_bar = self.f_p_bar(z, l_bar, v_bar)
            bar_data.append((tag, t_bar, l_bar, z, v_bar, p_bar))
            p_bar_err = abs(self.dp_dz(z, l_bar, v_bar) * z_err)
            bar_err.append(("L", t_bar_err, l_bar_err, z_err, v_bar_err, p_bar_err))

        upd_bar_data(tag=POINT_START, t_bar=0, l_bar=0, z=z_0, v_bar=0, t_bar_err=0, l_bar_err=0, z_err=0, v_bar_err=0)

        record.append((0, (0, self.psi_0, 0, p_bar_0)))

        z_i = z_0
        z_j = z_b
        n = 1
        delta_z = z_b - z_0
        t_bar_i, l_bar_i, v_bar_i = 0, 0, 0
        is_burn_out_contained = True

        """
        Instead of letting the integrator handle the heavy lifting, we
        partition Z and integrate upwards until either barrel exit or
        burnout has been achieved. This seek-and-validate inspired
        from bisection puts the group of value denoted by subscript
        i within or on the muzzle, with the propellant either still
        burning or right on the burnout point..
        """
        ztlv_record = [(z_0, (0, 0, 0))]
        p_max = 1e9  # 1GPa
        p_bar_max = p_max / p_scale

        def abort(x, ys, record):
            z = x
            t_bar, l_bar, v_bar = ys
            p_bar = self.f_p_bar(z, l_bar, v_bar)

            return l_bar > l_g_bar or p_bar > p_bar_max or v_bar < 0

        while z_i < z_b:  # terminates if burnout is achieved
            ztlv_record_i = []
            if z_j == z_i:
                raise ValueError("Numerical accuracy exhausted in search of exit/burnout point.")
            try:
                if z_j > z_b:
                    z_j = z_b
                z, (t_bar_j, l_bar_j, v_bar_j), _ = rkf78(
                    self.ode_z,
                    (t_bar_i, l_bar_i, v_bar_i),
                    z_i,
                    z_j,
                    rel_tol=tol,
                    abs_tol=tol**2,
                    abort_func=abort,
                    record=ztlv_record_i,
                )

                p_bar_j = self.f_p_bar(z_j, l_bar_j, v_bar_j)

            except ValueError as e:
                raise e

            if l_bar_j >= l_g_bar:
                if abs(l_bar_i - l_g_bar) / l_g_bar > tol or l_bar_i == 0:
                    n *= 2
                    z_j = z_i + delta_z / n
                else:
                    is_burn_out_contained = False
                    break  # l_bar_i is solved to within a tol of l_bar_g

            else:
                ztlv_record.extend(ztlv_record_i)
                if v_bar_j <= 0:

                    z, (t_bar, l_bar, v_bar) = ztlv_record[-1]

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
                    )

                if p_bar_j > p_bar_max:
                    raise ValueError(
                        "Nobel-Abel EoS is generally accurate enough below"
                        + " 600MPa. However, Unreasonably high pressure "
                        + "(>{:.0f} MPa) was encountered.".format(p_max / 1e6),
                    )  # in practice most press-related spikes are captured here

                t_bar_i, l_bar_i, v_bar_i = t_bar_j, l_bar_j, v_bar_j

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
        we wrote the SoE to be be piecewise continous from (0, z_b] and (z_b,
        +inf) it is necessary to do this to prevent the RKF integrator coming
        up with irreducible error estimates and driving the step size to 0
        around Z = z_b
        """
        if is_burn_out_contained:
            z_i = z_b + tol

        record.extend(
            (t_bar, (l_bar, self.f_psi_z(z), v_bar, self.f_p_bar(z, l_bar, v_bar)))
            for (z, (t_bar, l_bar, v_bar)) in ztlv_record
            if t_bar != 0
        )

        if is_burn_out_contained:
            logger.info("integrated to burnout point.")
        else:
            logger.warning("shot exited barrel before burnout.")

        """
        Subscript e indicate exit condition.
        At this point, since its guaranteed that point i will be further towards the chamber of the firearm than 
        point e, we do not have to worry about the dependency of the correct behaviour of this ODE
        on the positive direction of integration.
        """

        ltzv_record = []
        l_bar, (t_bar_e, z_e, v_bar_e), (t_bar_err, z_err, v_bar_err) = rkf78(
            self.ode_l, (t_bar_i, z_i, v_bar_i), l_bar_i, l_g_bar, rel_tol=tol, abs_tol=tol**2, record=ltzv_record
        )

        if l_bar != l_g_bar:
            if v_bar_e <= 0:
                raise ValueError(
                    "Squib load condition detected post burnout:"
                    + " Round stopped in bore at {:.0f} mm".format(l_bar * self.l_0 * 1e3)
                )

        record.extend(
            (t_bar, (l_bar, self.f_psi_z(z), v_bar, self.f_p_bar(z, l_bar, v_bar)))
            for (l_bar, (t_bar, z, v_bar)) in ltzv_record
        )

        upd_bar_data(
            tag=POINT_EXIT,
            t_bar=t_bar_e,
            l_bar=l_g_bar,
            z=z_e,
            v_bar=v_bar_e,
            t_bar_err=t_bar_err,
            l_bar_err=0,
            z_err=z_err,
            v_bar_err=v_bar_err,
        )

        t_bar_f = None
        if z_b > 1.0 and z_e >= 1.0:  # fracture point exist and is contained
            """
            Subscript f indicate fracture condition
            ODE w.r.t Z is integrated from Z_0 to 1, from onset of projectile movement to charge fracture
            """
            (_, (t_bar_f, l_bar_f, v_bar_f), (t_bar_err_f, l_bar_err_f, v_bar_err_f)) = rkf78(
                self.ode_z, (0, 0, 0), z_0, 1, rel_tol=tol, abs_tol=tol**2
            )

            upd_bar_data(
                tag=POINT_FRACTURE,
                t_bar=t_bar_f,
                l_bar=l_bar_f,
                z=1,
                v_bar=v_bar_f,
                t_bar_err=t_bar_err_f,
                l_bar_err=l_bar_err_f,
                z_err=0,
                v_bar_err=v_bar_err_f,
            )

        t_bar_b = None
        if is_burn_out_contained:
            """
            Subscript b indicated burnout condition
            ODE w.r.t Z is integrated from Z_0 to z_b, from onset of projectile movement to charge burnout.
            """

            (
                _,
                (t_bar_b, l_bar_b, v_bar_b),
                (t_bar_err_b, l_bar_err_b, v_bar_err_b),
            ) = rkf78(self.ode_z, (0, 0, 0), z_0, z_b, rel_tol=tol, abs_tol=tol**2)

            upd_bar_data(
                tag=POINT_BURNOUT,
                t_bar=t_bar_b,
                l_bar=l_bar_b,
                z=z_b,
                v_bar=v_bar_b,
                t_bar_err=t_bar_err_b,
                l_bar_err=l_bar_err_b,
                z_err=0,
                v_bar_err=v_bar_err_b,
            )

        """
        Subscript p indicate peak pressure
        In theory the time-domain equations should be flatter around the peak pressure point. As well as, not having 
        starting issues.

        we hereby simply solve p golden section searching it from origin to point e, i.e. inside the barrel.
        """

        def find_peak(f, tag):
            """
            tolerance is specified a bit differently for gold section search
            GSS tol is the length between the upper bound and lower bound
            of the maxima/minima, thus by our definition of tolerance (one-sided),
            we take the median value.
            """
            t_bar_tol = tol * min(t for t in (t_bar_e, t_bar_b, t_bar_f) if t is not None)

            t_bar_1, t_bar_2 = gss(f, 0, t_bar_e if t_bar_b is None else t_bar_b, x_tol=t_bar_tol, find_min=False)

            t_bar = 0.5 * (t_bar_1 + t_bar_2)

            _, (z, l_bar, v_bar), (z_err, l_bar_err, v_bar_err) = rkf78(
                self.ode_t, (z_0, 0, 0), 0, t_bar, rel_tol=tol, abs_tol=tol**2
            )
            t_bar_err = 0.5 * t_bar_tol

            upd_bar_data(
                tag=tag,
                t_bar=t_bar,
                l_bar=l_bar,
                z=z,
                v_bar=v_bar,
                t_bar_err=t_bar_err,
                l_bar_err=l_bar_err,
                z_err=z_err,
                v_bar_err=v_bar_err,
            )

        def g(t, tag) -> float:
            z, l_bar, v_bar = rkf78(self.ode_t, (z_0, 0, 0), 0, t, rel_tol=tol, abs_tol=tol**2)[1]

            p_bar = self.f_p_bar(z, l_bar, v_bar)
            if tag == POINT_PEAK_AVG:
                return p_bar
            else:
                ps, pb = self.to_ps_pb(l_bar * self.l_0, p_bar * p_scale)
                if tag == POINT_PEAK_SHOT:
                    return ps / p_scale
                elif tag == POINT_PEAK_BREECH:
                    return pb / p_scale
                else:
                    raise ValueError("tag not handled.")

        for i, peak in enumerate([POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH]):
            find_peak(lambda x: g(x, peak), peak)

        """
        populate data for output purposes
        """

        if dom == DOMAIN_TIME:
            (z_j, l_bar_j, v_bar_j, t_bar_j) = (z_0, 0, 0, 0)
            for j in range(step):
                t_bar_k = t_bar_e / (step + 1) * (j + 1)
                (_, (z_j, l_bar_j, v_bar_j), (z_err, l_bar_err, v_bar_err)) = rkf78(
                    self.ode_t, (z_j, l_bar_j, v_bar_j), t_bar_j, t_bar_k, rel_tol=tol, abs_tol=tol**2
                )
                t_bar_j = t_bar_k

                upd_bar_data(
                    tag=SAMPLE,
                    t_bar=t_bar_j,
                    l_bar=l_bar_j,
                    z=z_j,
                    v_bar=v_bar_j,
                    t_bar_err=0,
                    l_bar_err=l_bar_err,
                    z_err=z_err,
                    v_bar_err=v_bar_err,
                )

        elif dom == DOMAIN_LEN:
            """
            Due to two issues, i.e. 1.the length domain ODE cannot be integrated from the origin point, and 2.the
            correct behaviour can only be expected when starting from a point with active burning else dZ flat lines.
            we do another Z domain integration to seed the initial values to a point where ongoing burning is guaranteed.
            (in the case of gun barrel length >= burn length, the group of value by subscript i will not guarantee
            burning is still ongoing).
            """
            t_bar_j = 0.5 * t_bar_i
            z_j, l_bar_j, v_bar_j = rkf78(self.ode_t, (z_0, 0, 0), 0, t_bar_j, rel_tol=tol, abs_tol=tol**2)[1]

            for j in range(step):
                l_bar_k = l_g_bar / (step + 1) * (j + 1)

                (_, (t_bar_j, z_j, v_bar_j), (t_bar_err, z_err, v_bar_err)) = rkf78(
                    self.ode_l, (t_bar_j, z_j, v_bar_j), l_bar_j, l_bar_k, rel_tol=tol, abs_tol=tol**2
                )
                l_bar_j = l_bar_k

                upd_bar_data(
                    tag=SAMPLE,
                    t_bar=t_bar_j,
                    l_bar=l_bar_j,
                    z=z_j,
                    v_bar=v_bar_j,
                    t_bar_err=t_bar_err,
                    l_bar_err=0,
                    z_err=z_err,
                    v_bar_err=v_bar_err,
                )

        logger.info(f"sampled for {step} points.")

        """
        sort the data points
        """

        data = []
        error = []
        p_trace = []
        trace_steps = max(step, 1)

        for bar_dataLine, bar_errorLine in zip(bar_data, bar_err):
            dtag, t_bar, l_bar, z, v_bar, p_bar = bar_dataLine

            _, t_bar_err, l_bar_err, z_err, v_bar_err, p_bar_err = bar_errorLine

            t, t_err = t_bar * t_scale, t_bar_err * t_scale
            l, l_err = l_bar * self.l_0, l_bar_err * self.l_0
            psi, psi_err = self.f_psi_z(z), abs(self.f_sigma_z(z) * z_err)
            v, v_err = v_bar * self.v_j, v_bar_err * self.v_j
            p, p_err = p_bar * p_scale, p_bar_err * p_scale
            ps, pb = self.to_ps_pb(l, p)
            temp = self.get_temperature(psi, l, p)

            p_line = []
            for i in range(trace_steps):
                x = i / trace_steps * (l + self.l_c)
                p_x, _ = self.to_px_u(l, ps, pb, v, x)
                pp = PressureProbePoint(x, p_x)
                p_line.append(pp)

            p_line.append(PressureProbePoint(l + self.l_c, ps))
            p_trace.append(PressureTraceEntry(dtag, temp, p_line))

            table_entry = GunTableEntry(
                tag=dtag,
                time=t,
                travel=l,
                burnup=psi,
                velocity=v,
                breech_pressure=pb,
                avg_pressure=p,
                shot_pressure=ps,
                temperature=temp,
            )
            error_entry = GunErrorEntry(
                tag="L",
                time=t_err,
                travel=l_err,
                burnup=psi_err,
                velocity=v_err,
                breech_pressure=None,
                avg_pressure=p_err,
                shot_pressure=None,
                temperature=None,
            )
            data.append(table_entry)
            error.append(error_entry)

        """
        scale the records too
        """

        for t_bar, (l_bar, psi, v_bar, p_bar) in record:
            t = t_bar * t_scale
            if t in [tableEntry.time for tableEntry in data]:
                continue
            l = l_bar * self.l_0
            t = t_bar * t_scale
            p = p_bar * p_scale
            ps, pb = self.to_ps_pb(l, p)
            v = v_bar * self.v_j
            temp = self.get_temperature(psi, l, p)

            p_line = []
            for i in range(trace_steps):
                x = i / trace_steps * (l + self.l_c)
                p_x, _ = self.to_px_u(l, ps, pb, v, x)
                pp = PressureProbePoint(x, p_x)
                p_line.append(pp)

            p_line.append(PressureProbePoint(l + self.l_c, ps))
            p_trace.append(PressureTraceEntry(COMPUTE, temp, p_line))

            table_entry = GunTableEntry(
                tag=COMPUTE,
                time=t,
                travel=l,
                burnup=psi,
                velocity=v,
                breech_pressure=pb,
                avg_pressure=p,
                shot_pressure=ps,
                temperature=temp,
            )
            error_entry = GunErrorEntry("L")
            data.append(table_entry)
            error.append(error_entry)

        data, error, p_trace = zip(
            *(
                (tableEntry, errorEntry, pTrace)
                for tableEntry, errorEntry, pTrace in sorted(
                    zip(data, error, p_trace), key=lambda entries: entries[0].time
                )
            )
        )

        # calculate a pressure and flow velocity tracing.
        gun_result = GunResult(self, data, error, p_trace)

        if self.material is None:
            logger.warning("material is not specified, skipping structural calculation.")
        else:

            try:
                self.get_structural(gun_result, step, tol)

            except Exception:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                logger.error("exception occurred during structural calculation:")
                logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
                logger.info("structural calculation skipped.")

        return gun_result

    def to_ps_pb(self, l: float, p: float) -> tuple[float, float]:
        """
        Convert average chamber pressure at certain travel to
        shot base pressure, and breech face pressure

        l: travel of the projectile
        p: average pressure

        Ps: pressure at shot
        Pb: pressure at breech
        """
        labda_g = l / self.l_0
        labda_1_prime = self.labda_1 * (1 / self.chi_k + labda_g) / (1 + labda_g)
        labda_2_prime = self.labda_2 * (1 / self.chi_k + labda_g) / (1 + labda_g)

        factor_s = 1 + labda_2_prime * (self.omega / (self.phi_1 * self.m))  # factor_b = P/P_b = phi / phi_1

        factor_b = (self.phi_1 * self.m + labda_2_prime * self.omega) / (
            self.phi_1 * self.m + labda_1_prime * self.omega
        )

        return p / factor_s, p / factor_b

    def to_px_u(self, l: float, p_s: float, p_b: float, v: float, x: float) -> tuple[float, float]:
        """
        Convert the average chamber to pressure and gas flow speed
        at arbitrary point x for projectile travel of l and average pressure
        of p, **assuming the Lagrangian distribution**.

        Note that with the current state of research, only characteristic point
        values are available for other distributions, use to_ps_pb() instead for that.

        l: projectile travel
        p_s: pressure of shot
        p_b: pressure of breech
        x: probe point, start from the breech face.
        """
        l_1 = l
        l_c = self.l_0 / self.chi_k  # physical length of the chamber.
        a_1 = self.s
        a_0 = a_1 * self.chi_k
        r = self.chi_k * x if x < l_c else (x - l_c) + self.l_0
        k = (r / (self.l_0 + l)) ** 2
        p_x = p_s * k + p_b * (1 - k)

        if x < l_c:
            u = a_1 * x * v / (self.vol_0 + a_1 * l_1)
        else:
            u = (a_1 * x + (a_0 - a_1) * l_c) * v / (self.vol_0 + a_1 * l_1)

        return p_x, u

    def get_structural(self, gun_result: GunResult, step: int, tol: float):
        if not self.material:
            raise ValueError("Material must be supplied for structural calculation.")

        step = max(step, 1)

        # step 1. calculate the barrel mass
        r = 0.5 * self.caliber
        l_c = self.l_c
        l_g = self.l_g
        chi_k = self.chi_k
        sigma = self.material.yield_strength
        s = self.s
        r_b = r * chi_k**0.5
        x_probes = (
            [i / step * l_c for i in range(step)]
            + [l_c * (1 - tol)]
            + [i / step * l_g + l_c for i in range(step)]
            + [l_g + l_c]
        )
        p_probes = [0] * len(x_probes)

        for gun_table_entry in gun_result.table_data:
            l = gun_table_entry.travel
            v = gun_table_entry.velocity
            p_s = gun_table_entry.shot_pressure
            p_b = gun_table_entry.breech_pressure
            for i, x in enumerate(x_probes):
                if (x - l_c) <= l:
                    p_x, _ = self.to_px_u(l, p_s, p_b, v, x)
                    p_probes[i] = max(p_probes[i], p_x)
                else:
                    break

        # strength requirement given structural safety factor.
        for i in range(len(p_probes)):
            p_probes[i] *= self.ssf

        rho_probes = []
        v = 0
        if self.is_af:
            """
            m : r_a / r_i
            k : r_o / r_i
            n : p_vM_max / sigma

            1 < m < k

            The point of optimum autofrettage, or the minimum autofrettage
            necessary to contain the working pressure, is
            """
            i = step + 1
            x_c, p_c = x_probes[:i], p_probes[:i]  # c for chamber
            x_b, p_b = x_probes[i:], p_probes[i:]  # b for barrel

            v_c, rho_c = Gun.vrho_k(x_c, p_c, [s * chi_k for _ in x_c], sigma, tol)
            v_b, rho_b = Gun.vrho_k(
                x_b, p_b, [s for _ in x_b], sigma, tol, p_ref=max(p_c), k_max=rho_c[-1] * chi_k**0.5
            )

            v = v_c + v_b
            rho_probes = rho_c + rho_b

        else:
            """
            The yield criterion chosen here is the fourth strength
            theory (von Mises) as it is generally accepted to be the most
            applicable for this application.

            The limiting stress points circumferentially along the circumference
            of the barrel.

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
                        f"Limit to conventional construction ({sigma * 3*1e-6:.3f} MPa)" + " exceeded in section."
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

        hull = []
        for x, rho in zip(x_probes, rho_probes):
            if x < l_c:
                hull.append(OutlineEntry(x, r_b, rho * r_b))
            else:
                hull.append(OutlineEntry(x, r, rho * r))

        gun_result.outline = hull
        gun_result.tubeMass = v * self.material.rho

        logger.info("conducted structural calculation.")

    @staticmethod
    def vrho_k(
        x_s: list[float],  # probe location
        p_s: list[float],  # probe pressure
        s_s: list[float],  # probe cross-section area
        sigma_yield: float,  # yield strength
        tol: float,  # tolerance
        k_max: float = None,  # maximum autofrettage
        k_min: float = None,  # minimum autofrettage
        p_ref=None,  # manually specify maximum pressure
    ):

        if p_ref is not None:
            p_max = p_ref
        else:
            p_max = max(p_s)

        # optimal autofrettage for this pressure
        m_opt = exp(p_max / sigma_yield * 3**0.5 * 0.5)

        def sigma_min(m):
            return (sigma_yield * (1 - (1 + 2 * log(m)) / m**2) + 2 * p_max / m**2) * (3**0.5 * 0.5)

        if sigma_min(m_opt) > sigma_yield:
            """if the minimum junction stress at the optimal autofrettage
            fraction cannot be achieved down to material yield even as
            the thickness goes to infinity, raise an error and abort
            calculation"""
            raise ValueError(
                "Plastic-elastic junction stress exceeds material "
                + f"yield ({sigma_yield * 1e-6:.3f} MPa) for autofrettaged construction."
            )

        elif sigma_min(1) > sigma_yield:
            """if the minimum junction stress at an autofrettage fraction
            of 1 exceeds material yield, implies a certain amount of
            autofrettaging is required to contain the pressure"""
            m_min, _ = dekker(sigma_min, 1, m_opt, y=sigma_yield, y_rel_tol=tol)

            """safety, fudge code to ensure a valid minimum autofrettage
            fraction is found.
            """
            while sigma_min(m_min) > sigma_yield:
                m_min *= 1 + tol**2

        else:
            m_min = 1 + tol**2

        def f(m: float):
            rhos = []
            v = 0
            for p in p_s:
                sigma_max = Gun.sigma_von_mises(m, p, m, sigma_yield)
                """
                the limit as k -> +inf for the stress is:

                lim sigma_tresca =
                 k-> +inf
                  sigma * [1 - (1 + 2 ln(m))/m**2 ] + 2p/m**2
                """
                if sigma_max < sigma_yield:
                    rho = m
                else:
                    rho, _ = secant(
                        lambda k: Gun.sigma_von_mises(k, p, m, sigma_yield),
                        m,
                        m * 2,
                        y=sigma_yield,
                        x_min=m,
                        y_rel_tol=tol,
                    )

                rhos.append(rho)

            for i in range(len(x_s) - 1):
                x_0, x_1 = x_s[i], x_s[i + 1]
                s_0, s_1 = s_s[i], s_s[i + 1]
                rho_0, rho_1 = rhos[i], rhos[i + 1]
                dv = 0.5 * ((rho_0**2 - 1) * s_0 + (rho_1**2 - 1) * s_1) * (x_1 - x_0)
                v += dv
            return v, rhos

        m_max = m_opt
        if k_max is not None:
            """another constraint is the continuation criteria, i.e. the
            barrel base should not require a thickness jump from the front of
            breech, within reason.

            A maximum ratio corresponds to a minimum autofrettage ratio.
            manually setting f(m_opt) to -1 since this is always true,
            but the limitations to numerical accuracy can cause the result
            to float around ~ +/- 10 * tol

            Since in all use cases k_max is set using a section of the chamber
            with higher pressure ratings, it is almost guaranteed that the
            corresponding m_min must exist. In the case this is not true, the
            resulting error raised will cause the program to gracefully default
            to finding the minimum auto-frettaged mass.
            """
            try:
                m_k, _ = dekker(lambda m: f(m)[1][0] if m != m_opt else -1, m_min, m_opt, y=k_max, y_rel_tol=tol)
                m_min = max(m_k, m_min)
            except ValueError:
                pass

        if k_min is not None:
            try:
                m_k, _ = dekker(lambda m: f(m)[1][0] if m != m_opt else -1, m_min, m_opt, y=k_min, y_rel_tol=tol)
                m_max = min(m_k, m_max)
            except ValueError:
                pass

        m_best, _ = gss(lambda m: f(m)[0], m_min, m_max, y_rel_tol=tol, find_min=True)
        return f(m_best)

    @staticmethod
    def sigma_von_mises(k, p, m, sigma):
        """
        k   : probing point, radius ratio
        p   : working pressure (from within)
        m   : autofrettaged (plastically yielded region) rim radius over
            : barrel radius.

        Calculate the von Misses stress at point radius ratio k.
        When supplied with k = m, the result is for at plastic-elastic
        juncture. This is the limiting stress point for an auto-
        frettaged gun barrel under internal pressure loading.

        In general, for increasing m up until m=k, the ability of a barrel to
        tolerate stress increase.
        """
        sigma_tr = (
            sigma * (k / m) ** 2 * ((m / k) ** 2 - (1 - (m / k) ** 2 + 2 * log(m)) / (k**2 - 1))
            + 2 * p / (k**2 - 1) * (k / m) ** 2
        )
        return sigma_tr * 3**0.5 * 0.5  # convert Tresca to von Misses equivalent stress


if __name__ == "__main__":
    """standard 7 hole cylinder has d_0=e_1, hole dia = 0.5 * arc width
    d_0 = 4*2e_1+3*d_0 = 11 * e_1
    """
    compositions = Composition.read_file("prop/propellants.csv")

    M17 = compositions["M17"]

    from prop import MultPerfGeometry

    M17SHC = Propellant(M17, MultPerfGeometry.SEVEN_PERF_CYLINDER, 2, 2.5)

    lf = 0.5
    print("DELTA/rho:", lf)
    cm = 0.1
    test = Gun(
        caliber=0.05,
        shot_mass=1.0,
        propellant=M17SHC,
        web=6.66e-3,
        charge_mass=cm,
        chamber_volume=cm / M17SHC.rho_p / lf,
        start_pressure=30e6,
        length_gun=3.5,
        chambrage=1.5,
        drag_coefficient=0.05,
    )

    result = test.integrate(0, 1e-3, dom=DOMAIN_TIME, sol=SOL_LAGRANGE)
