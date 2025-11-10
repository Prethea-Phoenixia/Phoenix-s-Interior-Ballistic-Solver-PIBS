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
    GenericResult,
    OutlineEntry,
    PressureProbePoint,
    PressureTraceEntry,
)
from .num import dekker, gss, integrate, rkf
from .prop import Propellant
from .material import Material

logger = logging.getLogger(__name__)


@dataclass
class GunResult(GenericResult):
    gun: Gun
    table_data: list[GunTableEntry]


@dataclass
class GunTableEntry(GenericEntry):
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

        i, _ = integrate(lambda x: f(om, x), 0, 1, tol)

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
    i_u, _ = integrate(lambda x: g(omega, x), 0, 1, tol)
    i_l, _ = integrate(lambda x: f(omega, x), 0, 1, tol)
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
        self.w = charge_mass
        self.vol_0 = chamber_volume
        self.p_0 = start_pressure
        self.l_g = length_gun
        self.chi_k = chambrage  # ration of l_0 / l_chamber
        self.l_0 = self.vol_0 / self.s
        self.l_c = self.l_0 / self.chi_k
        self.delta = self.w / self.vol_0

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

        return max(p_bar, self.p_a_bar)

    def ode_t(self, _: float, zlv: tuple[float, float, float], __: float) -> tuple[float, float, float]:
        z, l_bar, v_bar = zlv
        p_bar = self.f_p_bar(z, l_bar, v_bar)

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

    def ode_l(self, l_bar: float, tzv: tuple[float, float, float], _: float) -> tuple[float, float, float]:
        """length domain ode of internal ballistics
        the 1/v_bar pose a starting problem that prevent us from using it from
        initial condition."""
        t_bar, z, v_bar = tzv

        p_bar = self.f_p_bar(z, l_bar, v_bar)
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

    def ode_z(self, z: float, tlv: tuple[float, float, float], _: float):
        t_bar, l_bar, v_bar = tlv
        p_bar = self.f_p_bar(z, l_bar, v_bar)

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
            return self.s * p * (l + l_psi) / (self.w * psi * r)
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
        ambient_p, ambient_gamma = max(ambient_p, 1), max(ambient_gamma, 1)

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
            lambda _z: self.propellant.f_psi_z(_z) - self.psi_0, 0, 1, x_tol=tol, y_rel_tol=tol, y_abs_tol=tol**2
        )

        record = []

        if any((step < 0, tol < 0)):
            raise ValueError("Invalid integration specification")

        if ambient_rho < 0:
            raise ValueError("Invalid ambient condition")

        if sol == SOL_LAGRANGE:
            labda_1, labda_2 = 1 / 2, 1 / 3
        elif sol == SOL_PIDDUCK:
            labda_1, labda_2 = pidduck(self.w / (self.phi_1 * self.m), self.theta + 1, tol)
        elif sol == SOL_MAMONTOV:
            labda_1, labda_2 = pidduck(self.w / (self.phi_1 * self.m), 1, tol)
        else:
            raise ValueError("Unknown Solution")

        self.labda_1, self.labda_2 = labda_1, labda_2

        labda = self.l_g / self.l_0
        cc = 1 - (1 - 1 / self.chi_k) * log(labda + 1) / labda  # chambrage correction factor

        self.phi = self.phi_1 + labda_2 * cc * self.w / self.m
        """
        见《枪炮内弹道学》（金，2014）p.70 式
        """

        self.b = (
            self.s**2
            * self.e_1**2
            / (self.f * self.phi * self.w * self.m * self.u_1**2)
            * (self.f * self.delta) ** (2 * (1 - self.n))
        )

        self.v_j = (2 * self.f * self.w / (self.theta * self.phi * self.m)) ** 0.5

        t_scale, p_scale = self.l_0 / self.v_j, self.f * self.delta

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

        self.append_bar_data(bar_data, tag=POINT_START, t_bar=0, l_bar=0, z=z_0, v_bar=0)

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
        p_max = 1e9
        p_bar_max = p_max / p_scale

        while z_i < z_b:  # terminates if burnout is achieved
            ztlv_record_i = []
            if z_j == z_i:
                raise ValueError("Numerical accuracy exhausted in search of exit/burnout point.")

            if z_j > z_b:
                z_j = z_b
            z_j, (t_bar_j, l_bar_j, v_bar_j) = rkf(
                self.ode_z,
                (t_bar_i, l_bar_i, v_bar_i),
                z_i,
                z_j,
                rel_tol=tol,
                abs_tol=tol**2,
                abort_func=lambda _x, _ys, _records: self.abort_z(
                    _x, _ys, _records, p_bar_max=p_bar_max, l_g_bar=l_g_bar
                ),
                record=ztlv_record_i,
            )

            p_bar_j = self.f_p_bar(z_j, l_bar_j, v_bar_j)

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
                    _, (t_bar, l_bar, v_bar) = ztlv_record[-1]
                    raise ValueError(
                        "Squib load condition detected: Shot stopped in bore.\n"
                        + "Shot is last calculated at {:.0f} mm at {:.0f} mm/s after {:.0f} ms".format(
                            l_bar * self.l_0 * 1e3, v_bar * self.v_j * 1e3, t_bar * t_scale * 1e3
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
                            p_bar_j * p_scale * 1e-6,
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
            (_t_bar, (_l_bar, self.f_psi_z(_z), _v_bar, self.f_p_bar(_z, _l_bar, _v_bar)))
            for (_z, (_t_bar, _l_bar, _v_bar)) in ztlv_record
            if _t_bar
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
        l_bar, (t_bar_e, z_e, v_bar_e) = rkf(
            self.ode_l, (t_bar_i, z_i, v_bar_i), l_bar_i, l_g_bar, rel_tol=tol, abs_tol=tol**2, record=ltzv_record
        )

        if l_bar != l_g_bar:
            if v_bar_e <= 0:
                raise ValueError(
                    "Squib load condition detected post burnout:"
                    + " Round stopped in bore at {:.0f} mm".format(l_bar * self.l_0 * 1e3)
                )

        record.extend(
            (_t_bar, (_l_bar, self.f_psi_z(_z), _v_bar, self.f_p_bar(_z, _l_bar, _v_bar)))
            for _l_bar, (_t_bar, _z, _v_bar) in ltzv_record
        )
        self.append_bar_data(bar_data, tag=POINT_EXIT, t_bar=t_bar_e, l_bar=l_g_bar, z=z_e, v_bar=v_bar_e)

        t_bar_f = None
        if z_b > 1.0 and z_e >= 1.0:  # fracture point exist and is contained
            """
            Subscript f indicate fracture condition
            ODE w.r.t Z is integrated from Z_0 to 1, from onset of projectile movement to charge fracture
            """
            t_bar_f, l_bar_f, v_bar_f = rkf(self.ode_z, (0, 0, 0), z_0, 1, rel_tol=tol, abs_tol=tol**2)[1]
            self.append_bar_data(bar_data, tag=POINT_FRACTURE, t_bar=t_bar_f, l_bar=l_bar_f, z=1, v_bar=v_bar_f)

        t_bar_b = None
        if is_burn_out_contained:
            """
            Subscript b indicated burnout condition
            ODE w.r.t Z is integrated from Z_0 to z_b, from onset of projectile movement to charge burnout.
            """

            t_bar_b, l_bar_b, v_bar_b = rkf(self.ode_z, (0, 0, 0), z_0, z_b, rel_tol=tol, abs_tol=tol**2)[1]
            self.append_bar_data(bar_data, tag=POINT_BURNOUT, t_bar=t_bar_b, l_bar=l_bar_b, z=z_b, v_bar=v_bar_b)

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
            t_bar_tol = tol * min(_t_bar for _t_bar in (t_bar_e, t_bar_b, t_bar_f) if _t_bar is not None)
            t_bar_p = 0.5 * sum(gss(f, 0, t_bar_e if t_bar_b is None else t_bar_b, x_tol=t_bar_tol, find_min=False))

            z_p, l_bar_p, v_bar_p = self.g(t_bar_p, tag, tol)[1]
            self.append_bar_data(bar_data, tag=tag, t_bar=t_bar_p, l_bar=l_bar_p, z=z_p, v_bar=v_bar_p)

        for i, peak in enumerate([POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH]):
            find_peak(lambda _t_bar: self.g(_t_bar, peak, z_0)[0], peak)

        """
        populate data for output purposes
        """

        if dom == DOMAIN_TIME:
            z_j, l_bar_j, v_bar_j, t_bar_j = z_0, 0, 0, 0
            for j in range(step):
                t_bar_k = t_bar_e / (step + 1) * (j + 1)
                z_j, l_bar_j, v_bar_j = rkf(
                    self.ode_t, (z_j, l_bar_j, v_bar_j), t_bar_j, t_bar_k, rel_tol=tol, abs_tol=tol**2
                )[1]
                t_bar_j = t_bar_k

                self.append_bar_data(bar_data, tag=SAMPLE, t_bar=t_bar_j, l_bar=l_bar_j, z=z_j, v_bar=v_bar_j)

        elif dom == DOMAIN_LEN:
            """
            Due to two issues, i.e. 1.the length domain ODE cannot be integrated from the origin point, and 2.the
            correct behaviour can only be expected when starting from a point with active burning else dZ flat lines.
            we do another Z domain integration to seed the initial values to a point where ongoing burning is guaranteed.
            (in the case of gun barrel length >= burn length, the group of value by subscript i will not guarantee
            burning is still ongoing).
            """
            t_bar_j = 0.5 * t_bar_i
            z_j, l_bar_j, v_bar_j = rkf(self.ode_t, (z_0, 0, 0), 0, t_bar_j, rel_tol=tol, abs_tol=tol**2)[1]

            for j in range(step):
                l_bar_k = l_g_bar / (step + 1) * (j + 1)

                t_bar_j, z_j, v_bar_j = rkf(
                    self.ode_l, (t_bar_j, z_j, v_bar_j), l_bar_j, l_bar_k, rel_tol=tol, abs_tol=tol**2
                )[1]
                l_bar_j = l_bar_k

                self.append_bar_data(bar_data, tag=SAMPLE, t_bar=t_bar_j, l_bar=l_bar_j, z=z_j, v_bar=v_bar_j)

        logger.info(f"sampled for {step} points.")

        """
        sort the data points
        """

        data = []
        p_trace = []
        trace_steps = max(step, 1)

        for bar_data_line in bar_data:
            dtag, t_bar, l_bar, z, v_bar, p_bar = bar_data_line

            t = t_bar * t_scale
            l = l_bar * self.l_0
            psi = self.f_psi_z(z)
            v = v_bar * self.v_j
            p = p_bar * p_scale
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

            data.append(table_entry)

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
            data.append(table_entry)

        data, p_trace = zip(
            *(
                (tableEntry, pTrace)
                for tableEntry, pTrace in sorted(zip(data, p_trace), key=lambda entries: entries[0].time)
            )
        )

        # calculate a pressure and flow velocity tracing.
        gun_result = GunResult(self, data, p_trace)

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

    def append_bar_data(
        self,
        bar_data: list[tuple[str, float, float, float, float, float]],
        tag: str,
        t_bar: float,
        l_bar: float,
        z: float,
        v_bar: float,
    ):
        bar_data.append((tag, t_bar, l_bar, z, v_bar, self.f_p_bar(z, l_bar, v_bar)))

    def g(self, t_bar: float, tag: str, tol: float) -> tuple[float, tuple[float, float, float]]:
        z, l_bar, v_bar = rkf(self.ode_t, (self.z_0, 0, 0), 0, t_bar, rel_tol=tol, abs_tol=tol**2)[1]
        p_bar = self.f_p_bar(z, l_bar, v_bar)
        if tag == POINT_PEAK_AVG:
            return p_bar, (z, l_bar, v_bar)
        else:
            ps_bar, pb_bar = self.to_ps_pb(l_bar * self.l_0, p_bar)
            if tag == POINT_PEAK_SHOT:
                return ps_bar, (z, l_bar, v_bar)
            elif tag == POINT_PEAK_BREECH:
                return pb_bar, (z, l_bar, v_bar)
            else:
                raise ValueError("tag not handled.")

    def abort_z(self, z: float, tlv: tuple[float, float, float], _, p_bar_max: float, l_g_bar: float):
        t_bar, l_bar, v_bar = tlv
        p_bar = self.f_p_bar(z, l_bar, v_bar)

        return l_bar > l_g_bar or p_bar > p_bar_max or v_bar < 0

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

        factor_s = 1 + labda_2_prime * (self.w / (self.phi_1 * self.m))  # factor_b = P/P_b = phi / phi_1
        factor_b = (self.phi_1 * self.m + labda_2_prime * self.w) / (self.phi_1 * self.m + labda_1_prime * self.w)

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
        r_c = r * chi_k**0.5
        x_probes = (
            [i / step * l_c for i in range(step)]
            + [l_c * (1 - tol)]
            + [i / step * l_g + l_c for i in range(step)]
            + [l_g + l_c]
        )
        p_probes = [0.0 for _ in range(len(x_probes))]

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

        i = step + 1
        x_c, p_c = x_probes[:i], p_probes[:i]  # c for chamber
        x_b, p_b = x_probes[i:], p_probes[i:]  # b for barrel

        if self.is_af:
            v_c, k_c, m_c = Gun.barrel_autofrettage(x_c, p_c, [s * chi_k for _ in x_c], sigma)
            v_b, k_b, m_b = Gun.barrel_autofrettage(x_b, p_b, [s for _ in x_b], sigma)
        else:
            v_c, k_c, m_c = Gun.barrel_monoblock(x_c, p_c, [s * chi_k for _ in x_c], sigma)
            v_b, k_b, m_b = Gun.barrel_monoblock(x_b, p_b, [s for _ in x_b], sigma)

        v = v_c + v_b
        k_probes = k_c + k_b
        m_probes = m_c + m_b

        hull = []
        for x, k, m in zip(x_probes, k_probes, m_probes):
            if x < l_c:
                hull.append(OutlineEntry(x, r_c, k * r_c, m * r_c))
            else:
                hull.append(OutlineEntry(x, r, k * r, m * r))

        gun_result.outline = hull
        gun_result.tubeMass = v * self.material.rho

        logger.info("conducted structural calculation.")

    @staticmethod
    def barrel_autofrettage(
        xs: list[float],  # probe location
        ps: list[float],  # probe pressure
        ss: list[float],  # probe cross-section area
        yield_strength: float,  # yield strength
    ) -> tuple[float, list[float], list[float]]:
        """
        m : r_a / r_i
        k : r_o / r_i
        n : p_vM_max / sigma

        1 < m < k

        m_opt = exp(p / sigma_y)
        P_y,i = sigma_y / 2 * (2 * ln(m) + 1 - (m/k)^2)
        P_y,o = sigma_y / 2 * (2 * ln(m) + k^2 - m^2)
        """
        v, ks, ms = 0.0, [], []

        for _p in ps:
            m_opt = exp(_p / yield_strength)  # Tresca
            k = m_opt  # full autofrettage
            ks.append(k)
            ms.append(m_opt)

        for i in range(len(xs) - 1):
            x_0, x_1 = xs[i], xs[i + 1]
            s_0, s_1 = ss[i], ss[i + 1]
            k_0, k_1 = ks[i], ks[i + 1]
            dv = 0.5 * ((k_0**2 - 1) * s_0 + (k_1**2 - 1) * s_1) * (x_1 - x_0)
            v += dv
        return v, ks, ms

    @staticmethod
    def barrel_monoblock(
        xs: list[float],  # probe location
        ps: list[float],  # probe pressure
        ss: list[float],  # probe cross-section area
        yield_strength: float,  # yield strength
    ) -> tuple[float, list[float], list[float]]:
        """
        P_tr, i = sigma_y (k^2 - 1)/(2 * k^2)
        P_tr, o = sigma_y (k^2 - 1)/2
        """
        v, ks, ms = 0.0, [], []

        for p in ps:
            if p > yield_strength * 0.5:
                raise ValueError(
                    f"Limit to conventional construction ({yield_strength * 0.5 * 1e-6:.3f} MPa)"
                    + " exceeded in section."
                )
            k = (1 - 2 * p / yield_strength) ** -0.5  # Tresca criteria
            ks.append(k)
            ms.append(1)

        for i in range(len(xs) - 1):
            x_0, x_1 = xs[i], xs[i + 1]
            s_0, s_1 = ss[i], ss[i + 1]
            k_0, k_1 = ks[i], ks[i + 1]
            dv = 0.5 * ((k_0**2 - 1) * s_0 + (k_1**2 - 1) * s_1) * (x_1 - x_0)
            v += dv

        return v, ks, ms


if __name__ == "__main__":
    pass
