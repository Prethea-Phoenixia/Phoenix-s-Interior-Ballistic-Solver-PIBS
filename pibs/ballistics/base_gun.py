from __future__ import annotations

import logging
import math

from .config import GunConfig, PropellantLoad, SolverConfig, StructuralConfig
from .num import dekker
from .prop import DelegatesPropellant


class BaseGun(DelegatesPropellant):

    def __init__(
        self,
        gcfg: GunConfig,
        load: PropellantLoad,
        solver: SolverConfig,
        structural: StructuralConfig | None = None,
        logger: logging.Logger | None = None,
    ):
        self.logger = logger if logger else logging.getLogger(__name__)
        super().__init__(propellant=load.propellant)

        self.gcfg = gcfg
        self.load = load
        self.solver = solver
        self.structural = structural

        self.caliber = gcfg.caliber
        self.e_1 = 0.5 * gcfg.web_thickness
        self.s = (0.5 * self.caliber) ** 2 * math.pi
        self.m = gcfg.shot_mass
        self.w = load.charge_mass
        self.vol_0 = gcfg.chamber_volume
        self.p_0 = load.start_pressure
        self.l_g = gcfg.barrel_length
        self.chi_k = gcfg.chambrage
        self.l_0 = self.vol_0 / self.s
        self.l_c = self.l_0 / self.chi_k
        self.delta = self.w / self.vol_0
        self.tol = self.solver.tolerance
        self.phi_1 = 1 / (1 - self.solver.drag_coefficient)

        self.psi_0 = (1 / self.delta - 1 / self.rho_p) / (self.f / self.p_0 + self.alpha - 1 / self.rho_p)
        if self.psi_0 <= 0:
            raise ValueError(
                "Initial burnup fraction is solved to be negative. This indicate an excessively high load density for start-pressure."
            )
        elif self.psi_0 >= 1:
            raise ValueError(
                "Initial burnup fraction is solved to be greater than unity. This indicate an excessively low loading density for start-pressure."
            )
        self.z_0 = dekker(self.propellant.f_psi_z, 0, self.propellant.z_b, y=self.psi_0, y_rel_tol=self.tol)

        self.phi_1 = 1 / (1 - self.solver.drag_coefficient)
        self.phi = self.phi_1 + self.w / (3 * self.m)

        self.v_j = (2 * self.f * self.w / (self.theta * self.phi * self.m)) ** 0.5

        self.b = (
            self.s**2
            * self.e_1**2
            / (self.f * self.phi * self.w * self.m * self.u_1**2)
            * (self.f * self.delta) ** (2 * (1 - self.n))
        )

    @staticmethod
    def barrel_monoblock(
        xs: list[float],
        ps: list[float],
        ss: list[float],
        yield_strength: float,
    ) -> tuple[float, list[float], list[float]]:
        v, ks, ms = 0.0, [], []

        for p in ps:
            if p > yield_strength * 0.5:
                raise ValueError(
                    f"Limit to conventional construction ({yield_strength * 0.5 * 1e-6:.3f} MPa)"
                    + " exceeded in section."
                )
            k = (1 - 2 * p / yield_strength) ** -0.5
            ks.append(k)
            ms.append(1)

        for i in range(len(xs) - 1):
            x_0, x_1 = xs[i], xs[i + 1]
            s_0, s_1 = ss[i], ss[i + 1]
            k_0, k_1 = ks[i], ks[i + 1]
            dv = 0.5 * ((k_0**2 - 1) * s_0 + (k_1**2 - 1) * s_1) * (x_1 - x_0)
            v += dv

        return v, ks, ms

    @staticmethod
    def barrel_autofrettage(
        xs: list[float],
        ps: list[float],
        ss: list[float],
        yield_strength: float,
    ) -> tuple[float, list[float], list[float]]:
        v, ks, ms = 0.0, [], []

        for _p in ps:
            m_opt = math.exp(_p / yield_strength)
            k = m_opt
            ks.append(k)
            ms.append(m_opt)

        for i in range(len(xs) - 1):
            x_0, x_1 = xs[i], xs[i + 1]
            s_0, s_1 = ss[i], ss[i + 1]
            k_0, k_1 = ks[i], ks[i + 1]
            dv = 0.5 * ((k_0**2 - 1) * s_0 + (k_1**2 - 1) * s_1) * (x_1 - x_0)
            v += dv
        return v, ks, ms
