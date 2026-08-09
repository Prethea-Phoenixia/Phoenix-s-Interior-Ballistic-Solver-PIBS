from __future__ import annotations

import json
import logging
import math
from dataclasses import asdict
from typing import TYPE_CHECKING

from . import JSONable
from .config import Environment, GunGeometry, PropellantLoad, Solver
from .num import dekker
from .prop import DelegatesPropellant

if TYPE_CHECKING:
    from .prop import Propellant


class BaseGun(DelegatesPropellant, JSONable):
    def __init__(
        self,
        geometry: GunGeometry,
        load: PropellantLoad,
        environment: Environment | None = None,
        solver: Solver | None = None,
        logger: logging.Logger | None = None,
    ):
        self.logger = logger if logger else logging.getLogger(__name__)
        super().__init__(propellant=load.propellant)

        self.geometry = geometry
        self.load = load
        self.environment = environment if environment else Environment()
        self.solver = solver if solver else Solver()

        self.caliber = geometry.caliber
        self.e_1 = 0.5 * geometry.web_thickness
        self.s = (0.5 * self.caliber) ** 2 * math.pi
        self.m = geometry.shot_mass
        self.w = load.charge_mass
        self.vol_0 = geometry.chamber_volume
        self.p_0 = load.start_pressure
        self.l_g = geometry.barrel_length
        self.chi_k = geometry.chambrage
        self.l_0 = self.vol_0 / self.s
        self.l_c = self.l_0 / self.chi_k
        self.delta = self.w / self.vol_0
        self.tol = self.solver.tolerance
        self.phi_1 = 1 / (1 - self.solver.drag_coefficient)

        self.ambient_density = self.environment.ambient_density
        self.ambient_pressure = self.environment.ambient_pressure
        self.ambient_adb_index = self.environment.adiabatic_index

        ambient_pressure, ambient_adb_index = max(self.ambient_pressure, 1), max(self.ambient_adb_index, 1)

        self.p_a_bar = ambient_pressure / (self.f * self.delta)
        self.c_a = (ambient_adb_index * ambient_pressure / self.ambient_density) ** 0.5 if self.ambient_density else 0
        self.k_1 = ambient_adb_index

        self.psi_0 = (1 / self.delta - 1 / self.rho_p) / (self.f / self.p_0 + self.alpha - 1 / self.rho_p)
        if self.psi_0 <= 0:
            raise ValueError(
                "Initial burnup fraction is solved to be negative. This indicate an excessively high load density for start-pressure."
            )
        elif self.psi_0 >= 1:
            raise ValueError(
                "Initial burnup fraction is solved to be greater than unity. This indicate an excessively low loading density for start-pressure."
            )
        self.z_0, _ = dekker(self.propellant.f_psi_z, 0, self.propellant.z_b, y=self.psi_0, y_rel_tol=self.tol)

        self.phi_1 = 1 / (1 - self.solver.drag_coefficient)
        self.phi = self.phi_1 + self.w / (3 * self.m)

        self.v_j = (2 * self.f * self.w / (self.theta * self.phi * self.m)) ** 0.5

        self.b = (
            self.s**2
            * self.e_1**2
            / (self.f * self.phi * self.w * self.m * self.u_1**2)
            * (self.f * self.delta) ** (2 * (1 - self.n))
        )

    def func_p_ad_bar(self, v_bar: float) -> float:
        if self.c_a and v_bar > 0:
            v_r = v_bar * self.v_j / self.c_a
            return (
                0.25 * self.k_1 * (self.k_1 + 1) * v_r**2
                + self.k_1 * v_r * (1 + (0.25 * (self.k_1 + 1)) ** 2 * v_r**2) ** 0.5
            ) * self.p_a_bar
        else:
            return 0.0

    def to_json(self) -> str:
        return json.dumps(
            {
                "geometry": asdict(self.geometry),
                "load": {
                    "charge_mass": self.w,
                    "start_pressure": self.p_0,
                    "propellant": json.loads(self.propellant.to_json()),
                },
                "environment": asdict(self.environment),
                "solver": asdict(self.solver),
            },
            ensure_ascii=False,
        )

    @classmethod
    def from_json(cls, json_dict: dict) -> BaseGun:
        from .prop import Propellant

        geometry = GunGeometry(**json_dict["geometry"])
        load_data = json_dict["load"]
        load = PropellantLoad(
            propellant=Propellant.from_json(load_data["propellant"]),
            charge_mass=load_data["charge_mass"],
            start_pressure=load_data["start_pressure"],
        )
        environment = Environment(**json_dict.get("environment", {}))
        solver = Solver(**json_dict.get("solver", {}))

        return cls(geometry=geometry, load=load, environment=environment, solver=solver)

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
