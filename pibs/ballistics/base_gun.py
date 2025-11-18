from __future__ import annotations
import math
import json


from . import JSONable
from .prop import Propellant, DelegatesPropellant
from .num import dekker


class BaseGun(DelegatesPropellant, JSONable):
    def __init__(
        self,
        caliber: float,
        shot_mass: float,
        propellant: Propellant,
        web: float,
        charge_mass: float,
        chamber_volume: float,
        start_pressure: float,
        length_gun: float,
        chambrage: float,
        tol: float,
        drag_coefficient: float,
        ambient_density: float,
        ambient_pressure: float,
        ambient_adb_index: float,
    ):

        super().__init__(propellant=propellant)
        self.caliber = caliber
        self.e_1 = 0.5 * web
        self.s = (0.5 * caliber) ** 2 * math.pi
        self.m = shot_mass
        self.w = charge_mass
        self.vol_0 = chamber_volume
        self.p_0 = start_pressure
        self.l_g = length_gun
        self.chi_k = chambrage
        self.l_0 = self.vol_0 / self.s
        self.l_c = self.l_0 / self.chi_k
        self.delta = self.w / self.vol_0
        self.tol = tol
        self.phi_1 = 1 / (1 - drag_coefficient)  # drag work coefficient
        self.ambient_density = ambient_density
        self.ambient_pressure = ambient_pressure
        self.ambient_adb_index = ambient_adb_index

        ambient_pressure, ambient_adb_index = max(ambient_pressure, 1), max(ambient_adb_index, 1)

        self.p_a_bar = ambient_pressure / (self.f * self.delta)
        self.c_a = (ambient_adb_index * ambient_pressure / ambient_density) ** 0.5 if ambient_density else 0
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
        self.z_0, _ = dekker(self.propellant.f_psi_z, 0, self.propellant.z_b, y=self.psi_0, y_rel_tol=tol)

        self.phi_1 = 1 / (1 - drag_coefficient)
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
                "caliber": self.caliber,
                "shot_mass": self.m,
                "propellant": json.loads(self.propellant.to_json()),
                "web": self.e_1 * 2,
                "charge_mass": self.w,
                "chamber_volume": self.vol_0,
                "start_pressure": self.p_0,
                "length_gun": self.l_g,
                "chambrage": self.chi_k,
                "tol": self.tol,
                "drag_coefficient": 1 - 1 / self.phi_1,
                "ambient_density": self.ambient_density,
                "ambient_pressure": self.ambient_pressure,
                "ambient_adb_index": self.ambient_adb_index,
            },
            ensure_ascii=False,
        )

    @classmethod
    def from_json(cls, json_dict: dict) -> BaseGun:
        deserialized_dict = {}
        for key, value in json_dict.items():
            if key == "propellant":
                value = Propellant.from_json(value)
            deserialized_dict[key] = value

        return cls(**deserialized_dict)

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
            m_opt = math.exp(_p / yield_strength)  # Tresca
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
