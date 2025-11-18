from __future__ import annotations

import json
import logging
from typing import Callable, Any
import math

from . import POINT_PEAK_AVG, Optimization_Targets, MIN_BARR_VOLUME, Points, MAX_ITER, MIN_PROJ_TRAVEL
from . import JSONable
from .num import gss

from .prop import Propellant, DelegatesPropellant


def probe_func(
    func: Callable[[float], Any], start: float, stop: float, tol: float, exceptions: tuple[Exception] = (ValueError,)
) -> float:
    delta = stop - start
    probe = new_probe = start
    while abs(2 * delta) > tol:
        try:
            func(new_probe)
            probe = new_probe
        except exceptions:
            delta *= 0.5
        finally:
            new_probe = probe + delta

    return probe


logger = logging.getLogger(__name__)


class Constrained(DelegatesPropellant, JSONable):
    def __init__(
        self,
        caliber: float,
        shot_mass: float,
        propellant: Propellant,
        start_pressure: float,
        drag_coefficient: float,
        design_pressure: float,
        design_velocity: float,
        tol: float,
        chambrage: float,
        min_web: float = 1e-6,
        max_length: float = 1e3,
        ambient_density: float = 1.204,
        ambient_pressure: float = 101.325e3,
        ambient_adb_index: float = 1.4,
        control: Points = POINT_PEAK_AVG,
    ):
        super().__init__(propellant=propellant)
        if any(
            (
                caliber <= 0,
                shot_mass <= 0,
                start_pressure <= 0,
                drag_coefficient < 0,
                drag_coefficient >= 1,
                chambrage < 1,
            )
        ):
            raise ValueError("Invalid parameters for constrained design")

        if any((design_pressure <= 0, design_velocity <= 0)):
            raise ValueError("Invalid design constraint")

        ambient_pressure = max(ambient_pressure, 1)
        ambient_adb_index = max(ambient_adb_index, 1)

        self.caliber = caliber

        self.s = (caliber / 2) ** 2 * math.pi
        self.m = shot_mass
        self.p_0 = start_pressure
        self.phi_1 = 1 / (1 - drag_coefficient)

        # design limits
        self.p_d = design_pressure
        self.v_d = design_velocity

        self.ambient_density = ambient_density
        self.ambient_pressure = ambient_pressure
        self.ambient_adb_index = ambient_adb_index
        self.control = control

        self.min_web = min_web
        self.max_length = max_length

        self.chi_k = chambrage
        self.tol = tol

    def to_json(self) -> str:
        return json.dumps(
            {
                "caliber": self.caliber,
                "shot_mass": self.m,
                "propellant": json.loads(self.propellant.to_json()),
                "start_pressure": self.p_0,
                "drag_coefficient": 1 - 1 / self.phi_1,
                "design_pressure": self.p_d,
                "design_velocity": self.v_d,
                "tol": self.tol,
                "min_web": self.min_web,
                "max_length": self.max_length,
                "ambient_density": self.ambient_density,
                "ambient_pressure": self.ambient_pressure,
                "ambient_adb_index": self.ambient_adb_index,
                "control": self.control,
                "chambrage": self.chi_k,
            },
            ensure_ascii=False,
        )

    @classmethod
    def from_json(cls, json_dict: dict) -> Constrained:
        deserialized_dict = {}
        for key, value in json_dict.items():
            if key == "propellant":
                value = Propellant.from_json(value)

            deserialized_dict[key] = value

        return cls(**deserialized_dict)

    def func_p_ad_bar(self, v_bar: float, c_a_bar: float, p_a_bar: float) -> float:
        if c_a_bar and v_bar > 0.0:
            v_r = v_bar / c_a_bar
            return (
                +0.25 * self.ambient_adb_index * (self.ambient_adb_index + 1) * v_r**2
                + self.ambient_adb_index * v_r * (1 + (0.25 * (self.ambient_adb_index + 1)) ** 2 * v_r**2) ** 0.5
            ) * p_a_bar

        else:
            return 0

    def solve(
        self,
        *_,
        load_fraction: float,
        charge_mass_ratio: float,
        length_gun: float | None = None,
        max_iteration: int = MAX_ITER,
        labda_1: float | None = None,
        labda_2: float | None = None,
        cc: float | None = None,
        it: int = 0,
        **__,
    ) -> tuple[float, float]: ...

    @staticmethod
    def validate_solve_inputs(solve):
        def wrapped_solve(self: Constrained, *args, load_fraction: float, charge_mass_ratio: float, **kwargs):
            if charge_mass_ratio <= 0:
                raise ValueError("Charge mass to projectile ratio must be positive")

            if load_fraction < self.minimum_load_fraction:
                raise ValueError(
                    "Design pressure cannot be achieved, in the limit of closed bomb operation, at the current load fraction."
                )
            if load_fraction >= 1:
                raise ValueError("Chamber is overfull (load fraction >= 1).")

            return solve(self, *args, load_fraction=load_fraction, charge_mass_ratio=charge_mass_ratio, **kwargs)

        return wrapped_solve

    @property
    def minimum_load_fraction(self) -> float:
        return (1 / (self.f / self.p_d + self.alpha)) / self.rho_p * (1 + self.tol)

    def get_f(self, charge_mass_ratio: float) -> Callable[[float], tuple[float, float, float]]:
        def _f(load_fraction: float) -> tuple[float, float, float]:
            e_1_delta, l_g_delta = self.solve(
                load_fraction=load_fraction, charge_mass_ratio=charge_mass_ratio, known_bore=False
            )
            l_0 = (self.m * charge_mass_ratio / (self.rho_p * load_fraction)) / self.s
            return e_1_delta, l_g_delta, l_g_delta + l_0

        return _f

    def maximum_load_fraction(self, charge_mass_ratio: float) -> float:
        return probe_func(
            self.get_f(charge_mass_ratio), start=self.minimum_load_fraction, stop=1 - self.tol, tol=self.tol
        )

    def find_min_v(
        self,
        charge_mass_ratio: float,
        opt_target: Optimization_Targets = MIN_BARR_VOLUME,
        **_,
    ) -> tuple[float, float, float]:
        """
        find the minimum volume solution.
        """
        logger.info("Optimizing under constraints.")
        """
        p = fΔ / (1 - αΔ)
        Δ = 1 / (f/p + α)
        """
        low = self.minimum_load_fraction
        logger.info(f"Min Δ/ρ = {low:.3%}.")

        high = self.maximum_load_fraction(charge_mass_ratio)
        logger.info(f"Max Δ/ρ = {high:.3%}.")

        if opt_target == MIN_PROJ_TRAVEL:
            _f_index = 1
        elif opt_target == MIN_BARR_VOLUME:
            _f_index = 2
        else:
            raise ValueError(f"Unknown target {opt_target}")

        _f = self.get_f(charge_mass_ratio)

        logger.info(f"Solution constrained to Δ/ρ : {low:.3%} - {high:.3%}")
        lf_low, lf_high = gss(
            lambda load_fraction: _f(load_fraction)[_f_index], low, high, x_tol=self.tol, find_min=True
        )
        lf = 0.5 * (lf_high + lf_low)
        e_1, l_g, _ = _f(lf)
        logger.info(f"Optimal Δ/ρ = {lf:.2f}")
        return lf, e_1, l_g
