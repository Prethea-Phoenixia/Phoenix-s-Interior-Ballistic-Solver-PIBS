from __future__ import annotations

import json
import logging
import math
from dataclasses import asdict
from typing import TYPE_CHECKING, Any, Callable

from . import MIN_BARR_VOLUME, MIN_PROJ_TRAVEL, JSONable, OptimizationTargets
from .config import DesignConstraint, GunGeometry, PropellantLoad, Solver
from .num import gss
from .prop import DelegatesPropellant

if TYPE_CHECKING:
    from .prop import Propellant


def probe_func(
    func: Callable[[float], Any],
    start: float,
    stop: float,
    tol: float,
    exceptions: tuple[type[Exception], ...] = (ValueError,),
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


class Constrained(DelegatesPropellant, JSONable):
    def __init__(
        self,
        geometry: GunGeometry,
        load: PropellantLoad,
        design: DesignConstraint,
        solver: Solver | None = None,
        logger: logging.Logger | None = None,
    ):
        self.logger = logger if logger else logging.getLogger(__name__)
        super().__init__(propellant=load.propellant)

        self.geometry = geometry
        self.load = load
        self.design = design
        self.solver = solver if solver else Solver()

        if any(
            (
                geometry.caliber <= 0,
                geometry.shot_mass <= 0,
                load.start_pressure <= 0,
                self.solver.drag_coefficient < 0,
                self.solver.drag_coefficient >= 1,
                geometry.chambrage < 1,
            )
        ):
            raise ValueError("Invalid parameters for constrained design")

        if any((design.design_pressure <= 0, design.design_velocity <= 0)):
            raise ValueError("Invalid design constraint")

        self.caliber = geometry.caliber

        self.s = (geometry.caliber / 2) ** 2 * math.pi
        self.m = geometry.shot_mass
        self.p_0 = load.start_pressure
        self.phi_1 = 1 / (1 - self.solver.drag_coefficient)

        self.p_d = design.design_pressure
        self.v_d = design.design_velocity

        self.min_web = design.min_web
        self.max_length = design.max_length

        self.chi_k = geometry.chambrage
        self.tol = self.solver.tolerance

    def to_json(self) -> str:
        return json.dumps(
            {
                "geometry": asdict(self.geometry),
                "load": {
                    "start_pressure": self.p_0,
                    "propellant": json.loads(self.propellant.to_json()),
                },
                "design": asdict(self.design),
                "solver": asdict(self.solver),
            },
            ensure_ascii=False,
        )

    @classmethod
    def from_json(cls, json_dict: dict) -> Constrained:
        from .prop import Propellant

        geometry = GunGeometry(**json_dict["geometry"])
        load_data = json_dict["load"]
        load = PropellantLoad(
            propellant=Propellant.from_json(load_data["propellant"]),
            start_pressure=load_data["start_pressure"],
        )
        design = DesignConstraint(**json_dict["design"])
        solver = Solver(**json_dict.get("solver", {}))

        return cls(geometry=geometry, load=load, design=design, solver=solver)

    def solve(
        self,
        *,
        load_fraction: float,
        charge_mass_ratio: float,
        length_gun: float | None = None,
        known_bore: bool = False,
        max_iterations: int | None = None,
        labda_1: float | None = None,
        labda_2: float | None = None,
        cc: float | None = None,
        it: int = 0,
    ) -> tuple[float, float]:
        raise NotImplementedError

    @staticmethod
    def validate_solve_inputs(solve):
        def wrapped_solve(
            self: Constrained,
            *,
            load_fraction: float,
            charge_mass_ratio: float,
            **kwargs,
        ):
            if charge_mass_ratio <= 0:
                raise ValueError("Charge mass to projectile ratio must be positive")

            if load_fraction < self.minimum_load_fraction:
                raise ValueError(
                    "Design pressure cannot be achieved, in the limit of closed bomb operation, at the current load fraction."
                )
            if load_fraction >= 1:
                raise ValueError("Chamber is overfull (load fraction >= 1).")

            return solve(self, load_fraction=load_fraction, charge_mass_ratio=charge_mass_ratio, **kwargs)

        return wrapped_solve

    @property
    def minimum_load_fraction(self) -> float:
        """
        p = fΔ / (1 - αΔ)
        Δ = 1 / (f/p + α)
        """
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

    def validate_load_fraction(self, charge_mass_ratio: float, proposed_lfs: list[float]) -> list[float]:
        # construct the load fraction ladder:
        f = self.get_f(charge_mass_ratio)
        proposed_lfs = [lf for lf in proposed_lfs if lf > self.minimum_load_fraction]
        for i in range(len(proposed_lfs)):
            try:
                f(proposed_lfs[i])
            except ValueError:
                return proposed_lfs[:i]

        return proposed_lfs

    def find_min_v(
        self,
        *,
        charge_mass_ratio: float,
        opt_target: OptimizationTargets = MIN_BARR_VOLUME,
        **_,
    ) -> tuple[float, float, float]:
        """
        find the minimum volume solution.
        """
        self.logger.info("Optimizing under constraints.")
        low = self.minimum_load_fraction
        self.logger.info(f"Min Δ/ρ = {low:.3%}.")

        high = self.maximum_load_fraction(charge_mass_ratio)
        self.logger.info(f"Max Δ/ρ = {high:.3%}.")

        if opt_target == MIN_PROJ_TRAVEL:
            _f_index = 1
        elif opt_target == MIN_BARR_VOLUME:
            _f_index = 2
        else:
            raise ValueError(f"Unknown target {opt_target}")

        _f: Callable[[float], tuple[float, float, float]] = self.get_f(charge_mass_ratio)

        self.logger.info(f"Solution constrained to Δ/ρ : {low:.3%} - {high:.3%}")
        lf_low, lf_high = gss(
            lambda load_fraction: _f(load_fraction)[_f_index], low, high, x_tol=self.tol, find_min=True
        )
        lf = 0.5 * (lf_high + lf_low)
        e_1 = _f(lf)[0]
        l_g = _f(lf)[1]
        self.logger.info(f"Optimal Δ/ρ = {lf :.2f}")
        return lf, e_1, l_g
