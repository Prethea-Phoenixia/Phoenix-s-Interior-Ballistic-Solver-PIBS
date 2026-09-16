from __future__ import annotations

import abc
import logging
import math
from typing import Callable

from . import OptimizationTarget
from .config import DesignConstraint, GunConfig, PropellantLoad, SolverConfig
from .num import gss
from .prop import DelegatesPropellant


class Constrained(DelegatesPropellant, abc.ABC):
    def __init__(
        self,
        gcfg: GunConfig,
        load: PropellantLoad,
        design: DesignConstraint,
        solver: SolverConfig,
        logger: logging.Logger | None = None,
    ):
        self.logger = logger if logger else logging.getLogger(__name__)
        super().__init__(propellant=load.propellant)

        self.gcfg = gcfg
        self.load = load
        self.design = design
        self.solver = solver

        if any(
            (
                gcfg.caliber <= 0,
                gcfg.shot_mass <= 0,
                load.start_pressure <= 0,
                self.solver.drag_coefficient < 0,
                self.solver.drag_coefficient >= 1,
                gcfg.chambrage < 1,
            )
        ):
            raise ValueError("Invalid parameters for constrained design")

        if any((design.design_pressure <= 0, design.design_velocity <= 0)):
            raise ValueError("Invalid design constraint")

        self.caliber = gcfg.caliber

        self.s = (gcfg.caliber / 2) ** 2 * math.pi
        self.m = gcfg.shot_mass
        self.p_0 = load.start_pressure
        self.phi_1 = 1 / (1 - self.solver.drag_coefficient)

        self.p_d = design.design_pressure
        self.v_d = design.design_velocity

        self.min_web = design.min_web
        self.max_length = design.max_length

        self.chi_k = gcfg.chambrage
        self.tol = self.solver.tolerance

    def solve(
        self,
        load_fraction: float,
        charge_mass_ratio: float,
        known_bore: bool,
        length_gun: float | None = None,
    ) -> tuple[float, float]:
        raise NotImplementedError()

    @staticmethod
    def validate_solve_inputs(solve):
        def wrapped_solve(
            self: Constrained,
            load_fraction: float,
            charge_mass_ratio: float,
            known_bore: bool,
            length_gun: float | None = None,
            **kwargs,
        ):
            if load_fraction < self.minimum_load_fraction:
                raise ValueError(
                    "Design pressure cannot be achieved, in the limit of closed bomb operation, at the current load fraction."
                )
            if load_fraction >= 1:
                raise ValueError("Chamber is overfull.")

            if known_bore:
                if length_gun is None:
                    raise ValueError("known_bore option requires length_gun")

            return solve(
                self,
                load_fraction=load_fraction,
                charge_mass_ratio=charge_mass_ratio,
                known_bore=known_bore,
                length_gun=length_gun,
                **kwargs,
            )

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
        f = self.get_f(charge_mass_ratio)
        start = self.minimum_load_fraction
        stop = 1 - self.tol
        tol = self.tol
        delta = stop - start
        probe = new_probe = start
        while abs(2 * delta) > tol:
            try:
                f(new_probe)
                probe = new_probe
            except ValueError:
                delta *= 0.5
            finally:
                new_probe = probe + delta

        return probe

    def validate_load_fraction(
        self, charge_mass_ratio: float, proposed_lfs: list[float]
    ) -> dict[float, tuple[float, float, float]]:
        # construct the load fraction ladder:
        f = self.get_f(charge_mass_ratio)
        proposed_lfs = [lf for lf in proposed_lfs if lf > self.minimum_load_fraction]
        results = {}
        for lf in proposed_lfs:
            try:
                results[lf] = f(lf)
            except ValueError:
                break

        return results

    def find_min_v(
        self,
        *,
        charge_mass_ratio: float,
        opt_target: OptimizationTarget = OptimizationTarget.MIN_BARR_VOLUME,
        **_,
    ) -> tuple[float, float, float]:
        """
        find the minimum volume solution.
        """
        self.logger.info("optimizing under constraints")
        low = self.minimum_load_fraction
        self.logger.info(f"min Δ/ρ = {low:.3%}")

        high = self.maximum_load_fraction(charge_mass_ratio)
        self.logger.info(f"max Δ/ρ = {high:.3%}")

        if opt_target == OptimizationTarget.MIN_PROJ_TRAVEL:
            _f_index = 1
        elif opt_target == OptimizationTarget.MIN_BARR_VOLUME:
            _f_index = 2
        else:
            raise ValueError(f"Unknown target {opt_target}")

        _f: Callable[[float], tuple[float, float, float]] = self.get_f(charge_mass_ratio)

        self.logger.info(f"Δ/ρ range: {low:.3%} - {high:.3%}")
        lf_low, lf_high = gss(
            lambda load_fraction: _f(load_fraction)[_f_index], low, high, x_tol=self.tol, find_min=True
        )

        lf = 0.5 * (lf_high + lf_low)
        e_1 = _f(lf)[0]
        l_g = _f(lf)[1]
        self.logger.info(f"optimal Δ/ρ = {lf:.2f}")
        return lf, e_1, l_g
