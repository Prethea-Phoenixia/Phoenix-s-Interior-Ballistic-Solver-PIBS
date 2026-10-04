from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from . import (
    VALID_PRESSURE_POINTS,
    VALID_SOLUTION_METHODS,
    Domain,
    GunType,
    OptimizationTarget,
    Point,
    SolutionMethod,
)

if TYPE_CHECKING:
    from .prop import Propellant


def _positive(name: str, value: float) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be > 0, got {value}")


def _non_negative(name: str, value: float) -> None:
    if value < 0:
        raise ValueError(f"{name} must be >= 0, got {value}")


def _unit_interval(name: str, value: float) -> None:
    if not 0 <= value <= 1:
        raise ValueError(f"{name} must be in [0, 1], got {value}")


def _one_of(name: str, value: str, options: tuple) -> None:
    if value not in options:
        raise ValueError(f"{name} must be one of {options}, got {value!r}")


@dataclass(frozen=True)
class Material:
    density: float
    yield_strength: float
    desc: str = ""


@dataclass
class SimConfig:
    """Complete simulation configuration. All values are in SI units.

    Single flat config consumed directly by the ballistics core
    (params.derive_params, simulate.simulate, constrained solvers).
    Constraint validation lives here in the domain layer; the UI only
    normalizes (empty -> 0.0, str -> float) before constructing this.
    """

    # Mode flags
    optimize: bool
    constrained: bool
    debug: bool
    lock_length: bool
    gun_type: GunType
    domain: Domain
    solution_method: SolutionMethod
    pressure_control_point: Point
    optimization_target: OptimizationTarget

    # Gun geometry
    caliber: float
    shot_mass: float
    barrel_length: float
    chamber_volume: float
    web_thickness: float
    chambrage: float

    # Nozzle (recoilless only)
    expansion_ratio: float
    efficiency: float

    # Propellant charge
    propellant: Propellant
    charge_mass: float
    start_pressure: float

    # Design constraints (constrained mode only)
    design_pressure: float
    design_velocity: float
    min_web: float
    max_length: float

    # Numerics
    max_iterations: int
    tolerance: float
    drag_coefficient: float

    # Structural (optional barrel strength calculation)
    structural_material: Material | None
    structural_safety_factor: float
    autofrettage: bool

    # Guide graph
    compute_guide: bool
    guide_min_cmr: float
    guide_max_cmr: float
    guide_step_cmr: float
    guide_step_lf: float

    # Sampling
    step: int

    # Runtime-injected (not from UI)
    logger: Any = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        _positive("caliber", self.caliber)
        _positive("shot_mass", self.shot_mass)
        if self.chambrage < 1:
            raise ValueError(f"chambrage must be >= 1, got {self.chambrage}")

        # geometry fields solved by the constrained designers are validated as
        # inputs only in the modes where they are actually consumed
        _non_negative("barrel_length", self.barrel_length)
        _non_negative("chamber_volume", self.chamber_volume)
        _non_negative("web_thickness", self.web_thickness)
        if not (self.constrained and not self.lock_length):
            _positive("barrel_length", self.barrel_length)
        if not (self.constrained and self.optimize):
            _positive("chamber_volume", self.chamber_volume)
        if not self.constrained:
            _positive("web_thickness", self.web_thickness)

        _positive("charge_mass", self.charge_mass)
        _positive("start_pressure", self.start_pressure)

        _non_negative("expansion_ratio", self.expansion_ratio)
        _unit_interval("efficiency", self.efficiency)
        if self.gun_type == GunType.RECOILLESS:
            if self.expansion_ratio < 1:
                raise ValueError(f"expansion_ratio must be >= 1, got {self.expansion_ratio}")
            _positive("efficiency", self.efficiency)

        _positive("tolerance", self.tolerance)
        _positive("max_iterations", self.max_iterations)
        if not 0 <= self.drag_coefficient < 1:
            raise ValueError(f"drag_coefficient must be in [0, 1), got {self.drag_coefficient}")
        _one_of("solution_method", self.solution_method, VALID_SOLUTION_METHODS)

        if self.constrained:
            _positive("design_pressure", self.design_pressure)
            if self.lock_length:  # velocity is not a design input in lock-length mode
                _non_negative("design_velocity", self.design_velocity)
            else:
                _positive("design_velocity", self.design_velocity)
            _positive("min_web", self.min_web)
            _positive("max_length", self.max_length)
            _one_of("pressure_control_point", self.pressure_control_point, VALID_PRESSURE_POINTS)

        if self.structural_material is not None:
            _positive("material density", self.structural_material.density)
            _positive("material yield strength", self.structural_material.yield_strength)
            _positive("structural_safety_factor", self.structural_safety_factor)
