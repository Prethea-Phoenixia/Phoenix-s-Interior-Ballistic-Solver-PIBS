from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from . import VALID_PRESSURE_POINTS, VALID_SOLUTION_METHODS

if TYPE_CHECKING:
    from .material import Material
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


@dataclass
class GunGeometry:
    caliber: float
    shot_mass: float
    barrel_length: float = 0.0
    chamber_volume: float = 0.0
    web_thickness: float = 0.0
    chambrage: float = 2.0

    def __post_init__(self):
        _positive("caliber", self.caliber)
        _positive("shot_mass", self.shot_mass)
        _non_negative("barrel_length", self.barrel_length)
        _non_negative("chamber_volume", self.chamber_volume)
        _non_negative("web_thickness", self.web_thickness)
        _positive("chambrage", self.chambrage)


@dataclass
class PropellantLoad:
    propellant: Propellant
    charge_mass: float = 0.0
    start_pressure: float = 0.0

    def __post_init__(self):
        _non_negative("charge_mass", self.charge_mass)
        _non_negative("start_pressure", self.start_pressure)


@dataclass
class Nozzle:
    expansion_ratio: float = 0.0
    efficiency: float = 0.92

    def __post_init__(self):
        _non_negative("expansion_ratio", self.expansion_ratio)
        _unit_interval("efficiency", self.efficiency)


@dataclass
class DesignConstraint:
    design_pressure: float = 0.0
    design_velocity: float = 0.0
    min_web: float = 1e-6
    max_length: float = 1e3
    pressure_control: str = "PEAK_BREECH_P"

    def __post_init__(self):
        _non_negative("design_pressure", self.design_pressure)
        _non_negative("design_velocity", self.design_velocity)
        _positive("min_web", self.min_web)
        _positive("max_length", self.max_length)
        _one_of("pressure_control", self.pressure_control, VALID_PRESSURE_POINTS)


@dataclass
class Solver:
    tolerance: float = 1e-5
    max_iterations: int = 10
    solution_method: str = "SOL_PIDDUCK"
    drag_coefficient: float = 0.0

    def __post_init__(self):
        _positive("tolerance", self.tolerance)
        _positive("max_iterations", self.max_iterations)
        _unit_interval("drag_coefficient", self.drag_coefficient)
        _one_of("solution_method", self.solution_method, VALID_SOLUTION_METHODS)


@dataclass
class Structural:
    material: Material | None = None
    safety_factor: float = 1.35
    autofrettage: bool = False

    def __post_init__(self):
        _positive("safety_factor", self.safety_factor)
