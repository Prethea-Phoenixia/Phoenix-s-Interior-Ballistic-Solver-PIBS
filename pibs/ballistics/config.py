from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .material import Material
    from .prop import Propellant


@dataclass
class GunGeometry:
    caliber: float
    shot_mass: float
    barrel_length: float = 0.0
    chamber_volume: float = 0.0
    web_thickness: float = 0.0
    chambrage: float = 2.0


@dataclass
class PropellantLoad:
    propellant: Propellant
    charge_mass: float = 0.0
    start_pressure: float = 0.0


@dataclass
class Nozzle:
    expansion_ratio: float = 0.0
    efficiency: float = 0.92


@dataclass
class DesignConstraint:
    design_pressure: float = 0.0
    design_velocity: float = 0.0
    min_web: float = 1e-6
    max_length: float = 1e3
    pressure_control: str = "PEAK_BREECH_P"


@dataclass
class Solver:
    tolerance: float = 1e-5
    max_iterations: int = 10
    solution_method: str = "SOL_PIDDUCK"
    drag_coefficient: float = 0.0


@dataclass
class Structural:
    material: Material | None = None
    safety_factor: float = 1.35
    autofrettage: bool = False
