from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .ballistics import Domain, GunType, OptimizationTarget, Point, SolutionMethod
from .ballistics.config import DesignConstraint, Material, PropellantLoad, SolverConfig, StructuralConfig
from .ballistics.gun import GunConfig
from .ballistics.prop import Propellant
from .ballistics.recoilless import NozzleConfig


@dataclass
class SimulationConfig:
    """Complete simulation configuration gathered from UI widgets.

    All values are in SI units unless noted. Ballistics config objects are
    constructed from this via dispatch.sim_config_to_ballistics().
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

    caliber: float
    shot_mass: float
    gun_length: float
    chamber_volume: float
    web: float
    chambrage: float
    nozzle_expansion: float
    nozzle_efficiency: float

    # Propellant
    propellant: Propellant
    charge_mass: float
    start_pressure: float

    # Constraints
    design_pressure: float
    design_velocity: float
    min_web: float
    max_length: float
    max_iterations: int
    tolerance: float

    # Material
    structural_material: Material | None
    structural_safety_factor: float
    autofrettage: bool

    # Guidegraph
    compute_guide: bool
    guide_min_cmr: float
    guide_max_cmr: float
    guide_step_cmr: float
    guide_step_lf: float

    # Sampling
    step: int
    drag_coefficient: float

    # Runtime-injected (not from UI)
    logger: Any = field(default=None, repr=False)

    @property
    def gcfg(self) -> GunConfig:
        return GunConfig(
            caliber=self.caliber,
            shot_mass=self.shot_mass,
            barrel_length=self.gun_length,
            chamber_volume=self.chamber_volume,
            web_thickness=self.web,
            chambrage=self.chambrage,
        )

    @property
    def load(self) -> PropellantLoad:
        return PropellantLoad(
            propellant=self.propellant,
            charge_mass=self.charge_mass,
            start_pressure=self.start_pressure,
        )

    @property
    def solver(self) -> SolverConfig:
        return SolverConfig(
            tolerance=self.tolerance,
            max_iterations=self.max_iterations,
            solution_method=self.solution_method,
            drag_coefficient=self.drag_coefficient,
        )

    @property
    def structural(self) -> StructuralConfig | None:
        return (
            StructuralConfig(
                material=self.structural_material,
                safety_factor=self.structural_safety_factor,
                autofrettage=self.autofrettage,
            )
            if self.structural_material
            else None
        )

    @property
    def design(self) -> DesignConstraint:
        return DesignConstraint(
            design_pressure=self.design_pressure,
            design_velocity=self.design_velocity,
            min_web=self.min_web,
            max_length=self.max_length,
            pressure_control=self.pressure_control_point,
        )

    @property
    def nozzle(self) -> NozzleConfig:
        return NozzleConfig(
            expansion_ratio=self.nozzle_expansion,
            efficiency=self.nozzle_efficiency,
        )
