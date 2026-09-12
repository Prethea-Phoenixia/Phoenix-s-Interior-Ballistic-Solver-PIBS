from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .ballistics.config import DesignConstraint, PropellantLoad, Solver, Structural
from .ballistics.gun import GunGeometry
from .ballistics.material import Material
from .ballistics.prop import Propellant
from .ballistics.recoilless import Nozzle


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
    gun_type: str
    domain: str
    solution_method: str
    pressure_control_point: str
    optimization_target: str

    # Geometry
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
    charge_mass_ratio: float
    load_fraction: float
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


def sim_config_to_ballistics(cfg: SimulationConfig):
    geometry = GunGeometry(
        caliber=cfg.caliber,
        shot_mass=cfg.shot_mass,
        barrel_length=cfg.gun_length,
        chamber_volume=cfg.chamber_volume,
        web_thickness=cfg.web,
        chambrage=cfg.chambrage,
    )
    load = PropellantLoad(
        propellant=cfg.propellant,
        charge_mass=cfg.charge_mass,
        start_pressure=cfg.start_pressure,
    )
    solver = Solver(
        tolerance=cfg.tolerance,
        max_iterations=cfg.max_iterations,
        solution_method=cfg.solution_method,
        drag_coefficient=cfg.drag_coefficient,
    )

    structural = (
        Structural(
            material=cfg.structural_material,
            safety_factor=cfg.structural_safety_factor,
            autofrettage=cfg.autofrettage,
        )
        if cfg.structural_material
        else None
    )

    design = DesignConstraint(
        design_pressure=cfg.design_pressure,
        design_velocity=cfg.design_velocity,
        min_web=cfg.min_web,
        max_length=cfg.max_length,
        pressure_control=cfg.pressure_control_point,
    )
    nozzle = Nozzle(
        expansion_ratio=cfg.nozzle_expansion,
        efficiency=cfg.nozzle_efficiency,
    )
    return geometry, load, solver, structural, design, nozzle
