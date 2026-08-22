from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .ballistics.material import Material
from .ballistics.prop import Propellant


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
    propellant: Propellant | None
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

    # Guide graph params
    guide_min_cmr: float
    guide_max_cmr: float
    guide_step_cmr: float
    guide_step_lf: float

    # Sampling
    step: int
    drag_coefficient: float

    # Runtime-injected (not from UI)
    logger: Any = field(default=None, repr=False)
