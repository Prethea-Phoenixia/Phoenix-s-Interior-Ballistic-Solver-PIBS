import logging
import sys
import traceback
from logging.handlers import QueueHandler

from .ballistics import (
    CONVENTIONAL,
    RECOILLESS,
    ConstrainedGun,
    ConstrainedRecoilless,
    DesignConstraint,
    Environment,
    GunGeometry,
    Gun,
    Nozzle,
    PropellantLoad,
    Recoilless,
    Solver,
    Structural,
)
from .config import SimulationConfig
from .guidegraph import guide_graph


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
    environment = Environment(
        ambient_pressure=cfg.ambient_pressure,
        ambient_density=cfg.ambient_density,
        adiabatic_index=cfg.ambient_adiabatic_index,
    )
    structural = Structural(
        material=cfg.structural_material,
        safety_factor=cfg.structural_safety_factor,
        autofrettage=cfg.autofrettage,
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
    return geometry, load, solver, environment, structural, design, nozzle


def calculate(job_queue, log_queue, cfg: SimulationConfig):
    logger = logging.getLogger(__name__)
    logger.addHandler(QueueHandler(log_queue))
    logger.setLevel(logging.INFO)
    logger.info("calculation started.")
    gun, gun_result = None, None
    try:
        cfg.logger = logger
        geo, load, solver, env, struct, design, nozzle = sim_config_to_ballistics(cfg)

        if cfg.constrained:
            if cfg.gun_type == CONVENTIONAL:
                constrained = ConstrainedGun(
                    geometry=geo, load=load, design=design, solver=solver, environment=env, logger=logger
                )
            elif cfg.gun_type == RECOILLESS:
                constrained = ConstrainedRecoilless(
                    geometry=geo, load=load, design=design, nozzle=nozzle, solver=solver, environment=env, logger=logger
                )
            else:
                raise ValueError("unknown gun type")

            if cfg.optimize:
                l_f, e_1, l_g = constrained.find_min_v(
                    charge_mass_ratio=cfg.charge_mass_ratio,
                    opt_target=cfg.optimization_target,
                )
                cfg.load_fraction = l_f
                cfg.chamber_volume = cfg.charge_mass / cfg.propellant.rho_p / cfg.load_fraction
            else:
                e_1, l_g = constrained.solve(
                    load_fraction=cfg.load_fraction,
                    charge_mass_ratio=cfg.charge_mass_ratio,
                    known_bore=cfg.lock_length,
                )

            geo.web_thickness = 2 * e_1
            if not cfg.lock_length:
                geo.barrel_length = l_g

        if cfg.gun_type == CONVENTIONAL:
            gun = Gun(geometry=geo, load=load, solver=solver, environment=env, logger=logger)
        elif cfg.gun_type == RECOILLESS:
            gun = Recoilless(geometry=geo, load=load, nozzle=nozzle, solver=solver, environment=env, logger=logger)
        else:
            raise ValueError("unknown gun type")

        gun_result = gun.integrate(step=cfg.step, dom=cfg.domain)

        if struct.material:
            gun.structure(gun_result, structural=struct, step=cfg.step)

        logger.info("calculation concluded.")

    except Exception:
        logger.error("exception while calculating:")
        exc_type, exc_value, exc_traceback = sys.exc_info()
        logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
    finally:
        job_queue.put((cfg, gun, gun_result))


def guide(guide_job_queue, log_queue, cfg: SimulationConfig):
    logger = logging.getLogger(__name__)
    logger.addHandler(QueueHandler(log_queue))
    logger.setLevel(logging.INFO)
    logger.info("guidance diagram calculation started")
    cfg.logger = logger

    geo, load, solver, env, struct, design, nozzle = sim_config_to_ballistics(cfg)

    guide_results = None
    try:
        if cfg.gun_type == CONVENTIONAL:
            target = ConstrainedGun(
                geometry=geo, load=load, design=design, solver=solver, environment=env, logger=logger
            )
        elif cfg.gun_type == RECOILLESS:
            target = ConstrainedRecoilless(
                geometry=geo, load=load, design=design, nozzle=nozzle, solver=solver, environment=env, logger=logger
            )
        else:
            raise ValueError("unknown gun type")

        guide_results = guide_graph(
            target=target,
            gun_type=cfg.gun_type,
            min_cmr=cfg.guide_min_cmr,
            max_cmr=cfg.guide_max_cmr,
            step_cmr=cfg.guide_step_cmr,
            step_lf=cfg.guide_step_lf,
            logger=logger,
        )
        logger.info("guidance diagram calculation concluded.")

    except Exception:
        guide_results = None
        logger.error("exception while calculating guidance diagram:")
        exc_type, exc_value, exc_traceback = sys.exc_info()
        logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))

    finally:
        guide_job_queue.put(guide_results)
