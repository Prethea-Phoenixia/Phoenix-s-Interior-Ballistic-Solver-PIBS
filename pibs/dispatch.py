import logging
import sys
import traceback
from logging.handlers import QueueHandler

from .ballistics import CONVENTIONAL, RECOILLESS
from .ballistics.cons_gun import ConstrainedGun
from .ballistics.cons_rcl import ConstrainedRecoilless
from .ballistics.gun import Gun
from .ballistics.recoilless import Recoilless
from .config import SimulationConfig, sim_config_to_ballistics
from .guidegraph import guide_graph


def calculate(job_queue, log_queue, cfg: SimulationConfig):
    logger = logging.getLogger(__name__)
    logger.addHandler(QueueHandler(log_queue))
    logger.setLevel(logging.INFO)
    logger.info("calculation started.")
    gun, gun_result, guide_results = None, None, None
    try:
        cfg.logger = logger
        geo, load, solver, struct, design, nozzle = sim_config_to_ballistics(cfg)

        """
        - constrained -> match: velocity & pressure
            - varies: web, gun length
        - constrained + optimize -> match: velocity & pressure
            - varies: chamber volume, web, gun length
        - constrained + lock_length -> match: pressure
            - varies: web
        """

        if cfg.constrained:
            if cfg.gun_type == CONVENTIONAL:
                constrained = ConstrainedGun(geometry=geo, load=load, design=design, solver=solver, logger=logger)
            elif cfg.gun_type == RECOILLESS:
                constrained = ConstrainedRecoilless(
                    geometry=geo, load=load, design=design, nozzle=nozzle, solver=solver, logger=logger
                )
            else:
                raise ValueError("unknown gun type")

            if cfg.optimize:  # constrained + optimize
                l_f, e_1, l_g = constrained.find_min_v(
                    charge_mass_ratio=cfg.charge_mass_ratio,
                    opt_target=cfg.optimization_target,
                )
                cfg.load_fraction = l_f
                cfg.chamber_volume = cfg.charge_mass / cfg.propellant.rho_p / cfg.load_fraction

            else:  # constrained and constrained + lock_length
                e_1, l_g = constrained.solve(
                    load_fraction=cfg.load_fraction,
                    charge_mass_ratio=cfg.charge_mass_ratio,
                    length_gun=geo.barrel_length,
                    known_bore=cfg.lock_length,
                )

            geo.web_thickness = 2 * e_1  # update web in all cases.
            if not cfg.lock_length:
                geo.barrel_length = l_g

        if cfg.gun_type == CONVENTIONAL:
            gun = Gun(geometry=geo, load=load, solver=solver, logger=logger)
        elif cfg.gun_type == RECOILLESS:
            gun = Recoilless(geometry=geo, load=load, nozzle=nozzle, solver=solver, logger=logger)
        else:
            raise ValueError("unknown gun type")

        gun_result = gun.integrate(step=cfg.step, dom=cfg.domain)

        if struct:
            gun.structure(gun_result, structural=struct, step=cfg.step)

        # Guide graph calculation
        if cfg.compute_guide:
            try:
                logger.info("guidance diagram calculation started")
                if cfg.gun_type == CONVENTIONAL:
                    target = ConstrainedGun(geometry=geo, load=load, design=design, solver=solver, logger=logger)
                elif cfg.gun_type == RECOILLESS:
                    target = ConstrainedRecoilless(
                        geometry=geo, load=load, design=design, nozzle=nozzle, solver=solver, logger=logger
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
                logger.error("exception while calculating guidance diagram:")
                exc_type, exc_value, exc_traceback = sys.exc_info()
                logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))

        logger.info("calculation concluded.")

    except Exception:
        logger.error("exception while calculating:")
        exc_type, exc_value, exc_traceback = sys.exc_info()
        logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
    finally:
        job_queue.put((cfg, gun, gun_result, guide_results))
