import logging
import sys
import traceback
from logging.handlers import QueueHandler

from .ballistics.config import SimConfig
from .ballistics.constrained import find_min_v, solve_constrained
from .ballistics.simulate import simulate
from .guidegraph import guide_graph


def calculate(job_queue, log_queue, cfg: SimConfig):
    logger = logging.getLogger(__name__)
    logger.addHandler(QueueHandler(log_queue))
    logger.setLevel(logging.INFO)
    result, guide_results = None, None
    try:
        cfg.logger = logger

        """
        - constrained -> match: velocity & pressure
            - varies: web, gun length
        - constrained + optimize -> match: velocity & pressure
            - varies: chamber volume, web, gun length
        - constrained + lock_length -> match: pressure
            - varies: web
        """

        if cfg.constrained:
            if cfg.optimize:  # constrained + optimize
                lf, e_1, l_g = find_min_v(
                    cfg,
                    charge_mass_ratio=cfg.charge_mass / cfg.shot_mass,
                    opt_target=cfg.optimization_target,
                    logger=logger,
                )
                cfg.chamber_volume = cfg.charge_mass / cfg.propellant.rho_p / lf

            else:  # constrained and constrained + lock_length
                e_1, l_g = solve_constrained(
                    cfg,
                    load_fraction=cfg.charge_mass / cfg.chamber_volume / cfg.propellant.rho_p,
                    charge_mass_ratio=cfg.charge_mass / cfg.shot_mass,
                    length_gun=cfg.barrel_length,
                    known_bore=cfg.lock_length,
                    logger=logger,
                )

            cfg.web_thickness = 2 * e_1  # update web in all cases.
            if not cfg.lock_length:
                cfg.barrel_length = l_g

        result = simulate(cfg, logger=logger)

        # guidegraph calculation
        if cfg.compute_guide and (cfg.constrained and not cfg.lock_length):
            guide_results = guide_graph(cfg=cfg, logger=logger)

    except Exception as e:
        if cfg.debug:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
        else:
            logger.error(str(e))
    finally:
        job_queue.put((cfg, result, guide_results))
