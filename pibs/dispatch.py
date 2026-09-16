import logging
import sys
import traceback
from logging.handlers import QueueHandler

from .ballistics import GunType
from .ballistics.constr_gun import ConstrainedGun
from .ballistics.constr_rcl import ConstrainedRecoilless
from .ballistics.gun import Gun
from .ballistics.recoilless import Recoilless
from .config import SimulationConfig
from .guidegraph import guide_graph


def calculate(job_queue, log_queue, cfg: SimulationConfig):
    logger = logging.getLogger(__name__)
    logger.addHandler(QueueHandler(log_queue))
    logger.setLevel(logging.INFO)
    logger.info("calculation started")
    gun, gun_result, guide_results = None, None, None
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
            if cfg.gun_type == GunType.CONVENTIONAL:
                constrained = ConstrainedGun(
                    gcfg=cfg.gcfg, load=cfg.load, design=cfg.design, solver=cfg.solver, logger=logger
                )
            elif cfg.gun_type == GunType.RECOILLESS:
                constrained = ConstrainedRecoilless(
                    gcfg=cfg.gcfg,
                    load=cfg.load,
                    design=cfg.design,
                    nozzle=cfg.nozzle,
                    solver=cfg.solver,
                    logger=logger,
                )
            else:
                raise ValueError("unknown gun type")

            if cfg.optimize:  # constrained + optimize
                lf, e_1, l_g = constrained.find_min_v(
                    charge_mass_ratio=cfg.charge_mass_ratio,
                    opt_target=cfg.optimization_target,
                )
                cfg.chamber_volume = cfg.charge_mass / cfg.propellant.rho_p / lf

            else:  # constrained and constrained + lock_length
                e_1, l_g = constrained.solve(
                    load_fraction=cfg.charge_mass / cfg.chamber_volume / cfg.propellant.rho_p,
                    charge_mass_ratio=cfg.charge_mass_ratio,
                    length_gun=cfg.gcfg.barrel_length,
                    known_bore=cfg.lock_length,
                )

            cfg.web = 2 * e_1  # update web in all cases.
            if not cfg.lock_length:
                cfg.gun_length = l_g

        if cfg.gun_type == GunType.CONVENTIONAL:
            gun = Gun(gcfg=cfg.gcfg, load=cfg.load, solver=cfg.solver, logger=logger)
        elif cfg.gun_type == GunType.RECOILLESS:
            gun = Recoilless(gcfg=cfg.gcfg, load=cfg.load, nozzle=cfg.nozzle, solver=cfg.solver, logger=logger)
        else:
            raise ValueError("unknown gun type")

        gun_result = gun.integrate(step=cfg.step, dom=cfg.domain)
        gun.structure(gun_result, structural=cfg.structural, step=cfg.step)

        # guide graph calculation
        if cfg.compute_guide and (cfg.constrained and not cfg.lock_length):
            try:

                if cfg.gun_type == GunType.CONVENTIONAL:
                    target = ConstrainedGun(
                        gcfg=cfg.gcfg, load=cfg.load, design=cfg.design, solver=cfg.solver, logger=logger
                    )
                elif cfg.gun_type == GunType.RECOILLESS:
                    target = ConstrainedRecoilless(
                        gcfg=cfg.gcfg,
                        load=cfg.load,
                        design=cfg.design,
                        nozzle=cfg.nozzle,
                        solver=cfg.solver,
                        logger=logger,
                    )
                else:
                    raise ValueError("unknown gun type")

                logger.info("guidance diagram calculation started")
                guide_results = guide_graph(
                    target=target,
                    gun_type=cfg.gun_type,
                    min_cmr=cfg.guide_min_cmr,
                    max_cmr=cfg.guide_max_cmr,
                    step_cmr=cfg.guide_step_cmr,
                    step_lf=cfg.guide_step_lf,
                    logger=logger,
                )
                logger.info("guidance diagram calculation concluded")
            except Exception:
                logger.error("exception while calculating guidance diagram:")
                exc_type, exc_value, exc_traceback = sys.exc_info()
                logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))

        logger.info("calculation concluded")

    except Exception:
        logger.error("exception while calculating")
        exc_type, exc_value, exc_traceback = sys.exc_info()
        logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
    finally:
        job_queue.put((cfg, gun, gun_result, guide_results))
