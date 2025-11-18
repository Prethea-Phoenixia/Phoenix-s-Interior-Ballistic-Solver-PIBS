import sys
import traceback
from logging.handlers import QueueHandler

from . import root_logger
from .ballistics import CONVENTIONAL, RECOILLESS, ConstrainedGun, ConstrainedRecoilless, Gun, Recoilless
from .guidegraph import guide_graph


def calculate(job_queue, log_queue, kwargs):
    root_logger.addHandler(QueueHandler(log_queue))
    root_logger.info("calculation started.")

    gun_type = kwargs["typ"]
    constrain = kwargs["con"]
    optimize = kwargs["opt"]
    lock = kwargs["lock"]

    gun, gun_result = None, None
    try:
        if constrain:
            if gun_type == CONVENTIONAL:
                constrained = ConstrainedGun(**kwargs)
            elif gun_type == RECOILLESS:
                constrained = ConstrainedRecoilless(**kwargs)
            else:
                raise ValueError("unknown gun type")

            if optimize:
                if gun_type == CONVENTIONAL or gun_type == RECOILLESS:
                    l_f, e_1, l_g = constrained.find_min_v(**kwargs)

                else:
                    raise ValueError("unknown gun type")

                kwargs.update({"load_fraction": l_f})
                chamber_volume = kwargs["charge_mass"] / kwargs["propellant"].rho_p / kwargs["load_fraction"]
                kwargs.update({"chamber_volume": chamber_volume})
            else:
                if gun_type == CONVENTIONAL or gun_type == RECOILLESS:
                    e_1, l_g = constrained.solve(**kwargs)

                else:
                    raise ValueError("unknown gun type")

            kwargs.update({"web": 2 * e_1})

            if not lock:
                kwargs.update({"length_gun": l_g})

        if gun_type == CONVENTIONAL:
            gun = Gun(**kwargs)
        elif gun_type == RECOILLESS:
            gun = Recoilless(**kwargs)
        else:
            raise ValueError("unknown gun type")

        gun_result = gun.integrate(**kwargs)

        if kwargs["structural_material"]:
            gun.structure(gun_result, **kwargs)

        root_logger.info("calculation concluded.")

    except Exception as e:
        gun, gun_result = None, None
        root_logger.error("exception while calculating:")
        if kwargs["deb"]:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            root_logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
        else:
            root_logger.error(str(e))
    finally:
        job_queue.put((kwargs, gun, gun_result))


def guide(guide_job_queue, log_queue, kwargs):
    root_logger.addHandler(QueueHandler(log_queue))
    root_logger.info("guidance diagram calculation started")

    guide_results = None
    try:
        guide_results = guide_graph(**kwargs)
        root_logger.info("guidance diagram calculation concluded.")

    except Exception as e:
        guide_results = None
        root_logger.error("exception while calculating guidance diagram:")
        if kwargs["deb"]:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            root_logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
        else:
            root_logger.error(str(e))

    finally:
        guide_job_queue.put(guide_results)
