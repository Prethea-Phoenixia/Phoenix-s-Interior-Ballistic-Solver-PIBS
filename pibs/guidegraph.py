from __future__ import annotations

import logging
import math
import multiprocessing
from dataclasses import dataclass, replace

import psutil
from tqdm import tqdm

from .ballistics import GunType, Point
from .ballistics.constrained import Constrained
from .ballistics.gun import Gun
from .ballistics.recoilless import Recoilless


def _pool_init():
    # Disable logging in child processes to avoid duplicate terminal output
    pibs_logger = logging.getLogger("pibs")
    pibs_logger.setLevel(logging.CRITICAL)
    pibs_logger.handlers.clear()
    pibs_logger.propagate = False


class TqdmLogger:
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.string_buffer = ""

    def write(self, msg: str) -> None:
        self.string_buffer += msg.strip("\r\n\t ")

    def flush(self) -> None:
        self.logger.info(self.string_buffer)
        self.string_buffer = ""


@dataclass
class GuideResultLine:
    load_density: float
    charge_mass: float
    half_web: float
    length_gun: float
    volume: float
    burnout: float


@dataclass
class GuideResults:
    lines: list[GuideResultLine]


def f(
    target: Constrained,
    gun_class,
    load_fraction: float,
    charge_mass_ratio: float,
    solve_cache: dict[tuple[float, float], tuple[float, float]],
    logger: logging.Logger,
) -> GuideResultLine | None:
    charge_mass = target.m * charge_mass_ratio
    load_density = load_fraction * target.propellant.rho_p
    try:
        half_web, length_gun = solve_cache[(load_fraction, charge_mass_ratio)]

        chamber_volume = charge_mass / load_density

        gcfg = replace(
            target.gcfg,
            web_thickness=2 * half_web,
            barrel_length=length_gun,
            chamber_volume=chamber_volume,
        )
        load = replace(target.load, charge_mass=charge_mass)

        nozzle = getattr(target, "nozzle", None)
        if nozzle is not None:
            gun = gun_class(
                gcfg=gcfg,
                load=load,
                nozzle=nozzle,
                solver=target.solver,
                logger=logger,
            )
        else:
            gun = gun_class(
                gcfg=gcfg,
                load=load,
                solver=target.solver,
                logger=logger,
            )

        gun_result = gun.integrate(step=0)

        try:
            burnout = gun_result.read_table_data(Point.BURNOUT).travel / length_gun
        except ValueError:
            burnout = 1 / gun_result.read_table_data(Point.EXIT).burnup

        volume = chamber_volume + length_gun * target.s

        return GuideResultLine(
            load_density=load_density,
            charge_mass=charge_mass,
            half_web=half_web,
            length_gun=length_gun,
            volume=volume,
            burnout=burnout,
        )

    except ValueError:
        return None


def guide_graph(
    *,
    target: Constrained,
    gun_type: GunType,
    min_cmr: float,
    max_cmr: float,
    step_cmr: float,
    step_lf: float,
    logger: logging.Logger | None = None,
) -> GuideResults:
    logger = logger if logger else logging.getLogger(__name__)
    tqdm_logger = TqdmLogger(logger)
    tqdm_kwargs = dict(
        file=tqdm_logger,
        ascii=True,
        ncols=20,
        bar_format="|{bar}|{percentage:5.0f}%",
        smoothing=0,
    )

    if gun_type == GunType.CONVENTIONAL:
        gun_class = Gun
    elif gun_type == GunType.RECOILLESS:
        gun_class = Recoilless
    else:
        raise ValueError("Unknown gun type")

    charge_mass_ratios = [i * step_cmr for i in range(math.ceil(min_cmr / step_cmr), math.ceil(max_cmr / step_cmr))]

    processes = psutil.cpu_count(logical=False) or 1
    logger.info(f"dispatching {processes} processes for max load fractions")

    with multiprocessing.Pool(processes=processes, initializer=_pool_init) as pool:
        proposed_lfs = [i * step_lf for i in range(math.ceil(1 / step_lf))]
        iterable = tuple((cmr, proposed_lfs) for cmr in charge_mass_ratios)
        validated_lfs = pool.starmap(
            func=target.validate_load_fraction,
            iterable=tqdm(iterable, **tqdm_kwargs),
        )

    parameters = []
    solve_cache = {}
    for charge_mass_ratio, lf_results in zip(charge_mass_ratios, validated_lfs):
        for load_fraction, (e_1, l_g, _) in lf_results.items():
            parameters.append((target, gun_class, load_fraction, charge_mass_ratio))
            solve_cache[(load_fraction, charge_mass_ratio)] = (e_1, l_g)

    logger.info(f"dispatching {processes} processes for guidance diagram")

    with multiprocessing.Pool(processes=processes, initializer=_pool_init) as pool:
        results = pool.starmap(
            func=f,
            iterable=tqdm(
                [(target, gun_class, lf, cmr, solve_cache, logger) for (target, gun_class, lf, cmr) in parameters],
                total=len(parameters),
                **tqdm_kwargs,
            ),
        )

    return GuideResults(lines=[result for result in results if result is not None])
