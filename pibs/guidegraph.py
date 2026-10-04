from __future__ import annotations

import logging
import math
import multiprocessing
from dataclasses import dataclass, replace

import psutil
from tqdm import tqdm

from .ballistics import Point
from .ballistics.config import SimConfig
from .ballistics.constrained import validate_load_fraction
from .ballistics.simulate import simulate


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
    cfg: SimConfig,
    load_fraction: float,
    charge_mass_ratio: float,
    solve_cache: dict[tuple[float, float], tuple[float, float]],
) -> GuideResultLine | None:
    charge_mass = cfg.shot_mass * charge_mass_ratio
    load_density = load_fraction * cfg.propellant.rho_p
    try:
        half_web, length_gun = solve_cache[(load_fraction, charge_mass_ratio)]

        chamber_volume = charge_mass / load_density

        sim_cfg = replace(
            cfg,
            web_thickness=2 * half_web,
            barrel_length=length_gun,
            chamber_volume=chamber_volume,
            charge_mass=charge_mass,
        )

        result = simulate(sim_cfg, step=0)

        try:
            burnout = result.read_table_data(Point.BURNOUT).travel / length_gun
        except ValueError:
            burnout = 1 / result.read_table_data(Point.EXIT).burnup

        volume = chamber_volume + length_gun * result.params.s

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


def guide_graph(*, cfg: SimConfig, logger: logging.Logger | None = None) -> GuideResults:
    logger = logger if logger else logging.getLogger(__name__)
    tqdm_logger = TqdmLogger(logger)
    tqdm_kwargs = dict(
        file=tqdm_logger,
        ascii=True,
        ncols=20,
        bar_format="|{bar}|{percentage:5.0f}%",
        smoothing=0,
    )

    # workers must not pickle the queue-connected logger
    worker_cfg = replace(cfg, logger=None)

    charge_mass_ratios = [
        i * cfg.guide_step_cmr
        for i in range(
            math.ceil(cfg.guide_min_cmr / cfg.guide_step_cmr), math.ceil(cfg.guide_max_cmr / cfg.guide_step_cmr)
        )
    ]

    processes = psutil.cpu_count(logical=False) or 1
    logger.info(f"dispatching {processes} processes for max load fractions")

    with multiprocessing.Pool(processes=processes, initializer=_pool_init) as pool:
        proposed_lfs = [i * cfg.guide_step_lf for i in range(math.ceil(1 / cfg.guide_step_lf))]
        iterable = tuple((worker_cfg, cmr, proposed_lfs) for cmr in charge_mass_ratios)
        validated_lfs = pool.starmap(
            func=validate_load_fraction,
            iterable=tqdm(iterable, **tqdm_kwargs),
        )

    parameters = []
    solve_cache = {}
    for charge_mass_ratio, lf_results in zip(charge_mass_ratios, validated_lfs):
        for load_fraction, (e_1, l_g, _) in lf_results.items():
            parameters.append((load_fraction, charge_mass_ratio))
            solve_cache[(load_fraction, charge_mass_ratio)] = (e_1, l_g)

    logger.info(f"dispatching {processes} processes for guidance diagram")

    with multiprocessing.Pool(processes=processes, initializer=_pool_init) as pool:
        results = pool.starmap(
            func=f,
            iterable=tqdm(
                [(worker_cfg, lf, cmr, solve_cache) for (lf, cmr) in parameters],
                total=len(parameters),
                **tqdm_kwargs,
            ),
        )

    return GuideResults(lines=[result for result in results if result is not None])
