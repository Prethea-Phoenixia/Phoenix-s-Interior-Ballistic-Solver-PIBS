from __future__ import annotations

import logging
import math
import multiprocessing
from dataclasses import replace

import psutil
from tqdm import tqdm

from .ballistics import CONVENTIONAL, POINT_BURNOUT, POINT_EXIT, RECOILLESS
from .ballistics.constrained import Constrained
from .ballistics.gun import Gun
from .ballistics.recoilless import Recoilless


class TqdmLogger:
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.string_buffer = ""

    def write(self, msg: str) -> None:
        self.string_buffer += msg.strip("\r\n\t ")

    def flush(self) -> None:
        self.logger.info(self.string_buffer)
        self.string_buffer = ""


def f(
    target: Constrained,
    gun_class,
    load_fraction: float,
    charge_mass_ratio: float,
) -> tuple[float, float, float | None, float | None, float | None, float | None]:

    charge_mass = target.m * charge_mass_ratio
    load_density = load_fraction * target.propellant.rho_p
    try:
        half_web, length_gun = target.solve(load_fraction=load_fraction, charge_mass_ratio=charge_mass_ratio)

        chamber_volume = charge_mass / load_density

        geo = replace(
            target.geometry,
            web_thickness=2 * half_web,
            barrel_length=length_gun,
            chamber_volume=chamber_volume,
        )
        load = replace(target.load, charge_mass=charge_mass)

        nozzle = getattr(target, "nozzle", None)
        if nozzle is not None:
            gun = gun_class(
                geometry=geo, load=load, nozzle=nozzle, solver=target.solver, environment=target.environment
            )
        else:
            gun = gun_class(geometry=geo, load=load, solver=target.solver, environment=target.environment)

        gun_result = gun.integrate(step=0)

        try:
            burnout = gun_result.read_table_data(POINT_BURNOUT).travel / length_gun
        except ValueError:
            burnout = 1 / gun_result.read_table_data(POINT_EXIT).burnup

        volume = chamber_volume + length_gun * target.s

    except ValueError:
        half_web, length_gun, volume, burnout = None, None, None, None

    return load_density, charge_mass, half_web, length_gun, volume, burnout


def guide_graph(
    *,
    target: Constrained,
    gun_type: str,
    min_cmr: float,
    max_cmr: float,
    step_cmr: float,
    step_lf: float,
    logger: logging.Logger | None = None,
):
    logger = logger if logger else logging.getLogger(__name__)
    tqdm_logger = TqdmLogger(logger)
    tqdm_kwargs = dict(
        file=tqdm_logger,
        ascii=False,
        miniters=1,
        ncols=40,
        bar_format="{percentage:5.1f}%",
        smoothing=0.3,
    )

    if gun_type == CONVENTIONAL:
        gun_class = Gun
    elif gun_type == RECOILLESS:
        gun_class = Recoilless
    else:
        raise ValueError("Unknown gun type")

    charge_mass_ratios = [i * step_cmr for i in range(math.ceil(min_cmr / step_cmr), math.ceil(max_cmr / step_cmr))]

    processes = psutil.cpu_count(logical=False) or 1
    logger.info(f"Dispatching {processes} processes for finding maximum load fractions.")

    with multiprocessing.Pool(processes=processes) as pool:
        proposed_lfs = [i * step_lf for i in range(math.ceil(1 / step_lf))]
        iterable = tuple((cmr, proposed_lfs) for cmr in charge_mass_ratios)
        validated_lfs = pool.starmap(
            func=target.validate_load_fraction,
            iterable=tqdm(iterable, **tqdm_kwargs),
        )

    parameters = []
    for charge_mass_ratio, lfs in zip(charge_mass_ratios, validated_lfs):
        for load_fraction in lfs:
            parameters.append((target, gun_class, load_fraction, charge_mass_ratio))

    logger.info(f"Dispatching {processes} processes for constructing guidance diagram.")

    with multiprocessing.Pool(processes=processes) as pool:
        results = pool.starmap(
            func=f,
            iterable=tqdm(parameters, total=len(parameters), **tqdm_kwargs),
        )

    return [result for result in results if result[2]]
