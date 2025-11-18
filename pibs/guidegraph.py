from __future__ import annotations

import logging
import multiprocessing
from itertools import repeat

import psutil

from .ballistics import CONVENTIONAL, POINT_BURNOUT, RECOILLESS, POINT_EXIT
from .ballistics.constrained_gun import ConstrainedGun
from .ballistics.constrained_recoilless import ConstrainedRecoilless
from .ballistics.gun import Gun
from .ballistics.recoilless import Recoilless

logger = logging.getLogger(__name__)


def starstarmap(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(apply_args_and_kwargs, args_for_starmap)


def apply_args_and_kwargs(fn, args, kwargs):
    return fn(*args, **kwargs)


def f(
    target: ConstrainedGun,
    load_fraction: float,
    charge_mass_ratio: float,
    **kwargs,
) -> tuple[float, float, float | None, float | None, float | None, float | None]:

    charge_mass = target.m * charge_mass_ratio
    load_density = load_fraction * target.propellant.rho_p
    try:
        half_web, length_gun = target.solve(load_fraction=load_fraction, charge_mass_ratio=charge_mass_ratio)

        chamber_volume = charge_mass / load_density

        if kwargs["typ"] == CONVENTIONAL:
            gun_class = Gun
        elif kwargs["typ"] == RECOILLESS:
            gun_class = Recoilless
        else:
            raise ValueError("Unknown gun type")

        gun = gun_class(
            **{
                **kwargs,
                **{
                    "web": 2 * half_web,
                    "charge_mass": charge_mass,
                    "chamber_volume": chamber_volume,
                    "length_gun": length_gun,
                },
            }
        )

        gun_result = gun.integrate(**{**kwargs, **{"step": 0}})

        try:
            burnout = gun_result.read_table_data(POINT_BURNOUT).travel / length_gun
        except ValueError:
            burnout = 1 / gun_result.read_table_data(POINT_EXIT).burnup

        volume = chamber_volume + length_gun * target.s  # convert to liters

    except ValueError:
        half_web, length_gun, volume, burnout = None, None, None, None

    return load_density, charge_mass, half_web, length_gun, volume, burnout


def guide_graph(*_, **kwargs):
    typ = kwargs["typ"]

    if typ == CONVENTIONAL:
        target = ConstrainedGun(**kwargs)
    elif typ == RECOILLESS:
        target = ConstrainedRecoilless(**kwargs)
    else:
        raise ValueError("unknown gun type")

    cmrs = []

    min_cmr, max_cmr, step_cmr, step_lf = (kwargs[k] for k in ("min_cmr", "max_cmr", "step_cmr", "step_lf"))

    charge_mass_ratio = min_cmr
    while charge_mass_ratio < max_cmr + 0.5 * step_cmr:
        cmrs.append(charge_mass_ratio)
        charge_mass_ratio += step_cmr

    processes = psutil.cpu_count(logical=False)
    logger.info(f"Dispatching {processes:} processes for finding maximum load fractions.")

    with multiprocessing.Pool(processes=processes) as pool:
        lfmaxs = pool.map(target.maximum_load_fraction, cmrs)

    parameters = []
    for charge_mass_ratio, max_lf in zip(cmrs, lfmaxs):
        load_fraction = target.minimum_load_fraction

        while load_fraction < max_lf + 0.5 * step_lf:
            load_fraction += step_lf

            kv = {k: v for k, v in kwargs.items()}
            kv.update({"load_fraction": load_fraction, "charge_mass_ratio": charge_mass_ratio})
            parameters.append(kv)

    logger.info(f"Dispatching {processes:} processes for constructing guidance diagram.")

    with multiprocessing.Pool(processes=processes) as pool:
        results = starstarmap(pool, f, repeat([target], len(parameters)), parameters)

    return [result for result in results if result[2]]
