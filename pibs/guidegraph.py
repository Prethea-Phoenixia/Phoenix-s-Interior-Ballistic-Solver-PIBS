import logging
import multiprocessing

from itertools import repeat
from math import inf
from typing import Optional

import psutil

from .ballistics import (
    CONVENTIONAL,
    POINT_BURNOUT,
    RECOILLESS,
    Constrained,
    ConstrainedRecoilless,
    Domains,
    Gun,
    GunTypes,
    Points,
    Propellant,
    Recoilless,
    Solutions,
)

logger = logging.getLogger(__name__)


def starstarmap(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(apply_args_and_kwargs, args_for_starmap)


def apply_args_and_kwargs(fn, args, kwargs):
    return fn(*args, **kwargs)


OptionalFloat = Optional[float]


def f(
    target: Constrained,
    load_fraction: float,
    charge_mass_ratio: float,
    caliber: float,
    propellant: Propellant,
    drag_coefficient: float,
    nozzle_expansion: float,
    nozzle_efficiency: float,
    tol: float,
    min_web: float,
    max_length: float,
    ambient_rho: float,
    ambient_p: float,
    ambient_gamma: float,
    control: Points,
    sol: Solutions,
    dom: Domains,
    typ: GunTypes,
    *_,
    **__,
):

    charge_mass = target.m * charge_mass_ratio
    load_density = load_fraction * target.propellant.rho_p
    try:
        half_web, length_gun = target.solve(
            load_fraction,
            charge_mass_ratio,
            tol=tol,
            minWeb=min_web,
            maxLength=max_length,
            sol=sol,
            ambientRho=ambient_rho,
            ambientP=ambient_p,
            ambientGamma=ambient_gamma,
            control=control,
        )

        chamber_volume = charge_mass / load_density

        if typ == CONVENTIONAL:
            gun = Gun(
                caliber=caliber,
                shot_mass=target.m,
                propellant=propellant,
                web=2 * half_web,
                charge_mass=charge_mass,
                chamber_volume=chamber_volume,
                start_pressure=target.p_0,
                length_gun=length_gun,
                chambrage=target.chi_k,
                drag_coefficient=drag_coefficient,
            )
        elif typ == RECOILLESS:
            gun = Recoilless(
                caliber=caliber,
                shot_mass=target.m,
                propellant=propellant,
                web=2 * half_web,
                charge_mass=charge_mass,
                chamber_volume=chamber_volume,
                start_pressure=target.p_0,
                length_gun=length_gun,
                chambrage=target.chi_k,
                drag_coefficient=drag_coefficient,
                nozzle_expansion=nozzle_expansion,
                nozzle_efficiency=nozzle_efficiency,
            )
        else:
            raise ValueError("unknown gun type")

        gun_result = gun.integrate(
            step=0, tol=tol, dom=dom, sol=sol, ambient_rho=ambient_rho, ambient_p=ambient_p, ambient_gamma=ambient_gamma
        )

        try:
            burnout = gun_result.read_table_data(POINT_BURNOUT).travel / length_gun
        except ValueError:
            burnout = inf

        tube_volume = length_gun * target.s
        volume = chamber_volume + tube_volume  # convert to liters

    except ValueError:
        half_web, length_gun, volume, burnout = None, None, None, None

    return load_density, charge_mass, half_web, length_gun, volume, burnout


def guide_graph(*_, **kwargs):
    typ = kwargs["typ"]

    if typ == CONVENTIONAL:
        target = Constrained(**kwargs)
    elif typ == RECOILLESS:
        target = ConstrainedRecoilless(**kwargs)
    else:
        raise ValueError("unknown gun type")

    load_fractions, charge_mass_ratios, charge_masses, load_densities = [], [], [], []

    min_cmr, max_cmr, step_cmr = (kwargs[k] for k in ("min_cmr", "max_cmr", "step_cmr"))
    charge_mass_ratio = min_cmr
    while charge_mass_ratio < max_cmr + 0.5 * step_cmr:
        charge_mass_ratios.append(charge_mass_ratio)
        charge_masses.append(charge_mass_ratio * target.m)
        charge_mass_ratio += step_cmr

    min_lf, max_lf, step_lf = (kwargs[k] for k in ("min_lf", "max_lf", "step_lf"))
    load_fraction = min_lf

    while load_fraction < max_lf + 0.5 * step_lf:
        load_fractions.append(load_fraction)
        load_densities.append(load_fraction * target.propellant.rho_p)
        load_fraction += step_lf

    parameters = []
    for charge_mass_ratio in charge_mass_ratios:
        for load_fraction in load_fractions:
            kv = {k: v for k, v in kwargs.items()}
            kv.update({"load_fraction": load_fraction, "charge_mass_ratio": charge_mass_ratio})
            parameters.append(kv)

    processes = psutil.cpu_count(logical=False)
    logger.info(f"dispatching {processes:}-process for constructing guidance diagram.")

    # parallel implementation
    with multiprocessing.Pool(processes=processes) as pool:
        results = starstarmap(pool, f, repeat([target], len(parameters)), parameters)

    shaped_results = [
        results[i * len(load_fractions) : (i + 1) * len(load_fractions)] for i in range(len(charge_mass_ratios))
    ]

    return shaped_results
