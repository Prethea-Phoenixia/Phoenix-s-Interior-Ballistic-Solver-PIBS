import logging
import multiprocessing
import os
from itertools import repeat
from math import inf
from typing import Optional

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
    loadFraction: float,
    chargeMassRatio: float,
    caliber: float,
    propellant: Propellant,
    dragCoefficient: float,
    nozzleExpansion: float,
    nozzleEfficiency: float,
    tol: float,
    minWeb: float,
    maxLength: float,
    ambientRho: float,
    ambientP: float,
    ambientGamma: float,
    control: Points,
    sol: Solutions,
    dom: Domains,
    typ: GunTypes,
    *_,
    **__,
):

    chargeMass = target.m * chargeMassRatio
    loadDensity = loadFraction * target.propellant.rho_p
    try:
        halfWeb, lengthGun = target.solve(
            loadFraction,
            chargeMassRatio,
            tol=tol,
            minWeb=minWeb,
            maxLength=maxLength,
            sol=sol,
            ambientRho=ambientRho,
            ambientP=ambientP,
            ambientGamma=ambientGamma,
            control=control,
        )

        chamberVolume = chargeMass / loadDensity

        if typ == CONVENTIONAL:
            gun = Gun(
                caliber=caliber,
                shotMass=target.m,
                propellant=propellant,
                grainSize=2 * halfWeb,
                chargeMass=chargeMass,
                chamberVolume=chamberVolume,
                startPressure=target.p_0,
                lengthGun=lengthGun,
                chambrage=target.chi_k,
                dragCoefficient=dragCoefficient,
            )
        elif typ == RECOILLESS:
            gun = Recoilless(
                caliber=caliber,
                shotMass=target.m,
                propellant=propellant,
                grainSize=2 * halfWeb,
                chargeMass=chargeMass,
                chamberVolume=chamberVolume,
                startPressure=target.p_0,
                lengthGun=lengthGun,
                chambrage=target.chi_k,
                dragCoefficient=dragCoefficient,
                nozzleExpansion=nozzleExpansion,
                nozzleEfficiency=nozzleEfficiency,
            )
        else:
            raise ValueError("unknown gun type")

        gunResult = gun.integrate(
            step=0, tol=tol, dom=dom, sol=sol, ambientRho=ambientRho, ambientP=ambientP, ambientGamma=ambientGamma
        )

        try:
            burnout = gunResult.readTableData(POINT_BURNOUT).travel / lengthGun
        except ValueError:
            burnout = inf

        tubeVolume = lengthGun * target.S
        volume = chamberVolume + tubeVolume  # convert to liters

    except ValueError:
        halfWeb, lengthGun, volume, burnout = None, None, None, None

    return loadDensity, chargeMass, halfWeb, lengthGun, volume, burnout


def guideGraph(*_, **kwargs):
    typ = kwargs["typ"]

    if typ == CONVENTIONAL:
        target = Constrained(**kwargs)
    elif typ == RECOILLESS:
        target = ConstrainedRecoilless(**kwargs)
    else:
        raise ValueError("unknown gun type")

    loadFractions, chargeMassRatios, chargeMasses, loadDensities = [], [], [], []

    minCMR, maxCMR, stepCMR = (kwargs[k] for k in ("minCMR", "maxCMR", "stepCMR"))
    chargeMassRatio = minCMR
    while chargeMassRatio < maxCMR + 0.5 * stepCMR:
        chargeMassRatios.append(chargeMassRatio)
        chargeMasses.append(chargeMassRatio * target.m)
        chargeMassRatio += stepCMR

    minLF, maxLF, stepLF = (kwargs[k] for k in ("minLF", "maxLF", "stepLF"))
    loadFraction = minLF

    while loadFraction < maxLF + 0.5 * stepLF:
        loadFractions.append(loadFraction)
        loadDensities.append(loadFraction * target.propellant.rho_p)
        loadFraction += stepLF

    parameters = []
    for chargeMassRatio in chargeMassRatios:
        for loadFraction in loadFractions:
            kv = {k: v for k, v in kwargs.items()}
            kv.update({"loadFraction": loadFraction, "chargeMassRatio": chargeMassRatio})
            parameters.append(kv)

    processes = os.cpu_count()
    logger.info(f"dispatching {processes:}-process for constructing guidance diagram.")

    # parallel implementation
    with multiprocessing.Pool(processes=processes) as pool:
        results = starstarmap(pool, f, repeat([target], len(parameters)), parameters)

    shapedResults = [
        results[i * len(loadFractions) : (i + 1) * len(loadFractions)] for i in range(len(chargeMassRatios))
    ]

    return shapedResults
