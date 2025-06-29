from math import inf
from typing import Optional
import matplotlib.pyplot as plt
import multiprocessing

from .ballistics import Constrained, Solutions, Points, Domains

from .ballistics import Propellant
from .ballistics import Gun
from .ballistics import DOMAIN_TIME
from .ballistics import SOL_LAGRANGE, SOL_PIDDUCK, SOL_MAMONTOV
from .ballistics import POINT_PEAK_AVG, POINT_PEAK_BREECH, POINT_PEAK_SHOT, POINT_BURNOUT

from labellines import labelLines

from itertools import repeat


def starstarmap(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(apply_args_and_kwargs, args_for_starmap)


def apply_args_and_kwargs(fn, args, kwargs):
    return fn(*args, **kwargs)


OptionalFloat = Optional[float]


def f(
    loadFraction: float,
    chargeMassRatio: float,
    target: Constrained,
    caliber: float,
    propellant: Propellant,
    dragCoefficient: float,
    tol: float,
    minWeb: float,
    maxLength: float,
    ambientRho: float,
    ambientP: float,
    ambientGamma: float,
    control: Points,
    sol: Solutions,
    dom: Domains,
):
    chargeMass = target.m * chargeMassRatio
    loadDensity = loadFraction * target.propellant.rho_p
    try:
        halfWeb, lengthGun = target.solve(
            loadFraction,
            chargeMassRatio,
            tol,
            minWeb=minWeb,
            maxLength=maxLength,
            sol=sol,
            ambientRho=ambientRho,
            ambientP=ambientP,
            ambientGamma=ambientGamma,
            control=control,
        )

        chamberVolume = chargeMass / loadDensity

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
        gunResult = gun.integrate(
            step=0, tol=tol, dom=dom, sol=sol, ambientRho=ambientRho, ambientP=ambientP, ambientGamma=ambientGamma
        )

        try:
            burnout = gunResult.readTableData(POINT_BURNOUT).travel / lengthGun
        except ValueError:
            burnout = 1.0

        tubeVolume = lengthGun * target.S
        volume = chamberVolume + tubeVolume  # convert to liters

    except ValueError:
        halfWeb, lengthGun, volume, burnout = None, None, None, None

    return loadDensity, chargeMass, halfWeb, lengthGun, burnout, volume


def guideGraph(
    caliber: float,
    shotMass: float,
    propellant: Propellant,
    startPressure: float,
    dragCoefficient: float,
    designPressure: float,
    designVelocity: float,
    chambrage: float,
    tol: float,
    minWeb: float = 1e-6,
    maxLength: float = 1e3,
    sol: Solutions = SOL_LAGRANGE,
    dom: Domains = DOMAIN_TIME,
    ambientRho: float = 1.204,
    ambientP: float = 101.325e3,
    ambientGamma: float = 1.4,
    control: Points = POINT_PEAK_AVG,
    minLF: float = 0.3,
    maxLF: float = 0.9,
    stepLF: float = 0.1,
    minCMR: float = 0.3,
    maxCMR: float = 0.9,
    stepCMR: float = 0.1,
    *_,
    **__,
):
    target = Constrained(
        tol=tol,
        caliber=caliber,
        shotMass=shotMass,
        propellant=propellant,
        startPressure=startPressure,
        dragCoefficient=dragCoefficient,
        designPressure=designPressure,
        designVelocity=designVelocity,
        chambrage=chambrage,
    )

    loadFractions = []
    chargeMassRatios = []

    chargeMasses = []
    loadDensities = []

    chargeMassRatio = minCMR
    while chargeMassRatio < maxCMR + 0.5 * stepCMR:
        chargeMassRatios.append(chargeMassRatio)
        chargeMasses.append(chargeMassRatio * target.m)
        chargeMassRatio += stepCMR

    loadFraction = minLF
    while loadFraction < maxLF + 0.5 * stepLF:
        loadFractions.append(loadFraction)
        loadDensities.append(loadFraction * target.propellant.rho_p)
        loadFraction += stepLF

    parameters = []
    for chargeMassRatio in chargeMassRatios:
        for loadFraction in loadFractions:
            parameters.append(
                {
                    "loadFraction": loadFraction,
                    "chargeMassRatio": chargeMassRatio,
                    "minWeb": minWeb,
                    "maxLength": maxLength,
                    "ambientRho": ambientRho,
                    "ambientP": ambientP,
                    "ambientGamma": ambientGamma,
                    "control": control,
                    "sol": sol,
                    "dom": dom,
                    "target": target,
                    "caliber": caliber,
                    "propellant": propellant,
                    "dragCoefficient": dragCoefficient,
                    "tol": tol,
                }
            )

    # parallel implementation
    with multiprocessing.Pool() as pool:
        results = starstarmap(pool, f, repeat([], len(parameters)), parameters)

    shapedResults = [
        results[i * len(loadFractions) : (i + 1) * len(loadFractions)] for i in range(len(chargeMassRatios))
    ]

    return shapedResults
