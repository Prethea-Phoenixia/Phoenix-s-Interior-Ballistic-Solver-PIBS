from pibs.ballistics import (
    CONVENTIONAL,
    POINT_PEAK_AVG,
    Composition,
    Propellant,
    SimpleGeometry,
)
from pibs.guidegraph import guide_graph

if __name__ == "__main__":
    compositions = Composition.read_file("../pibs/ballistics/resource/propellants.csv")
    pyroxylin = compositions["Pyroxylin"]

    guide_graph(
        caliber=125e-3,
        tol=1e-4,
        dragCoefficient=2e-2,
        control=POINT_PEAK_AVG,
        shotMass=5.67,
        startPressure=30e6,
        designPressure=392e6,
        designVelocity=1800,
        chambrage=1.25,
        propellant=Propellant(pyroxylin, SimpleGeometry.TUBE, 1, 100),
        minCMR=1,
        maxCMR=3,
        stepCMR=0.1,
        minLF=0.1,
        maxLF=0.9,
        stepLF=0.033,
        typ=CONVENTIONAL,
    )
