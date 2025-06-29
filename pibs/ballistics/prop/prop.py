"""
GEOMETRIES:
    lookup dictionary mapping geometrical descriptions
    to relevant geometrical object

"""

import csv
from enum import Enum
from math import pi


class MultPerfGeometry(Enum):
    """table 1-4 from ref[1] page 33"""

    def __init__(self, desc, A, B, C, rhoDiv, nHole):
        self.desc = desc
        self.A = A
        self.B = B
        self.C = C
        self.rhoDiv = rhoDiv  # rho divided by (e_1+d_0/2)
        self.nHole = nHole
        self.useAR = False

    def __str__(self):
        return self.desc

    SEVEN_PERF_CYLINDER = ("SEVEN_PERF_CYLINDER", 1, 7, 0, 0.2956, 7)

    SEVEN_PERF_ROSETTE = ("SEVEN_PERF_ROSETTE", 2, 8, 12 * 3**0.5 / pi, 0.1547, 7)
    FOURTEEN_PERF_ROSETTE = ("FOURTEEN_PERF_ROSETTE", 8 / 3, 47 / 3, 26 * 3**0.5 / pi, 0.1547, 14)
    NINETEEN_PERF_ROSETTE = ("NINETEEN_PERF_ROSETTE", 3, 21, 36 * 3**0.5 / pi, 0.1547, 19)  # rosette prism
    NINETEEN_PERF_CYLINDER = ("NINETEEN_PERF_CYLINDER", 1, 19, 0, 0.3559, 19)
    NINETEEN_PERF_HEXAGON = ("NINETEEN_PERF_HEXAGON", 18 / pi, 19, 18 * (3 * 3**0.5 - 1) / pi, 0.1864, 19)

    NINETEEN_PERF_ROUNDED_HEXAGON = (
        "NINETEEN_PERF_ROUNDED_HEXAGON",
        3**0.5 + 12 / pi,
        19,
        3 - 3**0.5 + 12 * (4 * 3**0.5 - 1) / pi,
        0.1977,
        19,
    )


class SimpleGeometry(Enum):
    """table 1-3 from ref[1] page 23"""

    def __init__(self, desc):
        self.desc = desc

    def __str__(self):
        return self.desc

    SPHERE = "SPHERE"
    ROD = "ROD"  # "Strip / Flake (Rect. Prism)"
    CYLINDER = "CYLINDER"
    TUBE = "TUBE"


GEOMETRIES = {i.desc: i for i in SimpleGeometry}
GEOMETRIES.update({i.desc: i for i in MultPerfGeometry})


class GrainComp:
    """
    detonation velocity and pressure exponent are related to:

    r_dot = u_1 * p^n

    tested velocity for n=1 roughly equals
    Nitrocellulose, 6e-7 - 9e-7 m/(MPa*s)
    Nitroglycerin, 7e-7 - 9e-7 m/(MPa*s)
    """

    allGrainCompDict = {}

    def __init__(
        self,
        name,
        desc,
        propellantForce,
        covolume,
        density,
        redAdbIndex,
        burnVel,
        pressureExp,
        flameTemp,
    ):
        self.name = name
        self.desc = desc
        self.f = propellantForce
        """
        Propellant force is related to the flame temperature
        by:
            f = R * T_v
            R = R_0/M
            R_0 = 8.314... J/(K*mol)
        where T_v is the temperature propellants develop when
        burning in an iso-volume chamber, ignoring losses.
        M is the molar mass of the gas developed. (kg/mol)
        """
        self.alpha = covolume
        self.rho_p = density
        self.theta = redAdbIndex
        self.u_1 = burnVel
        self.n = pressureExp
        self.T_v = flameTemp  # isochoric (const volume) adiabatic temperature
        self.R = self.f / self.T_v

    def __str__(self):
        return self.name

    @classmethod
    def readFile(cls, fileName):
        composition = []
        with open(fileName, newline="") as csvfile:
            propReader = csv.reader(
                csvfile,
                delimiter=",",
                quotechar='"',
                quoting=csv.QUOTE_MINIMAL,
            )

            skipFirstLine = True

            for prop in propReader:
                if skipFirstLine:
                    skipFirstLine = False
                    continue
                (
                    name,
                    desc,
                    adb,
                    density,
                    propellantForce,
                    covolume,
                    pressureExp,
                    burnRateCoe,
                    flameTemp,
                ) = prop

                redAdbIndex = float(adb) - 1

                newComp = cls(
                    name,
                    desc,
                    float(propellantForce),
                    float(covolume),
                    float(density),
                    float(redAdbIndex),
                    float(burnRateCoe),
                    float(pressureExp),
                    float(flameTemp),
                )

                composition.append(newComp)

        compoDict = {i.name: i for i in composition}
        cls.allGrainCompDict.update(compoDict)
        return {i.name: i for i in composition}

    def getLBR(self, p):
        """
        get linear burn rate, given a pressure supplied in Pa

        linear burn rate = (pressure)^self.n * self.u_1
        """
        return self.u_1 * p**self.n

    def getIsp(self, pRatio=None):
        # Tp = Tv/gamma
        if pRatio is None:
            # infinite pressure ratio, basically vacuum
            return (2 * self.f / self.theta) ** 0.5
        else:
            # pressure ratio is interpreted as chamber/exit
            return (2 * self.f / self.theta * (1 - pRatio ** (-self.theta / (self.theta + 1)))) ** 0.5

    @classmethod
    def check(cls):
        from tabulate import tabulate

        line = [
            [
                "Propellant",
                "Adb.Temp",
                "Force 10^3 ft-lb/lb",
                "Adb.Index",
                "Covolume in^3/lbs",
                "density",
                "b.r.coe mm/s/MPa",
                "b.r.exp",
            ]
        ]
        for comp in cls.allGrainCompDict.values():
            line.append(
                [
                    comp.name,
                    comp.T_v,
                    round(comp.f / 2.98907e3),
                    comp.theta + 1,
                    round(comp.alpha * 27680, 2),
                    comp.rho_p,
                    round(comp.u_1 * 1e6**comp.n * 1e3, 3),
                    comp.n,
                ]
            )

        line = [list(l) for l in zip(*line)]

        print(tabulate(line))


class Propellant:
    # assumed to be multi-holed propellants
    def __init__(self, composition, propGeom, R1, R2, fudge=0):
        """
        given e_1 as the maximum burn depth.
        R1: ratio w.r.t arc thickness
            interpreted as 1/alpha:
                Perf diameter to arc thickness ratio
                for perforated cylinders, d_0/(2*e_1)

                Single perf cylinder, d_0/(2*e_1)

                Secondary length to primary length ratio
                for rectangular rod shape (2*b)/(2*e_1)

        R2: ratio w.r.t arc thickness
            interpreted as 1/beta:
                Length to "effective diameter" ratio
                for perforated cylinders, 2*c/(2*e_1)

                tertiary length to primary length ratio
                for rectangular rod shapes, 2*c/(2*e_1)

        fudge: burn rate modifier, used to match experiment.


        Required attributes:
        .maxLF, geometrical volume fraction

        for geometries where burning is single phase
        (Z: 0->1, phi: 0->1)
        .Z_b = 1.0
        .chi
        .labda
        .mu

        for multi-perf propellants:
        (Z: 0->Z_b, phi: 0->1)
        .Z_b > 1.0
        .chi_s
        .labda_s

        """
        self.R1 = R1
        self.R2 = R2

        self.fudge = fudge

        if any((R1 < 0 if R1 is not None else False, R2 < 0 if R1 is not None else False)):
            raise ValueError("Geometry is impossible")

        self.composition = composition
        self.geometry = propGeom

        if self.geometry in MultPerfGeometry:
            arcThick, PR, LR = 1, R1, R2

            e_1 = 0.5 * arcThick
            d_0 = arcThick * PR

            if propGeom == MultPerfGeometry.SEVEN_PERF_CYLINDER:
                b, a = 3 * d_0 + 8 * e_1, 0

            elif propGeom == MultPerfGeometry.SEVEN_PERF_ROSETTE:
                b, a = d_0 + 4 * e_1, d_0 + 2 * e_1

            elif propGeom == MultPerfGeometry.FOURTEEN_PERF_ROSETTE:
                b, a = d_0 + 4 * e_1, d_0 + 2 * e_1

            elif propGeom == MultPerfGeometry.NINETEEN_PERF_ROSETTE:
                b, a = d_0 + 4 * e_1, d_0 + 2 * e_1

            elif propGeom == MultPerfGeometry.NINETEEN_PERF_CYLINDER:
                b, a = 5 * d_0 + 12 * e_1, 0

            elif propGeom == MultPerfGeometry.NINETEEN_PERF_HEXAGON:
                b, a = d_0 + 2 * e_1, d_0 + 2 * e_1

            elif propGeom == MultPerfGeometry.NINETEEN_PERF_ROUNDED_HEXAGON:
                b, a = d_0 + 2 * e_1, d_0 + 2 * e_1

            else:
                raise ValueError("unhandled propellant geometry {}".format(propGeom))

            A, B, C, n = propGeom.A, propGeom.B, propGeom.C, propGeom.nHole
            S_T = 0.25 * pi * (C * a**2 + A * b**2 - B * d_0**2)
            self.maxLF = S_T / (S_T + n * pi * 0.25 * d_0**2)

            D_0 = (C * a**2 + A * b**2 + (n - B) * d_0**2) ** 0.5
            # effective diameter, equals to diameter for perforated cylinders
            grainLength = LR * D_0
            # derive length based on "mean"/"effective" diameter
            c = 0.5 * grainLength

            rho = propGeom.rhoDiv * (e_1 + d_0 / 2)

            Pi_1 = (A * b + B * d_0) / (2 * c)
            Q_1 = (C * a**2 + A * b**2 - B * d_0**2) / (2 * c) ** 2

            """
            maximum load factor, the volume fraction of propellant
            within the geometry of the grain defined. In reality
            this will be a lot lower for any real weapons since
            this both cause undesirably high pressure spike at start
            of shot, and in addition is physically impossible given
            realistic grain packing behaviours
            """

            beta = e_1 / c

            # first phase of burning rate parameters, Z from 0 to 1
            self.chi = (Q_1 + 2 * Pi_1) / Q_1 * beta
            self.labda = (
                (n - 1 - 2 * Pi_1) / (Q_1 + 2 * Pi_1) * beta
            )  # deliberate misspell to prevent issues with python lambda keyword
            self.mu = -(n - 1) / (Q_1 + 2 * Pi_1) * beta**2

            self.Z_b = (e_1 + rho) / e_1  # second phase burning Z upper limit

            psi_s = self.chi * (1 + self.labda + self.mu)

            # second phase of burning rate parameters, Z from 1 to Z_b
            self.chi_s = (1 - psi_s * self.Z_b**2) / (self.Z_b - self.Z_b**2)
            self.labda_s = psi_s / self.chi_s - 1

        elif self.geometry in SimpleGeometry:
            self.Z_b = 1  # this will prevent the running of post fracture code

            if self.geometry == SimpleGeometry.SPHERE:
                self.maxLF = 1
                self.chi = 3
                self.labda = -1
                self.mu = 1 / 3

            elif self.geometry == SimpleGeometry.CYLINDER:
                self.maxLF = 1

                beta = 1 / R2

                self.chi = 2 + beta
                self.labda = -(1 + 2 * beta) / self.chi
                self.mu = beta / self.chi

            elif self.geometry == SimpleGeometry.TUBE:
                beta = 1 / R2

                self.maxLF = 4 * (R1 + 1) / (R1 + 2) ** 2

                self.chi = 1 + beta
                self.labda = -beta / (1 + beta)
                self.mu = 0

            elif self.geometry == SimpleGeometry.ROD:
                self.maxLF = 1
                beta, alpha = sorted((1 / R1, 1 / R2))  # ensure that alpha > beta, ascending order

                self.chi = 1 + alpha + beta
                self.labda = -(alpha + beta + alpha * beta) / self.chi
                self.mu = alpha * beta / self.chi

        else:
            raise ValueError("unhandled propellant geometry {}".format(propGeom))

    def __getattr__(self, attrName):
        if "composition" in vars(self) and not (attrName.startswith("__") and attrName.endswith("__")):
            try:
                if attrName == "u_1":
                    return getattr(self.composition, attrName) * (1 + self.fudge)
                else:
                    return getattr(self.composition, attrName)
            except AttributeError:
                raise AttributeError("%r object has no attribute %r" % (self.__class__.__name__, attrName))
        else:
            raise AttributeError

    def f_sigma_Z(self, Z):
        # is the first derivative of psi(Z)
        if Z <= 1.0:
            return self.chi * (1 + 2 * self.labda * Z + 3 * self.mu * Z**2)
        elif Z <= self.Z_b:
            return 1 + 2 * self.labda_s * Z
        else:
            return 0

    def f_ullim(self):
        if self.Z_b > 1:
            return (
                self.chi * (1 + 2 * self.labda + 3 * self.mu),
                1 + 2 * self.labda_s,
            )
        else:
            return 0

    def f_psi_Z(self, Z):
        if Z <= 1.0:
            return self.chi * Z * (1 + self.labda * Z + self.mu * Z**2)
        elif Z <= self.Z_b:
            return self.chi_s * Z * (1 + self.labda_s * Z)
        else:
            return 1.0


if __name__ == "__main__":
    compositions = GrainComp.readFile("propellants.csv")
    GrainComp.check()
