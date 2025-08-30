import difflib
from os import path

from ..num import secant
from .corner import balance
from .periodic import molarMasses


class Ingredient:
    allIngr = {}
    lastLine = 0
    lineIngr = {}

    def __init__(self, name, elements, Hf, rho, rho_u="lb/cu.in", Hf_u="cal/g", flag="", lineNo=None):
        if lineNo is not None:
            if lineNo > Ingredient.lastLine:
                Ingredient.lastLine = lineNo
            else:
                if lineNo not in Ingredient.lineIngr:
                    pass
                else:
                    raise ValueError("Line Number Collision")
        else:
            lineNo = Ingredient.lastLine + 1
            Ingredient.lastLine += 1

        self.lineNo = lineNo
        self.name = name
        self.flag = flag
        self.elements = elements

        if rho_u == "lb/cu.in":
            self.rho = rho * 27.680  # to g/cc
        elif rho_u == "g/cc":
            self.rho = rho
        else:
            raise ValueError("Unknown Unit For Density")

        A = 0
        for element, num in self.elements.items():
            A += molarMasses[element] * num

        if "C" in self.elements:
            self.Ci = self.elements["C"] / A
        else:
            self.Ci = 0
        if "H" in self.elements:
            self.Hi = self.elements["H"] / A
        else:
            self.Hi = 0
        if "O" in self.elements:
            self.Oi = self.elements["O"] / A
        else:
            self.Oi = 0
        if "N" in self.elements:
            self.Ni = self.elements["N"] / A
        else:
            self.Ni = 0

        self.A = A  # g/mol

        if Hf_u == "cal/g":
            self.Hf = Hf
        elif Hf_u == "cal/mol":
            self.Hf = Hf / A
        elif Hf_u == "J/g":
            self.Hf = Hf / 4.184
        elif Hf_u == "J/mol":
            self.Hf = Hf / (4.184 * A)
        elif Hf_u == "kJ/mol":
            self.Hf = Hf * 1000 / (4.184 * A)
        else:
            raise ValueError("Unknown Enthalpy Unit")

    @classmethod
    def readFile(cls, fileName=None):
        if fileName is None:
            fileName = path.join(path.dirname(__file__), "PEPCODED.DAF")
        # read data from PEP database
        with open(fileName, "r", encoding="ascii") as file:
            fileIngr = []
            for line in file:
                # print(line, end="")
                flag = line[0:2].strip()
                lineNo = int(line[2:9])

                if flag == "*":  # this line is a comment
                    pass
                elif flag == "+":
                    # multiline entry of ingredient name
                    if len(fileIngr) > 0:
                        newIngr = fileIngr[-1]
                        name = line[9:80].strip()
                        newIngr.name = newIngr.name + name
                else:
                    name = line[9:39].strip()
                    elements = dict()
                    for i in range(6):
                        nbr, element = (
                            int(line[39 + i * 5 : 42 + i * 5].strip()),
                            line[42 + i * 5 : 44 + i * 5].strip(),
                        )
                        if element != "":
                            elements.update({element: nbr})

                    Hf = int(line[69:74])
                    rho = float(line[74:80])

                    newIngr = Ingredient(name=name, elements=elements, Hf=Hf, rho=rho, flag=flag, lineNo=lineNo)

                    fileIngr.append(newIngr)

        for ingr in fileIngr:
            cls.allIngr.update({ingr.name: ingr})
            cls.lineIngr.update({ingr.lineNo: ingr})

    @classmethod
    def find(cls, name):
        closeIngrs = difflib.get_close_matches(
            name,
            list(cls.allIngr.keys()),
            n=3,
            cutoff=0.6,
        )
        n = 0
        thisIngr = None
        ingrs = []

        if len(closeIngrs) > 0:
            print("Found candidate:")
            for iname in closeIngrs:
                # if iname in cls.allIngr:
                ingr = cls.allIngr[iname]
                if ingr not in ingrs:
                    ingrs.append(ingr)
                if n == 0:
                    thisIngr = ingr
                n += 1

            for ingr in ingrs:
                print("-" + ingr.name)

            print("returning " + thisIngr.name + "\n")
            return thisIngr
        else:
            print('Unknown ingredient description "{:}"'.format(name) + "\n")
            return None

    @classmethod
    def getLine(cls, lineNo):
        if lineNo in cls.lineIngr:
            print("Returning line {:} : {:}".format(lineNo, cls.lineIngr[lineNo].name))
            return cls.lineIngr[lineNo]
        else:
            print("No such line as {:}\n".format(lineNo))
            return None

    @classmethod
    def nitrocellulose(cls, nitration):
        y = nitration * 100
        x = 162.14 * y / (1400.8 - 45 * y)
        elements = {"C": 6, "H": 10 - x, "O": 5 + 2 * x, "N": x}
        # see hunt SS 2.02

        return cls(
            name="Nitrocellulose ({:.2f}% N)".format(y),
            elements=elements,
            rho=0.0560,  # copied directly from database.
            Hf=-1417.029 + 6318.3 * nitration,  # fit from Tab.3 in R.S.Jessup & E.J.Prosen
        )


class Mixture:
    def __init__(self, name, compoDict, Delta=0.2, tol_z=1e-3, tol_b=1e-9, its=200):
        self.name = name
        self.Delta = Delta  # load density in g/cc
        self.tol_z = tol_z  # tolerance for zeta
        self.its = its
        self.tol_b = tol_b

        # Normalize the given composition such that the fractions sums to 1

        total = 0
        for ingr, fraction in compoDict.items():
            total += fraction

        self.compoDict = {ingr: fraction / total for ingr, fraction in compoDict.items()}

        # tally the relevant factors according to their mass fraction

        invRho = 0
        Ci, Hi, Ni, Oi = 0, 0, 0, 0
        Hf = 0

        for ingr, fraction in self.compoDict.items():
            if ingr.rho == 0:
                raise ValueError("{:} is not provided with density data".format(ingr.name))
            invRho += fraction / ingr.rho
            Ci += fraction * ingr.Ci  # mol/g
            Hi += fraction * ingr.Hi
            Oi += fraction * ingr.Oi
            Ni += fraction * ingr.Ni
            Hf += fraction * ingr.Hf

        self.rho = 1 / invRho
        self.Hf = Hf

        def f(T):
            zeta, _, _, _, _, _, _ = balance(self.Hf, T, Ci, Hi, Oi, Ni, V=1 / Delta, its=its, tol=tol_b)

            return zeta

        # zeta on the order of 0.5 per degree
        Tv = 0.5 * sum(secant(f, 3200, 3000, x_min=1600, x_max=4000, x_tol=tol_z))

        _, self.speciesList, n, E, self.b, self.p, self.f = balance(
            self.Hf, Tv, Ci, Hi, Oi, Ni, V=1 / Delta, its=its, tol=tol_b
        )
        # see Corner ss 3.4
        _, _, _, E1, _, _, _ = balance(self.Hf, Tv, Ci, Hi, Oi, Ni, V=1 / Delta, its=its, tol=tol_b)
        _, _, _, E2, _, _, _ = balance(self.Hf, 0.7 * Tv, Ci, Hi, Oi, Ni, V=1 / Delta, its=its, tol=tol_b)
        C_v = (E1 - E2) / (0.3 * Tv)
        # gas constant: 1.987 cal/(mol K)
        self.gamma = (n * 1.987 / C_v) + 1
        self.n = n

        self.Ci, self.Hi, self.Oi, self.Ni = Ci, Hi, Oi, Ni

        self.Tv = Tv
        self.Hf = Hf

    def balanceAt(self, T, verbose=True, its=200, tol=1e-9):
        Delta, speciesList, n, E, b, p, f = balance(
            self.Hf, T, self.Ci, self.Hi, self.Oi, self.Ni, V=1 / self.Delta, its=its, tol=tol
        )
        if verbose:
            print("Mixture: {:} At: {:}K".format(self.name, T))

            print(" @ Product  %mass  mol/g")
            print(
                *[
                    "{:>2} : {:^6} {:<6.1%} {:<6.4f}".format(i, name, mass, num)
                    for i, (name, mass, num) in enumerate(speciesList)
                ],
                sep="\n",
            )
            print("Average.Mol.Weight : {:>6.4g} g/mol Δ={:>6.1%}".format(1 / n, (1 / n - 1 / self.n) / (1 / self.n)))
            print("Covolume           : {:>6.4g} cc/g  Δ={:>6.1%}".format(b, (b - self.b) / self.b))

            C_v = E / (T - 300)
            gamma = (n * 1.987 / C_v) + 1
            print("Adiabatic index (T): {:>6.4g}".format(gamma))

        return speciesList

    def prettyPrint(self):
        C, H, O, N = (
            self.Ci * molarMasses["C"],
            self.Hi * molarMasses["H"],
            self.Oi * molarMasses["O"],
            self.Ni * molarMasses["N"],
        )
        print("Mixture: {:}".format(self.name))
        print("Specified Composition:---------------------------")
        for ingr, fraction in self.compoDict.items():
            print("--{:-<30}, {:>6.2%}".format(ingr.name, fraction))

        print("")
        print("Elemental Fractions:-----------------------------")
        print("C {:.2%} H {:.2%} N {:.2%} O {:.2%}".format(C, H, N, O))
        print("")
        print("Calculated Properties:---------------------------")
        print("Density            : {:>8.4g} g/cc".format(self.rho))
        print("Heat of Formation  : {:>8.0f} cal/g".format(self.Hf))
        print("Flame Temperature  : {:>8.6g} K (Isochoric Adiabatic)".format(self.Tv))
        print(" @ Product  %mass  mol/g")
        print(
            *[
                "{:>2} : {:^6} {:<6.1%} {:<6.4f}".format(i, name, mass, num)
                for i, (name, mass, num) in enumerate(self.speciesList)
            ],
            sep="\n",
        )
        print("Sp.Impetus (Force) : {:>8.5g} J/g".format(self.f))
        print("Covolume           : {:>8.4g} cc/g".format(self.b))
        # print(" @Temperature      : {:>6.0f} K".format(self.Tv))
        print(" @ Load Density    : {:>8.6g} g/cc".format(self.Delta))
        print(" @ Pressure        : {:>8.1f} MPa".format(self.p))
        print("avg Adb. index     : {:>8.5g}".format(self.gamma))
        print("")
        print("Estimated Burnrates:-----------------------------")
        b = -0.8340 + 8.3956e-4 * self.Tv
        print("after Ludwig Stiefel (ADA086093):")
        print(f"r <mm/s> = {b:.3f} (P <MPa>) ^ 0.8053")
        print(f"r <m/s> = {b * 1e-3 / 1e6**0.8053:.4g} (P <Pa>) ^ 0.8053")
        print("")


# def estU1(nitration, ng, cen, dbp, dnt, vas):
#     """
#     Estimate the burn velocity assuming a linear realtionship
#     between pressure and burn rate, i.e. u = u1 * P

#     nitration:
#         N, nitration % of nitrocellulose
#         硝化棉中氮的百分含量
#     ng:
#         H, percentage of nitroglycerin
#         硝化甘油的百分含量

#     cent:
#         E, percentage of centralite
#         中定剂的百分含量
#     dbp:
#         F, percentage of dibutyl phthalate
#         苯二甲酸二丁酯的百分含量
#     dnt:
#         G, percentage of dinitrotoluene
#         二硝基甲苯的百分含量
#     vaseline:
#         I, 凡士林的百分含量

#     The result is converted into m/s/Pa
#     """
#     N, H, E, F, G, I = nitration, ng, cen, dbp, dnt, vas

#     return (
#         10
#         ** (
#             0.7838
#             - 9
#             + 0.1366 * (N - 11.8)
#             + 0.008652 * H
#             - 0.02620 * E
#             - 0.02235 * F
#             - 0.007355 * G
#             - 0.03447 * I
#         )  # m^3/s-kg
#         / 9.8
#     )


# def estU8053(Tv):
#     return (-0.8340 + 8.3956e-4 * Tv) * 1e-3 / 1e6**0.8053

Ingredient.readFile()
