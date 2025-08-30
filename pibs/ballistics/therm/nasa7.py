"""
Jinpeng Zhai 翟锦鹏
2023 08 06
Contact: 914962409@qq.com

Python code to read the NASA Glenn Thermodynamic database,
generate fit for each specie and calcualte reaction constants
from thermalchemical data.

This author is not affiliated, nor is this code endorsed by any
official agency anywhere.

#The original NASA thermo data file format was documented in:
#
#S. Gordon and B.J. McBride, "Computer Program for Calculation of Complex
#Chemical Equilibrium Composition, Rocket Performance, Incident and Reflected
#Shocks and Chapman-Jouguet Detonations", NASA SP-273 (1971).
#
#The following information is quoted directly from that document.
#Eqs. 90, 91, and 92 give polynomial functions for the specific heat,
#enthalpy, and entropy.  The coefficients are given as data.
#
#Eq. 90: cp0/R=a1+a2T+a3*T^2+a4*T^3+a5*T^4
#Eq. 91: H0_T/(RT)=a1+a2/2*T+a3/3*T^2+a4/4*T^3+a5/5*T^4+a6/T
#Eq. 92: S0_T/R=a1*ln(T)+a2*T+a3/2*T^2+a4/3*T^3+a5/4*T^4+a7
#
#Appendix D: Thermo Data (Format and Listing)
#
#Card             Contents                              Format        Card
#Order                                                               Column
#1      THERMO                                           3A4         1 to 6
#2      Temperature ranges for 2 sets of coefficients:   3F10.3      1 to 30
#        lowest T, common T, and highest T
#3      Species name                                     3A4         1 to 12
#       Date                                             2A3         19 to 24
#       Atomic symbols and formula                       4(A2,F3.0)  25 to 44
#       Phase of species (S, L, or G for solid, liquid,  A1          45
#        or gas, respectively)
#       Temperature range                                2F10.3      46 to 65
#       Integer 1                                        I5          80
#4      Coefficients ai (i=1 to 5) in equations (90)     5(E15.8)    1 to 75
#        to (92) (for upper temperature interval)
#       Integer 2                                        I5          80
#5      Coefficients in equations (90) to (92) (a6, a7,  5(E15.8)    1 to 75
#        for upper temperature interval, and a1, a2, and
#        a3 for lower)
#       Integer 3                                        I5          80
#6      Coefficients in equations (90) to (92) (a4, a5,  4(E15.8)    1 to 60
#        a6, a7 for lower temperature interval)
#       Integer 4                                        I20         80
#(a)    Repeat cards numbered 1 to 4 in cc 80 for each
#       species
#(Final END (Indicates end of thermodynamic data)        3A4         1 to 3
# card)
#
#Gaseous species and condensed species with only one condensed phase can be in
#any order.  However, the sets for two or more condensed phases of the same
#species must be adjacent.  If there are more than two condensed phases of
#a species, their sets must be either in increasing or decreasing order
#according to their temperature intervals.
#
#The following example illustrates a thermo data file with a single species.
#Additional species would be added by duplicating lines 3 through 6:
#
#THERMO
# 200.0   1000.0   6000.0
#Br                J6/82BR  1    0     0    0G 200.000  6000.000                1
# 2.08851053E+00 7.12118611E-04-2.70003073E-07 4.14986299E-11-2.31188294E-15    2
# 1.28568767E+04 9.07351144E+00 2.48571711E+00 1.50647525E-04-5.37267333E-07    3
# 7.20921065E-10-2.50205558E-13 1.27092168E+04 6.86030804E+00 1.34535890E+04    4
#END
#
#Notes:
#1. The minimum and maximum temperatures given on the second line of the
# file are meant to be global min and max values, where each species
# provides a specific, presumably more restrictive, range.  The middle
# temperature on the second line is the common temperature for the high
# and low temperature ranges for all species in the file.
#2. Some files (such as the example, above) provide an additional value
# in the normally unused last field.  This is the enthalpy of formation,
# divided by the gas constant, and is redundant with Eq. 91.
#3. Although the use of columns 66 to 79 of the first line of species data
# is not specified, it is often used for the molecular weight.

The data, as well as the above description, is widely distributed on the internet.
One such source is at
>https://shepherd.caltech.edu/EDL/PublicResources/sdt/thermo.html

Note that in practice, the Phase of species only has to values: G for gas and C for
condensed phase? presumably.
"""

import difflib
import logging
from math import exp, log

from .periodic import molarMasses

# create logger
logger = logging.getLogger("THERMO")
logger.setLevel(logging.INFO)

# create console handler and set level to info
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)


class Specie:
    allFits = {}

    def __init__(
        self,
        cardNo,  # int
        specieName,  # string
        date,  # string, more like a reference
        formula,  # dict:{elementNameinCaps: count}
        phase,  # "G", "C", for gas and condensed phase (liquid + solid)
        fitTlow,  # float
        fitTmid,  # float
        fitThigh,  # float
        coes_l,  # 7 coefficients [a0...a7]
        coes_h,  # 7 coefficients [a0...a7]
        molWt=None,  # None, or float
        Hf__R=None,  # None, or float
    ):
        # publics
        self.cardNo = cardNo
        self.specieName = specieName
        self.date = date
        self.formula = formula
        self.phase = phase
        self.fitTlow = fitTlow
        self.fitThigh = fitThigh
        self.fitTmid = fitTmid

        # privates
        self._coe_l = coes_l
        self._coe_h = coes_h

        if molWt is None:
            molWt = 0
            for element, count in formula.items():
                molWt += molarMasses[element] * count
            logger.info(
                "{:}".format(specieName)
                + " is not provided with optional field molecular weight,"
                + " calculated as M = {:.5f}".format(molWt)
            )
        self.molWt = molWt

        if Hf__R is None:
            Hf__R = self(298.15, R=1, C="p")[1]
            logger.info(
                "{:}".format(specieName)
                + " is not provided with optional field standard enthalpy of formation"
                + ", calculated as ΔΗ⁰_f = {:.5f} * R".format(Hf__R),
            )
        self.Hf__R = Hf__R

    @classmethod
    def get(cls, specieName):
        """
        Fuzzy match the supplies specie name against available fits.
        If the supplies specie name is an exact match for available fits,
        the corresponding one will be returned.
        Otherwise, the fit corresponding to the first result returned by
        difflib.get_close_matches will be returned.
        If no fuzzy matched results can be found, a ValueErorr will be
        raised.
        """
        if specieName in cls.allFits:
            specie = cls.allFits[specieName]
            logger.info("Returning specie {:} [{:}].".format(specieName, specie.phase))
            return specie
        else:
            closeMatch = difflib.get_close_matches(
                specieName,
                list(cls.allFits.keys()),
                n=4,
                cutoff=0.6,
            )
            if len(closeMatch) > 1:
                logger.warning(
                    "Unknown specie name {:}, assuming {:} [{:}], did you mean: {:}.".format(
                        specieName, closeMatch[0], closeMatch[0].phase, ", ".join(closeMatch[1:])
                    )
                )
                return cls.allFits[closeMatch[0]]

            elif len(closeMatch) > 0:
                logger.warning(
                    "Unknown specie name {:}, assuming {:} [{:}].".format(
                        specieName, closeMatch[0], closeMatch[0].phase
                    )
                )
                return cls.allFits[closeMatch[0]]

            else:
                raise ValueError("No close match can be found for {:}".format(specieName))

    @classmethod
    def read(cls, curvePath):
        """
        Read a standard NASA thermo data file as documented in:
            S. Gordon and B.J. McBride, "Computer Program for Calculation of Complex
            Chemical Equilibrium Composition, Rocket Performance, Incident and Reflected
            Shocks and Chapman-Jouguet Detonations", NASA SP-273 (1971).
        """
        with open(curvePath, "r") as file:
            lines = file.readlines()
            """find the start of punched card file,
            marked by the word THERMO, allowing for
            comments to be written directly in the
            data file.

            the word END is checked to find the end
            of the punched card data region. If not
            found the file will be read in its
            entirety
            """
            startLine = None
            endLine = -1
            for lineNo, line in enumerate(lines):
                if line[0:6] == "THERMO":
                    startLine = lineNo
                if line[0:3] == "END":
                    endLine = lineNo

            if startLine is None:
                raise ValueError(
                    "Cannot read file as NASA 7 coefficient data:" + "Card data must start with line 'THERMO'"
                )

            lowT, comT, highT = None, None, None
            try:
                lowT, comT, highT = [
                    float(lines[startLine + 1][9 * i : 9 * (i + 1)]) for i in range(3)
                ]  # read the lowest, common and highest temperature ranges
            except ValueError:
                raise ValueError(
                    "Cannot read file as NASA 7 coefficient data:"
                    + "Temperature range must be specified on first"
                    + "data line."
                )

            lines = lines[startLine + 2 : endLine]  # crop lines to relevant data region
            linesCount = len(lines)

            if linesCount % 4 != 0:  # checks for data integrity
                logger.warning(
                    "read {:} data lines not divisible by 4.".format(linesCount)
                    + " Last datacard will not be read in full!"
                )
            else:
                logger.info("read {:} data lines from file.".format(linesCount))

            cards = len(lines) // 4

            for cardNo in range(cards):
                card = lines[cardNo * 4 : (cardNo + 1) * 4]
                pass

                # first line of card
                specieName = card[0][0:12].strip()
                date = card[0][18:24].strip()
                # formula = card[0][24:44]  # chemical formula, 4 elements
                formula = {}
                for i in range(4):
                    group = card[0][24 + 5 * i : 24 + 5 * (i + 1)]
                    element, count = str(group[:2]).replace(" ", ""), float(group[2:])
                    if element != "":
                        formula.update({element: count})

                phase = card[0][44]
                fitTlow = float(card[0][45:55])
                fitThigh = float(card[0][55:65])

                # the molecular weight is an optional value.
                molWt = None
                if card[0][65:79] != " " * 14:
                    try:
                        molWt = float(card[0][65:79])
                    except ValueError:
                        molWt = float(card[0][65:79].replace("-", "E-"))

                # second to fourth line: coefficients
                coes = [
                    (
                        float(card[j][15 * i : 15 * (i + 1)]) if i != 4 or j != 3 else card[j][15 * i : 15 * (i + 1)]
                    )  # spcial treatment for last optional entry
                    for j in range(1, 4)  # outer loop
                    for i in range(5)  # inner loop
                ]

                """
                last coefficient is optional enthalpy of formation divided
                by gas constant, which is technically redundant.
                """
                Hf__R = float(coes[-1])

                # coefficinets for the higher tempreature range
                coes_h = coes[:7]

                # coefficients for the lower temperature range
                coes_l = coes[7:-1]

                newFit = cls(
                    cardNo=cardNo,
                    specieName=specieName,
                    date=date,
                    formula=formula,
                    phase=phase,
                    fitTlow=fitTlow,
                    fitThigh=fitThigh,
                    fitTmid=comT,
                    molWt=molWt,
                    Hf__R=Hf__R,
                    coes_l=coes_l,
                    coes_h=coes_h,
                )
                cls.allFits.update({specieName: newFit})

        logger.info("Successfully logged {:} data cards".format(cardNo))

    @staticmethod
    def _R(u):
        """Pass through R supplied in value, and
        substitute numerical value when supplied in
        unit string."""
        if u == "J/(molK)":
            return 8.314
        elif u == "cal/(molK)":
            return 1.987

    @classmethod
    def dumpJSON(cls):
        import json

        allSpecies = {}

        for specie in cls.allFits.values():
            allSpecies.update(
                {
                    specie.specieName: {
                        "formula": specie.formula,
                        "phase": specie.phase,
                        "T_low": specie.fitTlow,
                        "T_mid": specie.fitTmid,
                        "T_high": specie.fitThigh,
                        "molWt": specie.molWt,
                        "Hf__R": specie.Hf__R,
                        "coes_l": specie._coe_l,
                        "coes_h": specie._coe_h,
                    }
                }
            )

        with open("nasa7.json", "w") as file:
            json.dump(allSpecies, file, indent="\t", ensure_ascii=False)

        # print(json.dumps(allSpecies, indent="\t", ensure_ascii=False))

    def __str__(self):
        return "Fit for " + self.specieName + " from {:.0f}K to {:.0f}K.".format(self.fitTlow, self.fitThigh)

    def __call__(self, T, R="J/(molK)", C="p"):  # 1.987 for calorie based unit
        """
        retrieve the Cp0, H0_T, and S0_T, with the unit implied by the supplied
        gas constant R.

                    R in  |  heat cap.(Cp0) | enthalpy (H0_T) | entropy (S0_T)
        ------------------+-----------------+-----------------+---------------
        8.314 ... J/mol K |  J/mol K        | J/mol           | J/mol K
        1.987 ... cal/molK| cal/mol K       | cal/mol         | cal/mol K

               T
        h(T) = ∫ Cp dT + h_ref = ΔΗ⁰_f @ 298.15 + (H⁰T - H⁰_298.15)
              T_ref

               T
        s(T) = ∫ 1/T Cp dT + s_ref
              T_ref

        By default the returned values are for constant pressure. The conversion
        to constant volume condition assumes ideal gas relations.

                      C as  |  heat cap.(Cp0) | enthalpy (H0_T) | entropy (S0_T)
        --------------------+-----------------+-----------------+---------------
        Constant pressure Cp|  as calculated  | as calculated   | as calculated
          Constant Volume Cv|       -R        |     -RT         |  - R log(T)
        """

        R = self._R(R)

        if T < self.fitTlow:
            raise ValueError("Temperature below fit minimum ({:.0f}K) for this specie.".format(self.fitTlow))

        elif T > self.fitThigh:
            raise ValueError("Temperature above fit maximum ({:.0f}K) for this specie.".format(self.fitThigh))

        if T < self.fitTmid:  # low range
            a1, a2, a3, a4, a5, a6, a7 = self._coe_l
        elif T < self.fitThigh:  # high range
            a1, a2, a3, a4, a5, a6, a7 = self._coe_h
        else:
            a1, a2, a3, a4, a5, a6, a7 = [None] * 7

        Cp0__R = a1 + a2 * T + a3 * T**2 + a4 * T**3 + a5 * T**4
        H0_T__RT = a1 + a2 / 2 * T + a3 / 3 * T**2 + a4 / 4 * T**3 + a5 / 5 * T**4 + a6 / T
        S0_T__R = a1 * log(T) + a2 * T + a3 / 2 * T**2 + a4 / 3 * T**3 + a5 / 4 * T**4 + a7
        if C == "p":
            return Cp0__R * R, H0_T__RT * R * T, S0_T__R * R

        elif C == "v":  # converted from Cp values assuming ideal gas relation!
            return (Cp0__R - 1) * R, (H0_T__RT - 1) * R * T, (S0_T__R - log(T)) * R
        else:
            raise ValueError("Unknown constant condition.")

    def getMMH(self, T, R="cal/(molK)", C="v", ref=300):
        """
        Mean molecular heat is defined as:
                    1     T
        MMH(T) = -------- ∫ [Cp or Cv] dT
                (T-T_ref) T_ref

        """
        if T == ref:
            return 0
        else:
            return (self(T, R, C)[1] - self(ref, R, C)[1]) / (T - ref)

    def Cv(self, T, R="J/(molK)"):
        """
        returns the constant volume specific heat
        """
        return self(T, R, C="v")[0]

    def Cp(self, T, R="J/(molK)"):
        """
        returns the constant pressure specific heat
        """
        return self(T, R, C="p")[0]

    def H(self, T, R="J/(molK)", C="p"):
        """
        returns the enthalpy
        """
        return self(T, R, C)[1]

    def S(self, T, R="J/(molK)", C="p"):
        """
        returns the entropy
        """
        return self(T, R, C)[2]

    def G(self, T, R="J/(molK)"):
        """
        Gibbs free energy, either in j or cal depending on
        unit of R.
        """
        C = "p"
        return self.H(T, R, C) - T * self.S(T, R, C)


class Reaction:
    def __init__(self, name, LHS, RHS):
        # check that the reaction is balanced
        self.name = name
        element_balance = {}
        for specie, count in LHS.items():  # Specie object, float
            for element, number in specie.formula.items():  # str, float
                if element in element_balance:
                    element_balance[element] += number * count
                else:
                    element_balance.update({element: number * count})
        for specie, count in LHS.items():  # Specie object, float
            for element, number in specie.formula.items():  # str, float
                if element in element_balance:
                    element_balance[element] -= number * count
                else:
                    raise ValueError("Unknwon element {:}".format(element) + "on the RHS.")
        if any((value != 0 for element, value in element_balance.items())):
            raise ValueError("Unbalanced reaction")
        else:
            self.LHS = LHS
            self.RHS = RHS
            logger.info("Reaction {:} checked as balanced".format(self))

    def __str__(self):
        reaction = self.name + ": "

        for i, (therm, count) in enumerate(self.LHS.items()):
            if i != 0:
                reaction += " + {:} {:}".format(count, therm.specieName)
            else:
                reaction += "{:} {:}".format(count, therm.specieName)

        reaction += " <-> "

        for i, (therm, count) in enumerate(self.RHS.items()):
            if i != 0:
                reaction += " + {:} {:}".format(count, therm.specieName)
            else:
                reaction += "{:} {:}".format(count, therm.specieName)

        return reaction

    def __call__(self, T, R="J/(molK)"):
        """
        Retrieve the equilibrium constant
        R in     | kp in |
        J/molK   |   1   |
        cal/molK |   1   |
        """
        Delta_G = 0  # Gibbs free energy for constant pressure case.
        for therm, count in self.LHS.items():
            Delta_G -= therm.G(T, R) * count
        for therm, count in self.RHS.items():
            Delta_G += therm.G(T, R) * count

        R = Specie._R(R)

        # return 10 ** (-Delta_G / (2.3026 * R * T))
        return exp(-Delta_G / (R * T))

    def Kp(self, T):
        """
        retrieve the constant pressure equilibrium constant
        """
        return self(T)

    def Hr(self, R="J/(molK)"):
        """
        retrieve the standard enthalpy of reaction:
        R in     |        DeltaH in
        J/molK   |   J/mol  (of reaction)
        cal/molK |   cal/mol (of reaction)
        """
        ref = 298.15
        C = "p"
        Delta_H = 0
        for therm, count in self.LHS.items():
            Delta_H -= therm.H(ref, R, C) * count
        for therm, count in self.RHS.items():
            Delta_H += therm.H(ref, R, C) * count

        return Delta_H


def main():
    # some simple examples
    Specie.read("/ballistics/resource/nasa7.dat")
    methane = Specie.get("CH4")
    ammonia = Specie.get("NH3")

    hydrogen = Specie.get("H2")
    nitrogen = Specie.get("N2")

    carbon = Specie.get("C(gr)")

    methaneDecomposition = Reaction(name="Methane decomposition", LHS={methane: 1}, RHS={carbon: 1, hydrogen: 2})

    ammoniaDecomposition = Reaction(
        name="Ammonia decomposition", RHS={ammonia: 1}, LHS={hydrogen: 1.5, nitrogen: 0.5}
    )  # this is k8

    for T in [2000, 2500, 3000, 3500, 4000]:
        print(T, nitrogen.Cp(T) / nitrogen.Cv(T))

    Specie.dumpJSON()


if __name__ == "__main__":
    main()
