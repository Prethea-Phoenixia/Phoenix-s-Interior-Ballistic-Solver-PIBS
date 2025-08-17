from pibs.ballistics.therm import Ingredient, Mixture


def main():
    NC1260 = Ingredient.getLine(683)
    RDX = Ingredient.getLine(847)

    EC = Ingredient(
        name="Ethyl Centralite",
        elements={"C": 17, "H": 20, "O": 1, "N": 2},
        rho=1.140,
        rho_u="g/cc",
        Hf=-391.5,
        Hf_u="J/g",
    )

    ATEC = Ingredient(
        name="Acetyl triethyl citrate",
        elements={"C": 14, "H": 22, "O": 8},
        rho=1.136,
        rho_u="g/cc",
        # Hf=-1257,
        Hf=-5459.6,
        Hf_u="J/g",
    )

    CAB = Ingredient(
        "Cellulose Acetate Butyrate",
        elements={"C": 15, "H": 22, "O": 8},
        Hf=-4933.76,
        Hf_u="J/g",
        rho=1.220,
        rho_u="g/cc",
    )

    BDNPA = Ingredient.getLine(189)
    BDNPF = Ingredient.getLine(190)

    XM39 = Mixture(
        "XM39",
        compoDict={RDX: 76, CAB: 12, NC1260: 4, ATEC: 7.6, EC: 0.4},
    )

    XM39.prettyPrint()

    M43 = Mixture(
        name="M43",
        compoDict={
            RDX: 76,
            CAB: 12,
            NC1260: 4,
            BDNPA: 7.6,
            # BDNPF: 7.6,
            EC: 0.4,
        },
        Delta=0.2,
    )
    M43.prettyPrint()

    NG = Ingredient.getLine(693)

    MeNENA = Ingredient(
        "Methyl-NENA",
        elements={"C": 3, "H": 7, "N": 3, "O": 5},
        Hf=-106.50,  # Burcat, 2015
        Hf_u="kJ/mol",
        rho=1.53,  # a.la ADA377866
        rho_u="g/cc",
    )

    EtNENA = Ingredient(
        "Ethyl-NENA",
        elements={"C": 4, "H": 9, "N": 3, "O": 5},
        Hf=-133.90,  # Burcat, 2015
        Hf_u="kJ/mol",
        rho=1.32,  # a.la ADA377866
        rho_u="g/cc",
    )

    ATKPRD20 = Mixture(
        name="ATK PRD20", compoDict={NC1260: 41.90, RDX: 25.71, MeNENA: 14.00, EtNENA: 10.00, NG: 7.69}, Delta=0.2
    )

    ATKPRD20.prettyPrint()

    ATKPRDS21 = Mixture(
        name="ATK PRD(S)21", compoDict={NC1260: 36.48, RDX: 30.33, MeNENA: 13.44, EtNENA: 9.57, NG: 9.46}, Delta=0.2
    )

    ATKPRDS21.prettyPrint()

    ATKPRDS22 = Mixture(
        name="ATK PRD(S)22", compoDict={NC1260: 31.11, RDX: 34.08, MeNENA: 12.57, EtNENA: 8.94, NG: 12.58}, Delta=0.2
    )

    ATKPRDS22.prettyPrint()

    # import matplotlib.pyplot as plt
    # from labellines import labelLines

    # fig, ax = plt.subplots(1, 1)
    #
    # speciesDict = {}
    #
    # print("HERE")
    #
    # Ts = list(range(1600, 4000, 100))
    # for T in Ts:
    #     speciesList = ATKPRDS22.balanceAt(T, False)
    #
    #     for entry in speciesList:
    #         specie, pop = entry[0], entry[1]
    #         if specie in speciesDict:
    #             speciesDict[specie].append(pop)
    #         else:
    #             speciesDict.update({specie: [pop]})
    #
    # for specie, pops in speciesDict.items():
    #     ax.plot(Ts, pops, label=specie)
    #
    # # ax.set_yscale("log")
    # ax.set_ylim((1e-6, 1))
    # labelLines(ax.get_lines())
    #
    # plt.show()

    NG = Ingredient.getLine(693)
    NC1316 = Ingredient.nitrocellulose(0.1316)
    NC1315 = Ingredient.nitrocellulose(0.1315)
    NC1311 = Ingredient.nitrocellulose(0.1311)
    DBP = Ingredient.find("DIBUTYL PHTHALATE")
    DNT = Ingredient.find("DINITRO TOLUENE")
    DPA = Ingredient.find("2 NITRO DIPHENYL AMINE")

    KNO3 = Ingredient.find("POTASSIUM NITRATE")
    K2SO4 = Ingredient.find("POTASSIUM SULFATE")

    ETC = Ingredient(
        name="Ethyl Centralite",
        elements={"C": 17, "H": 20, "O": 1, "N": 2},
        rho=1.140,
        rho_u="g/cc",
        Hf=-391.5,
        Hf_u="J/g",
    )
    MEC = Ingredient(
        name="Methyl Centralite",
        elements={"C": 15, "H": 16, "O": 1, "N": 2},
        rho=1.2,
        rho_u="g/cc",
        Hf=-73.1,
        Hf_u="kJ/mol",
    )

    WC846 = Mixture(name="WC846", compoDict={NC1316: 81.40, NG: 10.39, DBP: 5.61, DPA: 0.97, DNT: 0.06}, Delta=0.2)

    WC846.prettyPrint()

    WC870 = Mixture(
        name="WC870", compoDict={NC1311: 79.70, NG: 9.94, DBP: 5.68, DPA: 0.95, KNO3: 0.78, DNT: 0.69}, Delta=0.2
    )
    WC870.prettyPrint()

    IMR4227 = Mixture(name="IMR4227", compoDict={NC1315: 89.72, DNT: 6.74, DPA: 0.84, K2SO4: 0.73}, Delta=0.2)
    IMR4227.prettyPrint()

    IMR4350 = Mixture(name="IMR4350", compoDict={NC1315: 91.56, DNT: 4.80, DPA: 0.85, K2SO4: 0.73}, Delta=0.2)
    IMR4350.prettyPrint()

    IMR8138M = Mixture(name="IMR8138M", compoDict={NC1315: 92.96, ETC: 3.60, DPA: 0.89, K2SO4: 0.69}, Delta=0.2)
    IMR8138M.prettyPrint()

    IMR4895 = Mixture(name="IMR4895", compoDict={NC1315: 88.99, DNT: 7.55, DPA: 0.78, K2SO4: 0.93}, Delta=0.2)
    IMR4895.prettyPrint()

    IMR5010 = Mixture(name="IMR5010", compoDict={NC1315: 88.27, DNT: 8.76, K2SO4: 0.64, DPA: 0.55}, Delta=0.2)
    IMR5010.prettyPrint()

    CMR160 = Mixture(name="CMR160", compoDict={NC1315: 95.77, MEC: 2.43, DPA: 0.73, K2SO4: 0.84})
    CMR160.prettyPrint()
    #
    # P5BPfl = Mixture(name="P5BPFl", compoDict={NC1315: 95.77})


"""
    for T in (2855, 2782, 2878, 2927, 2816, 2874, 2869, 2776):
        print(f"{estU8053(T):.5g}")
"""

if __name__ == "__main__":
    main()
