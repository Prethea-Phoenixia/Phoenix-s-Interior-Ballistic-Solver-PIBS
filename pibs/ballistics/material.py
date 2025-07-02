from enum import Enum


def _30SIMN2MOVA_Y(T):
    return (
        -0.047 * (T / 293) ** 4 + 0.201 * (T / 293) ** 3 - 0.15 * (T / 293) ** 2 - 0.376 * (T / 293) + 1.375
    ) * 1195e6  # Circumferential Yield Strength in Pa


"""
def _30SIMN2MOVA_E(T):
    return (
        -0.025 * (T / 293) ** 2 + 0.04 * (T / 293) + 0.985
    ) * 213e9  # Youngs Modulus in Pa,
"""


class Material(Enum):
    def __init__(self, desc, ref, rho, Ys, Ts):
        self.desc = desc
        self.ref = ref
        self.rho = rho  # density, kg/m^3
        self.Ys = Ys  # Circumferential Yield Strength in Pa
        self.Ts = Ts  # Temperature condition of each entry.

    def __str__(self):
        return self.desc

    material_30SIMN2MOVA = (
        "30SiMn2MoVA",
        "Cao, Xu, (2017); Gu et al., (2018); Rp0.2",
        7801,
        [_30SIMN2MOVA_Y(293 + i * 100) for i in range(7)],
        [293 + i * 100 for i in range(7)],
    )

    material_38644Ti = (
        "38644Ti",
        "ADA196329; 0.1%",
        4820,
        [1056e6, 945e6],
        [294, 273 + 316],
    )

    material_A723 = (
        "ASTM A723 Grade 2",
        "ADA196329; 0.1%",
        7850,  # generic value
        [1050e6, 918e6],
        [294, 273 + 316],
    )

    def getTdict(self):
        return {f"{T:.0f} K": f"{T:.0f} K" for i, T in enumerate(self.Ts)}

    def createMaterialAtTemp(self, Tstr):
        for i, T in enumerate(self.Ts):
            if Tstr == f"{T:.0f} K":
                return MaterialAtTemp(
                    self.desc + "@{:.0f} K".format(self.Ts[i]),
                    self.rho,
                    self.Ys[i],
                )
        raise ValueError(f"Material has no condition at {Tstr:}")


class MaterialAtTemp:
    def __init__(self, desc, rho, Y):
        self.desc = desc
        self.Y = Y
        self.rho = rho


MATERIALS = {i.desc: i for i in Material}


if __name__ == "__main__":
    print(Material.material_30SIMN2MOVA.getTdict())
