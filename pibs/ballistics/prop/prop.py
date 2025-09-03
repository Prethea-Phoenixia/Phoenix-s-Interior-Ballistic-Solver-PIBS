"""
GEOMETRIES:
    lookup dictionary mapping geometrical descriptions
    to relevant geometrical object

"""

import csv
from abc import ABC, ABCMeta, abstractmethod
from enum import Enum, EnumType
from math import pi


class Geometry(ABC):
    all_geometries = list()

    def __init__(self):
        Geometry.all_geometries.append(self)

    @abstractmethod
    def get_form_function_coefficients(self, r1: float, r2: float) -> tuple[float, float, float, float, float, float]:
        """
        yield the
            chi, labda, mu, chi_s, labda_s, Z_b
        coefficients
        """
        ...

    @staticmethod
    def get_desc_geometry_dict() -> dict:
        desc_geometry_dict = {}
        for geometry in Geometry.all_geometries:
            desc_geometry_dict[geometry.desc] = geometry

        return desc_geometry_dict


class ABCEnum(ABCMeta, EnumType):
    pass


class SimpleGeometry(Geometry, Enum, metaclass=ABCEnum):
    """table 1-3 from ref[1] page 23"""

    def __init__(self, desc: str, alpha: int, beta: int):
        super().__init__()
        self.desc = desc
        self.alpha = alpha
        self.beta = beta

    def get_form_function_coefficients(self, r1: float, r2: float) -> tuple[float, float, float, float, float, float]:
        """
        alpha = e1 / b
        beta = e1 / c
        """

        if self.alpha >= 0:
            alpha = self.alpha
        else:
            alpha = 1.0 / r1

        if self.beta >= 0:
            beta = self.beta
        else:
            beta = 1.0 / r2

        if (alpha != 0 and not (1.0 >= alpha >= beta)) or (alpha == 0 and not (1.0 >= beta)):
            raise ValueError("Specified geometry is impossible.")

        chi = 1 + alpha + beta
        labda = -(alpha + beta + alpha * beta) / (1 + alpha + beta)
        mu = alpha * beta / (1.0 + alpha + beta)

        return chi, labda, mu, 0.0, 0.0, 1.0

    def __str__(self):
        return self.desc

    SPHERE = ("SPHERE", 1, 1)  # 立方体, alpha = beta = 1
    CYLINDER = ("CYLINDER", 1, -1)  # 方棍, alpha = 1, 1 > beta
    TUBE = ("TUBE", 0, -1)  # 管状, alpha = 0, 1 > beta
    STRIP = ("STRIP", -1, -1)  # 带状, 1 > alpha >= beta


class MultPerfGeometry(Geometry, Enum, metaclass=ABCEnum):
    """table 1-4 from ref[1] page 33"""

    def __init__(
        self,
        desc: str,
        A: float,
        B: float,
        C: float,
        rho_div: float,
        n_hole: int,
        a_factors: tuple[float, float],
        b_factors: tuple[float, float],
    ):

        super().__init__()
        self.desc = desc
        self.A, self.B, self.C = A, B, C
        self.rho_div, self.n_hole, self.a_factors, self.b_factors = rho_div, n_hole, a_factors, b_factors

    def __str__(self):
        return self.desc

    def get_form_function_coefficients(self, r1: float, r2: float) -> tuple[float, float, float, float, float, float]:
        """
        Parameters:
        -----------
        r1: float
            perforation diameter/web. d_0 / (2e_1)
        r2: float
            length/diameter (width) of grain.
        """
        web = 1
        e_1, d_0 = web / 2, web * r1
        a = sum(p * q for p, q in zip(self.a_factors, (d_0, e_1)))
        b = sum(p * q for p, q in zip(self.b_factors, (d_0, e_1)))

        A, B, C, n = self.A, self.B, self.C, self.n_hole
        c = 0.5 * r2 * (C * a**2 + A * b**2 + (n - B) * d_0**2) ** 0.5

        rho = self.rho_div * (e_1 + d_0 / 2)

        p = (A * b + B * d_0) / (2 * c)
        q = (C * a**2 + A * b**2 - B * d_0**2) / (2 * c) ** 2
        beta = e_1 / c

        # first phase of burning rate parameters, Z from 0 to 1
        chi, labda, mu = (q + 2 * p) / q * beta, (n - 1 - 2 * p) / (q + 2 * p) * beta, -(n - 1) / (q + 2 * p) * beta**2

        z_b = (e_1 + rho) / e_1  # second phase burning Z upper limit
        psi_s = chi * (1 + labda + mu)
        # second phase of burning rate parameters, Z from 1 to z_b
        chi_s = (1 - psi_s * z_b**2) / (z_b - z_b**2)
        labda_s = psi_s / chi_s - 1

        return chi, labda, mu, chi_s, labda_s, z_b

    # fmt: off
    SEVEN_PERF_CYLINDER = ("SEVEN_PERF_CYLINDER", 1, 7, 0, 0.2956, 7, (0, 0), (3, 8))
    SEVEN_PERF_ROSETTE = ("SEVEN_PERF_ROSETTE", 2, 8, 12 * 3**0.5 / pi, 0.1547, 7, (1, 2), (1, 4))
    FOURTEEN_PERF_ROSETTE = ("FOURTEEN_PERF_ROSETTE", 8 / 3, 47 / 3, 26 * 3**0.5 / pi, 0.1547, 14, (1, 2), (1, 4))
    NINETEEN_PERF_ROSETTE = ("NINETEEN_PERF_ROSETTE", 3, 21, 36 * 3**0.5 / pi, 0.1547, 19, (1, 2), (1, 4))
    NINETEEN_PERF_CYLINDER = ("NINETEEN_PERF_CYLINDER", 1, 19, 0, 0.3559, 19, (0, 0), (5, 12))
    NINETEEN_PERF_HEXAGON = ("NINETEEN_PERF_HEXAGON", 18 / pi, 19, 18 * (3 * 3**0.5 - 1) / pi, 0.1864, 19, (1, 2), (1, 2))
    NINETEEN_PERF_ROUNDED_HEXAGON = ("NINETEEN_PERF_ROUNDED_HEXAGON", 3**0.5 + 12 / pi, 19, 3 - 3**0.5 + 12 * (4 * 3**0.5 - 1) / pi, 0.1977, 19, (1, 2), (1, 2))
    # fmt: on


class Composition:

    all_compositions = []

    def __init__(
        self,
        name: str,
        desc: str,
        force: float,
        covolume: float,
        density: float,
        reduced_adiabatic_index: float,
        burn_rate_coefficient: float,
        pressure_exponent: float,
        adiabatic_flame_temperature: float = 0,
    ):
        self.name, self.desc = name, desc
        self.f = force
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
        self.theta = reduced_adiabatic_index
        self.u_1 = burn_rate_coefficient
        self.n = pressure_exponent
        self.T_v = adiabatic_flame_temperature  # isochoric (const volume) adiabatic temperature

    def __str__(self):
        return self.name

    @staticmethod
    def read_file(file_name: str):
        compositions = []
        with open(file_name, newline="", encoding="utf-8") as csvfile:
            sniffer = csv.Sniffer()
            content = csvfile.read()
            dialect = sniffer.sniff(content)
            has_header = sniffer.has_header(content)
            csvfile.seek(0)

            prop_reader = csv.reader(csvfile, dialect=dialect)
            for prop in prop_reader:
                if has_header:
                    has_header = False
                    continue
                name, desc, adb, density, force, covolume, pressure_exp, burn_rate_coef, flame_temp = prop

                reduced_adiabatic_index = float(adb) - 1

                new_comp = Composition(
                    name,
                    str(desc).replace("\\n", "\n"),
                    float(force),
                    float(covolume),
                    float(density),
                    float(reduced_adiabatic_index),
                    float(burn_rate_coef),
                    float(pressure_exp),
                    float(flame_temp) if flame_temp else None,
                )
                compositions.append(new_comp)

        Composition.all_compositions = compositions
        return Composition.get_name_composition_dict()

    def get_lbr(self, p):
        """
        get linear burn rate, given a pressure supplied in Pa

        linear burn rate = (pressure)^self.n * self.u_1
        """
        return self.u_1 * p**self.n

    def get_isp(self, p_ratio=None):
        # Tp = Tv/gamma
        if p_ratio:
            # pressure ratio is interpreted as chamber/exit
            return (2 * self.f / self.theta * (1 - p_ratio ** (-self.theta / (self.theta + 1)))) ** 0.5
        else:
            # infinite pressure ratio, basically vacuum
            return (2 * self.f / self.theta) ** 0.5

    @staticmethod
    def get_name_composition_dict() -> dict:
        name_composition_dict = {}
        for composition in Composition.all_compositions:
            name_composition_dict[composition.name] = composition

        return name_composition_dict


class Propellant:
    def __init__(
        self,
        composition: Composition,
        main_geom: Geometry,
        main_r1: float,
        main_r2: float,
        aux_geom: Geometry = SimpleGeometry.SPHERE,
        aux_r1: float = 1.0,
        aux_r2: float = 1.0,
        web_ratio: float = 1.0,  # e1_2 / e1_1
        mass_ratio: float = 0.0,  # w_2 / w_1
    ):
        """
        for propellant i = 1, 2:
            Z_i = e_i / e_0_i
            psi_i = ff_i(Z_i)
            dZ_i / dt = u_1 p ^ n_i / e_0_i

        In sum:
                        n
            psi = (1/w) Σ ff_i(Z_i) w_i
                       i=1
        """
        if any((main_r1 and main_r1 < 0, main_r2 and main_r2 < 0)):
            raise ValueError("Geometry is impossible")

        self.composition = composition
        self.geometry = main_geom

        self.main_params = main_geom.get_form_function_coefficients(r1=main_r1, r2=main_r2)
        self.aux_params = aux_geom.get_form_function_coefficients(r1=aux_r1, r2=aux_r2)

        self.chi, self.labda, self.mu, self.chi_s, self.labda_s, self.z_b = self.main_params
        _, _, _, _, _, z_b2 = self.aux_params

        self.web_ratio, self.mass_ratio = web_ratio, mass_ratio

        if self.mass_ratio > 0 and z_b2 > self.z_b / self.web_ratio:
            raise ValueError("Auxiliary grains must complete combustion in advance of the primary grains.")

    @property
    def f(self) -> float:
        return self.composition.f

    @property
    def alpha(self) -> float:
        return self.composition.alpha

    @property
    def rho_p(self) -> float:
        return self.composition.rho_p

    @property
    def theta(self) -> float:
        return self.composition.theta

    @property
    def u_1(self) -> float:
        return self.composition.u_1

    @property
    def n(self) -> float:
        return self.composition.n

    @property
    def T_v(self) -> float:
        return self.composition.T_v

    @staticmethod
    def f_sigma(z_i: float, params: tuple[float, float, float, float, float, float]) -> float:
        chi, labda, mu, chi_s, labda_s, z_b = params
        if z_i <= 1.0:
            return chi * (1 + 2 * labda * z_i + 3 * mu * z_i**2)
        elif z_i <= z_b:
            return chi_s * (1 + 2 * labda_s * z_i)
        else:
            return 0.0

    @staticmethod
    def f_psi(z_i: float, params: tuple[float, float, float, float, float, float]) -> float:
        chi, labda, mu, chi_s, labda_s, z_b = params
        if z_i <= 1.0:
            return chi * z_i * (1 + labda * z_i + mu * z_i**2)
        elif z_i <= z_b:
            return chi_s * z_i * (1 + labda_s * z_i)
        else:
            return 1.0

    def f_sigma_z(self, z: float) -> float:
        z_main, z_aux = z, z / self.web_ratio
        sigma_main, sigma_aux = self.f_sigma(z_main, self.main_params), self.f_sigma(z_aux, self.aux_params)

        return sigma_main / (self.mass_ratio + 1) + sigma_aux * self.mass_ratio / (self.mass_ratio + 1)

    def f_psi_z(self, z: float) -> float:
        z_main, z_aux = z, z / self.web_ratio
        psi_main, psi_aux = self.f_psi(z_main, self.main_params), self.f_psi(z_aux, self.aux_params)

        return psi_main / (self.mass_ratio + 1) + psi_aux * self.mass_ratio / (self.mass_ratio + 1)
