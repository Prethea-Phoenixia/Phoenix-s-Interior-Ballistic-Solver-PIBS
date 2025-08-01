from __future__ import annotations
from dataclasses import dataclass, asdict
from functools import wraps

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .gun import Gun
    from . import Points

from .prop import Propellant

from . import POINT_EXIT, POINT_PEAK_AVG


@dataclass
class HasGetRawLine:
    def get_raw_line(self):
        return list(asdict(self).values())


@dataclass
class GenericEntry(HasGetRawLine):
    tag: str
    time: float
    travel: float
    burnup: float
    velocity: float
    breechPressure: float
    avgPressure: float
    shotPressure: float
    temperature: float


@dataclass
class GenericErrorEntry(HasGetRawLine):
    tag: str
    time: float | None = None
    travel: float | None = None
    burnup: float | None = None
    velocity: float | None = None
    breechPressure: float | None = None
    avgPressure: float | None = None
    shotPressure: float | None = None
    temperature: float | None = None


@dataclass
class GenericResult:
    gun: Gun
    table_data: list[GenericEntry]
    error_data: list[GenericEntry]
    pressure_trace: list[PressureTraceEntry]

    tubeMass: float = None
    outline: list[OutlineEntry] = None
    thermalEfficiency: float = None
    ballisticEfficiency: float = None
    piezoEfficiency: float = None

    def read_table_data(self, tag: Points) -> GenericEntry:
        for tableEntry in self.table_data:
            if tableEntry.tag == tag:
                return tableEntry
        raise ValueError("no entry with tag")

    def get_raw_table_data(self):
        raw_lines = []
        for gunTableEntry in self.table_data:
            raw_lines.append(gunTableEntry.get_raw_line())

        return raw_lines

    def get_raw_error_data(self):
        raw_lines = []
        for gunErrorEntry in self.error_data:
            raw_lines.append(gunErrorEntry.get_raw_line())

        return raw_lines

    def get_raw_pressure_trace(self):
        raw_lines = []
        for pressureProbePoint in self.pressure_trace:
            raw_lines.append(pressureProbePoint.get_raw_line())

        return raw_lines

    def get_eff(self) -> tuple[float, float, float]:
        """
        te: thermal efficiency
        be: ballistic efficiency
        pe: piezoelectric efficiency
        """
        vg = self.read_table_data(POINT_EXIT).velocity
        p_max = self.read_table_data(POINT_PEAK_AVG).avgPressure

        te = (vg / self.gun.v_j) ** 2
        be = te / self.gun.phi
        pe = 0.5 * self.gun.phi * self.gun.m * vg**2 / (p_max * self.gun.s * self.gun.l_g)
        return te, be, pe


@dataclass
class PressureTraceEntry(HasGetRawLine):
    tag: str
    temperature: float
    pressure_trace: list[PressureProbePoint]


@dataclass
class PressureProbePoint:
    x: float
    p: float


@dataclass
class OutlineEntry(HasGetRawLine):
    x: float
    r_in: float
    r_ex: float


class DelegatesPropellant:
    def __init__(self, propellant: Propellant):
        self.propellant = propellant

    @property
    def f(self) -> float:
        return self.propellant.f

    @property
    def alpha(self) -> float:
        return self.propellant.alpha

    @property
    def rho_p(self) -> float:
        return self.propellant.rho_p

    @property
    def theta(self) -> float:
        return self.propellant.theta

    @property
    def u_1(self) -> float:
        return self.propellant.u_1

    @property
    def n(self) -> float:
        return self.propellant.n

    @property
    def T_v(self) -> float:
        return self.propellant.T_v

    @property
    def z_b(self) -> float:
        return self.propellant.z_b

    @wraps(Propellant.f_psi_z)
    def f_psi_z(self, *args, **kwargs):
        return self.propellant.f_psi_z(*args, **kwargs)

    @wraps(Propellant.f_sigma_z)
    def f_sigma_z(self, *args, **kwargs):
        return self.propellant.f_sigma_z(*args, **kwargs)
