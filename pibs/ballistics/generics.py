from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, TypeVar

T = TypeVar("T")
if TYPE_CHECKING:
    from .gun import Gun
    from . import Points

from . import POINT_EXIT, POINT_PEAK_AVG


@dataclass
class GenericEntry:
    tag: str
    time: float
    travel: float
    burnup: float
    velocity: float
    breech_pressure: float
    avg_pressure: float
    shot_pressure: float
    temperature: float


@dataclass
class GenericResult:
    gun: Gun
    table_data: list[GenericEntry]
    pressure_trace: list[PressureTraceEntry]

    tubeMass: float = None
    outline: list[OutlineEntry] = None
    thermal_efficiency: float = None
    ballistic_efficiency: float = None
    piezo_efficiency: float = None

    def read_table_data(self, tag: Points) -> GenericEntry:
        for tableEntry in self.table_data:
            if tableEntry.tag == tag:
                return tableEntry
        raise ValueError("no entry with tag")

    def get_eff(self) -> tuple[float, float, float]:
        """
        te: thermal efficiency
        be: ballistic efficiency
        pe: piezoelectric efficiency
        """
        vg = self.read_table_data(POINT_EXIT).velocity
        p_max = self.read_table_data(POINT_PEAK_AVG).avg_pressure
        te = (vg / self.gun.v_j) ** 2
        be = te / self.gun.phi
        pe = 0.5 * self.gun.phi * self.gun.m * vg**2 / (p_max * self.gun.s * self.gun.l_g)
        return te, be, pe


@dataclass
class PressureTraceEntry:
    tag: str
    temperature: float
    pressure_trace: list[PressureProbePoint]


@dataclass
class PressureProbePoint:
    x: float
    p: float


@dataclass
class OutlineEntry:
    x: float
    r_in: float
    r_ex: float
    r_pej: float  # plastic elastic junction
