from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, TypeVar

T = TypeVar("T")
if TYPE_CHECKING:
    from .gun import Gun

from . import Point


@dataclass
class GenericEntry:
    tag: Point
    time: float
    travel: float
    burnup: float
    velocity: float
    breech_pressure: float
    avg_pressure: float
    shot_pressure: float
    temperature: float | None


@dataclass
class GenericResult:
    gun: Gun
    table_data: list[GenericEntry]
    pressure_trace: list[PressureTraceEntry]

    tube_mass: float | None = None
    outline: list[OutlineEntry] | None = None

    def read_table_data(self, tag: Point) -> GenericEntry:
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
        vg = self.read_table_data(Point.EXIT).velocity
        p_max = self.read_table_data(Point.PEAK_AVG).avg_pressure
        te = (vg / self.gun.v_j) ** 2
        be = te / self.gun.phi
        pe = 0.5 * self.gun.phi * self.gun.m * vg**2 / (p_max * self.gun.s * self.gun.l_g)
        return te, be, pe


@dataclass
class PressureTraceEntry:
    tag: Point
    temperature: float | None
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
