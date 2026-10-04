from __future__ import annotations

from dataclasses import dataclass

from . import GunType, Point
from .params import Params


@dataclass
class TableEntry:
    tag: Point
    time: float
    travel: float
    burnup: float
    velocity: float
    breech_pressure: float
    avg_pressure: float
    shot_pressure: float
    temperature: float | None

    # recoilless-only fields (None for conventional guns)
    outflow_velocity: float | None = None
    stag_pressure: float | None = None
    outflow_fraction: float | None = None
    rel_stag_point: float | None = None


@dataclass
class SimResult:
    params: Params
    table_data: list[TableEntry]
    pressure_trace: list[PressureTraceEntry]

    tube_mass: float | None = None
    outline: list[OutlineEntry] | None = None

    @property
    def gun_type(self) -> GunType:
        return self.params.gun_type

    def read_table_data(self, tag: Point) -> TableEntry:
        for entry in self.table_data:
            if entry.tag == tag:
                return entry
        raise ValueError("no entry with tag")

    def get_eff(self) -> tuple[float, float, float]:
        """
        te: thermal efficiency
        be: ballistic efficiency
        pe: piezoelectric efficiency
        """
        params = self.params
        vg = self.read_table_data(Point.EXIT).velocity
        p_max = self.read_table_data(Point.PEAK_AVG).avg_pressure
        te = (vg / params.v_j) ** 2
        be = te / params.phi
        pe = 0.5 * params.phi * params.m * vg**2 / (p_max * params.s * params.l_g)
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
