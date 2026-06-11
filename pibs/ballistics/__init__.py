from __future__ import annotations

from enum import Enum
from typing import TypeVar


class StrEnum(str, Enum):
    def __str__(self) -> str:
        return self.value


class Domain(StrEnum):
    TIME = "DOMAIN_TIME"
    LEN = "DOMAIN_LEN"


class Point(StrEnum):
    START = "SHOT_START"
    PEAK_AVG = "PEAK_AVG_P"
    PEAK_BREECH = "PEAK_BREECH_P"
    PEAK_SHOT = "PEAK_SHOT_P"
    FRACTURE = "FRACTURE"
    BURNOUT = "BURNOUT"
    EXIT = "SHOT_EXIT"
    PEAK_STAG = "PEAK_STAG_P"
    SAMPLE = "SAMPLE"
    COMPUTE = "COMPUTE"


class SolutionMethod(StrEnum):
    LAGRANGE = "SOL_LAGRANGE"
    PIDDUCK = "SOL_PIDDUCK"
    MAMONTOV = "SOL_MAMONTOV"


class OptimizationTarget(StrEnum):
    MIN_BARR_VOLUME = "MIN_BARR_VOLUME"
    MIN_PROJ_TRAVEL = "MIN_PROJ_TRAVEL"


class GunType(StrEnum):
    CONVENTIONAL = "CONVENTIONAL"
    RECOILLESS = "RECOILLESS"


# Valid value collections
VALID_DOMAINS = tuple(d.value for d in Domain)
VALID_SOLUTION_METHODS = tuple(s.value for s in SolutionMethod)
VALID_PRESSURE_POINTS = (
    Point.PEAK_BREECH.value,
    Point.PEAK_SHOT.value,
    Point.PEAK_AVG.value,
    Point.PEAK_STAG.value,
)
GUN_PEAK_POINTS = (Point.PEAK_AVG.value, Point.PEAK_SHOT.value, Point.PEAK_BREECH.value)
RECOILLESS_PEAK_POINTS = (
    Point.PEAK_AVG.value,
    Point.PEAK_SHOT.value,
    Point.PEAK_BREECH.value,
    Point.PEAK_STAG.value,
)
VALID_OPT_TARGETS = tuple(o.value for o in OptimizationTarget)
VALID_GUN_TYPES = tuple(g.value for g in GunType)

MAX_ITER = 10


T = TypeVar("T")


class JSONable:
    def to_json(self) -> str:
        raise NotImplementedError

    @classmethod
    def from_json(cls: T, json_dict: dict) -> T:
        raise NotImplementedError
