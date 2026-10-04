from __future__ import annotations

from enum import Enum


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
VALID_SOLUTION_METHODS = tuple(s.value for s in SolutionMethod)
VALID_PRESSURE_POINTS = (
    Point.PEAK_BREECH.value,
    Point.PEAK_SHOT.value,
    Point.PEAK_AVG.value,
    Point.PEAK_STAG.value,
)
