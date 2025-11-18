from __future__ import annotations
from typing import Literal, Union, TypeVar

DOMAIN_TIME = "DOMAIN_TIME"
DOMAIN_LEN = "DOMAIN_LEN"

Domains = Union[Literal["DOMAIN_TIME", "DOMAIN_LEN"], str]

POINT_START = "SHOT_START"
POINT_PEAK_AVG = "PEAK_AVG_P"
POINT_PEAK_BREECH = "PEAK_BREECH_P"
POINT_PEAK_SHOT = "PEAK_SHOT_P"
POINT_FRACTURE = "FRACTURE"
POINT_BURNOUT = "BURNOUT"
POINT_EXIT = "SHOT_EXIT"
POINT_PEAK_STAG = "PEAK_STAG_P"
SAMPLE = "SAMPLE"
COMPUTE = "COMPUTE"

Points = Union[
    Literal[
        "SHOT_START",
        "PEAK_AVG_P",
        "PEAK_BREECH_P",
        "PEAK_SHOT_P",
        "FRACTURE",
        "BURNOUT",
        "SHOT_EXIT",
        "PEAK_STAG_P",
        "COMPUTE",
        "SAMPLE",
    ],
    str,
]

SOL_LAGRANGE = "SOL_LAGRANGE"
SOL_PIDDUCK = "SOL_PIDDUCK"
SOL_MAMONTOV = "SOL_MAMONTOV"


MIN_BARR_VOLUME = "MIN_BARR_VOLUME"  # minimum bore volume
MIN_PROJ_TRAVEL = "MIN_PROJ_TRAVEL"  # minimum barrel length
Optimization_Targets = Union[Literal["MIN_BARR_VOLUME", "MIN_PROJ_TRAVEL"], str]

Solutions = Union[Literal["SOL_LAGRANGE", "SOL_PIDDUCK", "SOL_MAMONTOV"], str]

CONVENTIONAL = "CONVENTIONAL"
RECOILLESS = "RECOILLESS"

GunTypes = Union[Literal["CONVENTIONAL", "RECOILLESS"], str]

# maximum iteration to correct for chambrage effects.
MAX_ITER = 10


T = TypeVar("T")


class JSONable:
    def to_json(self) -> str: ...

    @classmethod
    def from_json(cls: T, json_dict: dict) -> T: ...


from .gun import GenericEntry, GenericResult, Gun, OutlineEntry, PressureProbePoint, PressureTraceEntry
from .material import Material
from .constrained_gun import ConstrainedGun
from .constrained_recoilless import ConstrainedRecoilless
from .prop import Composition, Geometry, MultPerfGeometry, Propellant, SimpleGeometry
from .recoilless import Recoilless
