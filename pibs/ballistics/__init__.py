# from ._highlow import POINT_PEAK_BLEED, POINT_PEAK_HIGH, Highlow
# from ._optHighlow import ConstrainedHighlow
# from .enums import Domains, GunTypes, Points, Solutions
from __future__ import annotations

from typing import Union, Literal

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

Points = Union[
    Literal[
        "SHOT_START", "PEAK_AVG_P", "PEAK_BREECH_P", "PEAK_SHOT_P", "FRACTURE", "BURNOUT", "SHOT_EXIT", "PEAK_STAG_P"
    ],
    str,
]

SOL_LAGRANGE = "SOL_LAGRANGE"
SOL_PIDDUCK = "SOL_PIDDUCK"
SOL_MAMONTOV = "SOL_MAMONTOV"

Solutions = Union[Literal["SOL_LAGRANGE", "SOL_PIDDUCK", "SOL_MAMONTOV"], str]

CONVENTIONAL = "CONVENTIONAL"
RECOILLESS = "RECOILLESS"

GunTypes = Union[Literal["CONVENTIONAL", "RECOILLESS"], str]

from .gun import GenericEntry, GenericResult, Gun, OutlineEntry, PressureProbePoint, PressureTraceEntry
from .material import MATERIALS
from .optGun import Constrained
from .optRecoilless import ConstrainedRecoilless
from .prop import GEOMETRIES, GrainComp, MultPerfGeometry, Propellant, SimpleGeometry
from .recoilless import Recoilless
