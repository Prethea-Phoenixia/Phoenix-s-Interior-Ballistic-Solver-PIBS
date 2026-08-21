FLOAT_MIN = 1e-16

from .integrate import integrate
from .rkf import rkf45 as rkf
from .umf import dekker, gss
