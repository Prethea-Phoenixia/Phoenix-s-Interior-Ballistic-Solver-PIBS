import math
import sys
from typing import Callable

invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2


def gss(
    f: Callable[[float], float],
    a: float,
    b: float,
    x_tol: float = 0.0,
    find_min: bool = True,
) -> tuple[float, float]:
    """Golden-section search. improved from the example
    given on wikipedia. Reuse half the evaluations.

    Given a function f with a single local extremum in
    the interval [a,b], gss returns a subset interval
    [a,b] that contains the extremum with b-a <= x_tol
    (an absolute tolerance).

    a----c--d----b
    """

    a, b = (min(a, b), max(a, b))
    h = b - a
    if h <= x_tol:
        return a, b

    x_tol = max(abs(x_tol), sys.float_info.epsilon * max(abs(a), abs(b)))

    n = int(math.ceil(math.log(x_tol / h) / math.log(invphi))) - 1

    c, d = a + invphi2 * h, a + invphi * h
    yc = f(c)
    yd = f(d)

    for _ in range(n):
        h *= invphi
        if (yc < yd and find_min) or (yc > yd and not find_min):
            """
            from: a   c d b
            to:   a c d b
            """
            b, d = d, c
            yd = yc
            c = a + invphi2 * h
            yc = f(c)
        else:
            """
            from: a c d   b
            to:     a c d b
            """
            a, c = c, d
            yc = yd
            d = a + invphi * h
            yd = f(d)

    if (yc < yd and find_min) or (yc > yd and not find_min):
        return a, d
    else:
        return c, b
