import math
import sys
from typing import Callable

invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2


def gss(
    f: Callable[[float], float],
    a: float,
    b: float,
    x_tol: float = sys.float_info.epsilon,
    y_tol: float = sys.float_info.epsilon,
    find_min: bool = True,
) -> tuple[float, float]:
    """Golden-section search. improved from the example
    given on wikipedia. Reuse half the evaluations.

    Given a function f with a single local extremum in
    the interval [a,b], gss returns a subset interval
    [c,d] that contains the extremum with d-c <= relTol.

    a----c--d----b
    """

    a, b = (min(a, b), max(a, b))

    h = b - a
    if h <= x_tol:
        return a, b

    ya = f(a)
    yb = f(b)

    x_tol = max(x_tol, sys.float_info.epsilon)
    n = int(math.ceil(math.log(x_tol / h, 2) / math.log(invphi, 2))) - 1
    n = max(n, 1)  # at least one iteration should be ran

    c = a + invphi2 * h
    d = a + invphi * h
    yc = f(c)
    yd = f(d)

    k = 0
    while k < n:
        if (yc < yd and find_min) or (yc > yd and not find_min):
            # a---c---d  b
            b = d
            d = c
            yb = yd
            yd = yc
            h *= invphi
            c = a + invphi2 * h
            yc = f(c)

            if (abs(a - d) < x_tol) or (abs(ya - yd) < y_tol):
                break

        else:
            # a   c--d---b
            a = c
            c = d
            ya = yc
            yc = yd
            h *= invphi
            d = a + invphi * h
            yd = f(d)

            if (abs(c - b) < x_tol) or (abs(yc - yb) < y_tol):
                break

        k += 1

    if (yc < yd and find_min) or (yc > yd and not find_min):
        return a, d
    else:
        return c, b
