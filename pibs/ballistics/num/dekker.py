import math
import sys
from typing import Callable


def dekker(
    f: Callable[[float], float],
    x_0: float,
    x_1: float,
    y: float = 0.0,
    x_tol: float = 0.0,
    y_abs_tol: float = 0.0,
    y_rel_tol: float = 0.0,
    debug: bool = False,
) -> tuple[float, float]:

    x_tol = max(abs(x_tol), max(abs(x_0), abs(x_1), 1.0) * sys.float_info.epsilon)
    y_rel_tol = abs(y_rel_tol)
    y_abs_tol = max(abs(y_abs_tol), abs(y) * y_rel_tol, max(abs(y), 1.0) * sys.float_info.epsilon)

    fx_0 = f(x_0) - y
    fx_1 = f(x_1) - y

    if fx_0 * fx_1 > 0:
        raise ValueError(
            "Dekker method must be initiated by guesses bracketing root:\n"
            + "f({})-{}={}, f({})-{}={}".format(x_0, y, fx_0, x_1, y, fx_1)
        )

    elif fx_0 == 0:
        return x_0, x_0
    elif fx_1 == 0:
        return x_1, x_1

    if abs(fx_0) < abs(fx_1):
        b_j = x_0  # assign the better of the two initial guesses to b_j
        fb_j = fx_0

        b_i = a_j = x_1  # and the worse, the last guess of root b_i
        fb_i = fa_j = fx_1
    else:
        b_j = x_1
        fb_j = fx_1

        b_i = a_j = x_0
        fb_i = fa_j = fx_0

    record = []

    i = 0

    it = math.ceil(math.log2(abs(x_1 - x_0) / x_tol)) + 1
    for i in range(it):
        m = 0.5 * (a_j + b_j)
        if fb_i != fb_j:
            s = b_j - fb_j * (b_j - b_i) / (fb_j - fb_i)  # secant estimate
        else:
            s = m

        if min(b_j, m) < s < max(b_j, m):  # if secant estimate strictly between current estimate
            # and bisection estimate
            b_k = s  # assign the secant estimation to be the next estimate
        else:
            b_k = m

        fb_k = f(b_k) - y  # calculate new value of estimate

        if fa_j * fb_k < 0:  # if the contra-point is of different sign than current estimate
            a_k = a_j  # new contra-point is still the same
            fa_k = fa_j
        else:
            a_k = b_j  # otherwise, new contra-point should use the current est.
            fa_k = fb_j

        if abs(fa_k) < abs(fb_k):  # ensure b is still the best guess
            a_k, b_k = b_k, a_k
            fa_k, fb_k = fb_k, fa_k

        if debug:
            record.append((i, b_k, fb_k))

        if any(
            (abs(b_k - a_k) < x_tol, abs(fb_k) < y_abs_tol),
        ):
            if debug:
                print("{:>4}{:>24}{:>24}".format("I", "X", "FX"))
                record.sort(key=lambda v: v[1])
                for line in record:
                    print("{:>4}{:>24}{:>24}".format(*line))

            return b_k, a_k  # return the best, and the bracketing solution

        a_j = a_k
        fa_j = fa_k

        b_i, b_j = b_j, b_k
        fb_i, fb_j = fb_j, fb_k

    else:
        if debug:
            print("{:>4}{:>24}{:>24}".format("I", "X", "FX"))
            record.sort(key=lambda r: r[1])
            for line in record:
                print("{:>4}{:>24}{:>24}".format(*line))

        raise ValueError(
            "Dekker method called from {} to {}\n".format(x_0, x_1)
            + "Maximum iteration exceeded at it = {}/{}".format(i, it)
            + ",\nf({})-{}={}->\nf({})-{}={}".format(b_i, y, fb_i, b_j, y, fb_j)
        )
