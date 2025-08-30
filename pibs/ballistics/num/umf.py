"""
Jinpeng Zhai 翟锦鹏
2023 11 16
Contact: 914962409@qq.com

"Universal Mathematical Functions" that implements various optimization
and root-finding techniques. These have a common series of arguments as:

...(f, x_0, x_1, y, x_tol, y_abs_tol, y_rel_tol, it, debug)

basic parameters
    f:
        functions
    x_0:
        one argument
    x_1:
        another argument
    y:
        reference value, used as a target for zero finding, equivalent to
        lambda x: f(x) - y, y_abs_tol = y * y_rel_tol

tolerance specifications:
    x_tol:
        the change in x value that is considered acceptable
    y_abs_tol:
        the change in y value (or f(x)) that defines the upper limit
        of an acceptable solution.

    y_rel_tol:
        the relative change in y value (compared to a reference value,
        either the current y value, for optimization routines, or a target
        value, if given) that defines the upper limit of an acceptable
        solution.

    it:
        salvage mechanism, limit the maximum amount of iterations that will
        be accepted to prevent infinite runoff

    debug:
        boolean flag to print values, useful in debugging code.


The return signatures are:

    x_0, x_1:
        two best (first and second best) estimate of the true solution orn
        numerical optimum
"""

import math

invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2

FLOAT_MIN = 1e-16


def gss(
    f,
    a,
    b,
    x_tol=FLOAT_MIN,
    y_rel_tol=0,
    y_abs_tol=FLOAT_MIN,
    find_min=True,
    it=1e4,
    debug=False,
    f_report=None,
):
    """Golden-section search. improved from the example
    given on wikipedia. Reuse half the evaluations.

    Given a function f with a single local extremum in
    the interval [a,b], gss returns a subset interval
    [c,d] that contains the extremum with d-c <= relTol.

    a----c--d----b
    """
    record = []

    (a, b) = (min(a, b), max(a, b))

    h = b - a
    if h <= x_tol:
        return (a, b)

    ya = f(a)
    yb = f(b)

    # Required steps to achieve tolerance
    if x_tol != 0:
        n = int(math.ceil(math.log(x_tol / h, 2) / math.log(invphi, 2))) - 1
    else:
        n = math.inf
    n = min(n, it)
    n = max(n, 1)  # at least one iteration should be ran

    c = a + invphi2 * h
    d = a + invphi * h
    yc = f(c)
    yd = f(d)

    k = 0
    while k < n:
        if f_report is not None:
            f_report(k / n)

        if (yc < yd and find_min) or (yc > yd and not find_min):
            # a---c---d  b
            b = d
            d = c
            yb = yd
            yd = yc
            h *= invphi
            c = a + invphi2 * h
            yc = f(c)

            record.append((k, c, yc))

            if (
                (abs(a - d) <= x_tol)
                or (abs(ya - yd) < (y_rel_tol * min(abs(ya), abs(yd))))
                or (abs(ya - yd) <= y_abs_tol)
            ):
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

            record.append((k, d, yd))

            if (
                (abs(c - b) < x_tol)
                or (abs(yc - yb) < (y_rel_tol * min(abs(yc), abs(yb))))
                or (abs(yc - yd) < y_abs_tol)
            ):
                break

        k += 1

    if f_report is not None:
        f_report((k + 1) / n)

    if debug:
        print("GSS")
        print("{:>4}{:>24}{:>24}".format("I", "X", "FX"))
        record.sort(key=lambda line: line[1])
        for line in record:
            print("{:>4}{:>24}{:>24}".format(*line))

    if (yc < yd and find_min) or (yc > yd and not find_min):
        return a, d
    else:
        return c, b


# fmt: off
def secant(
    f, x_0, x_1, y=0, x_min=None, x_max=None, x_tol=FLOAT_MIN,
        y_rel_tol=0, y_abs_tol=FLOAT_MIN, it=100, debug=False):
    # fmt: on
    fx_0 = f(x_0) - y
    fx_1 = f(x_1) - y

    record = []

    if x_0 == x_1 or fx_0 == fx_1:
        errStr = "Impossible to calculate initial slope for secant search."
        errStr += "\nf({})-{}={}\nf({})-{}={}".format(x_0, y, fx_0, x_1, y, fx_1)
        raise ValueError(errStr)

    i = 0
    for i in range(it):
        x_2 = x_1 - fx_1 * (x_1 - x_0) / (fx_1 - fx_0)
        if x_min is not None and x_2 < x_min:
            x_2 = 0.9 * x_min + 0.1 * x_1
        if x_max is not None and x_2 > x_max:
            x_2 = 0.9 * x_max + 0.1 * x_1

        fx_2 = f(x_2) - y

        if debug:
            record.append((i, x_2, fx_2))

        if any(
            (
                abs(x_2 - x_1) < x_tol,
                abs(fx_2) < y_abs_tol,
                abs(fx_2) < (abs(y) * y_rel_tol),
            ),
        ):
            if debug:
                print("SECANT")
                print("{:>4}{:>24}{:>24}".format("I", "X", "FX"))
                record.sort(key=lambda line: line[1])
                for line in record:
                    print("{:>4}{:>24}{:>24}".format(*line))

            return x_2, x_1
        else:
            if fx_2 == fx_1:
                if debug:
                    print("SECANT")
                    print("{:>4}{:>24}{:>24}".format("I", "X", "FX"))
                    record.sort(key=lambda line: line[1])
                    for line in record:
                        print("{:>4}{:>24}{:>24}".format(*line))
                raise ValueError(
                    "Numerical plateau found at f({})-{}=f({})-{}={}".format(
                        x_1, y, x_2, y, fx_2
                    )
                )

            x_0, x_1, fx_0, fx_1 = x_1, x_2, fx_1, fx_2

    else:
        if debug:
            print("SECANT")
            print("{:>4}{:>24}{:>24}".format("I", "X", "FX"))
            record.sort(key=lambda line: line[1])
            for line in record:
                print("{:>4}{:>24}{:>24}".format(*line))

        raise ValueError(
            f"Secant method called from {x_min} to {x_max}\n"
            + f"Maximum iteration exceeded at it = {i}/{it}"
            + ",\n[0] f({})-{}={}->\n[1] f({})-{}={}->\n[2] f({})-{}={}".format(
                x_0, y, fx_0, x_1, y, fx_1, x_2, y, fx_2
            )
        )


# fmt: off
def dekker(
    f, x_0, x_1, y=0,
        x_tol=1e-16, y_rel_tol=0, y_abs_tol=1e-16, it=100, debug=False, f_report=None
):
    # fmt: on
    fx_0 = f(x_0) - y
    fx_1 = f(x_1) - y

    if fx_0 * fx_1 >= 0:
        raise ValueError(
            "Dekker method must be initiated by guesses bracketing root:\n"
            + "f({})-{}={}, f({})-{}={}".format(x_0, y, fx_0, x_1, y, fx_1)
        )

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

    y_abs_tol = max(y_abs_tol, abs(y) * y_rel_tol)

    i = 0
    for i in range(it):
        m = 0.5 * (a_j + b_j)
        if fb_i != fb_j:
            s = b_j - fb_j * (b_j - b_i) / (fb_j - fb_i)  # secant estimate
        else:
            s = m

        if (
            min(b_j, m) < s < max(b_j, m)
        ):  # if secant estimate strictly between current estimate
            # and bisection estimate
            b_k = s  # assign the secant estimation to be the next estimate
        else:
            b_k = m

        fb_k = f(b_k) - y  # calculate new value of estimate

        if (
            fa_j * fb_k < 0
        ):  # if the contrapoint is of different sign than current estimate
            a_k = a_j  # new contrapoint is still the same
            fa_k = fa_j
        else:
            a_k = b_j  # otherwise, new contrapoint should use the current est.
            fa_k = fb_j

        if abs(fa_k) < abs(fb_k):  # ensure b is still the best guess
            a_k, b_k = b_k, a_k
            fa_k, fb_k = fb_k, fa_k

        if f_report is not None:
            log_ini = math.log(max(abs(fx_0), abs(fx_1)))
            log_eps = math.log(abs(fb_k))
            log_fin = math.log(y_abs_tol)
            f_report(min((log_ini - log_eps) / (log_ini - log_fin), 1))

        if debug:
            record.append((i, b_k, fb_k))

        if any(
            (abs(b_k - a_k) < x_tol, abs(fb_k) < y_abs_tol),
        ):
            if debug:
                print("{:>4}{:>24}{:>24}".format("I", "X", "FX"))
                record.sort(key=lambda line: line[1])
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


def bisect(
    f, x_0, x_1, y=0, x_tol=1e-16, y_abs_tol=1e-16, y_rel_tol=0, debug=False
):
    """bisection method to numerically solve for zero
    two initial guesses must be of opposite sign.
    The root found is guaranteed to be within the range specified.
    """
    a, b = min(x_0, x_1), max(x_0, x_1)
    fa = f(a) - y
    fb = f(b) - y

    record = []

    if x_tol > 0:
        n = math.ceil(math.log((b - a) / x_tol, 2))
    else:
        n = math.inf

    if fa * fb >= 0:
        raise ValueError("Initial Guesses Must Be Of Opposite Sign")

    for i in range(n):
        if abs(fa - fb) < max(y_abs_tol, y * y_rel_tol):
            break

        c = 0.5 * (a + b)
        fc = f(c) - y

        record.append((i, c, fc))

        if fc * fa > 0:
            a = c
            fa = fc
        else:
            b = c
            fb = fc

    if debug:
        print("BISECT")
        print("{:>4}{:>24}{:>24}".format("I", "X", "FX"))
        record.sort(key=lambda r: r[1])
        for line in record:
            print("{:>4}{:>24}{:>24}".format(*line))

    return a, b
