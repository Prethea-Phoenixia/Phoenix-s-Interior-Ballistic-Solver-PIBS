from __future__ import annotations

import logging
import sys
import traceback
from math import inf
from typing import Callable, TypeVar

T = TypeVar("T")
logger = logging.getLogger(__name__)


def handle_record(record: list[tuple[float, T]]):
    output_string = "\nrecord:\n"
    for line in record:
        x, yval = line
        output_string += "{:^12.8g}|".format(x)
        for i, y in enumerate(yval):
            output_string += "{:^12.8g}|".format(y)
        output_string += "\n"

    logger.debug(output_string)


def rkf(
    order: int,
    d_func: Callable[[float, T, float], T],
    ini_val: T,
    x_0: float,
    x_1: float,
    rel_tol: float,
    abs_tol: float = 1e-16,
    min_tol: float = 1e-16,
    abort_func: Callable[[float, T, list[tuple[float, T]]], bool] = None,
    record: list[tuple[float, T]] = None,
    debug: bool = False,
    alphas: tuple[float, ...] = (),
    betas: tuple[tuple[float, ...], ...] = (),
    cs: tuple[float, ...] = (),
    c_hats: tuple[float, ...] = (),
) -> tuple[float, T, bool]:
    """
    use Runge Kutta Fehlberg of 7(8)th power to solve system of equation
    as defined by dFunc

    Arguments:
        d_func     : d/dx|x = dFunc(x, (y1, y2, y3...), dx)
        ini_val    : initial values for (y1, y2, y3...)
        x_0        : integration start point
        x_1        : integration end point
        rel_tol    : relative tolerance, per component
        abs_tol    : absolute tolerance, per component
        min_tol    : minimum tolerance, per component. This is added to the error
                    estimation, to encourage conservatism in the integrator, and to
                    guard against division by 0 if functional value tends to 0

        abort_func : optional, function that accepts arguments of
                    (x - current value of integrand, ys - current values of the SoE,
                    record - record of value up to that point)
                    and terminates the integrator on a boolean value of True

        min_tol    : optional, minimum magnitude of error
        record     : optional, if supplied will record all committed steps
        debug      : optional, enables additional debug printing when passed.

        order       : order of the runge-kutta method
        alphas      : constants of the Runge Kutta algorithm
        betas       : as above
        cs          : as above
        c_hats      : as above

    Returns:
        x_1, (y1, y2, y3...)|x = x_1, abort
    """
    if record is None:
        record = []
    x, y_this = x_0, ini_val

    beta = 0.84  # "safety" factor
    h = x_1 - x_0  # initial step size

    all_k = [list() for _ in range(13)]

    if h == 0:
        return x, y_this, False

    while (h > 0 and x < x_1) or (h < 0 and x > x_1):
        if (x + h) == x:
            break  # catch the error using the final lines
        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x  # this for handling the step size very close to x_1

        try:
            # initialize the next estimate, which is 7th order
            y_next = [y for y in y_this]
            # initialize the error estimate, which is 8th order
            y_next_hat = [y for y in y_this]

            for i, (bi, asi) in enumerate(zip(betas, alphas)):
                xi = x + asi * h  # x to use for calling dfunc
                yi = [y for y in y_this]  # initialize the current y vector
                for bij, kj in zip(bi[:i], all_k):
                    yi = [y + k * bij for y, k in zip(yi, kj)]

                # after the loop, yi is the new y we can call dFunc with.
                di = d_func(xi, yi, h)
                """
                ki   = h   *   di
                vector scalar  vector
                """
                ki = [h * d for d in di]
                all_k[i] = ki

                ci = cs[i]
                ci_hat = c_hats[i]

                # these two calculations propagate the values to each component
                y_next = [y + k * ci for y, k in zip(y_next, ki)]
                y_next_hat = [y + k * ci_hat for y, k in zip(y_next_hat, ki)]

            """
            truncation error is generated from the difference of the 7-th and 8-th
            order estimators.
            """
            truncation_errors = [y - y_hat for y, y_hat in zip(y_next, y_next_hat)]

        except (
            ValueError,  # catch complex numbers being supplied to functions, etc
            TypeError,  # catch complex numbers being used in comparisons, etc
            ZeroDivisionError,  # divide by zero in the equation
            OverflowError,  # numerical overflow, in practice very rare
        ):
            if debug:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                err_msg = "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
                logger.debug(f"Error encountered at x={x:.8g}")
                logger.debug(err_msg)

            h *= beta
            continue

        max_relative_error = 0.0  # initialize R
        for te, y1, y2 in zip(truncation_errors, y_this, y_next):
            ry = abs(te) / max((rel_tol * min(abs(y1), abs(y2))), abs_tol, min_tol)
            max_relative_error = max(max_relative_error, ry)

        if max_relative_error < 1:  # error is acceptable
            x, y_this = x + h, y_next

            if abort_func is not None and abort_func(x, y_this, record):
                # premature terminating cond. is met
                if debug:
                    handle_record(record=record)

                return x, y_this, True

            record.append((x, tuple(y_this)))

        delta = beta * abs(1 / max_relative_error) ** (1 / (order + 1)) if max_relative_error else inf
        h *= min(max(delta, 0.125), 2)

    if debug:
        logger.debug("exiting main loop normally")
        handle_record(record=record)

    if abs(x - x_1) > abs(x_1 - x_0) * rel_tol:
        raise ValueError(
            "Premature Termination of Integration due to vanishing step size," + " x at {}, h at {}.".format(x, h)
        )

    return x, y_this, False


def rkf78(
    d_func: Callable[[float, T, float], T],
    ini_val: T,
    x_0: float,
    x_1: float,
    rel_tol: float,
    abs_tol: float = 1e-16,
    min_tol: float = 1e-16,
    abort_func: Callable[[float, T, list[tuple[float, T]]], bool] = None,
    record: list[tuple[float, T]] = None,
    debug: bool = False,
) -> tuple[float, T, bool]:
    """
    use Runge Kutta Fehlberg of 7(8)th order to solve system of equation
    defined by dFunc

    Constants used for Runge-Kutta-Fehlberg 7(8), see Table X, pp.65 in ref:
    *Classical Fifth-, Sixth-, Seventh-
    and Eighth-Order Runge-Kutta Formulas With Stepsize Control,
    Erwin Fehlberg, George C. Marshall Spcae Flight Center
    Huntsville, Ala. NASA, Washington D.C., October 1968.*

    Arguments:
        d_func     : d/dx|x=x(y1, y2, y3....) = dFunc(x, y1, y2, y3..., dx)
        ini_val    : initial values for (y1, y2, y3...)
        x_0        : integration start point
        x_1        : integration end point
        rel_tol    : relative tolerance, per component
        abs_tol    : absolute tolerance, per component
        min_tol    : minimum tolerance, per component. This is added to the error
                    estimation, to encourage conservatism in the integrator, and to
                    guard against division by 0 if functional value tends to 0

        abort_func : optional, function that accepts arguments of
                    (x - current value of integrand, ys - current value of the SoE,
                    record - record of value up to that point)
                    and terminates the integrator on a boolean value of True

        min_tol    : optional, minimum magnitude of error
        record     : optional, if supplied will record all committed steps
        debug      : optional, enables additional debug printing when passed.

    Returns:
        x_1, (y1, y2, y3...)|x=x_1
    """

    betas = (
        (0,),
        (2 / 27,),
        (1 / 36, 1 / 12),
        (1 / 24, 0, 1 / 8),
        (5 / 12, 0, -25 / 16, 25 / 16),
        (1 / 20, 0, 0, 1 / 4, 1 / 5),
        (-25 / 108, 0, 0, 125 / 108, -65 / 27, 125 / 54),
        (31 / 300, 0, 0, 0, 61 / 225, -2 / 9, 13 / 900),
        (2, 0, 0, -53 / 6, 704 / 45, -107 / 9, 67 / 90, 3),
        (-91 / 108, 0, 0, 23 / 108, -976 / 135, 311 / 54, -19 / 60, 17 / 6, -1 / 12),
        (2383 / 4100, 0, 0, -341 / 164, 4496 / 1025, -301 / 82, 2133 / 4100, 45 / 82, 45 / 164, 18 / 41),
        (3 / 205, 0, 0, 0, 0, -6 / 41, -3 / 205, -3 / 41, 3 / 41, 6 / 41, 0),
        (-1777 / 4100, 0, 0, -341 / 164, 4496 / 1025, -289 / 82, 2193 / 4100, 51 / 82, 33 / 164, 12 / 41, 0, 1),
    )
    alphas = (0, 2 / 27, 1 / 9, 1 / 6, 5 / 12, 1 / 2, 5 / 6, 1 / 6, 2 / 3, 1 / 3, 1, 0, 1)
    cs = (41 / 840, 0, 0, 0, 0, 34 / 105, 9 / 35, 9 / 35, 9 / 280, 9 / 280, 41 / 840, 0, 0)
    c_hats = (0, 0, 0, 0, 0, 34 / 105, 9 / 35, 9 / 35, 9 / 280, 9 / 280, 0, 41 / 840, 41 / 840)

    return rkf(
        order=7,
        d_func=d_func,
        ini_val=ini_val,
        x_0=x_0,
        x_1=x_1,
        rel_tol=rel_tol,
        abs_tol=abs_tol,
        min_tol=min_tol,
        abort_func=abort_func,
        record=record,
        debug=debug,
        alphas=alphas,
        betas=betas,
        cs=cs,
        c_hats=c_hats,
    )


def rkf45(
    d_func: Callable[[float, T, float], T],
    ini_val: T,
    x_0: float,
    x_1: float,
    rel_tol: float,
    abs_tol: float = 1e-16,
    min_tol: float = 1e-16,
    abort_func: Callable[[float, T, list[tuple[float, T]]], bool] = None,
    record: list[tuple[float, T]] = None,
    debug: bool = False,
) -> tuple[float, T, bool]:
    """
    use Runge Kutta Fehlberg of 4(5)th order to solve system of equation
    as defined by dFunc

    Constants used for Runge-Kutta-Fehlberg 4(5), see Table II, pp.12 in ref:
    *Low Order Classical Runge-Kutta Formulas With Stepsize Control and Their Application
    to Some Heat Transfer Problems, Erwin Fehlberg, George C. Marshall Space Flight Center,
    Marshall, Alabama, NASA, Washington D.C., July, 1969*


    Arguments:
        d_func     : d/dx|x=x(y1, y2, y3....) = dFunc(x, y1, y2, y3..., dx)
        ini_val    : initial values for (y1, y2, y3...)
        x_0        : integration start point
        x_1        : integration end point
        rel_tol    : relative tolerance, per component
        abs_tol    : absolute tolerance, per component
        min_tol    : minimum tolerance, per component. This is added to the error
                    estimation, to encourage conservatism in the integrator, and to
                    guard against division by 0 if functional value tends to 0

        abort_func : optional, function that accepts arguments of
                    (x - current value of integrand, ys - current value of the SoE,
                    record - record of value up to that point)
                    and terminates the integrator on a boolean value of True

        min_tol    : optional, minimum magnitude of error
        record     : optional, if supplied will record all committed steps
        debug      : optional, enables additional debug printing when passed.

    Returns:
        x_1, (y1, y2, y3...)|x=x_1
    """
    alphas = (0, 2 / 9, 1 / 3, 3 / 4, 1, 5 / 6)
    cs = (1 / 9, 0, 9 / 20, 16 / 45, 1 / 12, 0)
    c_hats = (47 / 450, 0, 12 / 25, 32 / 225, 1 / 30, 6 / 25)
    betas = (
        (0,),
        (2 / 9,),
        (1 / 12, 1 / 4),
        (69 / 128, -243 / 128, 135 / 64),
        (-17 / 12, 27 / 4, -27 / 5, 16 / 15),
        (65 / 432, -5 / 16, 13 / 16, 4 / 27, 5 / 144),
    )
    return rkf(
        order=4,
        d_func=d_func,
        ini_val=ini_val,
        x_0=x_0,
        x_1=x_1,
        rel_tol=rel_tol,
        abs_tol=abs_tol,
        min_tol=min_tol,
        abort_func=abort_func,
        record=record,
        debug=debug,
        alphas=alphas,
        betas=betas,
        cs=cs,
        c_hats=c_hats,
    )


def rkf34(
    d_func: Callable[[float, T, float], T],
    ini_val: T,
    x_0: float,
    x_1: float,
    rel_tol: float,
    abs_tol: float = 1e-16,
    min_tol: float = 1e-16,
    abort_func: Callable[[float, T, list[tuple[float, T]]], bool] = None,
    record: list[tuple[float, T]] = None,
    debug: bool = False,
) -> tuple[float, T, bool]:
    """
    use Runge Kutta Fehlberg of 3(4)th order to solve system of equation
    as defined by dFunc

    Constants used for Runge-Kutta-Fehlberg 3(4), see Table VII, pp.22 in ref:
    *Low Order Classical Runge-Kutta Formulas With Stepsize Control and Their Application
    to Some Heat Transfer Problems, Erwin Fehlberg, George C. Marshall Space Flight Center,
    Marshall, Alabama, NASA, Washington D.C., July, 1969*


    Arguments:
        d_func     : d/dx|x=x(y1, y2, y3....) = dFunc(x, y1, y2, y3..., dx)
        ini_val    : initial values for (y1, y2, y3...)
        x_0        : integration start point
        x_1        : integration end point
        rel_tol    : relative tolerance, per component
        abs_tol    : absolute tolerance, per component
        min_tol    : minimum tolerance, per component. This is added to the error
                    estimation, to encourage conservatism in the integrator, and to
                    guard against division by 0 if functional value tends to 0

        abort_func : optional, function that accepts arguments of
                    (x - current value of integrand, ys - current value of the SoE,
                    record - record of value up to that point)
                    and terminates the integrator on a boolean value of True

        min_tol    : optional, minimum magnitude of error
        record     : optional, if supplied will record all committed steps
        debug      : optional, enables additional debug printing when passed.

    Returns:
        x_1, (y1, y2, y3...)|x=x_1
    """
    alphas = (0, 1 / 4, 4 / 9, 6 / 7, 1)
    cs = (1 / 6, 0, 27 / 52, 49 / 156, 0)
    c_hats = (43 / 288, 0, 243 / 416, 343 / 1872, 1 / 12)
    betas = ((0,), (1 / 4,), (4 / 81, 32 / 81), (57 / 98, -432 / 343, 1053 / 686), (1 / 6, 0, 27 / 52, 49 / 156))
    return rkf(
        order=3,
        d_func=d_func,
        ini_val=ini_val,
        x_0=x_0,
        x_1=x_1,
        rel_tol=rel_tol,
        abs_tol=abs_tol,
        min_tol=min_tol,
        abort_func=abort_func,
        record=record,
        debug=debug,
        alphas=alphas,
        betas=betas,
        cs=cs,
        c_hats=c_hats,
    )


def rkf23(
    d_func: Callable[[float, T, float], T],
    ini_val: T,
    x_0: float,
    x_1: float,
    rel_tol: float,
    abs_tol: float = 1e-16,
    min_tol: float = 1e-16,
    abort_func: Callable[[float, T, list[tuple[float, T]]], bool] = None,
    record: list[tuple[float, T]] = None,
    debug: bool = False,
) -> tuple[float, T, bool]:
    """
    use Runge Kutta Fehlberg of 2nd(3rd) order to solve system of equation
    defined by dFunc

    Constants used for Runge-Kutta-Fehlberg 2(3), see Table XII, pp.28 in ref:
    *Low Order Classical Runge-Kutta Formulas With Stepsize Control and Their Application
    to Some Heat Transfer Problems, Erwin Fehlberg, George C. Marshall Space Flight Center,
    Marshall, Alabama, NASA, Washington D.C., July, 1969*


    Arguments:
        d_func     : d/dx|x=x(y1, y2, y3....) = dFunc(x, y1, y2, y3..., dx)
        ini_val    : initial values for (y1, y2, y3...)
        x_0        : integration start point
        x_1        : integration end point
        rel_tol    : relative tolerance, per component
        abs_tol    : absolute tolerance, per component
        min_tol    : minimum tolerance, per component. This is added to the error
                    estimation, to encourage conservatism in the integrator, and to
                    guard against division by 0 if functional value tends to 0

        abort_func : optional, function that accepts arguments of
                    (x - current value of integrand, ys - current value of the SoE,
                    record - record of value up to that point)
                    and terminates the integrator on a boolean value of True

        min_tol    : optional, minimum magnitude of error
        record     : optional, if supplied will record all committed steps
        debug      : optional, enables additional debug printing when passed.

    Returns:
        x_1, (y1, y2, y3...)|x=x_1
    """
    alphas = (0, 1, 1 / 2)
    cs = (1 / 2, 1 / 2, 0)
    c_hats = (1 / 6, 1 / 6, 2 / 3)
    betas = ((0,), (1,), (1 / 4, 1 / 4))
    return rkf(
        order=2,
        d_func=d_func,
        ini_val=ini_val,
        x_0=x_0,
        x_1=x_1,
        rel_tol=rel_tol,
        abs_tol=abs_tol,
        min_tol=min_tol,
        abort_func=abort_func,
        record=record,
        debug=debug,
        alphas=alphas,
        betas=betas,
        cs=cs,
        c_hats=c_hats,
    )


def main():
    import time

    logging.basicConfig(encoding="utf-8", level=logging.DEBUG)

    def df(x: float, ys: tuple[float], dx: float):
        y = ys[0]
        return (7 * y**2 * x**3,)

    for order, function in zip(("7(8)th", "4(5)th", "3(4)th", "2(3)rd with 3 evals"), (rkf78, rkf45, rkf34, rkf23)):
        print(order)
        t_0 = time.time()
        v = (0,)
        for _ in range(100):
            _, v, _ = function(df, (3,), 2, 0, rel_tol=1e-4, abs_tol=1e-4, min_tol=1e-14, debug=False)
        t_1 = time.time()

        print(f"time: {t_1 - t_0}")
        print(f"computed value {v[0]}")
        true_val = -1 / (7 / 4 * 0**4 - 85 / 3)
        print(f"expected value {true_val}")
        print(f"absolute error {v[0]- true_val}")
        print(f"relative error {(v[0] - true_val)/true_val}")


if __name__ == "__main__":
    main()
