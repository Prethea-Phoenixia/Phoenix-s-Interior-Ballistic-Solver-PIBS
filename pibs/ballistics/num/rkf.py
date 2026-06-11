from __future__ import annotations

import logging
import sys
import traceback
from typing import Callable, Sequence, TypeVar

T = TypeVar("T", bound=Sequence[float])


def handle_record(record: list[tuple[float, T]], logger: logging.Logger):
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
    d_func: Callable[[float, T], T],
    ini_val: T,
    x_0: float,
    x_1: float,
    rel_tol: float,
    abs_tol: float = sys.float_info.epsilon,
    abort_func: Callable[[float, T, list[tuple[float, T]]], bool] | None = None,
    record: list[tuple[float, T]] | None = None,
    debug: bool = False,
    alphas: tuple[float, ...] = (),
    betas: tuple[tuple[float, ...], ...] = (),
    cs: tuple[float, ...] = (),
    c_hats: tuple[float, ...] = (),
    logger: logging.Logger | None = None,
) -> tuple[float, T, bool]:
    """
    drive an embedded Runge-Kutta-Fehlberg p(p+1) pair to solve a system of
    equations as defined by dFunc. The p-th order solution is committed at each
    accepted step, while the (p+1)-th order solution is used only for the local
    truncation error estimate that drives the adaptive step size control.

    This implementation uses pre-allocated arrays and in-place mutations for
    improved performance.

    Arguments:
        d_func     : d/dx|x = dFunc(x, (y1, y2, y3...))
        ini_val    : initial values for (y1, y2, y3...)
        x_0        : integration start point
        x_1        : integration end point
        rel_tol    : relative tolerance, per component
        abs_tol    : absolute tolerance, per component

        abort_func : optional, function that accepts arguments of
                    (x - current value of integrand, ys - current values of the SoE,
                    record - record of value up to that point)
                    and terminates the integrator on a boolean value of True

        record     : optional, if supplied will record all committed steps
        debug      : optional, enables additional debug printing when passed.

        order       : order p of the embedded p(p+1) pair
        alphas      : stage nodes of the Runge Kutta algorithm
        betas       : stage coefficients of the Runge Kutta algorithm
        cs          : weights of the p-th order (committed) solution
        c_hats      : weights of the (p+1)-th order (error estimate) solution
        logger      : optional, logger for debug output

    Returns:
        x_1, (y1, y2, y3...)|x = x_1, abort
    """

    if logger is None:
        logger = logging.getLogger(__name__)

    if record is None:
        record = []

    abs_tol = max(abs(abs_tol), sys.float_info.epsilon)
    x, y_this = x_0, list(ini_val)

    beta = 0.84  # "safety" factor
    h = x_1 - x_0  # initial step size

    n = len(y_this)
    n_stages = len(betas)

    # Pre-allocate arrays once
    all_k = [[0.0] * n for _ in range(n_stages)]
    y_next = [0.0] * n
    y_next_hat = [0.0] * n
    yi = [0.0] * n
    ki = [0.0] * n
    truncation_errors = [0.0] * n

    if h == 0:
        return x, y_this, False

    while (h > 0 and x < x_1) or (h < 0 and x > x_1):
        if (x + h) == x:
            break  # catch the error using the final lines
        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x  # this for handling the step size very close to x_1

        try:
            # Initialize the next estimates
            for j in range(n):
                y_next[j] = y_this[j]
                y_next_hat[j] = y_this[j]

            for i in range(n_stages):
                bi = betas[i]
                asi = alphas[i]
                xi = x + asi * h  # x to use for calling dfunc

                # Initialize yi from y_this
                for j in range(n):
                    yi[j] = y_this[j]

                # Add contributions from previous stages
                for j in range(n):
                    yi_j = yi[j]
                    for k_idx in range(i):
                        yi_j += all_k[k_idx][j] * bi[k_idx]
                    yi[j] = yi_j

                # after the loop, yi is the new y we can call dFunc with.
                di = d_func(xi, yi)
                """
                ki   = h   *   di
                vector scalar  vector
                """

                # Compute ki = h * di
                for j in range(n):
                    ki[j] = h * di[j]
                    all_k[i][j] = ki[j]

                ci = cs[i]
                ci_hat = c_hats[i]

                # these two calculations propagate the values to each component
                for j in range(n):
                    y_next[j] += ki[j] * ci
                    y_next_hat[j] += ki[j] * ci_hat

            """
            truncation error is generated from the difference of the p-th and
            (p+1)-th order estimators.
            """
            for j in range(n):
                truncation_errors[j] = y_next[j] - y_next_hat[j]

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

        max_relative_error = sys.float_info.epsilon  # initialize R
        for j in range(n):
            te = truncation_errors[j]
            y1 = y_this[j]
            y2 = y_next[j]
            ry = abs(te) / max((rel_tol * min(abs(y1), abs(y2))), abs_tol)
            if ry > max_relative_error:
                max_relative_error = ry

        if max_relative_error < 1:  # error is acceptable
            x = x + h
            for j in range(n):
                y_this[j] = y_next[j]

            if abort_func is not None and abort_func(x, y_this, record):
                if debug:
                    handle_record(record=record, logger=logger)

                return x, y_this, True

            record.append((x, tuple(y_this)))

        delta = beta * abs(1 / max_relative_error) ** (1 / (order + 1))
        h *= min(max(delta, 0.125), 2)

    if debug:
        logger.debug("exiting main loop normally")
        handle_record(record=record, logger=logger)

    if abs(x - x_1) > sys.float_info.epsilon * max(abs(x), abs(x_1)):
        raise ValueError(
            "Premature Termination of Integration due to vanishing step size," + " x at {}, h at {}.".format(x, h)
        )

    return x, y_this, False


def rkf45(
    d_func: Callable[[float, T], T],
    ini_val: T,
    x_0: float,
    x_1: float,
    rel_tol: float,
    abs_tol: float = sys.float_info.epsilon,
    abort_func: Callable[[float, T, list[tuple[float, T]]], bool] | None = None,
    record: list[tuple[float, T]] | None = None,
    debug: bool = False,
    logger: logging.Logger | None = None,
) -> tuple[float, T, bool]:
    """
    use Runge Kutta Fehlberg of 4(5)th order to solve system of equation
    as defined by dFunc

    Constants used for Runge-Kutta-Fehlberg 4(5), see Table II, pp.12 in ref:
    *Low Order Classical Runge-Kutta Formulas With Stepsize Control and Their Application
    to Some Heat Transfer Problems, Erwin Fehlberg, George C. Marshall Space Flight Center,
    Marshall, Alabama, NASA, Washington D.C., July, 1969*


    Arguments:
        d_func     : d/dx|x=x(y1, y2, y3....) = dFunc(x, (y1, y2, y3...))
        ini_val    : initial values for (y1, y2, y3...)
        x_0        : integration start point
        x_1        : integration end point
        rel_tol    : relative tolerance, per component
        abs_tol    : absolute tolerance, per component

        abort_func : optional, function that accepts arguments of
                    (x - current value of integrand, ys - current value of the SoE,
                    record - record of value up to that point)
                    and terminates the integrator on a boolean value of True

        record     : optional, if supplied will record all committed steps
        debug      : optional, enables additional debug printing when passed.
        logger     : optional, logger for debug output

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
        abort_func=abort_func,
        record=record,
        debug=debug,
        alphas=alphas,
        betas=betas,
        cs=cs,
        c_hats=c_hats,
        logger=logger,
    )


if __name__ == "__main__":
    import time

    logging.basicConfig(level=logging.DEBUG)

    def df(x: float, ys: tuple[float]):
        y = ys[0]
        return (7 * y**2 * x**3,)

    print("4(5)th")
    t_0 = time.time()
    v = (0,)
    for _ in range(100):
        _, v, _ = rkf45(df, (3.0,), 2, 0, rel_tol=1e-4, abs_tol=1e-4, debug=False)
    t_1 = time.time()

    print(f"time: {t_1 - t_0}")
    print(f"computed value {v[0]}")
    true_val = -1 / (7 / 4 * 0**4 - 85 / 3)
    print(f"expected value {true_val}")
    print(f"absolute error {v[0]- true_val}")
    print(f"relative error {(v[0] - true_val)/true_val}")
