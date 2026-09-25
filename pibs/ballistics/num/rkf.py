from __future__ import annotations

import logging
import sys
import traceback
import unittest
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
    abort_func: Callable[[float, T, list[tuple[float, T]]], bool] | None = None,
    record: list[tuple[float, T]] | None = None,
    alphas: tuple[float, ...] = (),
    betas: tuple[tuple[float, ...], ...] = (),
    cs: tuple[float, ...] = (),
    c_hats: tuple[float, ...] = (),
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

        abort_func : optional, function that accepts arguments of
                    (x - current value of integrand, ys - current values of the SoE,
                    record - record of value up to that point)
                    and terminates the integrator on a boolean value of True

        record     : optional, if supplied will record all committed steps
        order       : order p of the embedded p(p+1) pair
        alphas      : stage nodes of the Runge Kutta algorithm
        betas       : stage coefficients of the Runge Kutta algorithm
        cs          : weights of the p-th order (committed) solution
        c_hats      : weights of the (p+1)-th order (error estimate) solution

    Returns:
        x_1, (y1, y2, y3...)|x = x_1, abort
    """

    if record is None:
        record = []

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
            h_new = h * beta
            if h_new == h:
                # step size has underflowed to a denormal fixed point
                # (0.84 * h rounds back to h); give up and let the final
                # lines report premature termination
                break
            h = h_new
            continue

        max_relative_error = sys.float_info.epsilon  # initialize R
        for j in range(n):
            te = truncation_errors[j]
            y1: float = y_this[j]
            y2: float = y_next[j]
            # construct the error estimation:
            ry = abs(te) / (rel_tol * max(abs(y1), abs(y2), 1))
            # if ry > max_relative_error:
            #     max_relative_error = ry
            #

            max_relative_error = max(ry, max_relative_error)

        if max_relative_error < 1:  # error is acceptable
            x = x + h
            for j in range(n):
                y_this[j] = y_next[j]

            if abort_func is not None and abort_func(x, y_this, record):
                return x, y_this, True

            record.append((x, tuple(y_this)))

        delta = beta * abs(1 / max_relative_error) ** (1 / (order + 1))
        h *= min(max(delta, 0.125), 2)

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
    abort_func: Callable[[float, T, list[tuple[float, T]]], bool] | None = None,
    record: list[tuple[float, T]] | None = None,
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

        abort_func : optional, function that accepts arguments of
                    (x - current value of integrand, ys - current value of the SoE,
                    record - record of value up to that point)
                    and terminates the integrator on a boolean value of True

        record     : optional, if supplied will record all committed steps

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
        abort_func=abort_func,
        record=record,
        alphas=alphas,
        betas=betas,
        cs=cs,
        c_hats=c_hats,
    )


class TestRKF(unittest.TestCase):

    def test_func(self):

        def df(x: float, ys: tuple[float]) -> tuple[float]:
            return (7 * ys[0] ** 2 * x**3,)

        x_1, y_1, abort = rkf45(df, (3.0,), 2, 0, rel_tol=1e-4)

        true_val = -1 / (7 / 4 * 0**4 - 85 / 3)

        self.assertAlmostEqual(y_1[0], true_val, places=4)


if __name__ == "__main__":
    unittest.main()
