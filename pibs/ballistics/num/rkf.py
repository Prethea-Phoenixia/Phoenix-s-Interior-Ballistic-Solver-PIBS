from __future__ import annotations

import logging
import sys
from typing import Callable, Sequence, TypeVar

T = TypeVar("T", bound=Sequence[float])
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
    abs_tol: float = sys.float_info.epsilon,
    min_tol: float = sys.float_info.epsilon,
    abort_func: Callable[[float, T, list[tuple[float, T]]], bool] | None = None,
    record: list[tuple[float, T]] = None,
    debug: bool = False,
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

        record     : optional, if supplied will record all committed steps
        debug      : optional, enables additional debug printing when passed.

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
    x, y_this = x_0, ini_val

    beta = 0.84  # "safety" factor
    h = x_1 - x_0  # initial step size

    all_k = [list() for _ in range(len(betas))]

    if h == 0:
        return x, y_this, False

    while (h > 0 and x < x_1) or (h < 0 and x > x_1):
        if (x + h) == x:
            break  # catch the error using the final lines
        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x  # this for handling the step size very close to x_1

        try:
            # initialize the next estimate, which is p-th order
            y_next = [y for y in y_this]
            # initialize the error estimate, which is (p+1)-th order
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
            truncation error is generated from the difference of the p-th and
            (p+1)-th order estimators.
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

    if abs(x - x_1) > 8 * sys.float_info.epsilon * max(abs(x), abs(x_1)):
        raise ValueError(
            "Premature Termination of Integration due to vanishing step size," + " x at {}, h at {}.".format(x, h)
        )

    return x, y_this, False


def rkf45(
    d_func: Callable[[float, T, float], T],
    ini_val: T,
    x_0: float,
    x_1: float,
    rel_tol: float,
    abs_tol: float = sys.float_info.epsilon,
    min_tol: float = sys.float_info.epsilon,
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


if __name__ == "__main__":
    import time

    logging.basicConfig(level=logging.DEBUG)

    def df(x: float, ys: tuple[float], dx: float):
        y = ys[0]
        return (7 * y**2 * x**3,)

    print("4(5)th")
    t_0 = time.time()
    v = (0,)
    for _ in range(100):
        _, v, _ = rkf45(df, (3,), 2, 0, rel_tol=1e-4, abs_tol=1e-4, min_tol=1e-14, debug=False)
    t_1 = time.time()

    print(f"time: {t_1 - t_0}")
    print(f"computed value {v[0]}")
    true_val = -1 / (7 / 4 * 0**4 - 85 / 3)
    print(f"expected value {true_val}")
    print(f"absolute error {v[0]- true_val}")
    print(f"relative error {(v[0] - true_val)/true_val}")
