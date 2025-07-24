import inspect
import sys
import traceback

"""
Constants to be used for Runge-Kutta-Fehlberg 7(8), see:

Classical Fifth-, Sixth- Seventh- and Eighth-Order Runge-Kutta
Formulas With Stepsize Control
Erwin Fehlberg, George C. Marshall Spcae Flight Center
Huntsville, Ala.
NASA, Washington D.C., October 1968
"""


# fmt: off
_as = [0, 2/27, 1/9, 1/6, 5/12, 1/2, 5/6, 1/6, 2/3, 1/3, 1, 0, 1]


# Table X, pp.65 in ref.
bs = [
    [],
    [     2/27],
    [     1/26, 1/12],
    [     1/24,    0,    1/8], 
    [     5/12,    0, -25/16,    25/16],
    [     1/20,    0,      0,      1/4,       1/5],
    [  -25/108,    0,      0,  125/108,    -65/27,   125/54],
    [   31/300,    0,      0,        0,    61/225,    -2/9,    13/900],
    [        2,    0,      0,    -53/6,    704/45,  -107/9,     67/90,     3],
    [  -91/108,    0,      0,   23/108,  -976/135,  311/54,    -19/60,  17/6,  -1/12],
    [ 2383/4100,   0,      0, -341/164, 4496/1025, -301/82, 2133/4100, 45/82, 45/164, 18/41],
    [     3/205,   0,      0,        0,         0,   -6/41,    -3/205, -3/41,   3/41, 6/ 41, 0],
    [-1777/4100,   0,      0, -341/164, 4496/1025, -289/82, 2193/4100, 51/82, 33/164, 12/41, 0, 1]
]

cs     = [41/840,  0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 41/840,      0,      0]
cs_hat = [     0,  0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280,      0, 41/840, 41/840]
# fmt: on


def RKF78(
    dFunc,
    iniVal,
    x_0,
    x_1,
    relTol,
    absTol=1e-16,
    minTol=1e-16,
    adaptTo=True,
    abortFunc=None,
    record=None,
    debug=False,
):
    """
    use Runge Kutta Fehlberg of 7(8)th power to solve system of equation
    as defined by dFunc

    Arguments:
        dFunc       : d/dx|x=x(y1, y2, y3....) = dFunc(x, y1, y2, y3..., dx)
        iniVal      : initial values for (y1, y2, y3...)
        x_0         : integration start point
        x_1         : integration end point
        relTol      : relative tolerance, per component
        absTol      : absolute tolerance, per component
        minTol      : minimum tolerance, per component. This is added to the error
                    estimation, to encourage conservatism in the integrator, and to
                    guard against division by 0 if functional value tends to 0

        abortFunc   : optional, function that accepts arguments of
                    (x - current value of integrand, ys - current value of the SoE,
                    record - record of value up to that point)
                    and terminates the integrator on a boolean value of True

        adaptTo     : optional, values used to control error.
                    if True, adapt stepsize to error estimation in every component.
                    if [Boolean] * nbr. of components, adapt stepsize to error estimation in component that is true in adaptTo.
                    a value of False will cause an exception to be raised.

        minTol      : optional, minimum magnitude of error
        record      : optional, if supplied will record all committed steps


    Returns:
        (y1, y2, y3...)|x = x_1, (e1, e2, e3....)
        where e1, e2, e3...
        are the estimated maximum deviation (in absolute) for that individual
        component
    """
    # fmt: on
    if record is None:
        record = []
    x, y_this = x_0, iniVal

    beta = 0.84  # "safety" factor
    h = x_1 - x_0  # initial step size

    Rm = [0 for _ in iniVal]

    if adaptTo is True:
        adaptTo = [True] * len(iniVal)

    sig = inspect.signature(dFunc)
    params = len([param for param in sig.parameters.values() if param.kind == param.POSITIONAL_OR_KEYWORD])
    if debug:
        paramstr = [str(param) for param in sig.parameters.values() if param.kind == param.POSITIONAL_OR_KEYWORD]
        print("setup")

        for i, param in enumerate(paramstr[:-1]):
            print("{:_^12}|".format(param), end="")

        print("\n{:^12.8g}|".format(x_0), end="")
        for i, yval in enumerate(y_this):
            print("{:^12.8g}|".format(yval), end="")
        print()

    if adaptTo is False or ((params - 2) == len(adaptTo) == len(iniVal)):
        pass
    else:
        raise ValueError(
            "Argument number mismatch between dFunc, adapTo and iniVal.\n"
            + "dFunc(x, y_0...y_i, dx)\n"
            + "adaptTo = True or (boolean_0....boolean_i)\n"
            + "iniVal = (y_0.....y_i)"
        )

    allK = [None for _ in range(13)]

    if h == 0:
        return x, y_this, Rm

    while (h > 0 and x < x_1) or (h < 0 and x > x_1):
        if (x + h) == x:
            break  # catch the error using the final lines
        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x  # this for handling the step size very close to x_1

        try:
            """
            Wrap these in a error try-catch block, since an adaptive algorithm
            can very well generate estimates that brings the functions being
            integrated out of domain, and we do not want to meticulously write
            handling for every single cases (there are many).

            These cases typically involve division by zero, negative values being
            supplied to functions defined on positive domain, complex values
            being calculated for parameters, etc etc.,

            In these cases the stepsize will be reduced to salvage the calculation --
            very important for unsupervised automated usage! -- at least until
            the point where the step size approaches the limit of floating number
            precision, then we accept that the calculation has likely reached
            an asymptote, or explosion point, and that we should abort the integration.
            Code calling rkf7(8) should implement error handling in that case.
            """
            # initialize the next estimate, which is 7th order
            y_next = [y for y in y_this]
            # initialize the error esimate, which is 8th order
            y_next_hat = [y for y in y_this]

            for i, (bi, asi) in enumerate(zip(bs, _as)):
                """i ranging from 0 to 12
                bi retrieve the line of constant from Butcher Tableau
                """
                xi = x + asi * h  # x to use for calling dfunc
                yi = [y for y in y_this]  # initialize the current y vector
                for bij, Kj in zip(bi[:i], allK):
                    """Update the current y estimate using *all* values calculated
                    up to this point.
                    zip iterates up to the shortest list.
                    bij retrieve the "weight" for previous calculation
                    Kj retrieves previous calculation
                    """

                    """This is a fancy way of writing:
                    yi   = yi  +  bij  *  Kj
                    vector vector scalar  vector
                    without using fancy Python vector libraries like numpy
                    """
                    yi = [y + k * bij for y, k in zip(yi, Kj)]

                # after the loop, yi is the new y we can call dFunc with.
                di = dFunc(xi, *yi, h)
                """
                ki   = h   *   di
                vector scalar  vector
                """
                ki = [h * d for d in di]
                allK[i] = ki

                ci = cs[i]
                ci_hat = cs_hat[i]

                # these two calculations propagate the values to each component
                y_next = [y + k * ci for y, k in zip(y_next, ki)]
                y_next_hat = [y + k * ci_hat for y, k in zip(y_next_hat, ki)]

            """
            truncation error is generated from the difference of the 7-th and 8-th
            order estimators.
            """
            te = [y - y_hat for y, y_hat in zip(y_next, y_next_hat)]

        except (
            ValueError,  # catch complex numbers being supplied to funcitons, etc
            TypeError,  # catch complex numbers being used in comparisons, etc
            ZeroDivisionError,  # divide by zero in the equation
            OverflowError,  # numerical overflow, in practice very rare
        ) as e:
            # if debug:
            #     exc_type, exc_value, exc_traceback = sys.exc_info()
            #     errMsg = "".join(
            #         traceback.format_exception(
            #             exc_type, exc_value, exc_traceback
            #         )
            #     )
            #     print(f"Error encountered at x={x:.8g}")
            #     print(str(errMsg))
            h *= beta
            continue

        """
        Extrapolating global error from local truncation error.
        Using the entire range is considered more conservative (results in larger error)
        than scaling using the remaining, and also decouples error estimate from
        current position, which results in more consistent results.
        """
        Rs = [abs(e) * (x_1 - x_0) / h for e in te]

        """
        Construct a relative error specification, comparing the global extrapolated
        error to the smaller of current and next values.
        """
        R = 0  # initialize R
        for r, y1, y2, adapt in zip(Rs, y_this, y_next, adaptTo):
            """
            the generated relative error estimation for each component is to take
            the extrapolated global error (see above), divide by the least of
            current and future estimate (in absolute terms), or the absolute tolerance
            supplied, whichever is greater. This gives the "most lenient"
            interpretation of the two. This is intended to guard against cases
            where the integrated values approaches 0, such that a relative error
            will cause the step size to collapse to zero in pursuit of an error
            that's likely caused by floating point noise.

            Furthermore, a minimum tolerance term is added to prevent division
            by zero in the case of a both a supplied absolute tolerance of zero
            and an integrand at zero point.
            """
            Ry = abs(r) / (max((relTol * min(abs(y1), abs(y2))), absTol, minTol))
            if adapt:
                R = max(R, Ry)
            else:
                pass

        delta = 1  # initialize the prospective change in step size

        if R >= 1:  # error is greater than acceptable
            delta = beta * abs(1 / R) ** (1 / 8)

        else:  # error is acceptable
            y_this = y_next
            x += h
            Rm = [max(Rmi, Rsi) for Rmi, Rsi in zip(Rm, Rs)]

            if abortFunc is not None and abortFunc(x=x, ys=y_this, record=record):  # premature terminating cond. is met
                if debug:
                    print("exiting via debug")
                    print("record")
                    for line in record:
                        x, yval = line
                        print("{:^12.8g}|".format(x), end="")
                        for i, y in enumerate(yval):
                            print("{:^12.8g}|".format(y), end="")
                        print()

                return x, y_this, Rm

            record.append([x, [*y_this]])

            if R != 0:
                # adaptively modify the step size according to composite error estimate
                delta = beta * abs(1 / R) ** (1 / 7)

            else:  # R == 0
                """
                if this continues to be true, we are integrating a polynomial,
                in which case the error should be independent of the step size.
                Therefore we aggressively increase the step size to seek forward.
                """
                delta = 2

        """
        The step size should not change too rapidly, as this risks too rapidly
        diminishing the "accuracy reserve" on allowable step size during
        reduction in step size, or alternatively drive the integrator
        over important inflection points during growth in step size.
        """
        h *= min(max(delta, 0.125), 2)

    if debug:
        print("exiting main loop normally")
        print("record")
        for line in record:
            x, yval = line
            print("{:^12.8g}|".format(x), end="")
            for i, y in enumerate(yval):
                print("{:^12.8g}|".format(y), end="")
            print()

    if abs(x - x_1) > abs(x_1 - x_0) * relTol:
        raise ValueError(
            "Premature Termination of Integration due to vanishing step size," + " x at {}, h at {}.".format(x, h)
        )

    return x, y_this, Rm


def main():
    def df1(x, y, _):
        return (7 * y**2 * x**3,)

    _, v, e = RKF78(df1, (3,), 2, 0, relTol=1e-4, absTol=1e-4, minTol=1e-14, debug=True)

    print(v)
    print(e)

    print(e[0] / v[0])
    print("expected value")
    print(-1 / (7 / 4 * 0**4 - 85 / 3))


if __name__ == "__main__":
    main()
