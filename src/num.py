import math
import cmath


def sign(x):
    if x > 0:
        return 1
    elif x == 0:
        return 0
    elif x < 0:
        return -1


def intg(f, l, u, tol=1e-3):
    """
    To apply the quadrature procedure, first the problem is transformed on interval to:

    u              1                        given:
    ∫ f(x) dx -> a ∫ f(ax+b) dx             a = (u - l) / 2
    l             -1                        b = (u + l) / 2

    another transformation on the variable of integration eliminates the need to sample
    at either end points, which makes it possible to evaluate improper integrals if the
    asymptotes is at either end point.

    1                                        1
    ∫ f(u) du -- let u = 1.5v-0.5v**3 -> 1.5 ∫ f(1.5v-0.5v^3)*(1-v^2) dv
    -1                                      -1

    as the weight (1-v^2) is exactly 0 on both end points. We then sample evenly along
    v, take quadrature using the mid-point rule and doubling the number of nodes taken
    for each pass. This helps with suppressing harmonics if the function integrated is
    periodic. In addition, all of the previously calcualted quadratures can be reused
    in the next round, after dividing by half. This is especially important when funct-
    ion calls are expensive. Specifically, for pass k (k>=1 & integer) we consider 2^k-1
    points (besides the end points):

    v(i) = -1 + 2^(1-k) * i as i increments from 1 to 2^k-1 (inclusive).

                                   2^k-1
    let I(k) =  2^(1-k) * 1.5 * a * Σ f(1.5v-0.5v^3)*(1-v^2)
                                   i=1

                                     2^k+1
    then I(k+1) = 2^(-k) * 1.5 * a * Σ f(1.5v-0.5v^3)*(1-v^2) for every odd i + I(k)/2
                                     i=1

    as a rough approximation, the error is simply taken to be the change in estimated value
    between two successive evaluations:

    ΔI(k) = I(k) - I(k-1)

    if the quadrature procedure results in a converging result, then the error should dec-
    rease faster than the increment in the result, speaking in absolute terms. Although
    this is no-way guaranteed, it is convenient to take the increment as an upper bound
    on error. Therefore we check for three consecutive increments smaller than the specified
    tolerance before submitting the result as a good enough estimate for the integral.
    """

    tol = abs(tol)

    a = (u - l) / 2
    b = (u + l) / 2

    k = 1
    I = 0
    c = 0

    while c < 3:
        dI = 0
        for i in range(1, 2**k, 2):
            v = -1 + 2 ** (1 - k) * i
            u = 1.5 * v - 0.5 * v**3
            dI += f(a * u + b) * (1 - v**2)

        dI *= 1.5 * a * 2 ** (1 - k)
        d = abs(I * 0.5 - dI)
        I = I * 0.5 + dI
        k += 1

        if d < tol:
            c += 1
        else:
            c = 0

    return I, d


invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2


def gss(f, a, b, tol=1e-9, findMin=True):
    """Golden-section search. improved from the example
    given on wikipedia. Reuse half the evaluatins.

    Given a function f with a single local extremum in
    the interval [a,b], gss returns a subset interval
    [c,d] that contains the extremum with d-c <= tol.

    a----c--d----b
    """

    (a, b) = (min(a, b), max(a, b))
    h = b - a
    if h <= tol:
        return (a, b)

    # Required steps to achieve tolerance
    n = int(math.ceil(math.log(tol / h) / math.log(invphi)))

    c = a + invphi2 * h
    d = a + invphi * h
    yc = f(c)
    yd = f(d)

    for k in range(n - 1):
        if (yc < yd and findMin) or (yc > yd and not findMin):
            # a---c---d
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            yc = f(c)
        else:
            # c--d---b
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            yd = f(d)

    if (yc < yd and findMin) or (yc > yd and not findMin):
        return (a, d)
    else:
        return (c, b)


def RKF78(dFunc, iniVal, x_0, x_1, tol, absTol=1e-14, termAbv=None):
    """ """
    i = 0
    if termAbv is None:
        termAbv = tuple(None for _ in iniVal)
    y_this = iniVal
    x = x_0
    beta = 0.9  # "safety" factor
    h = x_1 - x_0  # initial step size
    Rm = tuple(0 for _ in iniVal)
    while (h > 0 and x < x_1) or (h < 0 and x > x_1):
        if (x + h) == x:
            break  # catch the error using the final lines
        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x
        try:
            allK = []
            K1 = tuple(k * h for k in dFunc(x, *y_this))
            allK.append(K1)

            K2 = tuple(
                k * h
                for k in dFunc(
                    x + 2 / 27 * h,
                    *(y + 2 / 27 * k1 for y, k1 in zip(y_this, *allK))
                )
            )
            allK.append(K2)

            K3 = tuple(
                k * h
                for k in dFunc(
                    x + 1 / 9 * h,
                    *(
                        y + 1 / 36 * k1 + 1 / 12 * k2
                        for y, k1, k2 in zip(y_this, *allK)
                    )
                )
            )
            allK.append(K3)

            K4 = tuple(
                k * h
                for k in dFunc(
                    x + 1 / 6 * h,
                    *(
                        y + 1 / 24 * k1 + 1 / 8 * k3
                        for y, k1, k2, k3 in zip(y_this, *allK)
                    )
                )
            )
            allK.append(K4)

            K5 = tuple(
                k * h
                for k in dFunc(
                    x + 5 / 12 * h,
                    *(
                        y + 5 / 12 * k1 - 25 / 16 * k3 + 25 / 16 * k4
                        for y, k1, k2, k3, k4 in zip(y_this, *allK)
                    )
                )
            )
            allK.append(K5)

            K6 = tuple(
                k * h
                for k in dFunc(
                    x + 1 / 2 * h,
                    *(
                        y + 1 / 20 * k1 + 1 / 4 * k4 + 1 / 5 * k5
                        for y, k1, k2, k3, k4, k5 in zip(y_this, *allK)
                    )
                )
            )
            allK.append(K6)

            K7 = tuple(
                k * h
                for k in dFunc(
                    x + 5 / 6 * h,
                    *(
                        y
                        - 25 / 108 * k1
                        + 125 / 108 * k4
                        - 65 / 27 * k5
                        + 125 / 54 * k6
                        for y, k1, k2, k3, k4, k5, k6 in zip(y_this, *allK)
                    )
                )
            )
            allK.append(K7)

            K8 = tuple(
                k * h
                for k in dFunc(
                    x + 1 / 6 * h,
                    *(
                        y
                        + 31 / 300 * k1
                        + 61 / 225 * k5
                        - 2 / 9 * k6
                        + 13 / 900 * k7
                        for y, k1, k2, k3, k4, k5, k6, k7 in zip(y_this, *allK)
                    )
                )
            )
            allK.append(K8)

            K9 = tuple(
                k * h
                for k in dFunc(
                    x + 2 / 3 * h,
                    *(
                        y
                        + 2 * k1
                        - 53 / 6 * k4
                        + 704 / 45 * k5
                        - 107 / 9 * k6
                        + 67 / 90 * k7
                        + 3 * k8
                        for y, k1, k2, k3, k4, k5, k6, k7, k8 in zip(
                            y_this, *allK
                        )
                    )
                )
            )
            allK.append(K9)

            K10 = tuple(
                k * h
                for k in dFunc(
                    x + 1 / 3 * h,
                    *(
                        y
                        - 91 / 108 * k1
                        + 23 / 108 * k4
                        - 976 / 135 * k5
                        + 311 / 54 * k6
                        - 19 / 60 * k7
                        + 17 / 6 * k8
                        - 1 / 12 * k9
                        for y, k1, k2, k3, k4, k5, k6, k7, k8, k9 in zip(
                            y_this, *allK
                        )
                    )
                )
            )
            allK.append(K10)

            K11 = tuple(
                k * h
                for k in dFunc(
                    x + h,
                    *(
                        y
                        + 2383 / 4100 * k1
                        - 341 / 164 * k4
                        + 4496 / 1025 * k5
                        - 301 / 82 * k6
                        + 2133 / 4100 * k7
                        + 45 / 82 * k8
                        + 45 / 164 * k9
                        + 18 / 41 * k10
                        for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 in zip(
                            y_this, *allK
                        )
                    )
                )
            )
            allK.append(K11)

            K12 = tuple(
                k * h
                for k in dFunc(
                    x,
                    *(
                        y
                        + 3 / 205 * k1
                        - 6 / 41 * k6
                        - 3 / 205 * k7
                        - 3 / 41 * k8
                        + 3 / 41 * k9
                        + 6 / 41 * k10
                        for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11 in zip(
                            y_this, *allK
                        )
                    )
                )
            )
            allK.append(K12)

            K13 = tuple(
                k * h
                for k in dFunc(
                    x + h,
                    *(
                        y
                        - 1777 / 4100 * k1
                        - 341 / 164 * k4
                        + 4496 / 1025 * k5
                        - 289 / 82 * k6
                        + 2193 / 4100 * k7
                        + 51 / 82 * k8
                        + 33 / 164 * k9
                        + 12 / 41 * k10
                        + k12
                        for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12 in zip(
                            y_this, *allK
                        )
                    )
                )
            )
            allK.append(K13)

            if any(isinstance(i, complex) for k in allK for i in k):
                raise TypeError

        except (
            TypeError,
            ZeroDivisionError,
        ):  # complex value has been encountered during calculation
            # or that through unfortuante chance we got a divide by zero

            h *= beta
            continue

        y_next = tuple(
            y
            + 41 / 840 * k1
            + 34 / 105 * k6
            + 9 / 35 * k7
            + 9 / 35 * k8
            + 9 / 280 * k9
            + 9 / 280 * k10
            + 41 / 840 * k11
            for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13 in zip(
                y_this, *allK
            )
        )
        """
        z_next = tuple(
            y
            + 34 / 105 * k6
            + 9 / 35 * k7
            + 9 / 35 * k8
            + 9 / 280 * k9
            + 9 / 280 * k10
            + 41 / 840 * k12
            + 41 / 840 * k13
            for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13 in zip(
                y_this, *allK
            )
        )
        """
        te = tuple(
            41 / 840 * (k1 + k11 - k12 - k13)
            for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13 in zip(
                y_this, *allK
            )
        )
        Rs = tuple(abs(e) / h for e in te)
        Rs = tuple(
            r / (absTol + tol * (abs(y) + abs(k1)))
            for r, y, k1 in zip(Rs, y_this, K1)
        )
        R = max(Rs)

        delta = 1

        if R >= tol:  # error is greater than acceptable
            delta = beta * abs(tol / R) ** (1 / 8)

        else:  # error is acceptable
            y_this = y_next
            x += h
            i += 1
            Rm = tuple(max(Rmi, Rsi) for Rmi, Rsi in zip(Rm, Rs))
            if any(
                cv > pv if pv is not None else False
                for cv, pv in zip(y_this, termAbv)
            ):  # premature terminating cond. is met
                return y_this
            if R != 0:  # sometimes the error can be estimated to be 0
                delta = beta * abs(tol / R) ** (1 / 7)

            else:
                """
                if this continues to be true, we are integrating a polynomial,
                in which case the error should be independent of the step size
                Therefore we aggressively increase the step size to seek forward.
                """
                delta = 2

        h *= min(max(delta, 0.3), 2)  # to ensure that this does not jump

    if abs(x - x_1) > tol:
        raise ValueError(
            "Premature Termination of Integration due to vanishing step size,"
            + " x at {}, h at {}.".format(x, h)
        )

    return y_this, Rm


def cubic(a, b, c, d):
    """
    returns the 3 roots of
    ax^3 + bx^2 + cx + d = 0
    assuming **real** coefficients.
    """
    if any(isinstance(i, complex) for i in (a, b, c, d)):
        raise ValueError("coefficients must be real")

    Delta = (
        18 * a * b * c * d
        - 4 * b**3 * d
        + b**2 * c**2
        - 4 * a * c**3
        - 27 * a**2 * d**2
    )
    """
    Δ>0: distinct real roots.
    Δ=0: repeating real roots.
    Δ<0: one real and 2 imaginary roots.
    """
    Delta_0 = b**2 - 3 * a * c
    Delta_1 = 2 * b**3 - 9 * a * b * c + 27 * a**2 * d

    C_1 = (0.5 * (Delta_1 + (Delta_1**2 - 4 * Delta_0**3) ** 0.5)) ** (
        1 / 3
    )
    C_2 = (0.5 * (Delta_1 - (Delta_1**2 - 4 * Delta_0**3) ** 0.5)) ** (
        1 / 3
    )

    xs = []
    if any(C != 0 for C in (C_1, C_2)):
        C = C_1 if C_1 != 0 else C_2
        epsilons = (
            1,
            complex(-0.5, 3**0.5 / 2),
            complex(-0.5, -(3**0.5) / 2),
        )
        for epsilon in epsilons:
            x = -1 / (3 * a) * (b + C * epsilon + Delta_0 / (C * epsilon))
            xs.append(x)
    else:
        for _ in range(3):
            xs.append(-b / (3 * a))

    if Delta >= 0:
        xs = list(z.real for z in xs)
    else:
        # one real and 2 imaginary roots.
        xs = list(
            z.real if abs(z.imag) == min(abs(z.imag) for z in xs) else z
            for z in xs
        )
    # put the first real solution at first.
    xs.sort(key=lambda z: 1 if isinstance(z, complex) else 0)
    return tuple(xs)


def bisect(f, x_0, x_1, tol, it=100):
    """bisection method to numerically solve for zero
    two initial guesses must be of opposite sign.
    The root found is guaranteed to be within the range specified.

    tol: tolerance of f(x) if x is root
    it: maximum iteration
    """
    a = x_0
    b = x_1

    if abs(f(a)) < tol:
        return a, f(a)
    elif abs(f(b)) < tol:
        return b, f(b)
    elif sign(f(a)) * sign(f(b)) > 0:
        raise ValueError("Initial guesses must be of opposite sign")

    for _ in range(it):
        c = (a + b) / 2

        if abs(f(c)) < tol:
            return (c, f(c))

        if sign(f(c)) == sign(f(a)):
            a = c
        else:
            b = c

    raise ValueError("Maximum iteration exceeded at ({},{})".format(c, f(c)))


def secant(f, x_0, x_1, x_min=None, x_max=None, tol=1e-6, it=100):
    """secant method that solves f(x) = 0 subjected to x in [x_min,x_max]"""
    if x_min is not None:
        if x_0 < x_min:
            x_0 = x_min
        if x_1 < x_min:
            x_1 = x_min
    if x_max is not None:
        if x_0 > x_max:
            x_0 = x_max
        if x_1 > x_max:
            x_1 = x_max

    fx_0 = f(x_0)
    fx_1 = f(x_1)

    if x_0 == x_1 or fx_0 == fx_1:
        raise ValueError("Initial guess must evaluate to different values")

    for _ in range(it):
        x_2 = x_1 - fx_1 * (x_1 - x_0) / (fx_1 - fx_0)
        if x_min is not None and x_2 < x_min:
            x_2 = x_min
        if x_max is not None and x_2 > x_max:
            x_2 = x_max
        x_0, x_1, fx_0, fx_1 = x_1, x_2, fx_1, f(x_2)
        if abs(fx_1) < tol or (fx_0 == fx_1):
            return x_1, fx_1

    raise ValueError("Maximum iteration exceeded at ({},{})".format(x_1, fx_1))


if __name__ == "__main__":
    from tabulate import tabulate

    print("verifying RKF78 is of 7th order")
    tols = [1e-3, 1e-4, 1e-5]
    devs = []
    for i in range(0, 12):
        dev_i = [i]
        for tol in tols:
            val = RKF78(
                lambda x, y: (x**i,),
                (0,),
                0.5,
                1.5,
                tol=tol,
            )

            act = (1.5 ** (i + 1) - (0.5) ** (i + 1)) / (i + 1)
            dev_i.append((act - val[0]) / act)

        devs.append(dev_i)

    print(tabulate(devs, headers=("order", *(str(t) for t in tols), "isSame?")))
