import math


def cubic(a, b, c, d):
    """
    returns the 3 roots of
    ax^3 + bx^2 + cx + d = 0
    assuming **real** coefficients.
    """
    if any(isinstance(i, complex) for i in (a, b, c, d)):
        raise ValueError("coefficients must be real")
    if a == 0:
        return quadratic(b, c, d)
    delta = 18 * a * b * c * d - 4 * b**3 * d + b**2 * c**2 - 4 * a * c**3 - 27 * a**2 * d**2
    """
    Δ>0: distinct real roots.
    Δ=0: repeating real roots.
    Δ<0: one real and 2 imaginary roots.
    """
    delta_0, delta_1 = b**2 - 3 * a * c, 2 * b**3 - 9 * a * b * c + 27 * a**2 * d

    c_1 = (0.5 * (delta_1 + (delta_1**2 - 4 * delta_0**3) ** 0.5)) ** (1 / 3)
    c_2 = (0.5 * (delta_1 - (delta_1**2 - 4 * delta_0**3) ** 0.5)) ** (1 / 3)

    xs = []
    if any(c != 0 for c in (c_1, c_2)):
        c = c_1 if c_1 != 0 else c_2
        epsilons = (
            1,
            complex(-0.5, 3**0.5 / 2),
            complex(-0.5, -(3**0.5) / 2),
        )
        for epsilon in epsilons:
            x = -1 / (3 * a) * (b + c * epsilon + delta_0 / (c * epsilon))
            xs.append(x)
    else:
        for _ in range(3):
            xs.append(-b / (3 * a))

    if delta >= 0:
        xs = list(z.real for z in xs)
    else:
        # one real and 2 imaginary roots.
        xs = list(z.real if abs(z.imag) == min(abs(z.imag) for z in xs) else z for z in xs)
    # put the first real solution at first.
    xs.sort(key=lambda z: 1 if isinstance(z, complex) else 0)
    return tuple(xs)


def quadratic(a, b, c):
    """
    solve the quadratic equation
    defined by:
    y = a*x**2 + b * x + c
    """

    delta = b**2 - 4 * a * c

    x_1 = 0.5 * (-b - delta**0.5) / a
    x_2 = 0.5 * (-b + delta**0.5) / a

    if delta > 0:
        return min(x_1, x_2), max(x_1, x_2)
    else:
        return x_1, x_2


def intg(f, l, u, tol=1e-3):
    """
    Integration, a.la the HP-34C. For more info see:
    "Handheld Calculator Evaluates Integrals", William M.Kahan
    Hewlett Packard Journal, August 1980 Volume 31, number 8.

    f: function, single variable.
    l: lower limit
    u: upper limit of integration
    tol: tolerance, see below

    To apply the quadrature procedure, first the problem is transformed on
    interval to:

    u              1                        given:
    ∫ f(x) dx -> a ∫ f(ax+b) dx             a = (u - l) / 2
    l             -1                        b = (u + l) / 2

    another transformation on the variable of integration eliminates the need
    to sample at either end points, which makes it possible to evaluate improper
    integrals if asymptotes are at either end point.

    1                                        1
    ∫ f(u) du -- let u = 1.5v-0.5v**3 -> 1.5 ∫ f(1.5v-0.5v^3)*(1-v^2) dv
    -1                                      -1

    as the weight (1-v^2) is exactly 0 on both end points. We then sample
    evenly along v, take quadrature using the mid-point rule and doubling
    the number of nodes taken for each pass. This helps with suppressing
    harmonics if the function integrated is periodic. In addition, all of
    the previously calcualted quadratures can be reused in the next round,
    after dividing by half. This is especially important when function calls
     are expensive. Specifically, for pass k (k>=1 & integer) we consider 2^k-1
    points (besides the end points):

    v(i) = -1 + 2^(1-k) * i as i increments from 1 to 2^k-1 (inclusive).

                                   2^k-1
    let I(k) =  2^(1-k) * 1.5 * a * Σ f(1.5v-0.5v^3)*(1-v^2)
                                   i=1

                                     2^k+1
    then I(k+1) = 2^(-k) * 1.5 * a * Σ f(1.5v-0.5v^3)*(1-v^2) for every odd i + I(k)/2
                                     i=1

    as a rough approximation, the error is simply taken to be the change in
    estimated value between two successive evaluations:

    ΔI(k) = I(k) - I(k-1)

    if the quadrature procedure results in a converging result, then the error
    should decrease faster than the increment in the result, speaking in
    absolute terms. Although this is no-way guaranteed, it is convenient to
    take the increment as an upper bound on error. Therefore we check for three
    consecutive increments smaller than the specified tolerance before
    submitting the result as a good enough estimate for the integral.
    """

    a = (u - l) / 2
    b = (u + l) / 2

    tol = abs(tol)  # ensure positive

    k = 1  # iteration counter
    i = 0  # integral counter
    c = 0  # trend counter, No. of iterations with reducing delta.
    d = math.inf  # delta, change per iteration

    while c < 3:
        di = 0  # change to integral
        for i in range(1, 2**k, 2):
            v = -1 + 2 ** (1 - k) * i
            u = 1.5 * v - 0.5 * v**3
            di += f(a * u + b) * (1 - v**2)

        di *= 1.5 * a * 2 ** (1 - k)
        i1 = i * 0.5 + di
        d = abs(i1 - i)
        i = i1
        k += 1

        if d < tol * (abs(i) + tol):
            c += 1
        else:
            c = 0

    return i, d


if __name__ == "__main__":
    print(cubic(1, 1, 2, 3))

    import sys
    import trace

    # create a Trace object, telling it what to ignore, and whether to
    # do tracing or line-counting or both.
    tracer = trace.Trace(ignoredirs=[sys.prefix, sys.exec_prefix], trace=0, count=1)

    # run the new command using the given tracer
    tracer.run("main()")

    # make a report, placing output in the current directory
    r = tracer.results()
    r.write_results(show_missing=True, coverdir=".")
