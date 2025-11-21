import math
from typing import Callable


def integrate(f: Callable[[float], float], l: float, u: float, tol: float = 1e-3) -> tuple[float, float]:
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
    harmonics if the function integrated is periodic. In addition, all
    the previously calculated quadratures can be reused in the next round,
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
    take the increment as an upper bound on error. Therefore, we check for three
    consecutive increments smaller than the specified tolerance before
    submitting the result as a good enough estimate for the integral.
    """

    a = (u - l) / 2
    b = (u + l) / 2

    tol = abs(tol)  # ensure positive

    it = 1  # iteration counter
    integral = 0  # integral counter
    count = 0  # trend counter, No. of iterations with reducing delta.
    delta = math.inf  # delta, change per iteration

    while count < 3:
        increment = 0  # change to integral
        for i in range(1, 2**it, 2):
            v = -1 + 2 ** (1 - it) * i
            u = 1.5 * v - 0.5 * v**3
            increment += f(a * u + b) * (1 - v**2)

        increment *= 1.5 * a * 2 ** (1 - it)
        next_integral = integral * 0.5 + increment

        delta = abs(next_integral - integral)
        integral = next_integral
        it += 1

        if delta < tol * (abs(integral) + tol):
            count += 1
        else:
            count = 0

    return integral, delta


if __name__ == "__main__":

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
