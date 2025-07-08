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
    Delta = 18 * a * b * c * d - 4 * b**3 * d + b**2 * c**2 - 4 * a * c**3 - 27 * a**2 * d**2
    """
    Δ>0: distinct real roots.
    Δ=0: repeating real roots.
    Δ<0: one real and 2 imaginary roots.
    """
    Delta_0 = b**2 - 3 * a * c
    Delta_1 = 2 * b**3 - 9 * a * b * c + 27 * a**2 * d

    C_1 = (0.5 * (Delta_1 + (Delta_1**2 - 4 * Delta_0**3) ** 0.5)) ** (1 / 3)
    C_2 = (0.5 * (Delta_1 - (Delta_1**2 - 4 * Delta_0**3) ** 0.5)) ** (1 / 3)

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

    Delta = b**2 - 4 * a * c

    x_1 = 0.5 * (-b - Delta**0.5) / a
    x_2 = 0.5 * (-b + Delta**0.5) / a

    if Delta > 0:
        return min(x_1, x_2), max(x_1, x_2)
    else:
        return x_1, x_2


def matMul(A, B):
    dimA = len(A), len(A[0])
    if any(len(row) != dimA[1] for row in A):
        raise ValueError("Matrix A is not consistent")
    dimB = len(B), len(B[0])
    if any(len(row) != dimB[1] for row in B):
        raise ValueError("Matrix B is not consistent")
    if dimA[1] != dimB[0]:
        raise ValueError("Dimension mistmatch for matrix A and B")

    R = [[0 for _ in range(dimB[1])] for _ in range(dimA[0])]
    BT = [*zip(*B)]

    i = 0
    for rowA in A:
        j = 0
        for columnB in BT:
            R[i][j] = sum(a * b for a, b in zip(rowA, columnB))
            j += 1
        i += 1
    return R


def solveMat(A, B):
    """
    Solve the linear system defined by Ax = B,
    where A is given in nested lists with the inner list representing the
    row entries, and B given in a flattened list representing the only column
    in the result vectory. A flattened list respresenting the x vector is
    returned.

    Specifically, we use Gauss-Jordanian elimination to calculate A^-1,
    and left multiply it such that A^-1*A*x = A^-1*B.

    """
    dim = len(A)

    if dim != len(B):
        raise ValueError("Dimension mismatch between A,x and B")

    if any(len(row) != dim for row in A):
        raise ValueError("Matrix A is not square")

    I = [[1 if i == j else 0 for i in range(dim)] for j in range(dim)]

    def swapRow(i, j):
        rowI = A[i], I[i]
        rowJ = A[j], I[j]

        A[i], I[i] = rowJ
        A[j], I[j] = rowI

    h = 0  # pivot row
    k = 0  # pivot column

    while h < dim and k < dim:
        # choose the largest possible absolute value as partial pivot
        imax = max(
            ((A[i][k], i) for i in range(h, dim)),
            key=lambda x: x[0],
        )[1]

        if A[imax][k] == 0:
            # no pivot in this column
            k += 1
        else:
            swapRow(h, imax)
            for i in range(h + 1, dim):
                f = A[i][k] / A[h][k]
                A[i][k] = 0  # fill the lower part of pivot column
                # do for all remaining elements in current row
                for j in range(k + 1, dim):
                    A[i][j] -= A[h][j] * f

                for j in range(0, dim):
                    # apply the same operation to the identity matrix.
                    I[i][j] -= I[h][j] * f

            h += 1
            k += 1

    for i in range(dim - 1, -1, -1):
        if A[i][i] != 0:
            for j in range(i):
                f = A[j][i] / A[i][i]
                A[j][i] = 0
                for k in range(0, dim):
                    I[j][k] -= I[i][k] * f

    # convert the leading entries to 1
    for i in range(dim):
        if A[i][i] != 0:
            f = 1 / A[i][i]
            for j in range(i, dim):
                A[i][j] *= f
            for j in range(0, dim):
                I[i][j] *= f

    # now the matrix I is converted into A^-1
    Ix = matMul(I, [[b] for b in B])
    result = [i[0] for i in Ix]

    return result


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
    I = 0  # integral counter
    c = 0  # trend counter, No. of iterations with reducing delta.
    d = math.inf  # delta, change per iteration

    while c < 3:
        dI = 0  # change to integral
        for i in range(1, 2**k, 2):
            v = -1 + 2 ** (1 - k) * i
            u = 1.5 * v - 0.5 * v**3
            dI += f(a * u + b) * (1 - v**2)

        dI *= 1.5 * a * 2 ** (1 - k)
        I1 = I * 0.5 + dI
        d = abs(I1 - I)
        I = I1
        k += 1

        if d < tol * (abs(I) + tol):
            c += 1
        else:
            c = 0

    return I, d


# def jarvis(points):  # TODO: write doc, fucking test this shit
#     def cosine(start, vertex, end):
#         if (start[0] == vertex[0]) and (start[1] == vertex[1]):
#             start = (vertex[0], vertex[1] - 1)
#
#         x_0, y_0 = start
#         x_1, y_1 = vertex
#         x_2, y_2 = end
#
#         v_x_0, v_y_0 = x_1 - x_0, y_1 - y_0
#         v_x_1, v_y_1 = x_2 - x_1, y_2 - y_1
#
#         l_0 = (v_x_0**2 + v_y_0**2) ** 0.5
#         l_1 = (v_x_1**2 + v_y_1**2) ** 0.5
#
#         v = max(min((v_x_0 * v_x_1 + v_y_0 * v_y_1) / (l_0 * l_1), 1), -1)
#
#         theta = math.acos(v)
#
#         mod = (l_0**2 + l_1**2 - 2 * l_0 * l_1 * math.cos(theta)) ** 0.5
#
#         return theta, mod
#
#     first = min(points)
#
#     start = first
#     prev = first
#     hull = [first]
#
#     for _ in range(len(points)):  # prevent infinite execution
#         theta_min, mod_max = math.inf, 0
#
#         for candidate in points:
#             if (candidate[0] == start[0]) and (candidate[1] == start[1]):
#                 continue
#
#             theta, mod = cosine(prev, start, candidate)
#             if theta < theta_min:
#                 end = candidate
#                 theta_min, mod_max = theta, mod
#
#             elif (theta == theta_min) and (mod > mod_max):
#                 end = candidate
#                 theta_min, mod_max = theta, mod
#
#         if (end[0] == first[0]) and (end[1] == first[1]):
#             return hull
#         else:
#             hull.append(end)
#             prev = start
#             start = end
#
#     raise ValueError(
#         "Something went wrong while finding the convex hull of the supplied "
#         + "points, current hull being {:}".format(hull),
#     )


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

    # # points = [[0, 0], [0, 1], [1, 1], [1, 0]]
    # points = []
    #
    # import random
    #
    # for _ in range(1000):
    #     points.append((random.random(), random.random()))
    #
    # import matplotlib.pyplot as plt
    #
    # plt.scatter(*zip(*points))
    # plt.plot(*zip(*jarvis(points)))
    #
    # plt.show()
    #
    # A = [[2, 1, -1], [-3, -1, 2], [-2, 1, 2]]
    # print(solveMat(A, [8, -11, -3]))
