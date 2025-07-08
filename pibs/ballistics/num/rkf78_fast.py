import inspect
import sys
import traceback

a2 = 2 / 27
a3 = 1 / 9
a4 = 1 / 6
a5 = 5 / 12
a6 = 1 / 2
a7 = 5 / 6
a8 = 1 / 6
a9 = 2 / 3
a10 = 1 / 3
a11 = 1
a12 = 0
a13 = 1


b21 = 2 / 27

b31 = 1 / 36
b32 = 1 / 12

b41 = 1 / 24
b43 = 1 / 8

b51 = 5 / 12
b53 = -25 / 16
b54 = 25 / 16

b61 = 1 / 20
b64 = 1 / 4
b65 = 1 / 5

b71 = -25 / 108
b74 = 125 / 108
b75 = -65 / 27
b76 = 125 / 54

b81 = 31 / 300
b85 = 61 / 225
b86 = -2 / 9
b87 = 13 / 900

b91 = 2
b94 = -53 / 6
b95 = 704 / 45
b96 = -107 / 9
b97 = 67 / 90
b98 = 3

b101 = -91 / 108
b104 = 23 / 108
b105 = -976 / 135
b106 = 311 / 54
b107 = -19 / 60
b108 = 17 / 6
b109 = -1 / 12

b111 = 2383 / 4100
b114 = -341 / 164
b115 = 4496 / 1025
b116 = -301 / 82
b117 = 2133 / 4100
b118 = 45 / 82
b119 = 45 / 164
b1110 = 18 / 41

b121 = 3 / 205
b126 = -6 / 41
b127 = -3 / 205
b128 = -3 / 41
b129 = 3 / 41
b1210 = 6 / 41

b131 = -1777 / 4100
b134 = -341 / 164
b135 = 4496 / 1025
b136 = -289 / 82
b137 = 2193 / 4100
b138 = 51 / 82
b139 = 33 / 164
b1310 = 12 / 41
b1312 = 1

c1 = 41 / 840
c6 = 34 / 105
c7 = 9 / 35
c8 = 9 / 35
c9 = 9 / 280
c10 = 9 / 280
c11 = 41 / 840

c_hat = -41 / 840

EPSILON = 1e-16


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
    if record is None:
        record = []
    n = 0
    y_this = iniVal
    x = x_0

    beta = 0.84

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
    elif abs(h) < EPSILON:  # does a linear extrapolation if the epsilon is small enough
        df = dFunc(x, *y_this, h)
        y_this = [y + df * h for y, df in zip(y_this, df)]
        return x_1, y_this, Rm

    while (h > EPSILON and x < x_1) or (h < -EPSILON and x > x_1):

        if (h > 0 and (x + h) > x_1) or (h < 0 and (x + h) < x_1):
            h = x_1 - x

        try:
            # fmt: off
            allK[0] = [*map((h).__mul__, dFunc(x, *y_this, h))]

            allK[1] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a2 * h,
                        *[y + b21 * k1 for y, k1 in zip(y_this, *allK[:1])], h
                    )
                )
            ]

            allK[2] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a3 * h,
                        *[
                            y + b31 * k1 + b32 * k2
                            for y, k1, k2 in zip(y_this, *allK[:2])
                        ], h
                    )
                )
            ]

            allK[3] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a4 * h,
                        *[
                            y + b41 * k1 + b43 * k3
                            for y, k1, k2, k3 in zip(y_this, *allK[:3])
                        ], h
                    )
                )
            ]

            allK[4] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a5 * h,
                        *[
                            y + b51 * k1 + b53 * k3 + b54 * k4
                            for y, k1, k2, k3, k4 in zip(y_this, *allK[:4])
                        ], h
                    )
                )
            ]

            allK[5] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a6 * h,
                        *[
                            y + b61 * k1 + b64 * k4 + b65 * k5
                            for y, k1, k2, k3, k4, k5 in zip(y_this, *allK[:5])
                        ], h
                    )
                )
            ]

            allK[6] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a7 * h,
                        *[
                            y + b71 * k1 + b74 * k4 + b75 * k5 + b76 * k6
                            for y, k1, k2, k3, k4, k5, k6 in zip(y_this, *allK[:6])
                        ], h
                    )
                )
            ]

            allK[7] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a8 * h,
                        *[
                            y + b81 * k1 + b85 * k5 + b86 * k6 + b87 * k7
                            for y, k1, k2, k3, k4, k5, k6, k7 in zip(
                                y_this, *allK[:7]
                            )
                        ], h
                    )
                )
            ]

            allK[8] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a9 * h,
                        *[
                            y + b91 * k1 + b94 * k4 + b95 * k5 + b96 * k6
                            + b97 * k7 + b98 * k8
                            for y, k1, k2, k3, k4, k5, k6, k7, k8 in
                            zip(y_this, *allK[:8])
                        ], h
                    )
                )
            ]

            allK[9] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a10 * h,
                        *[
                            y + b101 * k1 + b104 * k4 + b105 * k5 + b106 * k6
                            + b107 * k7 + b108 * k8 + b109 * k9
                            for y, k1, k2, k3, k4, k5, k6, k7, k8, k9
                            in zip(y_this, *allK[:9])
                        ], h
                    )
                )
            ]

            allK[10] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a11 * h,
                        *[
                            y + b111 * k1 + b114 * k4 + b115 * k5 + b116 * k6
                            + b117 * k7 + b118 * k8 + b119 * k9 + b1110 * k10
                            for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10
                            in zip(y_this, *allK[:10])
                        ], h
                    )
                )
            ]

            allK[11] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a12 * h,
                        *[
                            y + b121 * k1 + b126 * k6 + b127 * k7 + b128 * k8
                            + b129 * k9 + b1210 * k10
                            for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11
                            in zip(y_this, *allK[:11])
                        ], h
                    )
                )
            ]

            allK[12] = [
                *map(
                    (h).__mul__,
                    dFunc(
                        x + a13 * h,
                        *[
                            y + b131 * k1 + b134 * k4 + b135 * k5 + b136 * k6
                            + b137 * k7 + b138 * k8 + b139 * k9 + b1310 * k10
                            + b1312 * k12
                            for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12
                            in zip(y_this, *allK[:12])
                        ], h
                    )
                )
            ]

            y_next = [
                y + c1 * k1 + c6 * k6 + c7 * k7 + c8 * k8 + c9 * k9 + c10 * k10
                + c11 * k11
                for y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, _, _ in zip(
                    y_this, *allK
                )
            ]

            te = [
                c_hat * (k1 + k11 - k12 - k13)
                for y, k1, _, _, _, _, _, _, _, _, _, k11, k12, k13 in zip(
                    y_this, *allK
                )
            ]  # local truncation error, or difference per step

            # fmt: on

        except (
            ValueError,
            TypeError,
            ZeroDivisionError,
            OverflowError,
        ):
            if debug:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                errMsg = "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
                print(f"Error encountered at x={x:.8g}")
                print(str(errMsg))
            h *= beta
            continue

        Rs = [abs(e) * (x_1 - x_0) / h for e in te]

        R = max(
            abs(r)
            / (
                minTol
                + max(
                    (relTol * min(abs(y1), abs(y2))),
                    absTol,
                )
            )
            for r, y1, y2, adapt in zip(Rs, y_this, y_next, adaptTo)
            if adapt
        )

        delta = 1

        if R >= 1:  # error is greater than acceptable
            delta = beta * abs(1 / R) ** (1 / 8)

        else:  # error is acceptable
            y_this = y_next
            x += h
            Rm = [max(Rmi, Rsi) for Rmi, Rsi in zip(Rm, Rs)]

            if abortFunc is not None and abortFunc(x=x, ys=y_this, record=record):  # premature terminating cond. is met
                if debug:
                    print("exiting via abortFunc")
                    print("record")
                    for line in record:
                        xval, yvals = line
                        print("{:^12.8g}|".format(xval), end="")
                        for yval in yvals:
                            print("{:^12.8g}|".format(yval), end="")
                        print()

                return x, y_this, Rm

            record.append([x, [*y_this]])

            if R != 0:  # sometimes the error can be estimated to be 0
                delta = beta * abs(1 / R) ** (1 / 7)

            else:
                delta = 2

        h *= min(max(delta, 0.125), 2)
        n += 1

        # if debug and n % 100 == 0:
        #     print(x, *y_this, h)

    if debug:
        print("exiting main loop normally")
        print("record")
        for line in record:
            xval, yvals = line
            print("{:^12.8g}|".format(xval), end="")
            for yval in yvals:
                print("{:^12.8g}|".format(yval), end="")
            print()
        print(Rm)

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
