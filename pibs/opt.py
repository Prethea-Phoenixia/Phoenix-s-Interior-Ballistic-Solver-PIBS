from num import gss, RKF78, cubic, dekker
from prop import Propellant
from random import uniform
from math import pi, log
from gun import pidduck
from gun import SOL_LAGRANGE, SOL_PIDDUCK, SOL_MAMONTOV
from gun import POINT_PEAK_AVG, POINT_PEAK_BREECH, POINT_PEAK_SHOT

"""
Machine-accuracy factor, determines that, if a numerical method
is used within another, then how much more accurate should the
inner method be called as compared to the outter. A small value
is necessary to prevent numerical instability form causing sporadic
appearance of outlier, at the cost of increased computation times.
"""
N = 100


class Constrained:
    def __init__(
        self,
        caliber,
        shotMass,
        propellant,
        startPressure,
        dragCoefficient,
        designPressure,
        designVelocity,
        chambrage,
        **_,
    ):
        # constants for constrained designs

        if any(
            (
                caliber <= 0,
                shotMass <= 0,
                startPressure <= 0,
                dragCoefficient < 0,
                dragCoefficient >= 1,
            )
        ):
            raise ValueError("Invalid parameters for constrained design")

        if any(
            (
                designPressure <= 0,
                # designPressure <= startPressure,
                designVelocity <= 0,
            )
        ):
            raise ValueError("Invalid design constraint")

        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.propellant = propellant
        self.p_0 = startPressure
        self.phi_1 = 1 / (1 - dragCoefficient)

        # design limits
        self.p_d = designPressure
        self.v_d = designVelocity
        self.chi_k = chambrage

    def __getattr__(self, attrName):
        if "propellant" in vars(self) and not (
            attrName.startswith("__") and attrName.endswith("__")
        ):
            try:
                return getattr(self.propellant, attrName)
            except AttributeError:
                raise AttributeError(
                    "%r object has no attribute %r"
                    % (self.__class__.__name__, attrName)
                )
        else:
            raise AttributeError

    def solve(
        self,
        loadFraction,
        chargeMassRatio,
        tol,
        minWeb=1e-6,
        containBurnout=True,
        maxLength=1e3,
        labda_1=None,
        labda_2=None,
        cc=None,
        sol=SOL_LAGRANGE,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        control=POINT_PEAK_AVG,
        **_,
    ):
        # start = time.time()
        if any(
            (
                minWeb <= 0,
                tol <= 0,
                chargeMassRatio <= 0,
                loadFraction <= 0,
                loadFraction > 1,
                maxLength <= 0,
                ambientRho < 0,
                ambientP < 0,
                ambientGamma < 1,
            )
        ):
            raise ValueError(
                "Invalid parameters to solve constrained design problem"
            )
        """
        minWeb  : represents minimum possible grain size
        """
        m = self.m
        rho_p = self.rho_p
        theta = self.theta
        f = self.f
        chi = self.chi
        labda = self.labda
        mu = self.mu
        S = self.S
        maxLF = self.maxLF
        phi_1 = self.phi_1
        p_0 = self.p_0
        v_d = self.v_d
        p_d = self.p_d
        u_1 = self.u_1
        n = self.n
        alpha = self.alpha
        Z_b = self.Z_b
        chi_k = self.chi_k
        f_psi_Z = self.f_psi_Z

        if cc is None:
            cc = 1 / chi_k

        if loadFraction > maxLF:
            raise ValueError(
                "Specified Load Fraction Violates Geometrical Constraint"
            )

        omega = m * chargeMassRatio
        V_0 = omega / (rho_p * loadFraction)
        Delta = omega / V_0
        l_0 = V_0 / S
        gamma = theta + 1

        if any((labda_1 is None, labda_2 is None)):
            if sol == SOL_LAGRANGE:
                labda_1, labda_2 = 1 / 2, 1 / 3
            elif sol == SOL_PIDDUCK:
                labda_1, labda_2 = pidduck(omega / (phi_1 * m), gamma, tol)
            elif sol == SOL_MAMONTOV:
                labda_1, labda_2 = pidduck(omega / (phi_1 * m), 1, tol)
            else:
                raise ValueError("Unknown Solution")

        phi = phi_1 + labda_2 * omega / m * cc
        v_j = (2 * f * omega / (theta * phi * m)) ** 0.5
        v_bar_d = v_d / v_j

        if ambientRho != 0:
            c_a_bar = (ambientGamma * ambientP / ambientRho) ** 0.5 / v_j
            p_a_bar = ambientP / (f * Delta)
        else:
            c_a_bar = 0
            p_a_bar = 0

        gamma_1 = ambientGamma

        if v_j < v_d:
            raise ValueError(
                "Propellant load too low to achieve design velocity. "
                + " The 2nd ballistic limit for this loading conditions is"
                + " {:.4g} m/s.".format(v_j)
            )

        psi_0 = (1 / Delta - 1 / rho_p) / (f / p_0 + alpha - 1 / rho_p)
        Zs = cubic(a=chi * mu, b=chi * labda, c=chi, d=-psi_0)
        # pick a valid solution between 0 and 1
        Zs = tuple(
            Z
            for Z in Zs
            if not isinstance(Z, complex) and (Z > 0.0 and Z < 1.0)
        )  # evaluated from left to right, guards against complex >/< float
        if len(Zs) < 1:
            raise ValueError(
                "Propellant either could not develop enough pressure to "
                + "overcome start pressure, or has burnt to post fracture."
            )
        Z_0 = Zs[0]

        def _f_p_bar(Z, l_bar, v_bar):
            psi = f_psi_Z(Z)
            l_psi_bar = 1 - Delta / rho_p - Delta * (alpha - 1 / rho_p) * psi
            p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

            if control == POINT_PEAK_AVG:
                return p_bar

            else:
                prime = (1 / chi_k + l_bar) / (1 + l_bar)
                labda_1_prime = labda_1 * prime
                labda_2_prime = labda_2 * prime

                factor_s = 1 + labda_2_prime * (omega / (phi_1 * m))
                factor_b = (phi_1 * m + labda_2_prime * omega) / (
                    phi_1 * m + labda_1_prime * omega
                )

                if control == POINT_PEAK_SHOT:
                    return p_bar / factor_s
                elif control == POINT_PEAK_BREECH:
                    return p_bar / factor_b

        """
        step 1, find grain size that satisifies design pressure
        """
        p_bar_s = _f_p_bar(Z_0, 0, 0)
        p_bar_d = p_d / (f * Delta)  # convert to unitless

        if p_bar_d < p_bar_s:
            raise ValueError(
                "Interior ballistics of conventional gun precludes"
                + " a pressure lower than that at shot start."
            )
        l_bar_d = maxLength / l_0

        def abort_Z(x, ys, record):
            Z = x
            t_bar, l_bar, v_bar = ys
            p_bar = _f_p_bar(Z, l_bar, v_bar)

            o_x, o_ys = record[-1]

            oZ = o_x
            ot_bar, ol_bar, ov_bar = o_ys
            op_bar = _f_p_bar(oZ, ol_bar, ov_bar)

            return (p_bar < op_bar) or (p_bar > 2 * p_bar_d)

        def _f_p_e_1(e_1, tol=tol):
            """
            calculate either the peak pressure, given the arc thickness,
            or until the system develops 2x design pressure.
            """

            B = (
                S**2
                * e_1**2
                / (f * phi * omega * m * u_1**2)
                * (f * Delta) ** (2 * (1 - n))
            )

            def _ode_Z(Z, t_bar, l_bar, v_bar, _):
                """burnup domain ode of internal ballistics"""
                psi = f_psi_Z(Z)
                l_psi_bar = (
                    1 - Delta / rho_p - Delta * (alpha - 1 / rho_p) * psi
                )

                p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)
                if c_a_bar != 0 and v_bar > 0:
                    v_r = v_bar / c_a_bar
                    p_d_bar = (
                        +0.25 * gamma_1 * (gamma_1 + 1) * v_r**2
                        + gamma_1
                        * v_r
                        * (1 + (0.25 * (gamma_1 + 1)) ** 2 * v_r**2) ** 0.5
                    ) * p_a_bar

                else:
                    p_d_bar = 0

                if Z <= Z_b:
                    dt_bar = (2 * B / theta) ** 0.5 * p_bar**-n  # dt_bar/dZ
                    dl_bar = v_bar * dt_bar
                    dv_bar = 0.5 * theta * (p_bar - p_d_bar) * dt_bar

                else:
                    dt_bar = 0
                    dl_bar = 0
                    dv_bar = 0

                return [dt_bar, dl_bar, dv_bar]

            record = [[Z_0, [0, 0, 0]]]

            Z_j, (t_bar_j, l_bar_j, v_bar_j), e = RKF78(
                dFunc=_ode_Z,
                iniVal=(0, 0, 0),
                x_0=Z_0,
                x_1=Z_b,
                relTol=tol,
                absTol=tol**2,
                abortFunc=abort_Z,
                record=record,
            )

            p_bar_j = _f_p_bar(Z_j, l_bar_j, v_bar_j)

            if (
                p_bar_j >= 2 * p_bar_d
            ):  # case for abort due to excessive pressure
                return p_bar_j - p_bar_d, Z_j, t_bar_j, l_bar_j, v_bar_j

            # case for abort due to decreasing pressure

            def _f_p_Z(Z):
                i = record.index([v for v in record if v[0] <= Z][-1])
                x = record[i][0]
                ys = record[i][1]

                r = []
                _, (t_bar, l_bar, v_bar), _ = RKF78(
                    dFunc=_ode_Z,
                    iniVal=ys,
                    x_0=x,
                    x_1=Z,
                    relTol=tol,
                    absTol=tol**2,
                    record=r,
                )
                xs = [v[0] for v in record]
                record.extend(v for v in r if v[0] not in xs)
                record.sort()
                return _f_p_bar(Z, l_bar, v_bar)

            """
            find the peak pressure point.
            """

            if len(record) > 1:
                Z_i = record[-2][0]
            else:
                Z_i = Z_0

            Z_1, Z_2 = gss(_f_p_Z, Z_i, Z_j, y_rel_tol=tol, findMin=False)
            Z_p = 0.5 * (Z_1 + Z_2)

            if abs(Z_p - Z_b) < tol:
                Z_p = Z_b

            return _f_p_Z(Z_p) - p_bar_d, record[-1][0], *record[-1][-1]

        dp_bar_probe = _f_p_e_1(minWeb)[0]
        probeWeb = minWeb

        if dp_bar_probe < 0:
            raise ValueError(
                "Design pressure cannot be achieved by varying"
                + " web down to minimum"
            )

        while dp_bar_probe > 0:
            probeWeb *= 2
            dp_bar_probe = _f_p_e_1(probeWeb)[0]

        e_1, _ = dekker(
            lambda web: _f_p_e_1(web)[0],
            probeWeb,  # >0
            0.5 * probeWeb,  # ?0
            # x_tol=1e-14,
            y_abs_tol=p_bar_d * tol,
        )  # this is the e_1 that satisifies the pressure specification.

        (p_bar_dev, Z_i, t_bar_i, l_bar_i, v_bar_i) = _f_p_e_1(e_1)

        if abs(Z_i - Z_b) < tol:
            """
            fudging the starting Z value for velocity integration to prevent
            driving integrator to 0 at the transition point.
            """
            Z_i = Z_b + tol

        if abs(p_bar_dev) > tol * p_bar_d:
            raise ValueError(
                "Design pressure is not met, delta = {:.3g} Pa".format(
                    p_bar_dev * (f * Delta)
                )
            )

        if v_j * v_bar_i > v_d and containBurnout:
            raise ValueError(
                "Design velocity exceeded ({:.4g} m/s > {:.4g} m/s) before peak pressure.".format(
                    v_bar_i * v_j, v_bar_d * v_j
                )
            )

        """
        step 2, find the requisite muzzle length to achieve design velocity
        """

        B = (
            S**2
            * e_1**2
            / (f * phi * omega * m * u_1**2)
            * (f * Delta) ** (2 * (1 - n))
        )

        def _ode_v(v_bar, t_bar, Z, l_bar, _):
            psi = f_psi_Z(Z)

            l_psi_bar = 1 - Delta / rho_p - Delta * (alpha - 1 / rho_p) * psi
            p_bar = (psi - v_bar**2) / (l_bar + l_psi_bar)

            if c_a_bar != 0 and v_bar > 0:
                v_r = v_bar / c_a_bar
                p_d_bar = (
                    0.25 * gamma_1 * (gamma_1 + 1) * v_r**2
                    + gamma_1
                    * v_r
                    * (1 + (0.25 * (gamma_1 + 1)) ** 2 * v_r**2) ** 0.5
                ) * p_a_bar

            else:
                p_d_bar = 0

            dt_bar = 2 / (theta * (p_bar - p_d_bar))

            if Z <= Z_b:
                dZ = dt_bar * (0.5 * theta / B) ** 0.5 * p_bar**n
            else:
                dZ = 0

            dl_bar = v_bar * dt_bar

            return [dt_bar, dZ, dl_bar]

        def abort_v(x, ys, record):
            _, _, l_bar = ys
            return l_bar > l_bar_d

        try:
            vtzl_record = [[v_bar_i, (t_bar_i, Z_i, l_bar_i)]]
            (v_bar_g, (t_bar_g, Z_g, l_bar_g), _) = RKF78(
                dFunc=_ode_v,
                iniVal=(t_bar_i, Z_i, l_bar_i),
                x_0=v_bar_i,
                x_1=v_bar_d,
                relTol=tol,
                absTol=tol**2,
                abortFunc=abort_v,
                record=vtzl_record,
            )

        except ValueError:
            v_bar_m, (t_bar_m, Z_m, l_bar_m) = vtzl_record[-1]
            pmax = _f_p_bar(Z_m, l_bar_m, v_bar_m) * f * Delta
            vmax = v_bar_m * v_j
            lmax = l_bar_m * l_0
            raise ValueError(
                "Integration appears to be approaching asymptote, "
                + "last calculated to v = {:.4g} m/s, ".format(vmax)
                + "x = {:.4g} m, p = {:.4g} MPa. ".format(lmax, pmax * 1e-6)
                + "This indicates an excessive velocity target relative to pressure developed."
            )

        p_bar_g = _f_p_bar(Z_g, l_bar_g, v_bar_g)

        v_g = v_bar_g * v_j
        l_g = l_bar_g * l_0
        p_g = p_bar_g * f * Delta
        if l_bar_g > l_bar_d:
            raise ValueError(
                "Solution requires excessive tube length, last calculated to "
                + "v = {:.4g} m/s, x = {:.4g} m, ".format(v_g, l_g)
                + "p = {:.4g} MPa.".format(p_g * 1e-6)
            )

        if abs(v_bar_g - v_bar_d) > (tol * v_bar_d):
            raise ValueError(
                "Velocity specification is not met, last calculated to "
                + "v = {:.4g} m/s ({:+.3g} %), x = {:.4g} m, p = {:.4g} MPa".format(
                    v_g, (v_bar_g - v_bar_d) / v_bar_d * 1e2, l_g, p_g * 1e-6
                )
            )

        # calculate the averaged chambrage correction factor
        # implied by this solution
        cc_n = 1 - (1 - 1 / chi_k) * log(l_bar_g + 1) / l_bar_g

        # lengtime = time.time()

        if abs(cc_n - cc) > tol:
            # successive better approximations will eventually
            # result in value within tolerance.
            return self.solve(
                loadFraction=loadFraction,
                chargeMassRatio=chargeMassRatio,
                tol=tol,
                minWeb=minWeb,
                containBurnout=containBurnout,
                maxLength=maxLength,
                labda_1=labda_1,
                labda_2=labda_2,
                sol=sol,
                cc=cc_n,
                ambientRho=ambientRho,
                ambientP=ambientP,
                ambientGamma=ambientGamma,
                control=control,
            )
            # TODO: Maximum recursion depth exceeded in comparison is
            # occasionally thrown here. Investigate why.
        else:
            return e_1, l_bar_g * l_0

    def findMinV(
        self,
        chargeMassRatio,
        tol,
        minWeb,
        maxLength,
        sol=SOL_PIDDUCK,
        ambientRho=1.204,
        ambientP=101.325e3,
        ambientGamma=1.4,
        loadFraction=None,  # if a load fraction value is passed, it is used as a hint
        control=POINT_PEAK_AVG,
        **_,
    ):
        """
        find the minimum volume solution.
        """

        """
        Step 1, find a valid range of values for load fraction,
        using psuedo-bisection.

        high lf -> high web
        low lf -> low web
        """
        m = self.m
        omega = m * chargeMassRatio
        rho_p = self.rho_p
        S = self.S
        solve = self.solve
        phi_1 = self.phi_1
        theta = self.theta
        gamma = theta + 1

        if sol == SOL_LAGRANGE:
            labda_1, labda_2 = 1 / 2, 1 / 3
        elif sol == SOL_PIDDUCK:
            labda_1, labda_2 = pidduck(omega / (phi_1 * m), gamma, tol)
        elif sol == SOL_MAMONTOV:
            labda_1, labda_2 = pidduck(omega / (phi_1 * m), 1, tol)
        else:
            raise ValueError("Unknown Solution")

        def f(lf):
            V_0 = omega / (rho_p * lf)
            l_0 = V_0 / S

            e_1, l_g = solve(
                loadFraction=lf,
                chargeMassRatio=chargeMassRatio,
                tol=tol,
                minWeb=minWeb,
                containBurnout=False,
                maxLength=maxLength,
                labda_1=labda_1,
                labda_2=labda_2,
                sol=sol,
                ambientRho=ambientRho,
                ambientP=ambientP,
                ambientGamma=ambientGamma,
                control=control,
            )

            return e_1, (l_g + l_0), l_g

        records = []

        for i in range(N):
            startProbe = (
                loadFraction
                if (i == 1 and loadFraction is not None)
                else uniform(tol, 1 - tol)
            )
            try:
                _, lt_i, lg_i = f(startProbe)
                records.append((startProbe, lt_i))
                break
            except ValueError:
                pass

        if i == N - 1:
            raise ValueError(
                "Unable to find any valid load fraction"
                + " with {:d} random samples.".format(N)
            )

        low = tol
        probe = startProbe
        delta_low = low - probe

        new_low = probe + delta_low

        while abs(2 * delta_low) > tol:
            try:
                _, lt_i, lg_i = f(new_low)
                records.append((new_low, lt_i))
                probe = new_low
            except ValueError:
                delta_low *= 0.5
            finally:
                new_low = probe + delta_low

        low = probe

        high = 1 - tol
        probe = startProbe
        delta_high = high - probe

        new_high = probe + delta_high

        while abs(2 * delta_high) > tol and new_high < 1:
            try:
                _, lt_i, lg_i = f(new_high)
                records.append((new_high, lt_i))
                probe = new_high
            except ValueError:
                delta_high *= 0.5
            finally:
                new_high = probe + delta_high

        high = probe

        if abs(high - low) < tol:
            raise ValueError("No range of values satisfying constraint.")

        if len(records) > 2:
            records.sort(key=lambda x: x[0])
            for l, m, h in zip(records[:-2], records[1:-1], records[2:]):
                if l[1] > m[1] and h[1] > m[1]:
                    low = l[0]
                    high = h[0]

        delta = high - low

        low += delta * tol
        high -= delta * tol

        """
        Step 2, gss to min.

        It was found that at this step, setting the accuracy metric
        on the x-value (or the load fraction) gives more consistent
        result than requriing a relative tolerance on the function
        values.
        """
        lf_low, lf_high = gss(
            lambda lf: f(lf)[1],
            low,
            high,
            x_tol=tol,  # variable: load fraction
            findMin=True,
        )

        lf = 0.5 * (lf_high + lf_low)

        e_1, l_t, l_g = f(lf)

        return lf, e_1, l_g


if __name__ == "__main__":
    import cProfile

    from prop import GrainComp, SimpleGeometry

    compositions = GrainComp.readFile("data/propellants.csv")
    S22 = compositions["ATK PRD(S)22"]
    M1 = compositions["M1"]
    S22S = Propellant(S22, SimpleGeometry.SPHERE, 1, 2.5)
    M1S = Propellant(M1, SimpleGeometry.SPHERE, 1, 2.5)
    test = Constrained(
        caliber=50e-3,
        shotMass=1,
        propellant=M1S,
        startPressure=10e6,
        dragCoefficient=5e-2,
        designPressure=350e6,
        designVelocity=1200,
        chambrage=1.5,
    )

    datas = []
    """
    pr = cProfile.Profile()
    pr.enable()
    """
    for i in range(5):
        datas.append(
            test.findMinV(
                chargeMassRatio=1,
                tol=1e-3,
                minWeb=1e-6,
                maxLength=1000,
                control=POINT_PEAK_SHOT,
            )
        )
    for i in range(5):
        datas.append(
            test.findMinV(
                chargeMassRatio=1,
                tol=1e-3,
                minWeb=1e-6,
                maxLength=1000,
                control=POINT_PEAK_AVG,
            )
        )
    for i in range(5):
        datas.append(
            test.findMinV(
                chargeMassRatio=1,
                tol=1e-3,
                minWeb=1e-6,
                maxLength=1000,
                control=POINT_PEAK_BREECH,
            )
        )
    """
    pr.disable()
    pr.print_stats(sort="time")
    """
    from tabulate import tabulate

    means = [sum(x) / len(datas) for x in zip(*datas)]

    delta = []
    for line in datas:
        delta.append((v - m) / m for v, m in zip(line, means))

    print(tabulate(datas, headers=("load fract.", "web", "length")))
    print(*means)
    print(
        tabulate(
            delta, headers=("load fract.", "web", "length"), floatfmt=".3e"
        )
    )