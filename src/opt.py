from num import GSS, RKF78, cubic, secant
from prop import Propellant
from random import uniform
from math import pi


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
        **_,
    ):
        # constants for constrained designs

        if any(
            (
                caliber <= 0,
                shotMass <= 0,
                startPressure <= 0,
                dragCoefficient < 0,
                designPressure <= 0,
                designVelocity <= 0,
            )
        ):
            raise ValueError("Invalid parameters for constrained design")

        self.S = (caliber / 2) ** 2 * pi
        self.m = shotMass
        self.propellant = propellant
        self.p_0 = startPressure
        self.phi_1 = 1 + dragCoefficient

        # design limits
        self.p_d = designPressure
        self.v_d = designVelocity

    def __getattr__(self, attrName):
        try:
            return getattr(self.propellant, attrName)
        except:
            raise AttributeError("object has no '%s'" % attrName)

    def solve(
        self,
        loadFraction,
        chargeMassRatio,
        tol,
        minWeb=1e-6,
        containBurnout=True,
        maxLength=1e3,
        **_,
    ):
        if any(
            (
                minWeb <= 0,
                tol <= 0,
                chargeMassRatio <= 0,
                loadFraction <= 0,
                loadFraction > 1,
                maxLength <= 0,
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

        f_psi_Z = self.f_psi_Z
        f_sigma_Z = self.f_sigma_Z

        omega = m * chargeMassRatio
        V_0 = omega / (rho_p * maxLF * loadFraction)
        Delta = omega / V_0
        p_bar_0 = p_0 / (Delta * f)
        l_0 = V_0 / S
        phi = phi_1 + omega / (3 * m)
        """
        it is impossible to account for the chamberage effect given unspecified
        barrel length, in our formulation
        """
        v_j = (2 * f * omega / (theta * phi * m)) ** 0.5

        if v_j < v_d:
            raise ValueError(
                "Propellant load too low to achieve design velocity."
            )

        psi_0 = (1 / Delta - 1 / rho_p) / (f / p_0 + alpha - 1 / rho_p)

        Zs = cubic(
            a=chi * mu,
            b=chi * labda,
            c=chi,
            d=-psi_0,
        )
        # pick a valid solution between 0 and 1
        Zs = tuple(
            Z
            for Z in Zs
            if not isinstance(Z, complex) and (Z > 0.0 and Z < 1.0)
        )  # evaluated from left to right, guards against complex >/< float
        if len(Zs) < 1:
            raise ValueError(
                "Propellant either could not develop enough pressure to overcome"
                + " start pressure, or has burnt to post fracture."
            )
        Z_0 = Zs[0]

        def _fp_bar(Z, l_bar, v_bar):
            psi = f_psi_Z(Z)
            l_psi_bar = 1 - Delta / rho_p - Delta * (alpha - 1 / rho_p) * psi

            p_bar = (
                f * omega * psi - 0.5 * theta * phi * m * (v_bar * v_j) ** 2
            ) / (S * l_0 * (l_bar + l_psi_bar) * f * Delta)

            return p_bar

        """
        step 1, find grain size that satisifies design pressure
        """
        p_bar_d = p_d / (f * Delta)  # convert to unitless
        l_bar_d = maxLength / l_0

        def _fp_e_1(e_1, tol):
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

            # integrate this to end of burn

            def _ode_Z(Z, t_bar, l_bar, v_bar):
                """burnout domain ode of internal ballistics"""
                psi = f_psi_Z(Z)
                dpsi = f_sigma_Z(Z)  # dpsi/dZ

                l_psi_bar = (
                    1 - Delta / rho_p - Delta * (alpha - 1 / rho_p) * psi
                )
                p_bar = (
                    f * omega * psi - 0.5 * theta * phi * m * (v_bar * v_j) ** 2
                ) / (S * l_0 * (l_bar + l_psi_bar) * f * Delta)

                if Z <= Z_b:
                    dt_bar = (2 * B / theta) ** 0.5 * p_bar**-n  # dt_bar/dZ

                    dl_bar = v_bar * dt_bar

                    dv_bar = 0.5 * dt_bar * theta * p_bar

                else:
                    dt_bar = 0
                    dl_bar = 0
                    dv_bar = 0

                return [dt_bar, dl_bar, dv_bar]

            def abort(x, ys, o_x, o_ys):
                Z = x
                t_bar, l_bar, v_bar = ys

                p_bar = _fp_bar(Z, l_bar, v_bar)

                oZ = o_x
                ot_bar, ol_bar, ov_bar = o_ys

                op_bar = _fp_bar(oZ, ol_bar, ov_bar)

                return (
                    p_bar < op_bar or p_bar > 2 * p_bar_d
                )  # or l_bar > l_bar_d

            record = []

            (Z_j, (t_bar_j, l_bar_j, v_bar_j), (_, _, _)) = RKF78(
                dFunc=_ode_Z,
                iniVal=(0, 0, 0),
                x_0=Z_0,
                x_1=Z_b,
                relTol=0.1 * tol,
                absTol=0.1 * tol,
                abortFunc=abort,
                record=record,
            )

            p_bar_j = _fp_bar(Z_j, l_bar_j, v_bar_j)

            if len(record) > 1:
                Z_i, (t_bar_i, l_bar_i, v_bar_i) = record[-2]
            else:
                Z_i, (t_bar_i, l_bar_i, v_bar_i) = Z_0, (0, 0, 0)

            def _fp_Z(Z):
                _, (t_bar, l_bar, v_bar), (_, _, _) = RKF78(
                    dFunc=_ode_Z,
                    iniVal=(t_bar_i, l_bar_i, v_bar_i),
                    x_0=Z_i,
                    x_1=Z,
                    relTol=0.01 * tol,
                    absTol=0.01 * tol,
                )

                return _fp_bar(Z, l_bar, v_bar)

            """
            find the peak pressure point.
            The tolerance on Z must be guarded such that floating point
            error does not become rampant
            """
            Z_1, Z_2 = GSS(
                _fp_Z,
                Z_i,
                Z_j,
                yRelTol=0.1 * tol,
                findMin=False,
                xTol=0,
            )

            Z_p = 0.5 * (Z_1 + Z_2)

            return _fp_Z(Z_p) - p_bar_d, Z_j, v_bar_i

        """
        The two initial guesses are good enough for the majority of cases,
        guess one: 0.1mm, guess two: 1mm
        """
        dp_bar_probe = _fp_e_1(minWeb, tol)[0]
        probeWeb = minWeb

        if dp_bar_probe < 0:
            raise ValueError(
                "Design pressure cannot be achieved by varying web down to minimum"
            )

        while dp_bar_probe > 0:
            probeWeb *= 2
            dp_bar_probe = _fp_e_1(probeWeb, tol)[0]

        e_1, _ = secant(
            lambda x: _fp_e_1(x, tol)[0],
            probeWeb,  # >0
            0.5 * probeWeb,  # ?0
            tol=p_bar_d * tol,
            x_min=0.5 * probeWeb,  # <=0
        )  # this is the e_1 that satisifies the pressure specification.

        p_bar_dev, Z_i, v_bar_i = _fp_e_1(e_1, tol)

        if abs(p_bar_dev / p_bar_d) > tol:
            raise ValueError("Design pressure is not met")

        """
        step 2, find the requisite muzzle length to achieve design velocity
        """
        v_bar_d = v_d / v_j

        if v_bar_i > v_bar_d and containBurnout:
            raise ValueError("Design velocity exceeded before peak pressure")

        # TODO: find some way of making this cross constraint less troublesome.
        B = (
            S**2
            * e_1**2
            / (f * phi * omega * m * u_1**2)
            * (f * Delta) ** (2 * (1 - n))
        )

        def _ode_v(v_bar, t_bar, Z, l_bar):
            # p_bar = _fp_bar(Z, l_bar, v_bar)
            psi = f_psi_Z(Z)
            dpsi = f_sigma_Z(Z)  # dpsi/dZ

            l_psi_bar = 1 - Delta / rho_p - Delta * (alpha - 1 / rho_p) * psi
            # dp_bar/dt_bar
            p_bar = (
                f * omega * psi - 0.5 * theta * phi * m * (v_bar * v_j) ** 2
            ) / (S * l_0 * (l_bar + l_psi_bar) * f * Delta)

            if Z <= Z_b:
                dZ = (2 / (B * theta)) ** 0.5 * p_bar ** (n - 1)
            else:
                dZ = 0
            dl_bar = 2 * v_bar / (theta * p_bar)
            dt_bar = 2 / (theta * p_bar)
            return (dt_bar, dZ, dl_bar)

        """
        Integrating from 0 enforce consistency and improves numerical
        stability of the result when called with inputs that are in close
        proximity.
        """

        def abort(x, ys, o_x, o_ys):
            v_bar = x
            t_bar, Z, l_bar = ys
            return l_bar > l_bar_d

        (v_bar_g, (t_bar_g, Z_g, l_bar_g), (_, _, _)) = RKF78(
            dFunc=_ode_v,
            iniVal=(0, Z_0, 0),
            x_0=0,
            x_1=v_bar_d,
            relTol=tol,
            absTol=tol,
            abortFunc=abort,
        )
        if l_bar_g > l_bar_d:
            raise ValueError(
                "Solution requires excessive tube length ({:.3e} m)".format(
                    maxLength
                )
            )
        # print("p
        if abs(v_bar_g - v_bar_d) / (v_bar_d) > tol:
            raise ValueError("Velocity specification is not met")

        return e_1, l_bar_g * l_0

    def findMinV(self, chargeMassRatio, tol, minWeb, maxLength, **_):
        """
        find the minimum volume solution.
        """

        """
        Step 1, find a valid range of values for load fraction, 
        using psuedo-bisection.

        high lf -> high web
        low lf -> low web
        """
        omega = self.m * chargeMassRatio
        rho_p = self.rho_p
        maxLF = self.maxLF
        S = self.S
        solve = self.solve

        def f(lf, mW=minWeb):
            V_0 = omega / (rho_p * maxLF * lf)
            l_0 = V_0 / S

            e_1, l_g = solve(
                loadFraction=lf,
                chargeMassRatio=chargeMassRatio,
                tol=0.1 * tol,  # this is to ensure unimodality up to ~tol
                minWeb=mW,
                containBurnout=False,
                maxLength=maxLength,
            )
            # print(l_g + l_0)
            return e_1, (l_g + l_0), l_g

        records = []
        N = 33
        for i in range(N):
            startProbe = uniform(tol, 1 - tol)
            try:
                _, lt_i, lg_i = f(startProbe)
                records.append((startProbe, lt_i))
                break
            except ValueError:
                pass

        if i == N - 1:
            raise ValueError(
                "Unable to find any valid load fraction with {:d} random samples.".format(
                    N
                )
            )

        low = tol
        probe = startProbe
        delta_low = low - probe

        web_i = minWeb
        new_low = probe + delta_low

        while abs(2 * delta_low) > tol:
            try:
                web_i, lt_i, lg_i = f(new_low)
                records.append((new_low, lt_i))
                probe = new_low
            except ValueError as e:
                delta_low *= 0.5
            finally:
                new_low = probe + delta_low

        actMinWeb = web_i * (1 - tol)

        low = probe

        high = 1 - tol
        probe = startProbe
        delta_high = high - probe

        new_high = probe + delta_high
        while abs(2 * delta_high) > tol and new_high < 1:
            try:
                _, lt_i, lg_i = f(new_high, actMinWeb)
                records.append((new_high, lt_i))
                probe = new_high
            except ValueError as e:
                delta_high *= 0.5
            finally:
                new_high = probe + delta_high

        high = probe

        if abs(high - low) < tol:
            raise ValueError("No range of values satisfying constraint.")

        if len(records) > 3:
            records.sort(key=lambda x: x[0])
            for l, m, h in zip(records[:-2], records[1:-1], records[2:]):
                if l[1] > m[1] and h[1] > m[1]:
                    low = l[0]
                    high = h[0]

        delta = high - low

        # Edge values are some times only semi-stable, i.e. when calling
        # f() with the same value will spuriously raise value errors. Therefore
        # we conservatively shrink the range by tolerance to avoid this issue.

        low += delta * tol
        high -= delta * tol

        """
        Step 2, gss to min.
        """

        lf_low, lf_high = GSS(
            lambda x: f(x, actMinWeb)[1],
            low,
            high,
            yRelTol=tol,
            xTol=0,
            findMin=True,
        )

        lf = 0.5 * (lf_high + lf_low)

        e_1, l_t, l_g = f(lf, actMinWeb)

        return lf, e_1, l_g


if __name__ == "__main__":
    from prop import GrainComp, MultPerfGeometry, SimpleGeometry

    compositions = GrainComp.readFile("data/propellants.csv")
    S22 = compositions["ATK PRD(S)22"]
    S22S = Propellant(S22, SimpleGeometry.SPHERE, 1, 2.5)

    test = Constrained(
        caliber=50e-3,
        shotMass=1,
        propellant=S22S,
        startPressure=10e6,
        dragCoefficient=5e-2,
        designPressure=350e6,
        designVelocity=1500,
    )

    for i in range(10):
        print(test.findMinV(chargeMassRatio=1, tol=1e-3, minWeb=1e-6))
