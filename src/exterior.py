from atmos import atmosphere, R_e
from drag import KdCurve
from math import pi, sin, cos, atan2, acos

from num import RKF78, gss, bisect


def dec_to_dms(deg):
    sign = deg >= 0

    deg = abs(deg)
    h = int(deg)
    m = int((deg - h) * 60)
    s = int((deg - h - m / 60) * 3600)

    # return sign, h, m, s

    return "{:}{:03}°{:02}'{:02}\"".format("+" if sign else "-", h, m, s)


class Bullet:
    def __init__(self, name, mass, diam, form, Kd_curve):
        self.name = name
        self.W = mass
        self.D = diam
        self.i = form
        self.f_Kd = Kd_curve

        self.C = self.W / (self.i * self.D**2)  # ballistic coefficient

    def _ode_t(self, t, x, y, vx, vy):
        """
        t: time of flight
        x, y: cooridnate in 2-d geocentric coordinate system:
        vx, vy: component of velocity
        """

        v = (vx**2 + vy**2) ** 0.5
        r = (x**2 + y**2) ** 0.5
        h = r - R_e
        dx = vx
        dy = vy

        _, _, c, rho, g = self.f_env(h)

        M = v / c  # mach
        Kd = self.f_Kd.get(M)  # drag coefficient
        a = Kd * rho * v**2 / self.C
        dvx = -a * vx / v - g * x / r
        dvy = -a * vy / v - g * y / r

        return dx, dy, dvx, dvy

    def record_to_data(self, record, prettyprint=True):
        """convert a forward trajectory record into more sensible values:"""

        data = []
        for line in record:
            t, (x, y, vx, vy) = line
            r = (x**2 + y**2) ** 0.5
            h = r - R_e

            psi = -(atan2(y, x) - 0.5 * pi)
            """
            Geocentric angle measuring from shot start to shot spalsh
            shot
            start    shot
                +--_splash
                |  /
            R_e |θ/ R_e
                |/
                +
            """
            gr = psi * R_e
            v = (vx**2 + vy**2) ** 0.5

            """
            calculate the vector angle between r_vec and v_vec, which
            is the complementary angle of the shot-to-horizon
            """
            phi = 90 - acos((x * vx + y * vy) / (r * v)) * 180 / pi
            data.append((t, h, gr, v, phi))

        if prettyprint:
            from tabulate import tabulate

            print(
                tabulate(
                    data,
                    headers=("ToF", "H", "G.R.", "V.", "Angle"),
                )
            )

        return data

    def rangeTable(
        self,
        tol,
        vel,
        minR,
        maxR,
        deltaR,
        gunH=0,
        tgtH=0,
        env=atmosphere,
        t_max=1000,
        prettyprint=True,
        elev_min=-90,
        elev_max=90,
    ):
        R = minR
        Rs = []
        while R <= maxR:
            Rs.append(R)
            R += deltaR

        print(Rs)

        lTrajs, hTrajs = self.inverse(
            tol=tol,
            vel=vel,
            tgtR=Rs,
            gunH=gunH,
            tgtH=tgtH,
            env=env,
            t_max=t_max,
            elev_min=elev_min,
            elev_max=elev_max,
        )

        lTable, hTable = [], []

        x_0, y_0 = 0, R_e + gunH  # location of the gun.

        for R, lTraj in zip(Rs, lTrajs):
            if lTraj is not None:
                elev, t, (x, y, vx, vy) = lTraj
                theta = elev * pi / 180
                # geocentric height of the impact point
                r = (x**2 + y**2) ** 0.5
                # magnitude of velocity at impact
                v = (vx**2 + vy**2) ** 0.5
                # calculate the impact angle with relation to horizon at impact.
                phi = 90 - acos((x * vx + y * vy) / (r * v)) * 180 / pi

                """
                calculate the spatial angle between impact point and gun
                elevation, a.k.a bullet drop. This measure is probably point-
                less for larger guns.
                """
                l = (
                    (x - x_0) ** 2 + (y - y_0) ** 2
                ) ** 0.5  # line of sight distance to target,
                drop = acos(
                    ((x - x_0) * cos(theta) + (y - y_0) * sin(theta)) / l
                ) * (180 / pi)

                lTable.append(
                    (
                        R,
                        dec_to_dms(elev),
                        dec_to_dms(drop),
                        t,
                        v,
                        dec_to_dms(phi),
                    )
                )

        # repeat the above for the higher arcing trajectories.
        for R, hTraj in zip(Rs, hTrajs):
            if hTraj is not None:
                elev, t, (x, y, vx, vy) = hTraj
                r = (x**2 + y**2) ** 0.5
                v = (vx**2 + vy**2) ** 0.5
                phi = 90 - acos((x * vx + y * vy) / (r * v)) * 180 / pi
                hTable.append(
                    (
                        R,
                        dec_to_dms(elev),
                        None,
                        t,
                        v,
                        dec_to_dms(phi),
                    )
                )

        if prettyprint:
            from tabulate import tabulate

            headers = (
                "Ground\nRange m",
                "Launch\nElevation",
                "Drop\nAngle",
                "Time of\nFlight s",
                "Velocity\nm/s",
                "Impact\nAngle",
            )

            print("Low")
            print(tabulate(lTable, headers=headers))

            print("High")
            print(tabulate(hTable, headers=headers))

    def inverse(
        self,
        tol,
        vel,
        tgtR,
        gunH=0,
        tgtH=0,
        env=atmosphere,
        t_max=1000,
        elev_min=-90,
        elev_max=90,
    ):
        """
        Inverse calculation: given shot splash range, calculate in inverse the
        angle necessary to achieve said range.
        """
        elev_min = max(elev_min, -90)
        elev_max = min(elev_max, 90)

        def f_r(elev, r=0, DESCEND=True):
            try:
                record = self.forward(
                    tol=tol,
                    vel=vel,
                    elev=elev,
                    gunH=gunH,
                    tgtH=tgtH,
                    env=env,
                    t_max=t_max,
                    DESCEND=DESCEND,
                )
                t, (x, y, vx, vy) = record[-1]
            except ValueError as e:
                # this is when the supplied elevation cannot arc the bullet
                # higher than the target plane

                # print(elev, -r, DESCEND)
                return -r, None

            psi = -(atan2(y, x) - 0.5 * pi)
            gr = psi * R_e

            # print(elev, gr - r, DESCEND)

            return gr - r, record[-1]

        """we assume 0 deg 0 min 1 sec is the maximum practical precision in
        terms of elevation angle that we care about.
        """
        elev_opt = 0.5 * sum(
            gss(
                lambda ang: f_r(ang)[0],
                elev_min,
                elev_max,
                x_tol=3600**-1,
                findMin=False,
            )
        )
        """
        Find the cresting elevation that ascending solution achieves its
        maximum range. This is also the value which the elevation barely crests
        """
        elev_cre = 0.5 * sum(
            gss(
                lambda ang: f_r(ang, DESCEND=False)[0],
                elev_min,
                elev_max,
                x_tol=3600**-1,
                findMin=False,
            )
        )
        r_cre = f_r(elev_cre, DESCEND=False)[0]
        """the minimum range point is achieved either at minimum or maximum
        elevation specified. Therefore it suffice to compare the two:"""

        r_opt = f_r(elev_opt)[0]
        r_min = f_r(elev_min)[0]
        r_max = f_r(elev_max)[0]

        if isinstance(tgtR, int) or isinstance(tgtR, float):
            tgtR = [tgtR]

        # tgtR = [R for R in tgtR if R != 0]

        lTrajs = []
        hTrajs = []
        for R in tgtR:
            if r_min < R < r_cre:
                elev_i, elev_j = bisect(
                    lambda ang: f_r(ang, R, DESCEND=False)[0],
                    # elev_min,
                    elev_cre,
                    elev_opt,
                    x_tol=3600**-1,
                )
                l_elev = 0.5 * (elev_i + elev_j)
                _, rec = f_r(l_elev, R, DESCEND=False)
                lTrajs.append((l_elev, *rec))

            elif r_min < R < r_opt:
                elev_i, elev_j = bisect(
                    lambda ang: f_r(ang, R)[0],
                    elev_min,
                    elev_opt,
                    x_tol=3600**-1,
                )
                """
                if any((f_r(elev_j)[1] is None, f_r(elev_i)[1] is None)):
                    lTrajs.append(None)
                else:
                """
                l_elev = 0.5 * (elev_i + elev_j)
                _, rec = f_r(l_elev, R)
                lTrajs.append((l_elev, *rec))

            if r_max < R < r_opt:
                elev_i, elev_j = bisect(
                    lambda ang: f_r(ang, R)[0],
                    elev_opt,
                    elev_max,
                    x_tol=3600**-1,
                )
                """
                We reject any solution occuring within 1/3600 degree, in
                elevation, of the cresting point for descending solutions.
                """
                if any((f_r(elev_j)[1] is None, f_r(elev_i)[1] is None)):
                    hTrajs.append(None)
                else:
                    h_elev = 0.5 * (elev_i + elev_j)
                    rec = f_r(h_elev, R)[1]
                    hTrajs.append((h_elev, *rec))

            else:
                lTrajs.append(None)
                hTrajs.append(None)

        return lTrajs, hTrajs

    def forward(
        self,
        tol,
        vel,
        elev,
        gunH=0,
        tgtH=0,
        env=atmosphere,
        t_max=1000,
        DESCEND=True,
    ):
        """
        Forward calculation: given ballistic parameters, determine the
        trajector flown by the shot.
        """
        self.f_env = env

        """
         dv
        ---- = Kd(Mach) * rho  * V^2 / C
         dt

         C = M / (i D^2)
        """

        if (not DESCEND) and gunH > tgtH:
            raise ValueError(
                "No Ascending Solution Possible Given Gun Height"
                + " is Higher Than Target."
            )
        if gunH == tgtH:
            gunH += tol

        x_0, y_0 = 0, R_e + gunH
        theta = elev * pi / 180

        def abortTgt(x, ys, o_x, o_ys):
            x, y, vx, vy = ys
            h = (x**2 + y**2) ** 0.5 - (R_e + tgtH)

            return (
                (h < 0 and ((-x * vx + y * vy) < 0))
                if DESCEND  # abort the calculation on downward crossing of target plane
                else (
                    (h > 0 and ((-x * vx + y * vy) > 0))
                    or ((-x * vx + y * vy) < 0)
                )
            )

        vx_0 = vel * cos(theta)
        vy_0 = vel * sin(theta)

        record = [[0, [x_0, y_0, vx_0, vy_0]]]

        try:
            t_2, (x_2, y_2, _, _), _ = RKF78(
                self._ode_t,
                (x_0, y_0, vx_0, vy_0),
                0,
                t_max,
                relTol=tol,
                abortFunc=abortTgt,
                adaptTo=(False, False, True, True),
                record=record,
            )  # Coarse integration to maximum time to find approximate ToF
        except ValueError:
            print("forward")
            print(*record, sep="\n")  # debug code
            t, (x, y, vx, vy) = record[-1]
            raise ValueError(
                "Projectile Exited Environment Function Range\n"
                + "Last Calculated at {:.3f}s, {:.3f}m ASL".format(
                    t, (x**2 + y**2) ** 0.5 - R_e
                )
            )
        if t_2 == t_max:
            raise ValueError("Projectile Maximum Time-of-Flight t_max Exceeded")

        if len(record) > 1:
            t_1, (x_1, y_1, vx_1, vy_1) = record[-1]

        else:
            t_1, (x_1, y_1, vx_1, vy_1) = record[0]

        def f_tgt(t):
            _, (x, y, _, _), _ = RKF78(
                self._ode_t,
                (x_1, y_1, vx_1, vy_1),
                t_1,
                t,
                relTol=tol,
                adaptTo=(False, False, True, True),
            )  # fine integration from last point before impact
            return (x**2 + y**2) ** 0.5 - (R_e + tgtH)

        h_1 = (x_1**2 + y_1**2) ** 0.5 - (R_e + tgtH)
        h_2 = (x_2**2 + y_2**2) ** 0.5 - (R_e + tgtH)

        if (DESCEND and (h_1 < 0)) or (not DESCEND and (h_2 < 0)):
            """
            then we are probably cresting below the target plane.
            In case DESCNED:

                this--->
                 |    _++_
            h>0  | +1+    +-+
            -----|/----------\-------------
            h<0  /1   _+_    2\
                     1   2
                    /     \

            In case ASCEND:
                        <----this
                      _++_    |
            h>0    +2+    +-+ |
            ------/----------\|-------------
            h<0  /1   _+_    2\
                     1   2
                    /     \


            In this case:
            if we want to find the descending solution, then we need to try
            and raise the point t_1 to above the target plane, and determine
            time of impact between new time at peak, t_prime and t_2

            if we want to find the ascending solution, then we need to try
            and raise the point t_2 to above the target plane, and determine
            time of impact between t_1 and new time at peak, t_prime

            This is a very edge case scenario but nevertheless
            in the name of accuracy and rigouroness it needs
            be done.
            """

            t_prime = 0.5 * sum(
                gss(f_tgt, t_1, t_2, x_tol=max(t_2, 1) * tol, findMin=False)
            )
            h_prime = f_tgt(t_prime)
            if h_prime > 0:
                if DESCEND:
                    # the new peak barely crest the target plane.
                    t_t = 0.5 * sum(
                        bisect(f_tgt, t_prime, t_2, x_tol=max(t_2, 1) * tol)
                    )
                else:
                    t_t = 0.5 * sum(
                        bisect(f_tgt, t_1, t_prime, x_tol=max(t_2, 1) * tol)
                    )

            else:
                # even the new peak point found cannot crest the target plane.
                raise ValueError(
                    "Projectile Cresting Below Target at {:.3f} m".format(
                        h_prime
                    )
                )

        else:
            t_t = 0.5 * sum(bisect(f_tgt, t_1, t_2, x_tol=max(t_2, 1) * tol))

        _, _, _ = RKF78(
            self._ode_t,
            (x_1, y_1, vx_1, vy_1),
            t_1,
            t_t,
            relTol=tol,
            adaptTo=(False, False, True, True),
            record=record,
        )

        return record


if __name__ == "__main__":
    test = Bullet(
        "test", mass=9.0990629, diam=88e-3, Kd_curve=KdCurve["G8"], form=0.925
    )

    test.record_to_data(
        test.forward(tol=1e-3, vel=819.92, elev=1.419, tgtH=10, DESCEND=False)
    )
    # input()
    """
    print(
        *test.inverse(
            tol=1e-3,
            vel=819.92,
            tgtR=[0, 100, 200, 300, 400, 500],
            tgtH=10,
        ),
        sep="\n"
    )
    """
    test.rangeTable(
        tol=1e-6, vel=819.92, minR=0, maxR=5000, deltaR=250, tgtH=10
    )

    """
    test = Bullet(
        "M2 ball",
        mass=0.045942,
        form=0.86,
        diam=12.7e-3,
        Kd_curve=KdCurve["G5"],
    )

    test.record_to_data(test.forward(tol=1e-9, vel=856 * 0.99, elev=0.725 * 2))
    """
