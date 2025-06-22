from math import acos, atan2, cos, pi, sin
from multiprocessing import Pool
from random import uniform

from atmos import R_e, atmosphere
from drag import KdCurve
from num import RKF78, dekker, gss

USE_MP = True


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

    def _ode_t(self, t, x, y, vx, vy, _):
        """
        t: time of flight
        x, y: cooridnate in 2-d geocentric coordinate system:
        vx, vy: component of velocity
        """
        vsq = vx**2 + vy**2
        v = vsq**0.5
        r = (x**2 + y**2) ** 0.5
        h = r - R_e
        dx = vx
        dy = vy

        _, _, c, rho, g, _ = self.f_env(h)

        M = v / c  # mach
        Kd = self.f_Kd.Kd(M)  # drag coefficient
        a = Kd * rho * vsq / self.C
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

        kwargs = {
            "tol": tol,
            "vel": vel,
            "tgtR": Rs,
            "gunH": gunH,
            "tgtH": tgtH,
            "env": env,
            "t_max": t_max,
            "elev_min": elev_min,
            "elev_max": elev_max,
        }

        lTrajs, hTrajs = self.inverse(**kwargs)

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
                l = ((x - x_0) ** 2 + (y - y_0) ** 2) ** 0.5  # line of sight distance to target,
                drop = acos(((x - x_0) * cos(theta) + (y - y_0) * sin(theta)) / l) * (180 / pi)

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
            print(self.name + " @ {:.3f}m/s Emplaced at {:.1f}m ASL Target at {:}m ASL.".format(vel, gunH, tgtH))
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
        # package the context variables for easy calling of forwards.
        kwargs = {
            "tol": tol,
            "vel": vel,
            "gunH": gunH,
            "tgtH": tgtH,
            "env": env,
            "t_max": t_max,
        }
        elev_min = max(elev_min, -90)
        elev_max = min(elev_max, 90)

        N = 33

        for i in range(N):
            elev_probe = uniform(elev_min, elev_max)
            try:
                self.forward(**kwargs, elev=elev_probe, DESCEND=False)
                break
            except ValueError:
                pass

        else:
            # if i == N - 1:
            raise ValueError("No valid elevation can be found within specified in {:} samples".format(N))

        def bisect_serach(elev_tgt):
            try:
                self.forward(**kwargs, elev=elev_tgt)
                return elev_tgt
            except ValueError:
                pass

            delta = elev_tgt - elev_probe
            probe_i = elev_probe
            probe_j = probe_i + delta
            while abs(delta) > 3600**-1:
                try:
                    self.forward(**kwargs, elev=probe_j)
                    probe_i = probe_j
                except ValueError:
                    delta *= 0.5
                probe_j = probe_i + delta

            return probe_i

        elev_max = bisect_serach(elev_max)
        elev_min = bisect_serach(elev_min)

        def f_r(elev, r=0, DESCEND=True):
            try:
                record = self.forward(**kwargs, elev=elev, DESCEND=DESCEND)
                t, (x, y, vx, vy) = record[-1]
            except ValueError:
                # this is when the supplied elevation cannot arc the bullet
                # higher than the target plane
                return -r, None

            psi = -(atan2(y, x) - 0.5 * pi)
            gr = psi * R_e
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
        maximum range, or more technically a value less than 1/3600 deg above it.

        The corresponding range this elevation gives are the upper limit of the
        ascending solution, and the lower limit of the descending solution. Or
        in other words, the cresting elevation is the minimum elevation
        of all solutions.

        It is always desired to make this value small such that the exclusion zone
        where we cannot solve for either is minimal.

                  /- r_cre
               __+__
            _-/     \-_ des.soln
        ---+-----------+--------
         / asc.soln      \
        /                 \
        """

        r_opt_d, _ = f_r(elev_opt)

        # r_opt_a = f_r(elev_opt, DESCEND=False)[0]
        r_min_d = f_r(elev_min)[0]
        r_max_d = f_r(elev_max)[0]

        if gunH < tgtH:
            elev_cre_i, elev_cre_j = gss(
                lambda ang: f_r(ang, DESCEND=False)[0],
                elev_min,
                elev_max,
                x_tol=3600**-1,
                findMin=False,
            )
            elev_cre = 0.9 * elev_cre_i + 0.1 * elev_cre_j
            while f_r(elev_cre)[1] is None:
                elev_cre_i = elev_cre
                elev_cre = 0.9 * elev_cre_i + 0.1 * elev_cre_j

            elev_min = max(elev_min, elev_cre)
        r_min_d = f_r(elev_min, DESCEND=True)[0]
        r_min_a = f_r(elev_min, DESCEND=False)[0]

        r_max_a = f_r(elev_max, DESCEND=False)[0]

        if isinstance(tgtR, int) or isinstance(tgtR, float):
            tgtR = [tgtR]

        mpKwargs = {
            "forward": self.forward,
            "r_max_a": r_max_a,
            "r_max_d": r_max_d,
            "r_min_a": r_min_a,
            "r_min_d": r_min_d,
            "r_opt_d": r_opt_d,
            "elev_min": elev_min,
            "elev_max": elev_max,
            "elev_opt": elev_opt,
            "forward_kwargs": kwargs,
        }
        if USE_MP:
            with Pool() as p:
                lTrajs, hTrajs = zip(
                    *p.starmap(calc_wrapper, [(R, mpKwargs) for R in tgtR]),
                )

        else:
            lTrajs = []
            hTrajs = []

            for R in tgtR:
                if r_max_a < R < r_min_a:
                    elev_i, elev_j = dekker(
                        lambda ang: f_r(ang, R, DESCEND=False)[0],
                        elev_min,
                        elev_max,
                        x_tol=3600**-1,
                    )

                    l_elev = 0.5 * (elev_i + elev_j)
                    _, rec = f_r(l_elev, R, DESCEND=False)

                    lTrajs.append((l_elev, *rec))

                elif r_min_d < R < r_opt_d:
                    elev_i, elev_j = dekker(
                        lambda ang: f_r(ang, R)[0],
                        elev_min,
                        elev_opt,
                        x_tol=3600**-1,
                    )

                    l_elev = 0.5 * (elev_i + elev_j)
                    _, rec = f_r(l_elev, R)

                    lTrajs.append((l_elev, *rec))

                else:
                    lTrajs.append(None)

                if r_max_d < R < r_opt_d:
                    elev_i, elev_j = dekker(
                        lambda ang: f_r(ang, R)[0],
                        elev_opt,
                        elev_max,
                        x_tol=3600**-1,
                    )

                    h_elev = 0.5 * (elev_i + elev_j)
                    rec = f_r(h_elev, R)[1]
                    hTrajs.append((h_elev, *rec))

                else:
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

        if gunH == tgtH:
            gunH += tol

        if (not DESCEND) and gunH > tgtH:
            raise ValueError("No Ascending Solution Possible Given Gun Height" + " is Higher Than Target.")

        x_0, y_0 = 0, R_e + gunH
        theta = elev * pi / 180

        def abortTgt(x, ys, record):
            x, y, vx, vy = ys
            h = (x**2 + y**2) ** 0.5 - (R_e + tgtH)

            return (
                (h < 0 and ((-x * vx + y * vy) < 0))
                if DESCEND  # abort the calculation on downward crossing of target plane
                else ((h > 0 and ((-x * vx + y * vy) > 0)) or ((-x * vx + y * vy) < 0))
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
                + "Last Calculated at {:.3f}s, {:.3f}m ASL".format(t, (x**2 + y**2) ** 0.5 - R_e)
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

            t_prime = 0.5 * sum(gss(f_tgt, t_1, t_2, x_tol=max(t_2, 1) * tol, findMin=False))
            h_prime = f_tgt(t_prime)
            if h_prime > 0:
                if DESCEND:
                    # the new peak barely crest the target plane.
                    t_t, _ = dekker(f_tgt, t_prime, t_2, x_tol=max(t_2, 1) * tol)

                else:
                    t_t, _ = dekker(f_tgt, t_1, t_prime, x_tol=max(t_2, 1) * tol)

            else:
                # even the new peak point found cannot crest the target plane.
                raise ValueError("Projectile Cresting Below Target at {:.3f} m".format(h_prime))

        else:
            t_t, _ = dekker(f_tgt, t_1, t_2, x_tol=max(t_2, 1) * tol)

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


def calc_wrapper(R, kwargs):
    return calc_for_R(R=R, **kwargs)


def calc_for_R(
    forward,
    R,
    r_max_a,
    r_max_d,
    r_min_a,
    r_min_d,
    r_opt_d,
    elev_min,
    elev_max,
    elev_opt,
    forward_kwargs,
):
    def f_r(elev, r=0, DESCEND=True):
        try:
            record = forward(**forward_kwargs, elev=elev, DESCEND=DESCEND)
            t, (x, y, vx, vy) = record[-1]
        except ValueError:
            # this is when the supplied elevation cannot arc the bullet
            # higher than the target plane
            return -r, None

        psi = -(atan2(y, x) - 0.5 * pi)
        gr = psi * R_e

        return gr - r, record[-1]

    if r_max_a < R < r_min_a:
        elev_i, elev_j = dekker(
            lambda ang: f_r(ang, R, DESCEND=False)[0],
            elev_min,
            elev_max,
            x_tol=3600**-1,
        )

        l_elev = 0.5 * (elev_i + elev_j)
        _, rec = f_r(l_elev, R, DESCEND=False)

        lTraj = (l_elev, *rec)

    elif r_min_d < R < r_opt_d:
        elev_i, elev_j = dekker(
            lambda ang: f_r(ang, R)[0],
            elev_min,
            elev_opt,
            x_tol=3600**-1,
        )

        l_elev = 0.5 * (elev_i + elev_j)
        _, rec = f_r(l_elev, R)

        lTraj = (l_elev, *rec)

    else:
        lTraj = None

    if r_max_d < R < r_opt_d:
        elev_i, elev_j = dekker(
            lambda ang: f_r(ang, R)[0],
            elev_opt,
            elev_max,
            x_tol=3600**-1,
        )

        h_elev = 0.5 * (elev_i + elev_j)
        rec = f_r(h_elev, R)[1]
        hTraj = (h_elev, *rec)

    else:
        hTraj = None

    return lTraj, hTraj


if __name__ == "__main__":
    test = Bullet("test", mass=4.6, diam=27e-3, Kd_curve=KdCurve["M829"], form=1)
    # print(dec_to_dms(1 / 3600))
    """
    test.record_to_data(
        test.forward(
            tol=1e-3,
            vel=1819.92,
            elev=30,
            DESCEND=True,
        )
    )"""

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
    test.rangeTable(tol=1e-3, vel=5000, minR=0, maxR=200000, deltaR=2000, tgtH=100000)

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
