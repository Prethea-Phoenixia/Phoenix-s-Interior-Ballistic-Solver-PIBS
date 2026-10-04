import unittest

from pibs.ballistics.num.rkf import rkf45


class TestRKF(unittest.TestCase):
    def test_func(self):
        def df(x: float, ys: tuple[float]) -> tuple[float]:
            return (7 * ys[0] ** 2 * x**3,)

        x_1, y_1, abort = rkf45(df, (3.0,), 2, 0, rel_tol=1e-4)

        true_val = -1 / (7 / 4 * 0**4 - 85 / 3)

        self.assertAlmostEqual(y_1[0], true_val, places=4)


if __name__ == "__main__":
    unittest.main()
