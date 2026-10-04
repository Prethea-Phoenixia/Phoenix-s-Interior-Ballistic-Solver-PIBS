import unittest

from pibs.ballistics.num import gss


class TestGSS(unittest.TestCase):
    def test_quadratic(self):
        def f(x: float) -> float:
            return x**2

        x_0 = gss(f, -1, 1)
        self.assertAlmostEqual(x_0, 0)


if __name__ == "__main__":
    unittest.main()
