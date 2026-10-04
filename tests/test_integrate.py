import math
import unittest

from pibs.ballistics.num import intg


class TestIntg(unittest.TestCase):
    def test_cos(self):
        def f(x: float) -> float:
            return math.cos(x)

        v = intg(f, 0, 2 * math.pi, tol=1e-3)
        self.assertAlmostEqual(v, 0, places=3)

    def test_precision(self):
        def f(x: float) -> float:
            return math.exp(-(x**2))

        v = intg(f, 0, 1, tol=1e-3)
        self.assertAlmostEqual(v, math.erf(1) * math.sqrt(math.pi) / 2, places=3)


if __name__ == "__main__":
    unittest.main()
