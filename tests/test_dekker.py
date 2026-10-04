import unittest

from pibs.ballistics.num import dekker


class TestDekker(unittest.TestCase):
    def test_cubic(self):
        def f(x: float) -> float:
            return x**3

        self.assertAlmostEqual(dekker(f, -1, 1), 0)

    def test_quadratic(self):
        def f(x: float) -> float:
            return (x - 1) ** 2 - 1

        x_0 = dekker(f, -1, 1)
        self.assertAlmostEqual(x_0, 0)


if __name__ == "__main__":
    unittest.main()
