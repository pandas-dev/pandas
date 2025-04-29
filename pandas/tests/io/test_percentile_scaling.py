import unittest
from pandas.io.percentile_scaling import percentile_scaling

class TestPercentileScaling(unittest.TestCase):
    def test_scaling(self):
        data = [10, 20, 30, 40, 50]
        expected = [0.0, 25.0, 50.0, 75.0, 100.0]
        result = percentile_scaling(data)
        for r, e in zip(result, expected):
            self.assertAlmostEqual(r, e)

    def test_identical_values(self):
        with self.assertRaises(ValueError):
            percentile_scaling([5, 5, 5])

    def test_empty(self):
        with self.assertRaises(ValueError):
            percentile_scaling([])

if __name__ == "__main__":
    unittest.main()
