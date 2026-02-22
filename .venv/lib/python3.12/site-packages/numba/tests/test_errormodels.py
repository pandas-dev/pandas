"""
Test setting/overriding error models
"""

from numba import jit
import unittest


class TestErrorModel(unittest.TestCase):

    def test_div_by_zero_python(self):
        @jit   # python model is the default
        def model_python(val):
            return 1 / val

        with self.assertRaises(ZeroDivisionError):
            model_python(0)

    def test_div_by_zero_numpy(self):
        @jit(error_model='numpy')
        def model_numpy(val):
            return 1 / val

        self.assertEqual(model_numpy(0), float('inf'))


if __name__ == '__main__':
    unittest.main()
