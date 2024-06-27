from numba import njit
from numba.core import types
import unittest


def domax3(a, b, c):
    return max(a, b, c)


def domin3(a, b, c):
    return min(a, b, c)


class TestMaxMin(unittest.TestCase):
    def test_max3(self):
        pyfunc = domax3
        argtys = (types.int32, types.float32, types.double)
        cfunc = njit(argtys)(pyfunc)

        a = 1
        b = 2
        c = 3

        self.assertEqual(pyfunc(a, b, c), cfunc(a, b, c))

    def test_min3(self):
        pyfunc = domin3
        argtys = (types.int32, types.float32, types.double)
        cfunc = njit(argtys)(pyfunc)

        a = 1
        b = 2
        c = 3

        self.assertEqual(pyfunc(a, b, c), cfunc(a, b, c))


if __name__ == '__main__':
    unittest.main()
