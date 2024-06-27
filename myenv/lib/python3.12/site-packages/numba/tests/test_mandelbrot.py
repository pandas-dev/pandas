import unittest
from numba import njit
from numba.core import types


def is_in_mandelbrot(c):
    i = 0
    z = 0.0j
    for i in range(100):
        z = z ** 2 + c
        if (z.real * z.real + z.imag * z.imag) >= 4:
            return False
    return True


class TestMandelbrot(unittest.TestCase):

    def test_mandelbrot(self):
        pyfunc = is_in_mandelbrot
        cfunc = njit((types.complex64,))(pyfunc)

        points = [0+0j, 1+0j, 0+1j, 1+1j, 0.1+0.1j]
        for p in points:
            self.assertEqual(cfunc(p), pyfunc(p))


if __name__ == '__main__':
    unittest.main()
