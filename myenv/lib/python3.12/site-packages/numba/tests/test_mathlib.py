import itertools
import math
import sys
import unittest
import warnings

import numpy as np

from numba import njit, types
from numba.tests.support import TestCase
from numba.np import numpy_support


def sin(x):
    return math.sin(x)


def cos(x):
    return math.cos(x)


def tan(x):
    return math.tan(x)


def sinh(x):
    return math.sinh(x)


def cosh(x):
    return math.cosh(x)


def tanh(x):
    return math.tanh(x)


def asin(x):
    return math.asin(x)


def acos(x):
    return math.acos(x)


def atan(x):
    return math.atan(x)


def atan2(y, x):
    return math.atan2(y, x)


def asinh(x):
    return math.asinh(x)


def acosh(x):
    return math.acosh(x)


def atanh(x):
    return math.atanh(x)


def sqrt(x):
    return math.sqrt(x)


def npy_sqrt(x):
    return np.sqrt(x)


def exp(x):
    return math.exp(x)


def expm1(x):
    return math.expm1(x)


def log(x):
    return math.log(x)


def log1p(x):
    return math.log1p(x)


def log10(x):
    return math.log10(x)


def log2(x):
    return math.log2(x)


def floor(x):
    return math.floor(x)


def ceil(x):
    return math.ceil(x)


def trunc(x):
    return math.trunc(x)


def isnan(x):
    return math.isnan(x)


def isinf(x):
    return math.isinf(x)


def isfinite(x):
    return math.isfinite(x)


def hypot(x, y):
    return math.hypot(x, y)


def nextafter(x, y):
    return math.nextafter(x, y)


def degrees(x):
    return math.degrees(x)


def radians(x):
    return math.radians(x)


def erf(x):
    return math.erf(x)


def erfc(x):
    return math.erfc(x)


def gamma(x):
    return math.gamma(x)


def lgamma(x):
    return math.lgamma(x)


def pow(x, y):
    return math.pow(x, y)

def gcd(x, y):
    return math.gcd(x, y)

def copysign(x, y):
    return math.copysign(x, y)


def frexp(x):
    return math.frexp(x)


def ldexp(x, e):
    return math.ldexp(x, e)


def get_constants():
    return math.pi, math.e


class TestMathLib(TestCase):

    def test_constants(self):
        cfunc = njit(get_constants)
        self.assertPreciseEqual(cfunc(), cfunc.py_func())

    def run_unary(self, pyfunc, x_types, x_values, prec='exact', **kwargs):
        cfunc = njit(pyfunc)
        for tx, vx in zip(x_types, x_values):
            got = cfunc(vx)
            expected = pyfunc(vx)
            actual_prec = 'single' if tx is types.float32 else prec
            msg = 'for input %r' % (vx,)
            self.assertPreciseEqual(got, expected, prec=actual_prec, msg=msg,
                                    **kwargs)

    def run_binary(self, pyfunc, x_types, x_values, y_values, prec='exact'):
        cfunc = njit(pyfunc)
        for ty, x, y in zip(x_types, x_values, y_values):
            got = cfunc(x, y)
            expected = pyfunc(x, y)
            actual_prec = 'single' if ty is types.float32 else prec
            msg = 'for inputs (%r, %r)' % (x, y)
            self.assertPreciseEqual(got, expected, prec=actual_prec, msg=msg)

    def check_predicate_func(self, pyfunc):
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float32, types.float32,
                   types.float64, types.float64, types.float64]
        x_values = [0, 0, 0, 0, 0, 0,
                    float('inf'), 0.0, float('nan'),
                    float('inf'), 0.0, float('nan')]
        self.run_unary(pyfunc, x_types, x_values)

    def test_sin(self):
        pyfunc = sin
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [-2, -1, -2, 2, 1, 2, .1, .2]
        self.run_unary(pyfunc, x_types, x_values)

    @unittest.skipIf(sys.platform == 'win32',
                     "not exactly equal on win32 (issue #597)")
    def test_cos(self):
        pyfunc = cos
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [-2, -1, -2, 2, 1, 2, .1, .2]
        self.run_unary(pyfunc, x_types, x_values)

    def test_tan(self):
        pyfunc = tan
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [-2, -1, -2, 2, 1, 2, .1, .2]
        self.run_unary(pyfunc, x_types, x_values)

    def test_sqrt(self):
        pyfunc = sqrt
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [2, 1, 2, 2, 1, 2, .1, .2]
        self.run_unary(pyfunc, x_types, x_values)

    def test_npy_sqrt(self):
        pyfunc = npy_sqrt
        x_values = [2, 1, 2, 2, 1, 2, .1, .2]
        # XXX poor precision for int16 inputs
        x_types = [types.int16, types.uint16]
        self.run_unary(pyfunc, x_types, x_values, prec='single')
        x_types = [types.int32, types.int64,
                   types.uint32, types.uint64,
                   types.float32, types.float64]
        self.run_unary(pyfunc, x_types, x_values)

    def test_exp(self):
        pyfunc = exp
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [-2, -1, -2, 2, 1, 2, .1, .2]
        self.run_unary(pyfunc, x_types, x_values)

    def test_expm1(self):
        pyfunc = expm1
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [-2, -1, -2, 2, 1, 2, .1, .2]
        self.run_unary(pyfunc, x_types, x_values)

    def test_log(self):
        pyfunc = log
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 10, 100, 1000, 100000, 1000000, 0.1, 1.1]
        self.run_unary(pyfunc, x_types, x_values)

    def test_log1p(self):
        pyfunc = log1p
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 10, 100, 1000, 100000, 1000000, 0.1, 1.1]
        self.run_unary(pyfunc, x_types, x_values)

    def test_log10(self):
        pyfunc = log10
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 10, 100, 1000, 100000, 1000000, 0.1, 1.1]
        self.run_unary(pyfunc, x_types, x_values)

    def test_log2(self):
        pyfunc = log2
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 10, 100, 1000, 100000, 1000000, 0.1, 1.1]
        self.run_unary(pyfunc, x_types, x_values)

    def test_asin(self):
        pyfunc = asin
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 1, 1, 1, 1, 1, 1., 1.]
        self.run_unary(pyfunc, x_types, x_values)

    def test_acos(self):
        pyfunc = acos
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 1, 1, 1, 1, 1, 1., 1.]
        self.run_unary(pyfunc, x_types, x_values)

    def test_atan(self):
        pyfunc = atan
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [-2, -1, -2, 2, 1, 2, .1, .2]
        self.run_unary(pyfunc, x_types, x_values)

    def test_atan2(self):
        pyfunc = atan2
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [-2, -1, -2, 2, 1, 2, .1, .2]
        y_values = [x * 2 for x in x_values]
        self.run_binary(pyfunc, x_types, x_values, y_values)

    def test_asinh(self):
        pyfunc = asinh
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 1, 1, 1, 1, 1, 1., 1.]
        self.run_unary(pyfunc, x_types, x_values, prec='double')

    def test_acosh(self):
        pyfunc = acosh
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 1, 1, 1, 1, 1, 1., 1.]
        self.run_unary(pyfunc, x_types, x_values)

    def test_atanh(self):
        pyfunc = atanh
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [0, 0, 0, 0, 0, 0, 0.1, 0.1]
        self.run_unary(pyfunc, x_types, x_values, prec='double')

    def test_sinh(self):
        pyfunc = sinh
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 1, 1, 1, 1, 1, 1., 1.]
        self.run_unary(pyfunc, x_types, x_values)

    def test_cosh(self):
        pyfunc = cosh
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 1, 1, 1, 1, 1, 1., 1.]
        self.run_unary(pyfunc, x_types, x_values)

    def test_tanh(self):
        pyfunc = tanh
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [0, 0, 0, 0, 0, 0, 0.1, 0.1]
        self.run_unary(pyfunc, x_types, x_values)

    def test_floor(self):
        pyfunc = floor
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [0, 0, 0, 0, 0, 0, 0.1, 1.9]
        self.run_unary(pyfunc, x_types, x_values)

    def test_ceil(self):
        pyfunc = ceil
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [0, 0, 0, 0, 0, 0, 0.1, 1.9]
        self.run_unary(pyfunc, x_types, x_values)

    def test_trunc(self):
        pyfunc = trunc
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [0, 0, 0, 0, 0, 0, 0.1, 1.9]
        self.run_unary(pyfunc, x_types, x_values)

    def test_isnan(self):
        self.check_predicate_func(isnan)

    def test_isinf(self):
        self.check_predicate_func(isinf)

    def test_isfinite(self):
        self.check_predicate_func(isfinite)

    def test_hypot(self):
        pyfunc = hypot
        x_types = [types.int64, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 2, 3, 4, 5, 6, .21, .34]
        y_values = [x + 2 for x in x_values]
        # Issue #563: precision issues with math.hypot() under Windows.
        prec = 'single'
        self.run_binary(pyfunc, x_types, x_values, y_values, prec)
        # Check that values that overflow in naive implementations do not
        # in the numba impl

        def naive_hypot(x, y):
            return math.sqrt(x * x + y * y)

        cfunc = njit(pyfunc)
        for fltty in (types.float32, types.float64):
            dt = numpy_support.as_dtype(fltty).type
            val = dt(np.finfo(dt).max / 30.)
            nb_ans = cfunc(val, val)
            self.assertPreciseEqual(nb_ans, pyfunc(val, val), prec='single')
            self.assertTrue(np.isfinite(nb_ans))

            with warnings.catch_warnings():
                warnings.simplefilter("error", RuntimeWarning)
                self.assertRaisesRegex(RuntimeWarning,
                                        'overflow encountered in .*scalar',
                                        naive_hypot, val, val)

    def test_nextafter(self):
        pyfunc = nextafter
        x_types = [types.float32, types.float64,
                   types.int32, types.int64,
                   types.uint32, types.uint64]
        x_values = [0.0, .21, .34, 1005382.042, -25.328]
        y1_values = [x + 2 for x in x_values]
        y2_values = [x - 2 for x in x_values]

        self.run_binary(pyfunc, x_types, x_values, y1_values)
        self.run_binary(pyfunc, x_types, x_values, y2_values)

        # Test using pos/neg inf
        self.run_binary(pyfunc, x_types, [0.0, -.5, .5], [math.inf]*3)
        self.run_binary(pyfunc, x_types, [0.0, -.5, .5], [-math.inf]*3)

        # if both args to nextafter are equal, then it is returned unchanged.
        self.run_binary(pyfunc, x_types, x_values, x_values)

    def test_degrees(self):
        pyfunc = degrees
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 1, 1, 1, 1, 1, 1., 1.]
        self.run_unary(pyfunc, x_types, x_values)

    def test_radians(self):
        pyfunc = radians
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [1, 1, 1, 1, 1, 1, 1., 1.]
        self.run_unary(pyfunc, x_types, x_values)

    def test_erf(self):
        pyfunc = erf
        x_values = [1., 1., -1., -0.0, 0.0, 0.5, 5, float('inf')]
        x_types = [types.float32, types.float64] * (len(x_values) // 2)
        self.run_unary(pyfunc, x_types, x_values, prec='double', ulps=2)

    def test_erfc(self):
        pyfunc = erfc
        x_values = [1., 1., -1., -0.0, 0.0, 0.5, 5, float('inf')]
        x_types = [types.float32, types.float64] * (len(x_values) // 2)
        self.run_unary(pyfunc, x_types, x_values, prec='double', ulps=4)

    def test_gamma(self):
        pyfunc = gamma
        x_values = [1., -0.9, -0.5, 0.5]
        x_types = [types.float32, types.float64] * (len(x_values) // 2)
        self.run_unary(pyfunc, x_types, x_values, prec='double', ulps=3)
        x_values = [-0.1, 0.1, 2.5, 10.1, 50., float('inf')]
        x_types = [types.float64] * len(x_values)
        self.run_unary(pyfunc, x_types, x_values, prec='double', ulps=8)

    def test_lgamma(self):
        pyfunc = lgamma
        x_values = [1., -0.9, -0.1, 0.1, 200., 1e10, 1e30, float('inf')]
        x_types = [types.float32, types.float64] * (len(x_values) // 2)
        self.run_unary(pyfunc, x_types, x_values, prec='double')

    def test_pow(self):
        pyfunc = pow
        x_types = [types.int16, types.int32, types.int64,
                   types.uint16, types.uint32, types.uint64,
                   types.float32, types.float64]
        x_values = [-2, -1, -2, 2, 1, 2, .1, .2]
        y_values = [x * 2 for x in x_values]
        self.run_binary(pyfunc, x_types, x_values, y_values)

    def test_gcd(self):
        from itertools import product, repeat, chain
        pyfunc = gcd
        signed_args = product(
            sorted(types.signed_domain), *repeat((-2, -1, 0, 1, 2, 7, 10), 2)
        )
        unsigned_args = product(
            sorted(types.unsigned_domain), *repeat((0, 1, 2, 7, 9, 16), 2)
        )
        x_types, x_values, y_values = zip(*chain(signed_args, unsigned_args))
        self.run_binary(pyfunc, x_types, x_values, y_values)

    def test_copysign(self):
        pyfunc = copysign
        value_types = [types.float32, types.float64]
        values = [-2, -1, -0.0, 0.0, 1, 2, float('-inf'), float('inf'),
                  float('nan')]
        x_types, x_values, y_values = list(zip(
            *itertools.product(value_types, values, values)))
        self.run_binary(pyfunc, x_types, x_values, y_values)

    def test_frexp(self):
        pyfunc = frexp
        x_types = [types.float32, types.float64]
        x_values = [-2.5, -0.0, 0.0, 3.5,
                    float('-inf'), float('inf'), float('nan')]
        self.run_unary(pyfunc, x_types, x_values, prec='exact')

    def test_ldexp(self):
        pyfunc = ldexp
        cfunc = njit(pyfunc)
        for fltty in (types.float32, types.float64):
            for args in [(2.5, -2), (2.5, 1), (0.0, 0), (0.0, 1),
                         (-0.0, 0), (-0.0, 1),
                         (float('inf'), 0), (float('-inf'), 0),
                         (float('nan'), 0)]:
                msg = 'for input %r' % (args,)
                self.assertPreciseEqual(cfunc(*args), pyfunc(*args))


if __name__ == '__main__':
    unittest.main()
