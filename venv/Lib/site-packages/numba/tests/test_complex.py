import itertools
import math
import sys

from numba import jit, types
from numba.tests.support import TestCase, skip_if_py314
from .complex_usecases import *
import unittest

enable_pyobj_flags = {'forceobj': True}
no_pyobj_flags = {'nopython': True}


class BaseComplexTest(object):

    def basic_values(self):
        reals = [-0.0, +0.0, 1, -1, +1.5, -3.5,
                 float('-inf'), float('+inf')]
        if sys.platform != 'win32':
            reals += [float('nan')]
        return [complex(x, y) for x, y in itertools.product(reals, reals)]

    def more_values(self):
        reals = [-0.0, +0.0, 1, -1, -math.pi, +math.pi,
                 float('-inf'), float('+inf')]
        if sys.platform != 'win32':
            reals += [float('nan')]
        return [complex(x, y) for x, y in itertools.product(reals, reals)]

    def non_nan_values(self):
        reals = [-0.0, +0.0, 1, -1, -math.pi, +math.pi,
                 float('inf'), float('-inf')]
        return [complex(x, y) for x, y in itertools.product(reals, reals)]

    def run_unary(self, pyfunc, x_types, x_values, ulps=1, abs_tol=None,
                  flags=enable_pyobj_flags):
        for tx in x_types:
            cfunc = jit((tx,), **flags)(pyfunc)
            prec = 'single' if tx in (types.float32, types.complex64) else 'double'
            for vx in x_values:
                try:
                    expected = pyfunc(vx)
                except ValueError as e:
                    self.assertIn("math domain error", str(e))
                    continue
                got = cfunc(vx)
                msg = 'for input %r with prec %r' % (vx, prec)
                self.assertPreciseEqual(got, expected, prec=prec,
                                        ulps=ulps, abs_tol=abs_tol, msg=msg)

    def run_binary(self, pyfunc, value_types, values, ulps=1,
                   flags=enable_pyobj_flags):
        for tx, ty in value_types:
            cfunc = jit((tx, ty), **flags)(pyfunc)
            prec = ('single'
                    if set([tx, ty]) & set([types.float32, types.complex64])
                    else 'double')
            for vx, vy in values:
                try:
                    expected = pyfunc(vx, vy)
                except ValueError as e:
                    self.assertIn("math domain error", str(e))
                    continue
                except ZeroDivisionError:
                    continue
                got = cfunc(vx, vy)
                msg = 'for input %r with prec %r' % ((vx, vy), prec)
                self.assertPreciseEqual(got, expected, prec=prec,
                                        ulps=ulps, msg=msg)


class TestComplex(BaseComplexTest, TestCase):

    def test_real(self, flags=enable_pyobj_flags):
        self.run_unary(real_usecase, [types.complex64, types.complex128],
                       self.basic_values(), flags=flags)
        self.run_unary(real_usecase, [types.int8, types.int64],
                       [1, 0, -3], flags=flags)
        self.run_unary(real_usecase, [types.float32, types.float64],
                       [1.5, -0.5], flags=flags)

    def test_real_npm(self):
        self.test_real(flags=no_pyobj_flags)

    def test_imag(self, flags=enable_pyobj_flags):
        self.run_unary(imag_usecase, [types.complex64, types.complex128],
                       self.basic_values(), flags=flags)
        self.run_unary(imag_usecase, [types.int8, types.int64],
                       [1, 0, -3], flags=flags)
        self.run_unary(imag_usecase, [types.float32, types.float64],
                       [1.5, -0.5], flags=flags)

    def test_imag_npm(self):
        self.test_imag(flags=no_pyobj_flags)

    def test_conjugate(self, flags=enable_pyobj_flags):
        self.run_unary(conjugate_usecase, [types.complex64, types.complex128],
                       self.basic_values(), flags=flags)
        self.run_unary(conjugate_usecase, [types.int8, types.int64],
                       [1, 0, -3], flags=flags)
        self.run_unary(conjugate_usecase, [types.float32, types.float64],
                       [1.5, -0.5], flags=flags)

    def test_conjugate_npm(self):
        self.test_conjugate(flags=no_pyobj_flags)

    def test_div(self, flags=enable_pyobj_flags):
        """
        Test complex.__div__ implementation with non-trivial values.
        """
        # XXX Fold into test_operator?
        values = list(itertools.product(self.more_values(), self.more_values()))
        value_types = [(types.complex128, types.complex128),
                       (types.complex64, types.complex64)]
        self.run_binary(div_usecase, value_types, values, flags=flags)

    @skip_if_py314
    def test_div_npm(self):
        self.test_div(flags=no_pyobj_flags)


class TestCMath(BaseComplexTest, TestCase):
    """
    Tests for cmath module support.
    """

    def check_predicate_func(self, pyfunc, flags):
        self.run_unary(pyfunc, [types.complex128, types.complex64],
                       self.basic_values(), flags=flags)

    def check_unary_func(self, pyfunc, flags, ulps=1, abs_tol=None,
                         values=None):
        self.run_unary(pyfunc, [types.complex128],
                       values or self.more_values(), flags=flags, ulps=ulps,
                       abs_tol=abs_tol)
        # Avoid discontinuities around pi when in single precision.
        self.run_unary(pyfunc, [types.complex64],
                       values or self.basic_values(), flags=flags, ulps=ulps,
                       abs_tol=abs_tol)

    # Conversions

    def test_phase(self):
        self.check_unary_func(phase_usecase, enable_pyobj_flags)

    def test_phase_npm(self):
        self.check_unary_func(phase_usecase, no_pyobj_flags)

    def test_polar(self):
        self.check_unary_func(polar_usecase, enable_pyobj_flags)

    def test_polar_npm(self):
        self.check_unary_func(polar_usecase, no_pyobj_flags)

    def test_rect(self, flags=enable_pyobj_flags):
        def do_test(tp, seed_values):
            values = [(z.real, z.imag) for z in seed_values
                      if not math.isinf(z.imag) or z.real == 0]
            self.run_binary(rect_usecase, [(tp, tp)], values, flags=flags)
        do_test(types.float64, self.more_values())
        # Avoid discontinuities around pi when in single precision.
        do_test(types.float32, self.basic_values())

    def test_rect_npm(self):
        self.test_rect(flags=no_pyobj_flags)

    # Classification

    def test_isnan(self, flags=enable_pyobj_flags):
        self.check_predicate_func(isnan_usecase, enable_pyobj_flags)

    def test_isnan_npm(self):
        self.check_predicate_func(isnan_usecase, no_pyobj_flags)

    def test_isinf(self, flags=enable_pyobj_flags):
        self.check_predicate_func(isinf_usecase, enable_pyobj_flags)

    def test_isinf_npm(self):
        self.check_predicate_func(isinf_usecase, no_pyobj_flags)

    def test_isfinite(self, flags=enable_pyobj_flags):
        self.check_predicate_func(isfinite_usecase, enable_pyobj_flags)

    def test_isfinite_npm(self):
        self.check_predicate_func(isfinite_usecase, no_pyobj_flags)

    # Power and logarithms

    def test_exp(self):
        self.check_unary_func(exp_usecase, enable_pyobj_flags, ulps=2)

    def test_exp_npm(self):
        # Aggressive optimization fixes the following subnormal float problem.
        ## The two tests are failing due to subnormal float problems.
        ## We are seeing (6.9532198665326e-310+2.1221202807e-314j) != 0j
        self.check_unary_func(exp_usecase, no_pyobj_flags, ulps=2)

    def test_log(self):
        self.check_unary_func(log_usecase, enable_pyobj_flags)

    def test_log_npm(self):
        self.check_unary_func(log_usecase, no_pyobj_flags)

    def test_log_base(self, flags=enable_pyobj_flags):
        values = list(itertools.product(self.more_values(), self.more_values()))
        value_types = [(types.complex128, types.complex128),
                       (types.complex64, types.complex64)]
        self.run_binary(log_base_usecase, value_types, values, flags=flags,
                        ulps=3)

    @skip_if_py314
    def test_log_base_npm(self):
        self.test_log_base(flags=no_pyobj_flags)

    def test_log10(self):
        self.check_unary_func(log10_usecase, enable_pyobj_flags)

    def test_log10_npm(self):
        self.check_unary_func(log10_usecase, no_pyobj_flags)

    def test_sqrt(self):
        self.check_unary_func(sqrt_usecase, enable_pyobj_flags)

    def test_sqrt_npm(self):
        self.check_unary_func(sqrt_usecase, no_pyobj_flags)
        # issue #3499, check large complex128 values, values cross the float32
        # representation limit threshold
        values = [-10 ** i for i in range(36, 41)]
        self.run_unary(sqrt_usecase, [types.complex128],
                       values, flags=no_pyobj_flags)

    # Trigonometric functions

    def test_acos(self):
        self.check_unary_func(acos_usecase, enable_pyobj_flags, ulps=2)

    def test_acos_npm(self):
        self.check_unary_func(acos_usecase, no_pyobj_flags, ulps=2)

    def test_asin(self):
        self.check_unary_func(asin_usecase, enable_pyobj_flags, ulps=2)

    def test_asin_npm(self):
        self.check_unary_func(asin_usecase, no_pyobj_flags, ulps=2)

    def test_atan(self):
        self.check_unary_func(atan_usecase, enable_pyobj_flags, ulps=2,)

    def test_atan_npm(self):
        self.check_unary_func(atan_usecase, no_pyobj_flags, ulps=2,)

    def test_cos(self):
        self.check_unary_func(cos_usecase, enable_pyobj_flags, ulps=2)

    def test_cos_npm(self):
        self.check_unary_func(cos_usecase, no_pyobj_flags, ulps=2)

    def test_sin(self):
        # See test_sinh.
        self.check_unary_func(sin_usecase, enable_pyobj_flags, abs_tol='eps')

    def test_sin_npm(self):
        self.check_unary_func(sin_usecase, no_pyobj_flags, abs_tol='eps')

    def test_tan(self):
        self.check_unary_func(tan_usecase, enable_pyobj_flags, ulps=2)

    def test_tan_npm(self):
        self.check_unary_func(tan_usecase, enable_pyobj_flags, ulps=2)

    # Hyperbolic functions

    def test_acosh(self):
        self.check_unary_func(acosh_usecase, enable_pyobj_flags)

    def test_acosh_npm(self):
        self.check_unary_func(acosh_usecase, no_pyobj_flags)

    def test_asinh(self):
        self.check_unary_func(asinh_usecase, enable_pyobj_flags, ulps=2)

    def test_asinh_npm(self):
        self.check_unary_func(asinh_usecase, no_pyobj_flags, ulps=2)

    def test_atanh(self):
        self.check_unary_func(atanh_usecase, enable_pyobj_flags, ulps=2)

    def test_atanh_npm(self):
        self.check_unary_func(atanh_usecase, no_pyobj_flags, ulps=2)

    def test_cosh(self):
        self.check_unary_func(cosh_usecase, enable_pyobj_flags, ulps=2)

    def test_cosh_npm(self):
        self.check_unary_func(cosh_usecase, no_pyobj_flags, ulps=2)

    def test_sinh(self):
        self.check_unary_func(sinh_usecase, enable_pyobj_flags, abs_tol='eps')

    def test_sinh_npm(self):
        self.check_unary_func(sinh_usecase, no_pyobj_flags, abs_tol='eps')

    def test_tanh(self):
        self.check_unary_func(tanh_usecase, enable_pyobj_flags, ulps=2)

    def test_tanh_npm(self):
        self.check_unary_func(tanh_usecase, enable_pyobj_flags, ulps=2)


if __name__ == '__main__':
    unittest.main()
