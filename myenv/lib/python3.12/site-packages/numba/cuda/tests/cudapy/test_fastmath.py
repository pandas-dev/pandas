from typing import List
from dataclasses import dataclass, field
from numba import cuda, float32
from numba.cuda.compiler import compile_ptx_for_current_device, compile_ptx
from math import cos, sin, tan, exp, log, log10, log2, pow, tanh
from operator import truediv
import numpy as np
from numba.cuda.testing import (CUDATestCase, skip_on_cudasim,
                                skip_unless_cc_75)
import unittest


@dataclass
class FastMathCriterion:
    fast_expected: List[str] = field(default_factory=list)
    fast_unexpected: List[str] = field(default_factory=list)
    prec_expected: List[str] = field(default_factory=list)
    prec_unexpected: List[str] = field(default_factory=list)

    def check(self, test: CUDATestCase, fast: str, prec: str):
        test.assertTrue(all(i in fast for i in self.fast_expected))
        test.assertTrue(all(i not in fast for i in self.fast_unexpected))
        test.assertTrue(all(i in prec for i in self.prec_expected))
        test.assertTrue(all(i not in prec for i in self.prec_unexpected))


@skip_on_cudasim('Fastmath and PTX inspection not available on cudasim')
class TestFastMathOption(CUDATestCase):
    def _test_fast_math_common(self, pyfunc, sig, device, criterion):

        # Test jit code path
        fastver = cuda.jit(sig, device=device, fastmath=True)(pyfunc)
        precver = cuda.jit(sig, device=device)(pyfunc)

        criterion.check(
            self, fastver.inspect_asm(sig), precver.inspect_asm(sig)
        )

        # Test compile_ptx code path
        fastptx, _ = compile_ptx_for_current_device(
            pyfunc, sig, device=device, fastmath=True
        )
        precptx, _ = compile_ptx_for_current_device(
            pyfunc, sig, device=device
        )

        criterion.check(self, fastptx, precptx)

    def _test_fast_math_unary(self, op, criterion: FastMathCriterion):
        def kernel(r, x):
            r[0] = op(x)

        def device_function(x):
            return op(x)

        self._test_fast_math_common(
            kernel, (float32[::1], float32), device=False, criterion=criterion
        )
        self._test_fast_math_common(
            device_function, (float32,), device=True, criterion=criterion
        )

    def _test_fast_math_binary(self, op, criterion: FastMathCriterion):
        def kernel(r, x, y):
            r[0] = op(x, y)

        def device(x, y):
            return op(x, y)

        self._test_fast_math_common(
            kernel,
            (float32[::1], float32, float32), device=False, criterion=criterion
        )
        self._test_fast_math_common(
            device, (float32, float32), device=True, criterion=criterion
        )

    def test_cosf(self):
        self._test_fast_math_unary(
            cos,
            FastMathCriterion(
                fast_expected=['cos.approx.ftz.f32 '],
                prec_unexpected=['cos.approx.ftz.f32 ']
            )
        )

    def test_sinf(self):
        self._test_fast_math_unary(
            sin,
            FastMathCriterion(
                fast_expected=['sin.approx.ftz.f32 '],
                prec_unexpected=['sin.approx.ftz.f32 ']
            )
        )

    def test_tanf(self):
        self._test_fast_math_unary(
            tan,
            FastMathCriterion(fast_expected=[
                'sin.approx.ftz.f32 ',
                'cos.approx.ftz.f32 ',
                'div.approx.ftz.f32 '
            ], prec_unexpected=['sin.approx.ftz.f32 '])
        )

    @skip_unless_cc_75
    def test_tanhf(self):

        self._test_fast_math_unary(
            tanh,
            FastMathCriterion(
                fast_expected=['tanh.approx.f32 '],
                prec_unexpected=['tanh.approx.f32 ']
            )
        )

    def test_tanhf_compile_ptx(self):
        def tanh_kernel(r, x):
            r[0] = tanh(x)

        def tanh_common_test(cc, criterion):
            fastptx, _ = compile_ptx(tanh_kernel, (float32[::1], float32),
                                     fastmath=True, cc=cc)
            precptx, _ = compile_ptx(tanh_kernel, (float32[::1], float32),
                                     cc=cc)
            criterion.check(self, fastptx, precptx)

        tanh_common_test(cc=(7, 5), criterion=FastMathCriterion(
            fast_expected=['tanh.approx.f32 '],
            prec_unexpected=['tanh.approx.f32 ']
        ))

        tanh_common_test(cc=(7, 0),
                         criterion=FastMathCriterion(
            fast_expected=['ex2.approx.ftz.f32 ',
                           'rcp.approx.ftz.f32 '],
            prec_unexpected=['tanh.approx.f32 ']))

    def test_expf(self):
        self._test_fast_math_unary(
            exp,
            FastMathCriterion(
                fast_unexpected=['fma.rn.f32 '],
                prec_expected=['fma.rn.f32 ']
            )
        )

    def test_logf(self):
        # Look for constant used to convert from log base 2 to log base e
        self._test_fast_math_unary(
            log, FastMathCriterion(
                fast_expected=['lg2.approx.ftz.f32 ', '0f3F317218'],
                prec_unexpected=['lg2.approx.ftz.f32 '],
            )
        )

    def test_log10f(self):
        # Look for constant used to convert from log base 2 to log base 10
        self._test_fast_math_unary(
            log10, FastMathCriterion(
                fast_expected=['lg2.approx.ftz.f32 ', '0f3E9A209B'],
                prec_unexpected=['lg2.approx.ftz.f32 ']
            )
        )

    def test_log2f(self):
        self._test_fast_math_unary(
            log2, FastMathCriterion(
                fast_expected=['lg2.approx.ftz.f32 '],
                prec_unexpected=['lg2.approx.ftz.f32 ']
            )
        )

    def test_powf(self):
        self._test_fast_math_binary(
            pow, FastMathCriterion(
                fast_expected=['lg2.approx.ftz.f32 '],
                prec_unexpected=['lg2.approx.ftz.f32 '],
            )
        )

    def test_divf(self):
        self._test_fast_math_binary(
            truediv, FastMathCriterion(
                fast_expected=['div.approx.ftz.f32 '],
                fast_unexpected=['div.rn.f32'],
                prec_expected=['div.rn.f32'],
                prec_unexpected=['div.approx.ftz.f32 '],
            )
        )

    def test_divf_exception(self):
        # LTO optimizes away the exception status due to an oversight
        # in the way we generate it (it is not added to the used list).
        self.skip_if_lto("Exceptions not supported with LTO")

        def f10(r, x, y):
            r[0] = x / y

        sig = (float32[::1], float32, float32)
        fastver = cuda.jit(sig, fastmath=True, debug=True)(f10)
        precver = cuda.jit(sig, debug=True)(f10)
        nelem = 10
        ary = np.empty(nelem, dtype=np.float32)
        with self.assertRaises(ZeroDivisionError):
            precver[1, nelem](ary, 10.0, 0.0)

        try:
            fastver[1, nelem](ary, 10.0, 0.0)
        except ZeroDivisionError:
            self.fail("Divide in fastmath should not throw ZeroDivisionError")

    @unittest.expectedFailure
    def test_device_fastmath_propagation(self):
        # The fastmath option doesn't presently propagate to device functions
        # from their callees - arguably it should do, so this test is presently
        # an xfail.
        @cuda.jit("float32(float32, float32)", device=True)
        def foo(a, b):
            return a / b

        def bar(arr, val):
            i = cuda.grid(1)
            if i < arr.size:
                arr[i] = foo(i, val)

        sig = (float32[::1], float32)
        fastver = cuda.jit(sig, fastmath=True)(bar)
        precver = cuda.jit(sig)(bar)

        # Variants of the div instruction are further documented at:
        # https://docs.nvidia.com/cuda/parallel-thread-execution/index.html#floating-point-instructions-div

        # The fast version should use the "fast, approximate divide" variant
        self.assertIn('div.approx.f32', fastver.inspect_asm(sig))
        # The precise version should use the "IEEE 754 compliant rounding"
        # variant, and neither of the "approximate divide" variants.
        self.assertIn('div.rn.f32', precver.inspect_asm(sig))
        self.assertNotIn('div.approx.f32', precver.inspect_asm(sig))
        self.assertNotIn('div.full.f32', precver.inspect_asm(sig))


if __name__ == '__main__':
    unittest.main()
