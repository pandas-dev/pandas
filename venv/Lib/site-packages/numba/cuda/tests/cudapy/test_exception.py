import numpy as np

from numba import cuda
from numba.cuda.testing import unittest, xfail_unless_cudasim, CUDATestCase
from numba.core import config


class TestException(CUDATestCase):
    def setUp(self):
        super().setUp()
        # LTO optimizes away the exception status due to an oversight
        # in the way we generate it (it is not added to the used list).
        self.skip_if_lto("Exceptions not supported with LTO")

    def test_exception(self):
        def foo(ary):
            x = cuda.threadIdx.x
            if x == 2:
                # NOTE: indexing with a out-of-bounds constant can fail at
                # compile-time instead (because the getitem is rewritten as a
                # static_getitem)
                ary.shape[-x]

        unsafe_foo = cuda.jit(foo)
        safe_foo = cuda.jit(debug=True, opt=False)(foo)

        if not config.ENABLE_CUDASIM:
            # Simulator throws exceptions regardless of debug
            # setting
            unsafe_foo[1, 3](np.array([0, 1]))

        with self.assertRaises(IndexError) as cm:
            safe_foo[1, 3](np.array([0, 1]))
        self.assertIn("tuple index out of range", str(cm.exception))

    def test_user_raise(self):
        @cuda.jit(debug=True, opt=False)
        def foo(do_raise):
            if do_raise:
                raise ValueError

        foo[1, 1](False)
        with self.assertRaises(ValueError):
            foo[1, 1](True)

    def case_raise_causing_warp_diverge(self, with_debug_mode):
        """Testing issue #2655.

        Exception raising code can cause the compiler to miss location
        of unifying branch target and resulting in unexpected warp
        divergence.
        """
        with_opt_mode = not with_debug_mode

        @cuda.jit(debug=with_debug_mode, opt=with_opt_mode)
        def problematic(x, y):
            tid = cuda.threadIdx.x
            ntid = cuda.blockDim.x

            if tid > 12:
                for i in range(ntid):
                    y[i] += x[i] // y[i]

            cuda.syncthreads()
            if tid < 17:
                for i in range(ntid):
                    x[i] += x[i] // y[i]

        @cuda.jit
        def oracle(x, y):
            tid = cuda.threadIdx.x
            ntid = cuda.blockDim.x

            if tid > 12:
                for i in range(ntid):
                    if y[i] != 0:
                        y[i] += x[i] // y[i]

            cuda.syncthreads()
            if tid < 17:
                for i in range(ntid):
                    if y[i] != 0:
                        x[i] += x[i] // y[i]

        n = 32
        got_x = 1. / (np.arange(n) + 0.01)
        got_y = 1. / (np.arange(n) + 0.01)
        problematic[1, n](got_x, got_y)

        expect_x = 1. / (np.arange(n) + 0.01)
        expect_y = 1. / (np.arange(n) + 0.01)
        oracle[1, n](expect_x, expect_y)

        np.testing.assert_almost_equal(expect_x, got_x)
        np.testing.assert_almost_equal(expect_y, got_y)

    def test_raise_causing_warp_diverge(self):
        """Test case for issue #2655.
        """
        self.case_raise_causing_warp_diverge(with_debug_mode=False)

    # The following two cases relate to Issue #7806: Division by zero stops the
    # kernel. https://github.com/numba/numba/issues/7806.

    def test_no_zero_division_error(self):
        # When debug is False:
        # - Division by zero raises no exception
        # - Execution proceeds after a divide by zero
        @cuda.jit
        def f(r, x, y):
            r[0] = y[0] / x[0]
            r[1] = y[0]

        r = np.zeros(2)
        x = np.zeros(1)
        y = np.ones(1)

        f[1, 1](r, x, y)

        self.assertTrue(np.isinf(r[0]), 'Expected inf from div by zero')
        self.assertEqual(r[1], y[0], 'Expected execution to continue')

    def test_zero_division_error_in_debug(self):
        # When debug is True:
        # - Zero by division raises an exception
        # - Execution halts at the point of division by zero
        @cuda.jit(debug=True, opt=False)
        def f(r, x, y):
            r[0] = y[0] / x[0]
            r[1] = y[0]

        r = np.zeros(2)
        x = np.zeros(1)
        y = np.ones(1)

        # Simulator and device behaviour differs slightly in the exception
        # raised - in debug mode, the CUDA target uses the Python error model,
        # which gives a ZeroDivision error. The simulator uses NumPy with the
        # error mode for division by zero set to raise, which results in a
        # FloatingPointError instead.
        if config.ENABLE_CUDASIM:
            exc = FloatingPointError
        else:
            exc = ZeroDivisionError

        with self.assertRaises(exc):
            f[1, 1](r, x, y)

        self.assertEqual(r[0], 0, 'Expected result to be left unset')
        self.assertEqual(r[1], 0, 'Expected execution to stop')

    @xfail_unless_cudasim
    def test_raise_in_device_function(self):
        # This is an expected failure because reporting of exceptions raised in
        # device functions does not work correctly - see Issue #8036:
        # https://github.com/numba/numba/issues/8036
        msg = 'Device Function Error'

        @cuda.jit(device=True)
        def f():
            raise ValueError(msg)

        @cuda.jit(debug=True)
        def kernel():
            f()

        with self.assertRaises(ValueError) as raises:
            kernel[1, 1]()

        self.assertIn(msg, str(raises.exception))


if __name__ == '__main__':
    unittest.main()
