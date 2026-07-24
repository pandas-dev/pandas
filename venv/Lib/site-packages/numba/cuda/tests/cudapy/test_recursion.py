from numba import cuda
from numba.core.errors import TypingError
from numba.cuda.testing import CUDATestCase, skip_on_cudasim
import numpy as np
import unittest


class TestSelfRecursion(CUDATestCase):

    def setUp(self):
        # Avoid importing this module at the top level, as it triggers
        # compilation and can therefore fail
        from numba.cuda.tests.cudapy import recursion_usecases
        self.mod = recursion_usecases
        super().setUp()

    def check_fib(self, cfunc):
        @cuda.jit
        def kernel(r, x):
            r[0] = cfunc(x[0])

        x = np.asarray([10], dtype=np.int64)
        r = np.zeros_like(x)
        kernel[1, 1](r, x)

        actual = r[0]
        expected = 55
        self.assertPreciseEqual(actual, expected)

    def test_global_explicit_sig(self):
        self.check_fib(self.mod.fib1)

    def test_inner_explicit_sig(self):
        self.check_fib(self.mod.fib2)

    def test_global_implicit_sig(self):
        self.check_fib(self.mod.fib3)

    @skip_on_cudasim('Simulator does not compile')
    def test_runaway(self):
        with self.assertRaises(TypingError) as raises:
            cfunc = self.mod.runaway_self

            @cuda.jit('void()')
            def kernel():
                cfunc(1)

        self.assertIn("cannot type infer runaway recursion",
                      str(raises.exception))

    @unittest.skip('Needs insert_unresolved_ref support in target')
    def test_type_change(self):
        pfunc = self.mod.type_change_self.py_func
        cfunc = self.mod.type_change_self

        @cuda.jit
        def kernel(r, x, y):
            r[0] = cfunc(x[0], y[0])

        args = 13, 0.125
        x = np.asarray([args[0]], dtype=np.int64)
        y = np.asarray([args[1]], dtype=np.float64)
        r = np.zeros_like(x)

        kernel[1, 1](r, x, y)

        expected = pfunc(*args)
        actual = r[0]

        self.assertPreciseEqual(actual, expected)

    @unittest.expectedFailure
    def test_raise(self):
        # This is an expected failure because reporting of exceptions raised in
        # device functions does not work correctly - see Issue #8036:
        # https://github.com/numba/numba/issues/8036
        with self.assertRaises(ValueError) as raises:
            self.mod.raise_self_kernel[1, 1](3)

        self.assertEqual(str(raises.exception), "raise_self")

    @unittest.skip('Needs insert_unresolved_ref support in target')
    def test_optional_return(self):
        pfunc = self.mod.make_optional_return_case()
        cfunc = self.mod.make_optional_return_case(cuda.jit)

        @cuda.jit
        def kernel(r, x):
            res = cfunc(x[0])
            if res is None:
                res = 999
            r[0] = res

        def cpu_kernel(x):
            res = pfunc(x)
            if res is None:
                res = 999
            return res

        for arg in (0, 5, 10, 15):
            expected = cpu_kernel(arg)
            x = np.asarray([arg], dtype=np.int64)
            r = np.zeros_like(x)
            kernel[1, 1](r, x)
            actual = r[0]

            self.assertEqual(expected, actual)

    @skip_on_cudasim('Recursion handled because simulator does not compile')
    def test_growing_return_tuple(self):
        cfunc = self.mod.make_growing_tuple_case(cuda.jit)

        with self.assertRaises(TypingError) as raises:
            @cuda.jit('void()')
            def kernel():
                cfunc(100)

        self.assertIn(
            "Return type of recursive function does not converge",
            str(raises.exception),
        )


if __name__ == '__main__':
    unittest.main()
