import sys

import numpy as np

import unittest
from numba import jit, njit
from numba.core import types, errors, utils
from numba.tests.support import (captured_stdout, TestCase, EnableNRTStatsMixin)


def print_value(x):
    print(x)

def print_array_item(arr, i):
    print(arr[i].x)

def print_values(a, b, c):
    print(a, b, c)

def print_empty():
    print()

def print_string(x):
    print(x, "hop!", 3.5)

def print_vararg(a, b, c):
    print(a, b, *c)

def print_string_vararg(a, b, c):
    print(a, "hop!", b, *c)

def make_print_closure(x):
    def print_closure():
        return x
    return jit(nopython=True)(x)


class TestPrint(EnableNRTStatsMixin, TestCase):

    def check_values(self, typ, values):
        cfunc = njit((typ,))(print_value)
        for val in values:
            with captured_stdout():
                cfunc(val)
                self.assertEqual(sys.stdout.getvalue(), str(val) + '\n')

    def test_print_values(self):
        """
        Test printing a single argument value.
        """
        # Various scalars
        self.check_values(types.int32, (1, -234))
        self.check_values(types.int64, (1, -234,
                                        123456789876543210,
                                        -123456789876543210))
        self.check_values(types.uint64, (1, 234,
                                         123456789876543210, 2**63 + 123))
        self.check_values(types.boolean, (True, False))
        self.check_values(types.float64, (1.5, 100.0**10.0, float('nan')))
        self.check_values(types.complex64, (1+1j,))
        self.check_values(types.NPTimedelta('ms'), (np.timedelta64(100, 'ms'),))

        cfunc = njit((types.float32,))(print_value)
        with captured_stdout():
            cfunc(1.1)
            # Float32 will lose precision
            got = sys.stdout.getvalue()
            expect = '1.10000002384'
            self.assertTrue(got.startswith(expect))
            self.assertTrue(got.endswith('\n'))

        # Test array
        arraytype = types.Array(types.int32, 1, 'C')
        cfunc = njit((arraytype,))(print_value)
        with captured_stdout():
            cfunc(np.arange(10, dtype=np.int32))
            self.assertEqual(sys.stdout.getvalue(),
                             '[0 1 2 3 4 5 6 7 8 9]\n')

    @unittest.skip("Issue with intermittent NRT leak, see #9355.")
    def test_print_nrt_type(self):
        # NOTE: this check is extracted from the above as it started
        # intermittently leaking since the merge of #9330 (compile_isolated
        # removal patch). It's not clear why this happens, see #9355 for
        # thoughts/details. This test is skipped until it is resolved.

        # NRT-enabled type
        with self.assertNoNRTLeak():
            x = [1, 3, 5, 7]
            with self.assertRefCount(x):
                self.check_values(types.List(types.intp, reflected=True), (x,))

    def test_print_array_item(self):
        """
        Test printing a Numpy character sequence
        """
        dtype = np.dtype([('x', 'S4')])
        arr = np.frombuffer(bytearray(range(1, 9)), dtype=dtype)

        pyfunc = print_array_item
        cfunc = jit(nopython=True)(pyfunc)
        for i in range(len(arr)):
            with captured_stdout():
                cfunc(arr, i)
                self.assertEqual(sys.stdout.getvalue(), str(arr[i]['x']) + '\n')

    def test_print_multiple_values(self):
        pyfunc = print_values
        cfunc = njit((types.intp,) * 3)(pyfunc)
        with captured_stdout():
            cfunc(1, 2, 3)
            self.assertEqual(sys.stdout.getvalue(), '1 2 3\n')

    def test_print_nogil(self):
        pyfunc = print_values
        cfunc = jit(nopython=True, nogil=True)(pyfunc)
        with captured_stdout():
            cfunc(1, 2, 3)
            self.assertEqual(sys.stdout.getvalue(), '1 2 3\n')

    def test_print_empty(self):
        pyfunc = print_empty
        cfunc = njit((),)(pyfunc)
        with captured_stdout():
            cfunc()
            self.assertEqual(sys.stdout.getvalue(), '\n')

    def test_print_strings(self):
        pyfunc = print_string
        cfunc = njit((types.intp,))(pyfunc)
        with captured_stdout():
            cfunc(1)
            self.assertEqual(sys.stdout.getvalue(), '1 hop! 3.5\n')

    def test_print_vararg(self):
        # Test *args support for print().  This is desired since
        # print() can use a dedicated IR node.
        pyfunc = print_vararg
        cfunc = jit(nopython=True)(pyfunc)
        with captured_stdout():
            cfunc(1, (2, 3), (4, 5j))
            self.assertEqual(sys.stdout.getvalue(), '1 (2, 3) 4 5j\n')

        pyfunc = print_string_vararg
        cfunc = jit(nopython=True)(pyfunc)
        with captured_stdout():
            cfunc(1, (2, 3), (4, 5j))
            self.assertEqual(sys.stdout.getvalue(), '1 hop! (2, 3) 4 5j\n')

    def test_inner_fn_print(self):
        @jit(nopython=True)
        def foo(x):
            print(x)

        @jit(nopython=True)
        def bar(x):
            foo(x)
            foo('hello')

        # Printing an array requires the Env.
        # We need to make sure the inner function can obtain the Env.
        x = np.arange(5)
        with captured_stdout():
            bar(x)
            self.assertEqual(sys.stdout.getvalue(), '[0 1 2 3 4]\nhello\n')

    def test_print_w_kwarg_raises(self):
        @jit(nopython=True)
        def print_kwarg():
            print('x', flush=True)

        with self.assertRaises(errors.UnsupportedError) as raises:
            print_kwarg()
        expected = ("Numba's print() function implementation does not support "
                    "keyword arguments.")
        self.assertIn(raises.exception.msg, expected)

    def test_print_no_truncation(self):
        ''' See: https://github.com/numba/numba/issues/3811
        '''
        @jit(nopython=True)
        def foo():
            print(''.join(['a'] * 10000))
        with captured_stdout():
            foo()
            self.assertEqual(sys.stdout.getvalue(), ''.join(['a'] * 10000) + '\n')

if __name__ == '__main__':
    unittest.main()
