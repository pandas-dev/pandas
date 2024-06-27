from ctypes import *
import sys
import threading

import numpy as np


from numba import jit, njit
from numba.core import types, errors
from numba.core.typing import ctypes_utils
from numba.tests.support import MemoryLeakMixin, tag, TestCase
from numba.tests.ctypes_usecases import *
import unittest


class TestCTypesTypes(TestCase):

    def _conversion_tests(self, check):
        check(c_double, types.float64)
        check(c_int, types.intc)
        check(c_uint16, types.uint16)
        check(c_size_t, types.size_t)
        check(c_ssize_t, types.ssize_t)

        check(c_void_p, types.voidptr)
        check(POINTER(c_float), types.CPointer(types.float32))
        check(POINTER(POINTER(c_float)),
              types.CPointer(types.CPointer(types.float32)))

        check(None, types.void)

    def test_from_ctypes(self):
        """
        Test converting a ctypes type to a Numba type.
        """
        def check(cty, ty):
            got = ctypes_utils.from_ctypes(cty)
            self.assertEqual(got, ty)

        self._conversion_tests(check)

        # An unsupported type
        with self.assertRaises(TypeError) as raises:
            ctypes_utils.from_ctypes(c_wchar_p)
        self.assertIn("Unsupported ctypes type", str(raises.exception))

    def test_to_ctypes(self):
        """
        Test converting a Numba type to a ctypes type.
        """
        def check(cty, ty):
            got = ctypes_utils.to_ctypes(ty)
            self.assertEqual(got, cty)

        self._conversion_tests(check)

        # An unsupported type
        with self.assertRaises(TypeError) as raises:
            ctypes_utils.to_ctypes(types.ellipsis)
        self.assertIn("Cannot convert Numba type '...' to ctypes type",
                      str(raises.exception))


class TestCTypesUseCases(MemoryLeakMixin, TestCase):

    def test_c_sin(self):
        pyfunc = use_c_sin
        cfunc = njit((types.double,))(pyfunc)
        x = 3.14
        self.assertEqual(pyfunc(x), cfunc(x))

    def test_two_funcs(self):
        # Check that two constant functions don't get mixed up.
        pyfunc = use_two_funcs
        cfunc = njit((types.double,))(pyfunc)
        x = 3.14
        self.assertEqual(pyfunc(x), cfunc(x))

    @unittest.skipUnless(is_windows, "Windows-specific test")
    def test_stdcall(self):
        # Just check that it doesn't crash
        cfunc = njit((types.uintc,))(use_c_sleep)

        cfunc(1)

    def test_ctype_wrapping(self):
        pyfunc = use_ctype_wrapping
        cfunc = njit((types.double,))(pyfunc)
        x = 3.14
        self.assertEqual(pyfunc(x), cfunc(x))

    def test_ctype_voidptr(self):
        pyfunc = use_c_pointer
        # pyfunc will segfault if called
        cfunc = njit((types.int32,))(pyfunc)
        x = 123
        self.assertEqual(cfunc(x), x + 1)

    def test_function_pointer(self):
        pyfunc = use_func_pointer
        cfunc = jit(nopython=True)(pyfunc)
        for (fa, fb, x) in [
            (c_sin, c_cos, 1.0),
            (c_sin, c_cos, -1.0),
            (c_cos, c_sin, 1.0),
            (c_cos, c_sin, -1.0)]:
            expected = pyfunc(fa, fb, x)
            got = cfunc(fa, fb, x)
            self.assertEqual(got, expected)
        # A single specialization was compiled for all calls
        self.assertEqual(len(cfunc.overloads), 1, cfunc.overloads)

    def test_untyped_function(self):
        with self.assertRaises(TypeError) as raises:
            njit((types.double,))(use_c_untyped)
        self.assertIn("ctypes function '_numba_test_exp' doesn't define its argument types",
                      str(raises.exception))

    def test_python_call_back(self):
        mydct = {'what': 1232121}

        def call_me_maybe(arr):
            return mydct[arr[0].decode('ascii')]

        # Create a callback into the python interpreter
        py_call_back = CFUNCTYPE(c_int, py_object)(call_me_maybe)

        def pyfunc(a):
            what = py_call_back(a)
            return what

        cfunc = jit(nopython=True, nogil=True)(pyfunc)
        arr = np.array(["what"], dtype='S10')
        self.assertEqual(pyfunc(arr), cfunc(arr))

    def test_python_call_back_threaded(self):
        def pyfunc(a, repeat):
            out = 0
            for _ in range(repeat):
                out += py_call_back(a)
            return out

        cfunc = jit(nopython=True, nogil=True)(pyfunc)

        arr = np.array(["what"], dtype='S10')
        repeat = 1000

        expected = pyfunc(arr, repeat)
        outputs = []

        # Warm up
        cfunc(arr, repeat)

        # Test the function in multiple threads to exercise the
        # GIL ensure/release code

        def run(func, arr, repeat):
            outputs.append(func(arr, repeat))

        threads = [threading.Thread(target=run, args=(cfunc, arr, repeat))
                   for _ in range(10)]

        # Start threads
        for th in threads:
            th.start()

        # End threads
        for th in threads:
            th.join()

        # Check results
        for got in outputs:
            self.assertEqual(expected, got)

    def test_passing_array_ctypes_data(self):
        """
        Test the ".ctypes.data" attribute of an array can be passed
        as a "void *" parameter.
        """
        def pyfunc(arr):
            return c_take_array_ptr(arr.ctypes.data)

        cfunc = jit(nopython=True, nogil=True)(pyfunc)

        arr = np.arange(5)

        expected = pyfunc(arr)
        got = cfunc(arr)

        self.assertEqual(expected, got)

    def check_array_ctypes(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)

        arr = np.linspace(0, 10, 5)
        expected = arr ** 2.0
        got = cfunc(arr)
        self.assertPreciseEqual(expected, got)
        return cfunc

    def test_passing_array_ctypes_voidptr(self):
        """
        Test the ".ctypes" attribute of an array can be passed
        as a "void *" parameter.
        """
        self.check_array_ctypes(use_c_vsquare)

    def test_passing_array_ctypes_voidptr_pass_ptr(self):
        """
        Test the ".ctypes" attribute of an array can be passed
        as a pointer parameter of the right type.
        """
        cfunc = self.check_array_ctypes(use_c_vcube)

        # Non-compatible pointers are not accepted (here float32* vs. float64*)
        with self.assertRaises(errors.TypingError) as raises:
            cfunc(np.float32([0.0]))

        self.assertIn("No implementation of function ExternalFunctionPointer",
                      str(raises.exception))

    def test_storing_voidptr_to_int_array(self):
        # Make C callback that returns a void*
        cproto = CFUNCTYPE(c_void_p)

        @cproto
        def get_voidstar():
            return 0xdeadbeef

        # Make python functions that use the C callback
        def pyfunc(a):
            ptr = get_voidstar()
            a[0] = ptr
            return ptr

        # Compile it
        cfunc = njit((types.uintp[::1],))(pyfunc)

        # Setup inputs
        arr_got = np.zeros(1, dtype=np.uintp)
        arr_expect = arr_got.copy()

        # Run functions
        ret_got = cfunc(arr_got)
        ret_expect = pyfunc(arr_expect)

        # Check
        self.assertEqual(ret_expect, 0xdeadbeef)
        self.assertPreciseEqual(ret_got, ret_expect)
        self.assertPreciseEqual(arr_got, arr_expect)


if __name__ == '__main__':
    unittest.main()

