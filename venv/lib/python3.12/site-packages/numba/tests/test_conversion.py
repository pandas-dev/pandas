import gc
import itertools

import numpy as np

import unittest
from numba import jit, njit
from numba.core import types
from numba.tests.support import TestCase
from numba.np import numpy_support


def identity(x):
    return x

def addition(x, y):
    return x + y

def equality(x, y):
    return x == y

def foobar(x, y, z):
    return x


class TestConversion(TestCase):
    """
    Testing Python to Native conversion
    """

    def test_complex_identity(self):
        pyfunc = identity
        cfunc = njit(types.complex64(types.complex64))(pyfunc)

        xs = [1.0j, (1+1j), (-1-1j), (1+0j)]
        for x in xs:
            self.assertEqual(cfunc(x), x)
        for x in np.complex64(xs):
            self.assertEqual(cfunc(x), x)

        cfunc = njit(types.complex128(types.complex128))(pyfunc)

        xs = [1.0j, (1+1j), (-1-1j), (1+0j)]
        for x in xs:
            self.assertEqual(cfunc(x), x)
        for x in np.complex128(xs):
            self.assertEqual(cfunc(x), x)

    def test_complex_addition(self):
        pyfunc = addition
        cfunc = njit(types.complex64(types.complex64, types.complex64))(pyfunc)

        xs = [1.0j, (1+1j), (-1-1j), (1+0j)]
        for x in xs:
            y = x
            self.assertEqual(cfunc(x, y), x + y)
        for x in np.complex64(xs):
            y = x
            self.assertEqual(cfunc(x, y), x + y)


        cfunc = njit(types.complex128(types.complex128,
                                      types.complex128))(pyfunc)

        xs = [1.0j, (1+1j), (-1-1j), (1+0j)]
        for x in xs:
            y = x
            self.assertEqual(cfunc(x, y), x + y)
        for x in np.complex128(xs):
            y = x
            self.assertEqual(cfunc(x, y), x + y)

    def test_boolean_as_int(self):
        pyfunc = equality
        cfunc = njit((types.boolean, types.intp))(pyfunc)

        xs = True, False
        ys = -1, 0, 1

        for xs, ys in itertools.product(xs, ys):
            self.assertEqual(pyfunc(xs, ys), cfunc(xs, ys))

    def test_boolean_as_float(self):
        pyfunc = equality
        cfunc = njit((types.boolean, types.float64))(pyfunc)

        xs = True, False
        ys = -1, 0, 1

        for xs, ys in itertools.product(xs, ys):
            self.assertEqual(pyfunc(xs, ys), cfunc(xs, ys))

    def test_boolean_eq_boolean(self):
        pyfunc = equality
        cfunc = njit((types.boolean, types.boolean))(pyfunc)

        xs = True, False
        ys = True, False

        for xs, ys in itertools.product(xs, ys):
            self.assertEqual(pyfunc(xs, ys), cfunc(xs, ys))

    # test when a function parameters are jitted as unsigned types
    # the function is called with negative parameters the Python error
    # that it generates is correctly handled -- a Python error is returned to the user
    # For more info, see the comment in Include/longobject.h for _PyArray_AsByteArray
    # which PyLong_AsUnsignedLongLong calls
    def test_negative_to_unsigned(self):
        def f(x):
            return x
        with self.assertRaises(OverflowError):
            jit('uintp(uintp)', nopython=True)(f)(-5)

    # test the switch logic in callwraper.py:build_wrapper() works for more than one argument
    # and where the error occurs 
    def test_multiple_args_negative_to_unsigned(self): 
        pyfunc = foobar
        cfunc = njit(types.uint64(types.uint64, types.uint64,
                                  types.uint64),)(pyfunc)

        test_fail_args = ((-1, 0, 1), (0, -1, 1), (0, 1, -1))
        with self.assertRaises(OverflowError):
            for a, b, c in test_fail_args:
                cfunc(a, b, c)

    # test switch logic of callwraper.py:build_wrapper() with records as function parameters
    def test_multiple_args_records(self): 
        pyfunc = foobar

        mystruct_dt = np.dtype([('p', np.float64),
                           ('row', np.float64),
                           ('col', np.float64)])
        mystruct = numpy_support.from_dtype(mystruct_dt)

        cfunc = njit(mystruct[:](mystruct[:], types.uint64,
                                 types.uint64),)(pyfunc)

        st1 = np.recarray(3, dtype=mystruct_dt)

        st1.p = np.arange(st1.size) + 1
        st1.row = np.arange(st1.size) + 1
        st1.col = np.arange(st1.size) + 1

        with self.assertRefCount(st1):
            test_fail_args = ((st1, -1, 1), (st1, 1, -1))

            for a, b, c in test_fail_args:
                with self.assertRaises(OverflowError):
                    cfunc(a, b, c)

            del test_fail_args, a, b, c
            gc.collect()

    # test switch logic of callwraper.py:build_wrapper() with no function parameters
    def test_with_no_parameters(self):
        def f():
            pass 
        self.assertEqual(f(), jit('()', nopython=True)(f)())

    def check_argument_cleanup(self, typ, obj):
        """
        Check that argument cleanup doesn't leak references.
        """
        def f(x, y):
            pass

        def _objects(obj):
            objs = [obj]
            if isinstance(obj, tuple):
                for v in obj:
                    objs += _objects(v)
            return objs

        objects = _objects(obj)

        cfunc = njit((typ, types.uint32))(f)
        with self.assertRefCount(*objects):
            cfunc(obj, 1)
        with self.assertRefCount(*objects):
            with self.assertRaises(OverflowError):
                cfunc(obj, -1)

        cfunc = njit((types.uint32, typ))(f)
        with self.assertRefCount(*objects):
            cfunc(1, obj)
        with self.assertRefCount(*objects):
            with self.assertRaises(OverflowError):
                cfunc(-1, obj)

    def test_cleanup_buffer(self):
        mem = memoryview(bytearray(b"xyz"))
        self.check_argument_cleanup(types.MemoryView(types.byte, 1, 'C'), mem)

    def test_cleanup_record(self):
        dtype = np.dtype([('x', np.float64), ('y', np.float64)])
        recarr = np.zeros(1, dtype=dtype)
        self.check_argument_cleanup(numpy_support.from_dtype(dtype), recarr[0])

    def test_cleanup_tuple(self):
        mem = memoryview(bytearray(b"xyz"))
        tp = types.UniTuple(types.MemoryView(types.byte, 1, 'C'), 2)
        self.check_argument_cleanup(tp, (mem, mem))

    def test_cleanup_optional(self):
        mem = memoryview(bytearray(b"xyz"))
        tp = types.Optional(types.MemoryView(types.byte, 1, 'C'))
        self.check_argument_cleanup(tp, mem)

    def test_stringliteral_to_unicode(self):
        # See issue #6907, explicit signature on bar() takes a unicode_type but
        # the call to bar() in foo() is with a StringLiteral

        @jit(types.void(types.unicode_type), nopython=True)
        def bar(string):
            pass

        @jit(types.void(), nopython=True)
        def foo2():
            bar("literal string")


if __name__ == '__main__':
    unittest.main()
