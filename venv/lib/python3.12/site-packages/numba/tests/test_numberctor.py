import numpy as np

from numba import jit, njit
from numba.core import types

from numba.tests.support import TestCase, tag
import unittest


def dobool(a):
    return bool(a)


def doint(a):
    return int(a)


def dofloat(a):
    return float(a)


def docomplex(a):
    return complex(a)


def docomplex2(a, b):
    return complex(a, b)


def complex_calc(a):
    z = complex(a)
    return z.real ** 2 + z.imag ** 2


def complex_calc2(a, b):
    z = complex(a, b)
    return z.real ** 2 + z.imag ** 2


def converter(tp):
    def f(a):
        return tp(a)
    return f


def real_np_types():
    for tp_name in ('int8', 'int16', 'int32', 'int64',
                    'uint8', 'uint16', 'uint32', 'uint64',
                    'intc', 'uintc', 'intp', 'uintp',
                    'float32', 'float64', 'bool_'):
        yield tp_name

def complex_np_types():
    for tp_name in ('complex64', 'complex128'):
        yield tp_name


class TestScalarNumberCtor(TestCase):
    """
    Test <number class>(some scalar)
    """

    def check_int_constructor(self, pyfunc):
        x_types = [
            types.boolean, types.int32, types.int64, types.float32, types.float64
        ]
        x_values = [1, 0, 1000, 12.2, 23.4]

        for ty, x in zip(x_types, x_values):
            cfunc = njit((ty,))(pyfunc)
            self.assertPreciseEqual(pyfunc(x), cfunc(x))

    def test_bool(self):
        self.check_int_constructor(dobool)

    def test_int(self):
        self.check_int_constructor(doint)

    def test_float(self):
        pyfunc = dofloat

        x_types = [
            types.int32, types.int64, types.float32, types.float64
        ]
        x_values = [1, 1000, 12.2, 23.4]

        for ty, x in zip(x_types, x_values):
            cfunc = njit((ty,))(pyfunc)
            self.assertPreciseEqual(pyfunc(x), cfunc(x),
                prec='single' if ty is types.float32 else 'exact')

    def test_complex(self):
        pyfunc = docomplex

        x_types = [
            types.int32, types.int64, types.float32, types.float64,
            types.complex64, types.complex128,
        ]
        x_values = [1, 1000, 12.2, 23.4, 1.5-5j, 1-4.75j]

        for ty, x in zip(x_types, x_values):
            cfunc = njit((ty,))(pyfunc)
            got = cfunc(x)
            expected = pyfunc(x)
            self.assertPreciseEqual(pyfunc(x), cfunc(x),
                prec='single' if ty is types.float32 else 'exact')

        # Check that complex(float32) really creates a complex64,
        # by checking the accuracy of computations.
        pyfunc = complex_calc
        x = 1.0 + 2**-50
        cfunc = njit((types.float32,))(pyfunc)
        self.assertPreciseEqual(cfunc(x), 1.0)
        # Control (complex128)
        cfunc = njit((types.float64,))(pyfunc)
        self.assertGreater(cfunc(x), 1.0)

    def test_complex2(self):
        pyfunc = docomplex2

        x_types = [
            types.int32, types.int64, types.float32, types.float64
        ]
        x_values = [1, 1000, 12.2, 23.4]
        y_values = [x - 3 for x in x_values]

        for ty, x, y in zip(x_types, x_values, y_values):
            cfunc = njit((ty, ty))(pyfunc)
            self.assertPreciseEqual(pyfunc(x, y), cfunc(x, y),
                prec='single' if ty is types.float32 else 'exact')

        # Check that complex(float32, float32) really creates a complex64,
        # by checking the accuracy of computations.
        pyfunc = complex_calc2
        x = 1.0 + 2**-50
        cfunc = njit((types.float32, types.float32))(pyfunc)
        self.assertPreciseEqual(cfunc(x, x), 2.0)
        # Control (complex128)
        cfunc = njit((types.float64, types.float32))(pyfunc)
        self.assertGreater(cfunc(x, x), 2.0)

    def check_type_converter(self, tp, np_type, values):
        pyfunc = converter(tp)
        cfunc = jit(nopython=True)(pyfunc)
        if issubclass(np_type, np.integer):
            # Converting from a Python int to a small Numpy int on 32-bit
            # builds can raise "OverflowError: Python int too large to
            # convert to C long".  Work around by going through a large
            # Numpy int first.
            np_converter = lambda x: np_type(np.int64(x))
        else:
            np_converter = np_type
        dtype = np.dtype(np_type)
        for val in values:
            if dtype.kind == 'u' and isinstance(val, float) and val < 0.0:
                # Converting negative float to unsigned int yields undefined
                # behaviour (and concretely different on ARM vs. x86)
                continue
            expected = np_converter(val)
            got = cfunc(val)
            self.assertPreciseEqual(got, expected,
                                    msg="for type %s with arg %s" % (np_type, val))

    def check_number_types(self, tp_factory):
        values = [0, 1, -1, 100003, 10000000000007, -100003, -10000000000007,
                  1.5, -3.5]
        for tp_name in real_np_types():
            np_type = getattr(np, tp_name)
            tp = tp_factory(tp_name)
            self.check_type_converter(tp, np_type, values)
        values.append(1.5+3j)
        for tp_name in complex_np_types():
            np_type = getattr(np, tp_name)
            tp = tp_factory(tp_name)
            self.check_type_converter(tp, np_type, values)

    def test_numba_types(self):
        """
        Test explicit casting to Numba number types.
        """
        def tp_factory(tp_name):
            return getattr(types, tp_name)
        self.check_number_types(tp_factory)

    def test_numpy_types(self):
        """
        Test explicit casting to Numpy number types.
        """
        def tp_factory(tp_name):
            return getattr(np, tp_name)
        self.check_number_types(tp_factory)


class TestArrayNumberCtor(TestCase):
    """
    Test <number class>(some sequence)
    """

    def check_type_constructor(self, np_type, values):
        pyfunc = converter(np_type)
        cfunc = jit(nopython=True)(pyfunc)
        for val in values:
            expected = np_type(val)
            got = cfunc(val)
            self.assertPreciseEqual(got, expected)

    def test_1d(self):
        values = [
            (1.0, 2.5),
            (1, 2.5),
            [1.0, 2.5],
            (),
            ]
        for tp_name in real_np_types():
            np_type = getattr(np, tp_name)
            self.check_type_constructor(np_type, values)
        values = [
            (1j, 2.5),
            [1.0, 2.5],
            ]
        for tp_name in complex_np_types():
            np_type = getattr(np, tp_name)
            self.check_type_constructor(np_type, values)

    def test_2d(self):
        values = [
            ((1.0, 2.5), (3.5, 4)),
            [(1.0, 2.5), (3.5, 4.0)],
            ([1.0, 2.5], [3.5, 4.0]),
            [(), ()],
            ]
        for tp_name in real_np_types():
            np_type = getattr(np, tp_name)
            self.check_type_constructor(np_type, values)
        for tp_name in complex_np_types():
            np_type = getattr(np, tp_name)
            self.check_type_constructor(np_type, values)


if __name__ == '__main__':
    unittest.main()
