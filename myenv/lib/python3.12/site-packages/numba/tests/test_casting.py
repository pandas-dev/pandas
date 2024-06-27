import numpy as np
from numba.core.errors import TypingError
from numba import njit
from numba.core import types
import struct
import unittest


def float_to_int(x):
    return types.int32(x)


def int_to_float(x):
    return types.float64(x) / 2


def float_to_unsigned(x):
    return types.uint32(x)


def float_to_complex(x):
    return types.complex128(x)


def numpy_scalar_cast_error():
    np.int32(np.zeros((4,)))

class TestCasting(unittest.TestCase):
    def test_float_to_int(self):
        pyfunc = float_to_int
        cfunc = njit((types.float32,))(pyfunc)

        self.assertEqual(cfunc.nopython_signatures[0].return_type, types.int32)
        self.assertEqual(cfunc(12.3), pyfunc(12.3))
        self.assertEqual(cfunc(12.3), int(12.3))
        self.assertEqual(cfunc(-12.3), pyfunc(-12.3))
        self.assertEqual(cfunc(-12.3), int(-12.3))

    def test_int_to_float(self):
        pyfunc = int_to_float
        cfunc = njit((types.int64,))(pyfunc)

        self.assertEqual(cfunc.nopython_signatures[0].return_type,
                         types.float64)
        self.assertEqual(cfunc(321), pyfunc(321))
        self.assertEqual(cfunc(321), 321. / 2)

    def test_float_to_unsigned(self):
        pyfunc = float_to_unsigned
        cfunc = njit((types.float32,))(pyfunc)

        self.assertEqual(cfunc.nopython_signatures[0].return_type, types.uint32)
        self.assertEqual(cfunc(3.21), pyfunc(3.21))
        self.assertEqual(cfunc(3.21), struct.unpack('I', struct.pack('i',
                                                                      3))[0])

    def test_float_to_complex(self):
        pyfunc = float_to_complex
        cfunc = njit((types.float64,))(pyfunc)
        self.assertEqual(cfunc.nopython_signatures[0].return_type,
                         types.complex128)
        self.assertEqual(cfunc(-3.21), pyfunc(-3.21))
        self.assertEqual(cfunc(-3.21), -3.21 + 0j)

    def test_array_to_array(self):
        """Make sure this compiles.

        Cast C to A array
        """
        @njit("f8(f8[:])")
        def inner(x):
            return x[0]

        inner.disable_compile()

        @njit("f8(f8[::1])")
        def driver(x):
            return inner(x)

        x = np.array([1234], dtype=np.float64)
        self.assertEqual(driver(x), x[0])
        self.assertEqual(len(inner.overloads), 1)

    def test_0darrayT_to_T(self):
        @njit
        def inner(x):
            return x.dtype.type(x)

        inputs = [
            (np.bool_, True),
            (np.float32, 12.3),
            (np.float64, 12.3),
            (np.int64, 12),
            (np.complex64, 2j+3),
            (np.complex128, 2j+3),
            (np.timedelta64, np.timedelta64(3, 'h')),
            (np.datetime64, np.datetime64('2016-01-01')),
            ('<U3', 'ABC'),
        ]

        for (T, inp) in inputs:
            x = np.array(inp, dtype=T)
            self.assertEqual(inner(x), x[()])

    def test_array_to_scalar(self):
        """
        Ensure that a TypingError exception is raised if
        user tries to convert numpy array to scalar
        """

        with self.assertRaises(TypingError) as raises:
            njit(())(numpy_scalar_cast_error)

        self.assertIn("Casting array(float64, 1d, C) to int32 directly is unsupported.",
                      str(raises.exception))

    def test_optional_to_optional(self):
        """
        Test error due mishandling of Optional to Optional casting

        Related issue: https://github.com/numba/numba/issues/1718
        """
        # Attempt to cast optional(intp) to optional(float64)
        opt_int = types.Optional(types.intp)
        opt_flt = types.Optional(types.float64)
        sig = opt_flt(opt_int)

        @njit(sig)
        def foo(a):
            return a

        self.assertEqual(foo(2), 2)
        self.assertIsNone(foo(None))


if __name__ == '__main__':
    unittest.main()
