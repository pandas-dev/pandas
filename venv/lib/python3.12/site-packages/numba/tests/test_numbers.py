# (Some) Tests for targets/numbers.py

import numpy as np

from numba import njit
from numba.core import types
from numba.core.errors import TypingError
from numba.tests.support import TestCase


def gen_view(a,b):
    def impl(x):
        return a(x).view(b)
    return impl


class TestViewIntFloat(TestCase):
    """ This tests the 'view' method on NumPy scalars. """

    def do_testing(self, inputs, dtypes):
        for value, initial_type, expected in inputs:
            for target_type, result in zip(dtypes, expected):
                view = njit(gen_view(initial_type, target_type))
                if not np.isnan(result):
                    # check against predefined value
                    self.assertEqual(view(value), target_type(result))
                    # check against numpy
                    self.assertEqual(view(value),
                                     view.py_func(value))
                else:
                    # check that our implementation results in nan
                    self.assertTrue(np.isnan(view(value)))
                    # check that numpy results in nan
                    self.assertTrue(np.isnan(view.py_func(value)))

    def test_8_bits(self):
        dtypes = (np.uint8, np.int8)
        #        Value  Initial Type   Expected answers using dtypes
        inputs = ((1,   np.uint8,    (1, 1)),
                  (-1,  np.int8,     (255, -1)))
        self.do_testing(inputs, dtypes)

    def test_32_bits(self):
        dtypes = (np.uint32, np.int32, np.float32)
        #        Value  Initial Type   Expected answers using dtypes
        inputs = ((1,   np.uint32,    (1, 1, 1.401298464324817e-45)),
                  (-1,  np.int32,     (4294967295, -1, np.nan)),
                  (1.0, np.float32,   (1065353216, 1065353216, 1.0)))
        self.do_testing(inputs, dtypes)

    def test_64_bits(self):
        dtypes = (np.uint64, np.int64, np.float64)
        #        Value  Initial Type   Expected answers using dtypes
        inputs = ((1,   np.uint64,    (1, 1, 5e-324)),
                  (-1,  np.int64,     (18446744073709551615, -1, np.nan)),
                  (1.0, np.float64,   (4607182418800017408,
                                       4607182418800017408,
                                       1.0))
                  )
        self.do_testing(inputs, dtypes)

    def test_python_scalar_exception(self):
        intty = getattr(np, 'int{}'.format(types.intp.bitwidth))

        @njit
        def myview():
            a = 1
            a.view(intty)

        with self.assertRaises(TypingError) as e:
            myview()
        self.assertIn("'view' can only be called on NumPy dtypes, "
                      "try wrapping the variable 'a' with 'np.<dtype>()'",
                      str(e.exception))

    def do_testing_exceptions(self, pair):
        with self.assertRaises(TypingError) as e:
            view = njit(gen_view(pair[0], pair[1]))
            view(1)
        self.assertIn("Changing the dtype of a 0d array is only supported "
                      "if the itemsize is unchanged",
                      str(e.exception))

    def test_exceptions32(self):
        for pair in ((np.int32, np.int8), (np.int8, np.int32)):
            self.do_testing_exceptions(pair)

    def test_exceptions64(self):
        for pair in ((np.int32, np.int64), (np.int64, np.int32)):
            self.do_testing_exceptions(pair)
