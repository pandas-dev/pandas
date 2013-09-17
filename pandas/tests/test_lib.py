import unittest
from datetime import datetime

import pandas.lib as lib
import numpy as np


class TestLib(unittest.TestCase):
    def test_maybe_convert_objects_uint64(self):
        # GH4471 - array with objects too big for int64
        arr = np.array([2 ** 63 + 1], dtype=object)
        result = lib.maybe_convert_objects(arr)
        expected = np.array([2 ** 63 + 1], dtype='uint64')
        self.assertEqual(result.dtype, np.dtype('uint64'))
        np.testing.assert_array_equal(result, expected)

        arr2 = np.array([5, 2, 3, 4, 5, 1, 2, 3, 22, 1000, 2**63 + 5,
                         2 ** 63 + 1000], dtype=object)
        result = lib.maybe_convert_objects(arr2)
        expected = arr2.copy().astype('uint64')
        self.assertEqual(result.dtype, np.dtype('uint64'))
        np.testing.assert_array_equal(result, expected)

    def test_maybe_convert_objects_uint64_unconvertible(self):
        # can't convert because negative number
        neg = np.array([-5, 2 ** 63 + 5, 3], dtype=object)
        neg2 = np.array([2 ** 63 + 100, -3], dtype=object)
        # can't convert because of datetime
        dt = np.array([datetime(2011, 5, 3), 2 ** 63 + 2], dtype=object)
        # can't convert because of complex
        cmplx = np.array([2 ** 63 + 5, 1+3j, 22], dtype=object)
        # can't convert b/c of float
        flt = np.array([3.25, 1, 3, 2 ** 63 +4], dtype=object)
        # can't convert b/c of nan
        null = np.array([5, 2, 2 ** 63 + 2, np.nan], dtype=object)
        null2 = np.array([np.nan, 2 ** 63 + 2], dtype=object)
        for arr in (neg, neg2, dt, cmplx, flt, null, null2):
            result = lib.maybe_convert_objects(arr.copy())
            self.assertEqual(result.dtype, np.object_)
            np.testing.assert_array_equal(result, arr)

