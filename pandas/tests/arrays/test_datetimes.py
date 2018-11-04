"""
Tests for DatetimeArray
"""
import operator

import numpy as np

import pandas as pd
from pandas.core.arrays import DatetimeArrayMixin as DatetimeArray
import pandas.util.testing as tm


class TestDatetimeArrayComparisons(object):
    # TODO: merge this into tests/arithmetic/test_datetime64 once it is
    #  sufficiently robust

    def test_cmp_dt64_arraylike_tznaive(self):
        # arbitrary tz-naive DatetimeIndex
        dti = pd.date_range('2016-01-1', freq='MS', periods=9, tz=None)
        arr = DatetimeArray(dti)
        assert arr.freq == dti.freq
        assert arr.tz == dti.tz

        right = dti

        expected = np.ones(len(arr), dtype=bool)

        for op in [operator.eq, operator.le, operator.ge]:
            result = op(arr, arr)
            tm.assert_numpy_array_equal(result, expected)
            for other in [right, np.array(right)]:
                # TODO: add list and tuple, and object-dtype once those
                #  are fixed in the constructor
                result = op(arr, other)
                tm.assert_numpy_array_equal(result, expected)

                result = op(other, arr)
                tm.assert_numpy_array_equal(result, expected)

        # !=, <, >
        expected = np.zeros(len(dti), dtype=bool)
        tm.assert_numpy_array_equal(arr != arr, expected)

        for op in [operator.ne, operator.lt, operator.gt]:
            result = op(arr, arr)
            tm.assert_numpy_array_equal(result, expected)
            for other in [right, np.array(right)]:
                # TODO: add list and tuple, and object-dtype once those
                #  are fixed in the constructor
                result = op(arr, other)
                tm.assert_numpy_array_equal(result, expected)

                result = op(other, arr)
                tm.assert_numpy_array_equal(result, expected)
