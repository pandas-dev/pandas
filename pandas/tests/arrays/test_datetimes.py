# -*- coding: utf-8 -*-
"""
Tests for DatetimeArray
"""
import operator

import numpy as np

import pandas as pd
from pandas.core.arrays import DatetimeArrayMixin as DatetimeArray
import pandas.util.testing as tm


class TestDatetimeArray(object):
    def test_repr(self):
        dti = pd.date_range('1994-07-01', periods=9, freq='W', tz='US/Central')
        arr = DatetimeArray(dti)

        # non-truncated
        expected = (
            "<DatetimeArrayMixin>\n"
            "['1994-07-03 00:00:00-05:00', '1994-07-10 00:00:00-05:00']\n"
            "Length: 2, dtype: datetime64[ns, US/Central], freq: W-SUN")
        result = repr(arr[:2])
        assert result == expected

        # truncated
        expected = (
            "<DatetimeArrayMixin>\n"
            "["
            "'1994-07-03 00:00:00-05:00', "
            "'1994-07-10 00:00:00-05:00', "
            "'1994-07-17 00:00:00-05:00', "
            "'...', "
            "'1994-08-14 00:00:00-05:00', "
            "'1994-08-21 00:00:00-05:00', "
            "'1994-08-28 00:00:00-05:00'"
            "]\n"
            "Length: 9, dtype: datetime64[ns, US/Central], freq: W-SUN")
        result = repr(arr)
        assert result == expected


class TestDatetimeArrayComparisons(object):
    # TODO: merge this into tests/arithmetic/test_datetime64 once it is
    #  sufficiently robust

    def test_cmp_dt64_arraylike_tznaive(self, all_compare_operators):
        # arbitrary tz-naive DatetimeIndex
        opname = all_compare_operators.strip('_')
        op = getattr(operator, opname)

        dti = pd.date_range('2016-01-1', freq='MS', periods=9, tz=None)
        arr = DatetimeArray(dti)
        assert arr.freq == dti.freq
        assert arr.tz == dti.tz

        right = dti

        expected = np.ones(len(arr), dtype=bool)
        if opname in ['ne', 'gt', 'lt']:
            # for these the comparisons should be all-False
            expected = ~expected

        result = op(arr, arr)
        tm.assert_numpy_array_equal(result, expected)
        for other in [right, np.array(right)]:
            # TODO: add list and tuple, and object-dtype once those
            #  are fixed in the constructor
            result = op(arr, other)
            tm.assert_numpy_array_equal(result, expected)

            result = op(other, arr)
            tm.assert_numpy_array_equal(result, expected)
