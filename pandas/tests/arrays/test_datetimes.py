# -*- coding: utf-8 -*-
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
    def test_box(self):
        df = pd.DataFrame({'year': [2015, 2016], 'month': [2, 3], 'day': [4, 5]})
        res = pd.to_datetime(df, box=False)
        assert isinstance(res, np.ndarray) == True
        res = pd.to_datetime(df, box=True)
        assert isinstance(res, np.ndarray) == False

    def test_utc(self):
        df = pd.DataFrame({'year': [2015, 2016], 'month': [2, 3], 'day': [4, 5]})
        res = pd.to_datetime(df, utc=True)
        assert str(res[0].tz) == 'UTC'
        res = pd.to_datetime(df, utc=False)
        assert str(res[0].tz) != 'UTC'
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
