# -*- coding: utf-8 -*-
"""
Tests for DatetimeArray
"""
import operator

import numpy as np
import pytest

from pandas.core.dtypes.dtypes import DatetimeTZDtype

import pandas as pd
from pandas.core.arrays import DatetimeArrayMixin as DatetimeArray
import pandas.util.testing as tm


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


class TestDatetimeArray(object):
    def test_astype_to_same(self):
        arr = DatetimeArray._from_sequence(['2000'], tz='US/Central')
        result = arr.astype(DatetimeTZDtype(tz="US/Central"), copy=False)
        assert result is arr

    def test_tz_setter_raises(self):
        arr = DatetimeArray._from_sequence(['2000'], tz='US/Central')
        with pytest.raises(AttributeError, match='tz_localize'):
            arr.tz = 'UTC'
