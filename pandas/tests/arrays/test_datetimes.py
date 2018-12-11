# -*- coding: utf-8 -*-
"""
Tests for DatetimeArray
"""
import operator

import numpy as np
import pytest

from pandas.errors import IncompatibleTimeZoneError

from pandas.core.dtypes.dtypes import DatetimeTZDtype

import pandas as pd
from pandas.core.arrays import DatetimeArrayMixin as DatetimeArray
import pandas.util.testing as tm


class TestDatetimeArrayConstructor(object):
    def test_mismatched_timezone_raises(self):
        a = DatetimeArray(np.array(['2000-01-01T06:00:00'], dtype='M8[ns]'),
                          dtype=DatetimeTZDtype(tz='US/Central'))
        dtype = DatetimeTZDtype(tz='US/Eastern')
        with pytest.raises(ValueError, match='Timezones'):
            DatetimeArray(a, dtype=dtype)

    def test_non_array_raises(self):
        with pytest.raises(ValueError, match='list'):
            DatetimeArray([1, 2, 3])

    def test_other_type_raises(self):
        with pytest.raises(ValueError,
                           match="The dtype of 'values' is incorrect"):
            DatetimeArray(np.array([1, 2, 3], dtype='bool'))

    def test_incorrect_dtype_raises(self):
        with pytest.raises(ValueError, match="Unexpected value for 'dtype'."):
            DatetimeArray(np.array([1, 2, 3], dtype='i8'), dtype='category')

    def test_freq_infer_raises(self):
        with pytest.raises(ValueError, match='Frequency inference'):
            DatetimeArray(np.array([1, 2, 3]), freq="infer")

    def test_copy(self):
        data = np.array([1, 2, 3], dtype='M8[ns]')
        arr = DatetimeArray(data, copy=False)
        assert arr._data is data

        arr = DatetimeArray(data, copy=True)
        assert arr._data is not data


class TestSetitem(object):
    def test_set_different_tz_raises(self):
        data = np.array([1, 2, 3], dtype='M8[ns]')
        arr = DatetimeArray(data, copy=False,
                            dtype=DatetimeTZDtype(tz="US/Central"))
        with pytest.raises(IncompatibleTimeZoneError, match="None"):
            arr[0] = pd.Timestamp('2000')

        with pytest.raises(IncompatibleTimeZoneError, match="US/Central"):
            arr[0] = pd.Timestamp('2000', tz="US/Eastern")


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
