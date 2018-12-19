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
from pandas.core.arrays.datetimes import sequence_to_dt64ns
import pandas.util.testing as tm


class TestDatetimeArrayConstructor(object):
    def test_mismatched_timezone_raises(self):
        a = DatetimeArray(np.array(['2000-01-01T06:00:00'], dtype='M8[ns]'),
                          dtype=DatetimeTZDtype(tz='US/Central'))
        dtype = DatetimeTZDtype(tz='US/Eastern')
        with pytest.raises(TypeError, match='do not match'):
            DatetimeArray(a, dtype=dtype)

    def test_non_array_raises(self):
        with pytest.raises(ValueError, match='list'):
            DatetimeArray([1, 2, 3])

    def test_other_type_raises(self):
        with pytest.raises(ValueError,
                           match="The dtype of 'values' is incorrect.*bool"):
            DatetimeArray(np.array([1, 2, 3], dtype='bool'))

    def test_incorrect_dtype_raises(self):
        with pytest.raises(ValueError, match="Unexpected value for 'dtype'."):
            DatetimeArray(np.array([1, 2, 3], dtype='i8'), dtype='category')

    def test_freq_infer_raises(self):
        with pytest.raises(ValueError, match='Frequency inference'):
            DatetimeArray(np.array([1, 2, 3], dtype='i8'), freq="infer")

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
        with pytest.raises(ValueError, match="None"):
            arr[0] = pd.Timestamp('2000')

        with pytest.raises(ValueError, match="US/Central"):
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


class TestDatetimeArray(object):
    def test_astype_to_same(self):
        arr = DatetimeArray._from_sequence(['2000'], tz='US/Central')
        result = arr.astype(DatetimeTZDtype(tz="US/Central"), copy=False)
        assert result is arr

    def test_tz_setter_raises(self):
        arr = DatetimeArray._from_sequence(['2000'], tz='US/Central')
        with pytest.raises(AttributeError, match='tz_localize'):
            arr.tz = 'UTC'


class TestSequenceToDT64NS(object):

    def test_tz_dtype_mismatch_raises(self):
        arr = DatetimeArray._from_sequence(['2000'], tz='US/Central')
        with pytest.raises(TypeError, match='do not match'):
            sequence_to_dt64ns(arr, dtype=DatetimeTZDtype(tz="UTC"))

    def test_tz_dtype_matches(self):
        arr = DatetimeArray._from_sequence(['2000'], tz='US/Central')
        result, _, _ = sequence_to_dt64ns(
            arr, dtype=DatetimeTZDtype(tz="US/Central"))
        tm.assert_extension_array_equal(arr, result)


class TestReductions(object):

    @pytest.mark.parametrize("tz", [None, "US/Central"])
    def test_min_max(self, tz):
        arr = DatetimeArray._from_sequence([
            '2000-01-03',
            '2000-01-03',
            'NaT',
            '2000-01-02',
            '2000-01-05',
            '2000-01-04',
        ], tz=tz)

        result = arr.min()
        expected = pd.Timestamp('2000-01-02', tz=tz)
        assert result == expected

        result = arr.max()
        expected = pd.Timestamp('2000-01-05', tz=tz)
        assert result == expected

        result = arr.min(skipna=False)
        assert result is pd.NaT

        result = arr.max(skipna=False)
        assert result is pd.NaT

    @pytest.mark.parametrize("tz", [None, "US/Central"])
    @pytest.mark.parametrize('skipna', [True, False])
    def test_min_max_empty(self, skipna, tz):
        arr = DatetimeArray._from_sequence([], tz=tz)
        result = arr.min(skipna=skipna)
        assert result is pd.NaT

        result = arr.max(skipna=skipna)
        assert result is pd.NaT
