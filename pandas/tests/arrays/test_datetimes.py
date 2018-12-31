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
        arr = DatetimeArray(np.array(['2000-01-01T06:00:00'], dtype='M8[ns]'),
                            dtype=DatetimeTZDtype(tz='US/Central'))
        dtype = DatetimeTZDtype(tz='US/Eastern')
        with pytest.raises(TypeError, match='data is already tz-aware'):
            DatetimeArray(arr, dtype=dtype)

    def test_incorrect_dtype_raises(self):
        with pytest.raises(ValueError, match="Unexpected value for 'dtype'."):
            DatetimeArray(np.array([1, 2, 3], dtype='i8'), dtype='category')

    def test_copy(self):
        data = np.array([1, 2, 3], dtype='M8[ns]')
        arr = DatetimeArray(data, copy=False)
        assert arr._data is data

        arr = DatetimeArray(data, copy=True)
        assert arr._data is not data


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

    @pytest.mark.parametrize("dtype", [
        int, np.int32, np.int64, 'uint32', 'uint64',
    ])
    def test_astype_int(self, dtype):
        arr = DatetimeArray._from_sequence([pd.Timestamp('2000'),
                                            pd.Timestamp('2001')])
        result = arr.astype(dtype)

        if np.dtype(dtype).kind == 'u':
            expected_dtype = np.dtype('uint64')
        else:
            expected_dtype = np.dtype('int64')
        expected = arr.astype(expected_dtype)

        assert result.dtype == expected_dtype
        tm.assert_numpy_array_equal(result, expected)

    def test_tz_setter_raises(self):
        arr = DatetimeArray._from_sequence(['2000'], tz='US/Central')
        with pytest.raises(AttributeError, match='tz_localize'):
            arr.tz = 'UTC'

    def test_setitem_different_tz_raises(self):
        data = np.array([1, 2, 3], dtype='M8[ns]')
        arr = DatetimeArray(data, copy=False,
                            dtype=DatetimeTZDtype(tz="US/Central"))
        with pytest.raises(ValueError, match="None"):
            arr[0] = pd.Timestamp('2000')

        with pytest.raises(ValueError, match="US/Central"):
            arr[0] = pd.Timestamp('2000', tz="US/Eastern")

    def test_setitem_clears_freq(self):
        a = DatetimeArray(pd.date_range('2000', periods=2, freq='D',
                                        tz='US/Central'))
        a[0] = pd.Timestamp("2000", tz="US/Central")
        assert a.freq is None

    def test_repeat_preserves_tz(self):
        dti = pd.date_range('2000', periods=2, freq='D', tz='US/Central')
        arr = DatetimeArray(dti)

        repeated = arr.repeat([1, 1])

        # preserves tz and values, but not freq
        expected = DatetimeArray(arr.asi8, freq=None, tz=arr.tz)
        tm.assert_equal(repeated, expected)

    def test_value_counts_preserves_tz(self):
        dti = pd.date_range('2000', periods=2, freq='D', tz='US/Central')
        arr = DatetimeArray(dti).repeat([4, 3])

        result = arr.value_counts()

        # Note: not tm.assert_index_equal, since `freq`s do not match
        assert result.index.equals(dti)

        arr[-2] = pd.NaT
        result = arr.value_counts()
        expected = pd.Series([1, 4, 2],
                             index=[pd.NaT, dti[0], dti[1]])
        tm.assert_series_equal(result, expected)


class TestSequenceToDT64NS(object):

    def test_tz_dtype_mismatch_raises(self):
        arr = DatetimeArray._from_sequence(['2000'], tz='US/Central')
        with pytest.raises(TypeError, match='data is already tz-aware'):
            sequence_to_dt64ns(arr, dtype=DatetimeTZDtype(tz="UTC"))

    def test_tz_dtype_matches(self):
        arr = DatetimeArray._from_sequence(['2000'], tz='US/Central')
        result, _, _ = sequence_to_dt64ns(
            arr, dtype=DatetimeTZDtype(tz="US/Central"))
        tm.assert_numpy_array_equal(arr._data, result)


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
