# -*- coding: utf-8 -*-

import numpy as np
import pytest

import pandas as pd
from pandas.core.arrays import TimedeltaArrayMixin as TimedeltaArray
import pandas.util.testing as tm


class TestTimedeltaArrayConstructor(object):
    def test_other_type_raises(self):
        with pytest.raises(TypeError,
                           match="dtype bool cannot be converted"):
            TimedeltaArray(np.array([1, 2, 3], dtype='bool'))

    def test_incorrect_dtype_raises(self):
        # TODO: why TypeError for 'category' but ValueError for i8?
        with pytest.raises(TypeError,
                           match='data type "category" not understood'):
            TimedeltaArray(np.array([1, 2, 3], dtype='i8'), dtype='category')

        with pytest.raises(ValueError,
                           match=r"Only timedelta64\[ns\] dtype is valid"):
            TimedeltaArray(np.array([1, 2, 3], dtype='i8'),
                           dtype=np.dtype(int))

    def test_copy(self):
        data = np.array([1, 2, 3], dtype='m8[ns]')
        arr = TimedeltaArray(data, copy=False)
        assert arr._data is data

        arr = TimedeltaArray(data, copy=True)
        assert arr._data is not data
        assert arr._data.base is not data


class TestTimedeltaArray(object):
    def test_from_sequence_dtype(self):
        msg = r"Only timedelta64\[ns\] dtype is valid"
        with pytest.raises(ValueError, match=msg):
            TimedeltaArray._from_sequence([], dtype=object)
        with pytest.raises(ValueError, match=msg):
            TimedeltaArray([], dtype=object)

    def test_abs(self):
        vals = np.array([-3600 * 10**9, 'NaT', 7200 * 10**9], dtype='m8[ns]')
        arr = TimedeltaArray(vals)

        evals = np.array([3600 * 10**9, 'NaT', 7200 * 10**9], dtype='m8[ns]')
        expected = TimedeltaArray(evals)

        result = abs(arr)
        tm.assert_timedelta_array_equal(result, expected)

    def test_neg(self):
        vals = np.array([-3600 * 10**9, 'NaT', 7200 * 10**9], dtype='m8[ns]')
        arr = TimedeltaArray(vals)

        evals = np.array([3600 * 10**9, 'NaT', -7200 * 10**9], dtype='m8[ns]')
        expected = TimedeltaArray(evals)

        result = -arr
        tm.assert_timedelta_array_equal(result, expected)

    def test_neg_freq(self):
        tdi = pd.timedelta_range('2 Days', periods=4, freq='H')
        arr = TimedeltaArray(tdi, freq=tdi.freq)

        expected = TimedeltaArray(-tdi._data, freq=-tdi.freq)

        result = -arr
        tm.assert_timedelta_array_equal(result, expected)

    @pytest.mark.parametrize("dtype", [
        int, np.int32, np.int64, 'uint32', 'uint64',
    ])
    def test_astype_int(self, dtype):
        arr = TimedeltaArray._from_sequence([pd.Timedelta('1H'),
                                             pd.Timedelta('2H')])
        result = arr.astype(dtype)

        if np.dtype(dtype).kind == 'u':
            expected_dtype = np.dtype('uint64')
        else:
            expected_dtype = np.dtype('int64')
        expected = arr.astype(expected_dtype)

        assert result.dtype == expected_dtype
        tm.assert_numpy_array_equal(result, expected)

    def test_setitem_clears_freq(self):
        a = TimedeltaArray(pd.timedelta_range('1H', periods=2, freq='H'))
        a[0] = pd.Timedelta("1H")
        assert a.freq is None
