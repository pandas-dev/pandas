# -*- coding: utf-8 -*-

import numpy as np
import pytest

import pandas as pd
from pandas.core.arrays import TimedeltaArrayMixin as TimedeltaArray
import pandas.util.testing as tm


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
