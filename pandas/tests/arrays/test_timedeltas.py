# -*- coding: utf-8 -*-

import numpy as np

import pandas as pd
from pandas.core.arrays import TimedeltaArrayMixin as TimedeltaArray
import pandas.util.testing as tm


class TestTimedeltaArray(object):

    def test_abs(self):
        vals = np.array([-3600, 'NaT', 7200], dtype='M8[s]')
        arr = TimedeltaArray(vals)

        evals = np.array([3600, 'NaT', 7200], dtype='M8[s]')
        expected = TimedeltaArray(evals)

        result = abs(arr)
        tm.assert_timedelta_array_equal(result, expected)

    def test_neg(self):
        vals = np.array([-3600, 'NaT', 7200], dtype='M8[s]')
        arr = TimedeltaArray(vals)

        evals = np.array([3600, 'NaT', -7200], dtype='M8[s]')
        expected = TimedeltaArray(evals)

        result = -arr
        tm.assert_timedelta_array_equal(result, expected)

    def test_neg_freq(self):
        tdi = pd.timedelta_range('2 Days', periods=4, freq='H')
        arr = TimedeltaArray(tdi)

        expected = TimedeltaArray(-tdi._data, freq=-dti.freq)

        result = -arr
        tm.assert_timedelta_array_equal(result, expected)
