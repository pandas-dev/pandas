""" test Series-boxed versions of the scalar Timedelta """
from datetime import timedelta

import pytest

import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import Timedelta, timedelta_range, Series, to_timedelta


class TestSeriesTimedeltas(object):
    def test_apply_to_timedelta(self):
        timedelta_NaT = to_timedelta('NaT')

        list_of_valid_strings = ['00:00:01', '00:00:02']
        a = to_timedelta(list_of_valid_strings)
        b = Series(list_of_valid_strings).apply(to_timedelta)
        # Can't compare until apply on a Series gives the correct dtype
        # assert_series_equal(a, b)

        list_of_strings = ['00:00:01', np.nan, pd.NaT, timedelta_NaT]

        # TODO: unused?
        a = to_timedelta(list_of_strings)  # noqa
        b = Series(list_of_strings).apply(to_timedelta)  # noqa
        # Can't compare until apply on a Series gives the correct dtype
        # assert_series_equal(a, b)

    def test_components(self):
        rng = timedelta_range('1 days, 10:11:12', periods=2, freq='s')
        rng.components

        # with nat
        s = Series(rng)
        s[1] = np.nan

        result = s.dt.components
        assert not result.iloc[0].isna().all()
        assert result.iloc[1].isna().all()

    def test_timedelta_arithmetic(self):
        data = Series(['nat', '32 days'], dtype='timedelta64[ns]')
        deltas = [timedelta(days=1), Timedelta(1, unit='D')]
        for delta in deltas:
            result_method = data.add(delta)
            result_operator = data + delta
            expected = Series(['nat', '33 days'], dtype='timedelta64[ns]')
            tm.assert_series_equal(result_operator, expected)
            tm.assert_series_equal(result_method, expected)

            result_method = data.sub(delta)
            result_operator = data - delta
            expected = Series(['nat', '31 days'], dtype='timedelta64[ns]')
            tm.assert_series_equal(result_operator, expected)
            tm.assert_series_equal(result_method, expected)
            # GH 9396
            result_method = data.div(delta)
            result_operator = data / delta
            expected = Series([np.nan, 32.], dtype='float64')
            tm.assert_series_equal(result_operator, expected)
            tm.assert_series_equal(result_method, expected)

    def test_overflow(self):
        # GH 9442
        s = Series(pd.date_range('20130101', periods=100000, freq='H'))
        s[0] += Timedelta('1s 1ms')

        # mean
        result = (s - s.min()).mean()
        expected = Timedelta((pd.DatetimeIndex((s - s.min())).asi8 / len(s)
                              ).sum())

        # the computation is converted to float so
        # might be some loss of precision
        assert np.allclose(result.value / 1000, expected.value / 1000)

        # sum
        pytest.raises(ValueError, lambda: (s - s.min()).sum())
        s1 = s[0:10000]
        pytest.raises(ValueError, lambda: (s1 - s1.min()).sum())
        s2 = s[0:1000]
        result = (s2 - s2.min()).sum()
