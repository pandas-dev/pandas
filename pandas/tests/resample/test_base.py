# pylint: disable=E1101

from datetime import timedelta

import numpy as np
import pytest

from pandas.compat import range, zip

import pandas as pd
from pandas import DataFrame, Series
from pandas.core.groupby.groupby import DataError
from pandas.core.indexes.datetimes import date_range
from pandas.core.indexes.period import PeriodIndex, period_range
from pandas.core.indexes.timedeltas import TimedeltaIndex, timedelta_range
from pandas.core.resample import TimeGrouper
import pandas.util.testing as tm
from pandas.util.testing import (
    assert_almost_equal, assert_frame_equal, assert_index_equal,
    assert_series_equal)


@pytest.fixture(params=[date_range, period_range, timedelta_range])
def resample_fixture(request, date_range_fixture, period_range_fixture,
                     timedelta_range_fixture):
    if request.param == date_range:
        return date_range_fixture
    if request.param == period_range:
        return period_range_fixture
    if request.param == timedelta_range:
        return timedelta_range_fixture


class TestBase(object):
    """base class for resampling testing"""

    @pytest.mark.parametrize('obj_type', ['series', 'frame'])
    @pytest.mark.parametrize('freq', ['2D', '1H'])
    def test_asfreq(self, freq, obj_type, resample_fixture):
        if resample_fixture._index_factory == period_range:
            pytest.skip('test overridden in test_period_index.py')

        obj = getattr(resample_fixture, obj_type)
        result = obj.resample(freq).asfreq()
        new_index = resample_fixture.create_index(
            obj.index[0], obj.index[-1], freq=freq)
        expected = obj.reindex(new_index)
        assert_almost_equal(result, expected)

    def test_asfreq_fill_value(self, resample_fixture):
        # test for fill value during resampling, issue 3715

        if resample_fixture._index_factory == period_range:
            pytest.skip('test overridden in test_period_index.py')

        s = resample_fixture.series

        result = s.resample('1H').asfreq()
        new_index = resample_fixture.create_index(
            s.index[0], s.index[-1], freq='1H')
        expected = s.reindex(new_index)
        assert_series_equal(result, expected)

        frame = s.to_frame('value')
        frame.iloc[1] = None
        result = frame.resample('1H').asfreq(fill_value=4.0)
        new_index = resample_fixture.create_index(frame.index[0],
                                                  frame.index[-1], freq='1H')
        expected = frame.reindex(new_index, fill_value=4.0)
        assert_frame_equal(result, expected)

    def test_resample_interpolate(self, series):
        # # 12925
        df = series.to_frame('value')
        assert_frame_equal(
            df.resample('1T').asfreq().interpolate(),
            df.resample('1T').interpolate())

    def test_raises_on_non_datetimelike_index(self):
        # this is a non datetimelike index
        xp = DataFrame()
        pytest.raises(TypeError, lambda: xp.resample('A').mean())

    @pytest.mark.parametrize('freq', ['M', 'D', 'H'])
    def test_resample_empty_series(self, freq, resample_method, series):
        # GH12771 & GH12868

        if resample_method == 'ohlc':
            pytest.skip('need to test for ohlc from GH13083')

        s = series[:0]
        result = getattr(s.resample(freq), resample_method)()

        expected = s.copy()
        expected.index = s.index._shallow_copy(freq=freq)
        assert_index_equal(result.index, expected.index)
        assert result.index.freq == expected.index.freq
        assert_series_equal(result, expected, check_dtype=False)

    @pytest.mark.parametrize('freq', ['M', 'D', 'H'])
    def test_resample_empty_dataframe(self, freq, resample_method, series):
        # GH13212
        index = series.index[:0]
        f = DataFrame(index=index)

        # count retains dimensions too
        result = getattr(f.resample(freq), resample_method)()
        if resample_method != 'size':
            expected = f.copy()
        else:
            # GH14962
            expected = Series([])

        expected.index = f.index._shallow_copy(freq=freq)
        assert_index_equal(result.index, expected.index)
        assert result.index.freq == expected.index.freq
        assert_almost_equal(result, expected, check_dtype=False)

        # test size for GH13212 (currently stays as df)

    @pytest.mark.parametrize("tm_index", tm.all_timeseries_index_generator(0))
    @pytest.mark.parametrize(
        "dtype",
        [np.float, np.int, np.object, 'datetime64[ns]'])
    def test_resample_empty_dtypes(self, tm_index, dtype, resample_method):

        # Empty series were sometimes causing a segfault (for the functions
        # with Cython bounds-checking disabled) or an IndexError.  We just run
        # them to ensure they no longer do.  (GH #10228)
        empty_series = Series([], tm_index, dtype)
        try:
            getattr(empty_series.resample('d'), resample_method)()
        except DataError:
            # Ignore these since some combinations are invalid
            # (ex: doing mean with dtype of np.object)
            pass

    def test_resample_loffset_arg_type(self, resample_fixture):
        # GH 13218, 15002
        df = resample_fixture.series.to_frame('value')
        expected_means = [df.values[i:i + 2].mean()
                          for i in range(0, len(df.values), 2)]
        expected_index = resample_fixture.create_index(
            df.index[0], periods=len(df.index) / 2, freq='2D')

        # loffset coerces PeriodIndex to DateTimeIndex
        if isinstance(expected_index, PeriodIndex):
            expected_index = expected_index.to_timestamp()

        expected_index += timedelta(hours=2)
        expected = DataFrame({'value': expected_means}, index=expected_index)

        for arg in ['mean', {'value': 'mean'}, ['mean']]:

            result_agg = df.resample('2D', loffset='2H').agg(arg)

            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                result_how = df.resample('2D', how=arg, loffset='2H')

            if isinstance(arg, list):
                expected.columns = pd.MultiIndex.from_tuples([('value',
                                                               'mean')])

            # GH 13022, 7687 - TODO: fix resample w/ TimedeltaIndex
            if isinstance(expected.index, TimedeltaIndex):
                with pytest.raises(AssertionError):
                    assert_frame_equal(result_agg, expected)
                    assert_frame_equal(result_how, expected)
            else:
                assert_frame_equal(result_agg, expected)
                assert_frame_equal(result_how, expected)

    def test_apply_to_empty_series(self, series):
        # GH 14313
        series = series[:0]

        for freq in ['M', 'D', 'H']:
            result = series.resample(freq).apply(lambda x: 1)
            expected = series.resample(freq).apply(np.sum)

            assert_series_equal(result, expected, check_dtype=False)

    def test_resampler_is_iterable(self, series):
        # GH 15314
        freq = 'H'
        tg = TimeGrouper(freq, convention='start')
        grouped = series.groupby(tg)
        resampled = series.resample(freq)
        for (rk, rv), (gk, gv) in zip(resampled, grouped):
            assert rk == gk
            assert_series_equal(rv, gv)

    def test_resample_quantile(self, series):
        # GH 15023
        s = series
        q = 0.75
        freq = 'H'
        result = s.resample(freq).quantile(q)
        expected = s.resample(freq).agg(lambda x: x.quantile(q))
        tm.assert_series_equal(result, expected)
