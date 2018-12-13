# pylint: disable=E1101

from datetime import datetime, timedelta

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


@pytest.fixture
def create_index(_index_factory):
    def _create_index(*args, **kwargs):
        """ return the _index_factory created using the args, kwargs """
        return _index_factory(*args, **kwargs)
    return _create_index


class Base(object):
    """
    base class for resampling testing, calling
    .create_series() generates a series of each index type
    """

    @pytest.mark.parametrize('freq', ['2D', '1H'])
    def test_asfreq(self, series_and_frame, freq, create_index):
        obj = series_and_frame

        result = obj.resample(freq).asfreq()
        new_index = create_index(obj.index[0], obj.index[-1], freq=freq)
        expected = obj.reindex(new_index)
        assert_almost_equal(result, expected)

    def test_asfreq_fill_value(self, create_index):
        # test for fill value during resampling, issue 3715

        s = self.create_series()

        result = s.resample('1H').asfreq()
        new_index = create_index(s.index[0], s.index[-1], freq='1H')
        expected = s.reindex(new_index)
        assert_series_equal(result, expected)

        frame = s.to_frame('value')
        frame.iloc[1] = None
        result = frame.resample('1H').asfreq(fill_value=4.0)
        new_index = create_index(frame.index[0],
                                 frame.index[-1], freq='1H')
        expected = frame.reindex(new_index, fill_value=4.0)
        assert_frame_equal(result, expected)

    def test_resample_interpolate(self):
        # # 12925
        df = self.create_series().to_frame('value')
        assert_frame_equal(
            df.resample('1T').asfreq().interpolate(),
            df.resample('1T').interpolate())

    def test_raises_on_non_datetimelike_index(self):
        # this is a non datetimelike index
        xp = DataFrame()
        pytest.raises(TypeError, lambda: xp.resample('A').mean())

    @pytest.mark.parametrize('freq', ['M', 'D', 'H'])
    def test_resample_empty_series(self, freq, resample_method):
        # GH12771 & GH12868

        if resample_method == 'ohlc':
            pytest.skip('need to test for ohlc from GH13083')

        s = self.create_series()[:0]
        result = getattr(s.resample(freq), resample_method)()

        expected = s.copy()
        expected.index = s.index._shallow_copy(freq=freq)
        assert_index_equal(result.index, expected.index)
        assert result.index.freq == expected.index.freq
        assert_series_equal(result, expected, check_dtype=False)

    @pytest.mark.parametrize('freq', ['M', 'D', 'H'])
    def test_resample_empty_dataframe(self, freq, resample_method):
        # GH13212
        index = self.create_series().index[:0]
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

    @pytest.mark.parametrize("index", tm.all_timeseries_index_generator(0))
    @pytest.mark.parametrize(
        "dtype",
        [np.float, np.int, np.object, 'datetime64[ns]'])
    def test_resample_empty_dtypes(self, index, dtype, resample_method):

        # Empty series were sometimes causing a segfault (for the functions
        # with Cython bounds-checking disabled) or an IndexError.  We just run
        # them to ensure they no longer do.  (GH #10228)
        empty_series = Series([], index, dtype)
        try:
            getattr(empty_series.resample('d'), resample_method)()
        except DataError:
            # Ignore these since some combinations are invalid
            # (ex: doing mean with dtype of np.object)
            pass

    def test_resample_loffset_arg_type(self, create_index):
        # GH 13218, 15002
        df = self.create_series().to_frame('value')
        expected_means = [df.values[i:i + 2].mean()
                          for i in range(0, len(df.values), 2)]
        expected_index = create_index(df.index[0],
                                      periods=len(df.index) / 2,
                                      freq='2D')

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

    def test_apply_to_empty_series(self):
        # GH 14313
        series = self.create_series()[:0]

        for freq in ['M', 'D', 'H']:
            result = series.resample(freq).apply(lambda x: 1)
            expected = series.resample(freq).apply(np.sum)

            assert_series_equal(result, expected, check_dtype=False)

    def test_resampler_is_iterable(self):
        # GH 15314
        series = self.create_series()
        freq = 'H'
        tg = TimeGrouper(freq, convention='start')
        grouped = series.groupby(tg)
        resampled = series.resample(freq)
        for (rk, rv), (gk, gv) in zip(resampled, grouped):
            assert rk == gk
            assert_series_equal(rv, gv)

    def test_resample_quantile(self):
        # GH 15023
        s = self.create_series()
        q = 0.75
        freq = 'H'
        result = s.resample(freq).quantile(q)
        expected = s.resample(freq).agg(lambda x: x.quantile(q))
        tm.assert_series_equal(result, expected)


class TestDatetimeIndex(Base):
    @pytest.fixture
    def _index_factory(self):
        return date_range

    @pytest.fixture
    def _series_name(self):
        return 'dti'

    def create_series(self):
        i = date_range(datetime(2005, 1, 1),
                       datetime(2005, 1, 10), freq='D')

        return Series(np.arange(len(i)), index=i, name='dti')


class TestPeriodIndex(Base):
    @pytest.fixture
    def _index_factory(self):
        return period_range

    @pytest.fixture
    def _series_name(self):
        return 'pi'

    def create_series(self):
        # TODO: replace calls to .create_series() by injecting the series
        # fixture
        i = period_range(datetime(2005, 1, 1),
                         datetime(2005, 1, 10), freq='D')

        return Series(np.arange(len(i)), index=i, name='pi')

    @pytest.mark.skip()
    def test_asfreq(self):
        pass

    @pytest.mark.skip()
    def test_asfreq_fill_value(self):
        pass


class TestTimedeltaIndex(Base):
    @pytest.fixture
    def _index_factory(self):
        return timedelta_range

    @pytest.fixture
    def _index_start(self):
        return '1 day'

    @pytest.fixture
    def _index_end(self):
        return '10 day'

    @pytest.fixture
    def _series_name(self):
        return 'tdi'

    def create_series(self):
        i = timedelta_range('1 day',
                            '10 day', freq='D')

        return Series(np.arange(len(i)), index=i, name='tdi')
