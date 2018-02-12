import pytest

from datetime import datetime, timedelta
from numpy.random import randn
import numpy as np

import pandas as pd
from pandas import (Series, DataFrame, bdate_range,
                    Timestamp, Index)
import pandas.core.window as rwindow
from pandas.errors import UnsupportedFunctionCall
import pandas.util.testing as tm
import pandas.util._test_decorators as td
from pandas.compat import range

N, K = 100, 10


def assert_equal(left, right):
    if isinstance(left, Series):
        tm.assert_series_equal(left, right)
    else:
        tm.assert_frame_equal(left, right)


class Base(object):

    _nan_locs = np.arange(20, 40)
    _inf_locs = np.array([])

    def _create_data(self):
        arr = randn(N)
        arr[self._nan_locs] = np.NaN

        self.arr = arr
        self.rng = bdate_range(datetime(2009, 1, 1), periods=N)
        self.series = Series(arr.copy(), index=self.rng)
        self.frame = DataFrame(randn(N, K), index=self.rng,
                               columns=np.arange(K))


class TestRolling(Base):

    def setup_method(self, method):
        self._create_data()

    def test_doc_string(self):

        df = DataFrame({'B': [0, 1, 2, np.nan, 4]})
        df
        df.rolling(2).sum()
        df.rolling(2, min_periods=1).sum()

    def test_constructor(self):
        # GH 12669

        for o in [self.series, self.frame]:
            c = o.rolling

            # valid
            c(window=2)
            c(window=2, min_periods=1)
            c(window=2, min_periods=1, center=True)
            c(window=2, min_periods=1, center=False)

            # GH 13383
            c(0)
            with pytest.raises(ValueError):
                c(-1)

            # not valid
            for w in [2., 'foo', np.array([2])]:
                with pytest.raises(ValueError):
                    c(window=w)
                with pytest.raises(ValueError):
                    c(window=2, min_periods=w)
                with pytest.raises(ValueError):
                    c(window=2, min_periods=1, center=w)

    @td.skip_if_no_scipy
    def test_constructor_with_win_type(self):
        # GH 13383
        for o in [self.series, self.frame]:
            c = o.rolling
            c(0, win_type='boxcar')
            with pytest.raises(ValueError):
                c(-1, win_type='boxcar')

    def test_constructor_with_timedelta_window(self):
        # GH 15440
        n = 10
        df = DataFrame({'value': np.arange(n)},
                       index=pd.date_range('2015-12-24', periods=n, freq="D"))
        expected_data = np.append([0., 1.], np.arange(3., 27., 3))
        for window in [timedelta(days=3), pd.Timedelta(days=3)]:
            result = df.rolling(window=window).sum()
            expected = DataFrame({'value': expected_data},
                                 index=pd.date_range('2015-12-24', periods=n,
                                                     freq="D"))
            tm.assert_frame_equal(result, expected)
            expected = df.rolling('3D').sum()
            tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        'window', [timedelta(days=3), pd.Timedelta(days=3), '3D'])
    def test_constructor_with_timedelta_window_and_minperiods(self, window):
        # GH 15305
        n = 10
        df = DataFrame({'value': np.arange(n)},
                       index=pd.date_range('2017-08-08', periods=n, freq="D"))
        expected = DataFrame(
            {'value': np.append([np.NaN, 1.], np.arange(3., 27., 3))},
            index=pd.date_range('2017-08-08', periods=n, freq="D"))
        result_roll_sum = df.rolling(window=window, min_periods=2).sum()
        result_roll_generic = df.rolling(window=window,
                                         min_periods=2).apply(sum)
        tm.assert_frame_equal(result_roll_sum, expected)
        tm.assert_frame_equal(result_roll_generic, expected)

    def test_numpy_compat(self):
        # see gh-12811
        r = rwindow.Rolling(Series([2, 4, 6]), window=2)

        msg = "numpy operations are not valid with window objects"

        for func in ('std', 'mean', 'sum', 'max', 'min', 'var'):
            tm.assert_raises_regex(UnsupportedFunctionCall, msg,
                                   getattr(r, func), 1, 2, 3)
            tm.assert_raises_regex(UnsupportedFunctionCall, msg,
                                   getattr(r, func), dtype=np.float64)

    def test_closed(self):
        df = DataFrame({'A': [0, 1, 2, 3, 4]})
        # closed only allowed for datetimelike
        with pytest.raises(ValueError):
            df.rolling(window=3, closed='neither')

    @pytest.mark.parametrize('roller', ['1s', 1])
    def tests_empty_df_rolling(self, roller):
        # GH 15819 Verifies that datetime and integer rolling windows can be
        # applied to empty DataFrames
        expected = DataFrame()
        result = DataFrame().rolling(roller).sum()
        tm.assert_frame_equal(result, expected)

        # Verifies that datetime and integer rolling windows can be applied to
        # empty DataFrames with datetime index
        expected = DataFrame(index=pd.DatetimeIndex([]))
        result = DataFrame(index=pd.DatetimeIndex([])).rolling(roller).sum()
        tm.assert_frame_equal(result, expected)

    def test_missing_minp_zero(self):
        # https://github.com/pandas-dev/pandas/pull/18921
        # minp=0
        x = pd.Series([np.nan])
        result = x.rolling(1, min_periods=0).sum()
        expected = pd.Series([0.0])
        tm.assert_series_equal(result, expected)

        # minp=1
        result = x.rolling(1, min_periods=1).sum()
        expected = pd.Series([np.nan])
        tm.assert_series_equal(result, expected)

    def test_missing_minp_zero_variable(self):
        # https://github.com/pandas-dev/pandas/pull/18921
        x = pd.Series([np.nan] * 4,
                      index=pd.DatetimeIndex(['2017-01-01', '2017-01-04',
                                              '2017-01-06', '2017-01-07']))
        result = x.rolling(pd.Timedelta("2d"), min_periods=0).sum()
        expected = pd.Series(0.0, index=x.index)
        tm.assert_series_equal(result, expected)

    def test_multi_index_names(self):

        # GH 16789, 16825
        cols = pd.MultiIndex.from_product([['A', 'B'], ['C', 'D', 'E']],
                                          names=['1', '2'])
        df = DataFrame(np.ones((10, 6)), columns=cols)
        result = df.rolling(3).cov()

        tm.assert_index_equal(result.columns, df.columns)
        assert result.index.names == [None, '1', '2']


class TestRollingTS(object):

    # rolling time-series friendly
    # xref GH13327

    def setup_method(self, method):

        self.regular = DataFrame({'A': pd.date_range('20130101',
                                                     periods=5,
                                                     freq='s'),
                                  'B': range(5)}).set_index('A')

        self.ragged = DataFrame({'B': range(5)})
        self.ragged.index = [Timestamp('20130101 09:00:00'),
                             Timestamp('20130101 09:00:02'),
                             Timestamp('20130101 09:00:03'),
                             Timestamp('20130101 09:00:05'),
                             Timestamp('20130101 09:00:06')]

    def test_doc_string(self):

        df = DataFrame({'B': [0, 1, 2, np.nan, 4]},
                       index=[Timestamp('20130101 09:00:00'),
                              Timestamp('20130101 09:00:02'),
                              Timestamp('20130101 09:00:03'),
                              Timestamp('20130101 09:00:05'),
                              Timestamp('20130101 09:00:06')])
        df
        df.rolling('2s').sum()

    def test_valid(self):

        df = self.regular

        # not a valid freq
        with pytest.raises(ValueError):
            df.rolling(window='foobar')

        # not a datetimelike index
        with pytest.raises(ValueError):
            df.reset_index().rolling(window='foobar')

        # non-fixed freqs
        for freq in ['2MS', pd.offsets.MonthBegin(2)]:
            with pytest.raises(ValueError):
                df.rolling(window=freq)

        for freq in ['1D', pd.offsets.Day(2), '2ms']:
            df.rolling(window=freq)

        # non-integer min_periods
        for minp in [1.0, 'foo', np.array([1, 2, 3])]:
            with pytest.raises(ValueError):
                df.rolling(window='1D', min_periods=minp)

        # center is not implemented
        with pytest.raises(NotImplementedError):
            df.rolling(window='1D', center=True)

    def test_on(self):

        df = self.regular

        # not a valid column
        with pytest.raises(ValueError):
            df.rolling(window='2s', on='foobar')

        # column is valid
        df = df.copy()
        df['C'] = pd.date_range('20130101', periods=len(df))
        df.rolling(window='2d', on='C').sum()

        # invalid columns
        with pytest.raises(ValueError):
            df.rolling(window='2d', on='B')

        # ok even though on non-selected
        df.rolling(window='2d', on='C').B.sum()

    def test_monotonic_on(self):

        # on/index must be monotonic
        df = DataFrame({'A': pd.date_range('20130101',
                                           periods=5,
                                           freq='s'),
                        'B': range(5)})

        assert df.A.is_monotonic
        df.rolling('2s', on='A').sum()

        df = df.set_index('A')
        assert df.index.is_monotonic
        df.rolling('2s').sum()

        # non-monotonic
        df.index = reversed(df.index.tolist())
        assert not df.index.is_monotonic

        with pytest.raises(ValueError):
            df.rolling('2s').sum()

        df = df.reset_index()
        with pytest.raises(ValueError):
            df.rolling('2s', on='A').sum()

    def test_frame_on(self):

        df = DataFrame({'B': range(5),
                        'C': pd.date_range('20130101 09:00:00',
                                           periods=5,
                                           freq='3s')})

        df['A'] = [Timestamp('20130101 09:00:00'),
                   Timestamp('20130101 09:00:02'),
                   Timestamp('20130101 09:00:03'),
                   Timestamp('20130101 09:00:05'),
                   Timestamp('20130101 09:00:06')]

        # we are doing simulating using 'on'
        expected = (df.set_index('A')
                    .rolling('2s')
                    .B
                    .sum()
                    .reset_index(drop=True)
                    )

        result = (df.rolling('2s', on='A')
                    .B
                    .sum()
                  )
        tm.assert_series_equal(result, expected)

        # test as a frame
        # we should be ignoring the 'on' as an aggregation column
        # note that the expected is setting, computing, and resetting
        # so the columns need to be switched compared
        # to the actual result where they are ordered as in the
        # original
        expected = (df.set_index('A')
                      .rolling('2s')[['B']]
                      .sum()
                      .reset_index()[['B', 'A']]
                    )

        result = (df.rolling('2s', on='A')[['B']]
                    .sum()
                  )
        tm.assert_frame_equal(result, expected)

    def test_frame_on2(self):

        # using multiple aggregation columns
        df = DataFrame({'A': [0, 1, 2, 3, 4],
                        'B': [0, 1, 2, np.nan, 4],
                        'C': Index([Timestamp('20130101 09:00:00'),
                                    Timestamp('20130101 09:00:02'),
                                    Timestamp('20130101 09:00:03'),
                                    Timestamp('20130101 09:00:05'),
                                    Timestamp('20130101 09:00:06')])},
                       columns=['A', 'C', 'B'])

        expected1 = DataFrame({'A': [0., 1, 3, 3, 7],
                               'B': [0, 1, 3, np.nan, 4],
                               'C': df['C']},
                              columns=['A', 'C', 'B'])

        result = df.rolling('2s', on='C').sum()
        expected = expected1
        tm.assert_frame_equal(result, expected)

        expected = Series([0, 1, 3, np.nan, 4], name='B')
        result = df.rolling('2s', on='C').B.sum()
        tm.assert_series_equal(result, expected)

        expected = expected1[['A', 'B', 'C']]
        result = df.rolling('2s', on='C')[['A', 'B', 'C']].sum()
        tm.assert_frame_equal(result, expected)

    def test_basic_regular(self):

        df = self.regular.copy()

        df.index = pd.date_range('20130101', periods=5, freq='D')
        expected = df.rolling(window=1, min_periods=1).sum()
        result = df.rolling(window='1D').sum()
        tm.assert_frame_equal(result, expected)

        df.index = pd.date_range('20130101', periods=5, freq='2D')
        expected = df.rolling(window=1, min_periods=1).sum()
        result = df.rolling(window='2D', min_periods=1).sum()
        tm.assert_frame_equal(result, expected)

        expected = df.rolling(window=1, min_periods=1).sum()
        result = df.rolling(window='2D', min_periods=1).sum()
        tm.assert_frame_equal(result, expected)

        expected = df.rolling(window=1).sum()
        result = df.rolling(window='2D').sum()
        tm.assert_frame_equal(result, expected)

    def test_min_periods(self):

        # compare for min_periods
        df = self.regular

        # these slightly different
        expected = df.rolling(2, min_periods=1).sum()
        result = df.rolling('2s').sum()
        tm.assert_frame_equal(result, expected)

        expected = df.rolling(2, min_periods=1).sum()
        result = df.rolling('2s', min_periods=1).sum()
        tm.assert_frame_equal(result, expected)

    def test_closed(self):

        # xref GH13965

        df = DataFrame({'A': [1] * 5},
                       index=[Timestamp('20130101 09:00:01'),
                              Timestamp('20130101 09:00:02'),
                              Timestamp('20130101 09:00:03'),
                              Timestamp('20130101 09:00:04'),
                              Timestamp('20130101 09:00:06')])

        # closed must be 'right', 'left', 'both', 'neither'
        with pytest.raises(ValueError):
            self.regular.rolling(window='2s', closed="blabla")

        expected = df.copy()
        expected["A"] = [1.0, 2, 2, 2, 1]
        result = df.rolling('2s', closed='right').sum()
        tm.assert_frame_equal(result, expected)

        # default should be 'right'
        result = df.rolling('2s').sum()
        tm.assert_frame_equal(result, expected)

        expected = df.copy()
        expected["A"] = [1.0, 2, 3, 3, 2]
        result = df.rolling('2s', closed='both').sum()
        tm.assert_frame_equal(result, expected)

        expected = df.copy()
        expected["A"] = [np.nan, 1.0, 2, 2, 1]
        result = df.rolling('2s', closed='left').sum()
        tm.assert_frame_equal(result, expected)

        expected = df.copy()
        expected["A"] = [np.nan, 1.0, 1, 1, np.nan]
        result = df.rolling('2s', closed='neither').sum()
        tm.assert_frame_equal(result, expected)

    def test_ragged_sum(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 3, 3, 7]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=2).sum()
        expected = df.copy()
        expected['B'] = [np.nan, np.nan, 3, np.nan, 7]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='3s', min_periods=1).sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 3, 5, 7]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='3s').sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 3, 5, 7]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='4s', min_periods=1).sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 3, 6, 9]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='4s', min_periods=3).sum()
        expected = df.copy()
        expected['B'] = [np.nan, np.nan, 3, 6, 9]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 3, 6, 10]
        tm.assert_frame_equal(result, expected)

    def test_ragged_mean(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).mean()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).mean()
        expected = df.copy()
        expected['B'] = [0.0, 1, 1.5, 3.0, 3.5]
        tm.assert_frame_equal(result, expected)

    def test_ragged_median(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).median()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).median()
        expected = df.copy()
        expected['B'] = [0.0, 1, 1.5, 3.0, 3.5]
        tm.assert_frame_equal(result, expected)

    def test_ragged_quantile(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).quantile(0.5)
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).quantile(0.5)
        expected = df.copy()
        expected['B'] = [0.0, 1, 1.5, 3.0, 3.5]
        tm.assert_frame_equal(result, expected)

    def test_ragged_std(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).std(ddof=0)
        expected = df.copy()
        expected['B'] = [0.0] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='1s', min_periods=1).std(ddof=1)
        expected = df.copy()
        expected['B'] = [np.nan] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='3s', min_periods=1).std(ddof=0)
        expected = df.copy()
        expected['B'] = [0.0] + [0.5] * 4
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).std(ddof=1)
        expected = df.copy()
        expected['B'] = [np.nan, 0.707107, 1.0, 1.0, 1.290994]
        tm.assert_frame_equal(result, expected)

    def test_ragged_var(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).var(ddof=0)
        expected = df.copy()
        expected['B'] = [0.0] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='1s', min_periods=1).var(ddof=1)
        expected = df.copy()
        expected['B'] = [np.nan] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='3s', min_periods=1).var(ddof=0)
        expected = df.copy()
        expected['B'] = [0.0] + [0.25] * 4
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).var(ddof=1)
        expected = df.copy()
        expected['B'] = [np.nan, 0.5, 1.0, 1.0, 1 + 2 / 3.]
        tm.assert_frame_equal(result, expected)

    def test_ragged_skew(self):

        df = self.ragged
        result = df.rolling(window='3s', min_periods=1).skew()
        expected = df.copy()
        expected['B'] = [np.nan] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).skew()
        expected = df.copy()
        expected['B'] = [np.nan] * 2 + [0.0, 0.0, 0.0]
        tm.assert_frame_equal(result, expected)

    def test_ragged_kurt(self):

        df = self.ragged
        result = df.rolling(window='3s', min_periods=1).kurt()
        expected = df.copy()
        expected['B'] = [np.nan] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).kurt()
        expected = df.copy()
        expected['B'] = [np.nan] * 4 + [-1.2]
        tm.assert_frame_equal(result, expected)

    def test_ragged_count(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).count()
        expected = df.copy()
        expected['B'] = [1.0, 1, 1, 1, 1]
        tm.assert_frame_equal(result, expected)

        df = self.ragged
        result = df.rolling(window='1s').count()
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).count()
        expected = df.copy()
        expected['B'] = [1.0, 1, 2, 1, 2]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=2).count()
        expected = df.copy()
        expected['B'] = [np.nan, np.nan, 2, np.nan, 2]
        tm.assert_frame_equal(result, expected)

    def test_regular_min(self):

        df = DataFrame({'A': pd.date_range('20130101',
                                           periods=5,
                                           freq='s'),
                        'B': [0.0, 1, 2, 3, 4]}).set_index('A')
        result = df.rolling('1s').min()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        df = DataFrame({'A': pd.date_range('20130101',
                                           periods=5,
                                           freq='s'),
                        'B': [5, 4, 3, 4, 5]}).set_index('A')

        tm.assert_frame_equal(result, expected)
        result = df.rolling('2s').min()
        expected = df.copy()
        expected['B'] = [5.0, 4, 3, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling('5s').min()
        expected = df.copy()
        expected['B'] = [5.0, 4, 3, 3, 3]
        tm.assert_frame_equal(result, expected)

    def test_ragged_min(self):

        df = self.ragged

        result = df.rolling(window='1s', min_periods=1).min()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).min()
        expected = df.copy()
        expected['B'] = [0.0, 1, 1, 3, 3]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).min()
        expected = df.copy()
        expected['B'] = [0.0, 0, 0, 1, 1]
        tm.assert_frame_equal(result, expected)

    def test_perf_min(self):

        N = 10000

        dfp = DataFrame({'B': np.random.randn(N)},
                        index=pd.date_range('20130101',
                                            periods=N,
                                            freq='s'))
        expected = dfp.rolling(2, min_periods=1).min()
        result = dfp.rolling('2s').min()
        assert ((result - expected) < 0.01).all().bool()

        expected = dfp.rolling(200, min_periods=1).min()
        result = dfp.rolling('200s').min()
        assert ((result - expected) < 0.01).all().bool()

    def test_ragged_max(self):

        df = self.ragged

        result = df.rolling(window='1s', min_periods=1).max()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).max()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).max()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

    def test_ragged_apply(self):

        df = self.ragged

        f = lambda x: 1
        result = df.rolling(window='1s', min_periods=1).apply(f)
        expected = df.copy()
        expected['B'] = 1.
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).apply(f)
        expected = df.copy()
        expected['B'] = 1.
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).apply(f)
        expected = df.copy()
        expected['B'] = 1.
        tm.assert_frame_equal(result, expected)

    def test_all(self):

        # simple comparison of integer vs time-based windowing
        df = self.regular * 2
        er = df.rolling(window=1)
        r = df.rolling(window='1s')

        for f in ['sum', 'mean', 'count', 'median', 'std',
                  'var', 'kurt', 'skew', 'min', 'max']:

            result = getattr(r, f)()
            expected = getattr(er, f)()
            tm.assert_frame_equal(result, expected)

        result = r.quantile(0.5)
        expected = er.quantile(0.5)
        tm.assert_frame_equal(result, expected)

        result = r.apply(lambda x: 1)
        expected = er.apply(lambda x: 1)
        tm.assert_frame_equal(result, expected)

    def test_all2(self):

        # more sophisticated comparison of integer vs.
        # time-based windowing
        df = DataFrame({'B': np.arange(50)},
                       index=pd.date_range('20130101',
                                           periods=50, freq='H')
                       )
        # in-range data
        dft = df.between_time("09:00", "16:00")

        r = dft.rolling(window='5H')

        for f in ['sum', 'mean', 'count', 'median', 'std',
                  'var', 'kurt', 'skew', 'min', 'max']:

            result = getattr(r, f)()

            # we need to roll the days separately
            # to compare with a time-based roll
            # finally groupby-apply will return a multi-index
            # so we need to drop the day
            def agg_by_day(x):
                x = x.between_time("09:00", "16:00")
                return getattr(x.rolling(5, min_periods=1), f)()
            expected = df.groupby(df.index.day).apply(
                agg_by_day).reset_index(level=0, drop=True)

            tm.assert_frame_equal(result, expected)

    def test_groupby_monotonic(self):

        # GH 15130
        # we don't need to validate monotonicity when grouping

        data = [
            ['David', '1/1/2015', 100], ['David', '1/5/2015', 500],
            ['David', '5/30/2015', 50], ['David', '7/25/2015', 50],
            ['Ryan', '1/4/2014', 100], ['Ryan', '1/19/2015', 500],
            ['Ryan', '3/31/2016', 50], ['Joe', '7/1/2015', 100],
            ['Joe', '9/9/2015', 500], ['Joe', '10/15/2015', 50]]

        df = DataFrame(data=data, columns=['name', 'date', 'amount'])
        df['date'] = pd.to_datetime(df['date'])

        expected = df.set_index('date').groupby('name').apply(
            lambda x: x.rolling('180D')['amount'].sum())
        result = df.groupby('name').rolling('180D', on='date')['amount'].sum()
        tm.assert_series_equal(result, expected)

    def test_non_monotonic(self):
        # GH 13966 (similar to #15130, closed by #15175)

        dates = pd.date_range(start='2016-01-01 09:30:00',
                              periods=20, freq='s')
        df = DataFrame({'A': [1] * 20 + [2] * 12 + [3] * 8,
                        'B': np.concatenate((dates, dates)),
                        'C': np.arange(40)})

        result = df.groupby('A').rolling('4s', on='B').C.mean()
        expected = df.set_index('B').groupby('A').apply(
            lambda x: x.rolling('4s')['C'].mean())
        tm.assert_series_equal(result, expected)

        df2 = df.sort_values('B')
        result = df2.groupby('A').rolling('4s', on='B').C.mean()
        tm.assert_series_equal(result, expected)

    def test_rolling_cov_offset(self):
        # GH16058

        idx = pd.date_range('2017-01-01', periods=24, freq='1h')
        ss = Series(np.arange(len(idx)), index=idx)

        result = ss.rolling('2h').cov()
        expected = Series([np.nan] + [0.5] * (len(idx) - 1), index=idx)
        tm.assert_series_equal(result, expected)

        expected2 = ss.rolling(2, min_periods=1).cov()
        tm.assert_series_equal(result, expected2)

        result = ss.rolling('3h').cov()
        expected = Series([np.nan, 0.5] + [1.0] * (len(idx) - 2), index=idx)
        tm.assert_series_equal(result, expected)

        expected2 = ss.rolling(3, min_periods=1).cov()
        tm.assert_series_equal(result, expected2)
