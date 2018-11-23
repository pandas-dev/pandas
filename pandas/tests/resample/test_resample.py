# pylint: disable=E1101

from datetime import datetime, timedelta

import numpy as np
import pytest

from pandas.compat import OrderedDict, range, zip
from pandas.errors import AbstractMethodError

import pandas as pd
from pandas import DataFrame, Series
from pandas.core.groupby.groupby import DataError
from pandas.core.indexes.datetimes import date_range
from pandas.core.indexes.period import PeriodIndex, period_range
from pandas.core.indexes.timedeltas import TimedeltaIndex
from pandas.core.resample import DatetimeIndex, TimeGrouper
import pandas.util.testing as tm
from pandas.util.testing import (
    assert_almost_equal, assert_frame_equal, assert_index_equal,
    assert_series_equal)

from pandas.tseries.offsets import BDay

bday = BDay()

# The various methods we support
downsample_methods = ['min', 'max', 'first', 'last', 'sum', 'mean', 'sem',
                      'median', 'prod', 'var', 'ohlc', 'quantile']
upsample_methods = ['count', 'size']
series_methods = ['nunique']
resample_methods = downsample_methods + upsample_methods + series_methods


def _simple_ts(start, end, freq='D'):
    rng = date_range(start, end, freq=freq)
    return Series(np.random.randn(len(rng)), index=rng)


def _simple_pts(start, end, freq='D'):
    rng = period_range(start, end, freq=freq)
    return Series(np.random.randn(len(rng)), index=rng)


class TestResampleAPI(object):

    def setup_method(self, method):
        dti = DatetimeIndex(start=datetime(2005, 1, 1),
                            end=datetime(2005, 1, 10), freq='Min')

        self.series = Series(np.random.rand(len(dti)), dti)
        self.frame = DataFrame(
            {'A': self.series, 'B': self.series, 'C': np.arange(len(dti))})

    def test_str(self):

        r = self.series.resample('H')
        assert ('DatetimeIndexResampler [freq=<Hour>, axis=0, closed=left, '
                'label=left, convention=start, base=0]' in str(r))

    def test_api(self):

        r = self.series.resample('H')
        result = r.mean()
        assert isinstance(result, Series)
        assert len(result) == 217

        r = self.series.to_frame().resample('H')
        result = r.mean()
        assert isinstance(result, DataFrame)
        assert len(result) == 217

    def test_groupby_resample_api(self):

        # GH 12448
        # .groupby(...).resample(...) hitting warnings
        # when appropriate
        df = DataFrame({'date': pd.date_range(start='2016-01-01',
                                              periods=4,
                                              freq='W'),
                        'group': [1, 1, 2, 2],
                        'val': [5, 6, 7, 8]}).set_index('date')

        # replication step
        i = pd.date_range('2016-01-03', periods=8).tolist() + \
            pd.date_range('2016-01-17', periods=8).tolist()
        index = pd.MultiIndex.from_arrays([[1] * 8 + [2] * 8, i],
                                          names=['group', 'date'])
        expected = DataFrame({'val': [5] * 7 + [6] + [7] * 7 + [8]},
                             index=index)
        result = df.groupby('group').apply(
            lambda x: x.resample('1D').ffill())[['val']]
        assert_frame_equal(result, expected)

    def test_groupby_resample_on_api(self):

        # GH 15021
        # .groupby(...).resample(on=...) results in an unexpected
        # keyword warning.
        df = DataFrame({'key': ['A', 'B'] * 5,
                        'dates': pd.date_range('2016-01-01', periods=10),
                        'values': np.random.randn(10)})

        expected = df.set_index('dates').groupby('key').resample('D').mean()

        result = df.groupby('key').resample('D', on='dates').mean()
        assert_frame_equal(result, expected)

    def test_pipe(self):
        # GH17905

        # series
        r = self.series.resample('H')
        expected = r.max() - r.mean()
        result = r.pipe(lambda x: x.max() - x.mean())
        tm.assert_series_equal(result, expected)

        # dataframe
        r = self.frame.resample('H')
        expected = r.max() - r.mean()
        result = r.pipe(lambda x: x.max() - x.mean())
        tm.assert_frame_equal(result, expected)

    def test_getitem(self):

        r = self.frame.resample('H')
        tm.assert_index_equal(r._selected_obj.columns, self.frame.columns)

        r = self.frame.resample('H')['B']
        assert r._selected_obj.name == self.frame.columns[1]

        # technically this is allowed
        r = self.frame.resample('H')['A', 'B']
        tm.assert_index_equal(r._selected_obj.columns,
                              self.frame.columns[[0, 1]])

        r = self.frame.resample('H')['A', 'B']
        tm.assert_index_equal(r._selected_obj.columns,
                              self.frame.columns[[0, 1]])

    def test_select_bad_cols(self):

        g = self.frame.resample('H')
        pytest.raises(KeyError, g.__getitem__, ['D'])

        pytest.raises(KeyError, g.__getitem__, ['A', 'D'])
        with pytest.raises(KeyError, match='^[^A]+$'):
            # A should not be referenced as a bad column...
            # will have to rethink regex if you change message!
            g[['A', 'D']]

    def test_attribute_access(self):

        r = self.frame.resample('H')
        tm.assert_series_equal(r.A.sum(), r['A'].sum())

    def test_api_compat_before_use(self):

        # make sure that we are setting the binner
        # on these attributes
        for attr in ['groups', 'ngroups', 'indices']:
            rng = pd.date_range('1/1/2012', periods=100, freq='S')
            ts = Series(np.arange(len(rng)), index=rng)
            rs = ts.resample('30s')

            # before use
            getattr(rs, attr)

            # after grouper is initialized is ok
            rs.mean()
            getattr(rs, attr)

    def tests_skip_nuisance(self):

        df = self.frame
        df['D'] = 'foo'
        r = df.resample('H')
        result = r[['A', 'B']].sum()
        expected = pd.concat([r.A.sum(), r.B.sum()], axis=1)
        assert_frame_equal(result, expected)

        expected = r[['A', 'B', 'C']].sum()
        result = r.sum()
        assert_frame_equal(result, expected)

    def test_downsample_but_actually_upsampling(self):

        # this is reindex / asfreq
        rng = pd.date_range('1/1/2012', periods=100, freq='S')
        ts = Series(np.arange(len(rng), dtype='int64'), index=rng)
        result = ts.resample('20s').asfreq()
        expected = Series([0, 20, 40, 60, 80],
                          index=pd.date_range('2012-01-01 00:00:00',
                                              freq='20s',
                                              periods=5))
        assert_series_equal(result, expected)

    def test_combined_up_downsampling_of_irregular(self):

        # since we are reallydoing an operation like this
        # ts2.resample('2s').mean().ffill()
        # preserve these semantics

        rng = pd.date_range('1/1/2012', periods=100, freq='S')
        ts = Series(np.arange(len(rng)), index=rng)
        ts2 = ts.iloc[[0, 1, 2, 3, 5, 7, 11, 15, 16, 25, 30]]

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            result = ts2.resample('2s', how='mean', fill_method='ffill')
        expected = ts2.resample('2s').mean().ffill()
        assert_series_equal(result, expected)

    def test_transform(self):

        r = self.series.resample('20min')
        expected = self.series.groupby(
            pd.Grouper(freq='20min')).transform('mean')
        result = r.transform('mean')
        assert_series_equal(result, expected)

    def test_fillna(self):

        # need to upsample here
        rng = pd.date_range('1/1/2012', periods=10, freq='2S')
        ts = Series(np.arange(len(rng), dtype='int64'), index=rng)
        r = ts.resample('s')

        expected = r.ffill()
        result = r.fillna(method='ffill')
        assert_series_equal(result, expected)

        expected = r.bfill()
        result = r.fillna(method='bfill')
        assert_series_equal(result, expected)

        with pytest.raises(ValueError):
            r.fillna(0)

    def test_apply_without_aggregation(self):

        # both resample and groupby should work w/o aggregation
        r = self.series.resample('20min')
        g = self.series.groupby(pd.Grouper(freq='20min'))

        for t in [g, r]:
            result = t.apply(lambda x: x)
            assert_series_equal(result, self.series)

    def test_agg_consistency(self):

        # make sure that we are consistent across
        # similar aggregations with and w/o selection list
        df = DataFrame(np.random.randn(1000, 3),
                       index=pd.date_range('1/1/2012', freq='S', periods=1000),
                       columns=['A', 'B', 'C'])

        r = df.resample('3T')

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            expected = r[['A', 'B', 'C']].agg({'r1': 'mean', 'r2': 'sum'})
            result = r.agg({'r1': 'mean', 'r2': 'sum'})
        assert_frame_equal(result, expected)

    # TODO: once GH 14008 is fixed, move these tests into
    # `Base` test class
    def test_agg(self):
        # test with all three Resampler apis and TimeGrouper

        np.random.seed(1234)
        index = date_range(datetime(2005, 1, 1),
                           datetime(2005, 1, 10), freq='D')
        index.name = 'date'
        df = DataFrame(np.random.rand(10, 2), columns=list('AB'), index=index)
        df_col = df.reset_index()
        df_mult = df_col.copy()
        df_mult.index = pd.MultiIndex.from_arrays([range(10), df.index],
                                                  names=['index', 'date'])
        r = df.resample('2D')
        cases = [
            r,
            df_col.resample('2D', on='date'),
            df_mult.resample('2D', level='date'),
            df.groupby(pd.Grouper(freq='2D'))
        ]

        a_mean = r['A'].mean()
        a_std = r['A'].std()
        a_sum = r['A'].sum()
        b_mean = r['B'].mean()
        b_std = r['B'].std()
        b_sum = r['B'].sum()

        expected = pd.concat([a_mean, a_std, b_mean, b_std], axis=1)
        expected.columns = pd.MultiIndex.from_product([['A', 'B'],
                                                       ['mean', 'std']])
        for t in cases:
            result = t.aggregate([np.mean, np.std])
            assert_frame_equal(result, expected)

        expected = pd.concat([a_mean, b_std], axis=1)
        for t in cases:
            result = t.aggregate({'A': np.mean,
                                  'B': np.std})
            assert_frame_equal(result, expected, check_like=True)

        expected = pd.concat([a_mean, a_std], axis=1)
        expected.columns = pd.MultiIndex.from_tuples([('A', 'mean'),
                                                      ('A', 'std')])
        for t in cases:
            result = t.aggregate({'A': ['mean', 'std']})
            assert_frame_equal(result, expected)

        expected = pd.concat([a_mean, a_sum], axis=1)
        expected.columns = ['mean', 'sum']
        for t in cases:
            result = t['A'].aggregate(['mean', 'sum'])
        assert_frame_equal(result, expected)

        expected = pd.concat([a_mean, a_sum], axis=1)
        expected.columns = pd.MultiIndex.from_tuples([('A', 'mean'),
                                                      ('A', 'sum')])
        for t in cases:
            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                result = t.aggregate({'A': {'mean': 'mean', 'sum': 'sum'}})
            assert_frame_equal(result, expected, check_like=True)

        expected = pd.concat([a_mean, a_sum, b_mean, b_sum], axis=1)
        expected.columns = pd.MultiIndex.from_tuples([('A', 'mean'),
                                                      ('A', 'sum'),
                                                      ('B', 'mean2'),
                                                      ('B', 'sum2')])
        for t in cases:
            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                result = t.aggregate({'A': {'mean': 'mean', 'sum': 'sum'},
                                      'B': {'mean2': 'mean', 'sum2': 'sum'}})
            assert_frame_equal(result, expected, check_like=True)

        expected = pd.concat([a_mean, a_std, b_mean, b_std], axis=1)
        expected.columns = pd.MultiIndex.from_tuples([('A', 'mean'),
                                                      ('A', 'std'),
                                                      ('B', 'mean'),
                                                      ('B', 'std')])
        for t in cases:
            result = t.aggregate({'A': ['mean', 'std'],
                                  'B': ['mean', 'std']})
            assert_frame_equal(result, expected, check_like=True)

        expected = pd.concat([a_mean, a_sum, b_mean, b_sum], axis=1)
        expected.columns = pd.MultiIndex.from_tuples([('r1', 'A', 'mean'),
                                                      ('r1', 'A', 'sum'),
                                                      ('r2', 'B', 'mean'),
                                                      ('r2', 'B', 'sum')])

    def test_agg_misc(self):
        # test with all three Resampler apis and TimeGrouper

        np.random.seed(1234)
        index = date_range(datetime(2005, 1, 1),
                           datetime(2005, 1, 10), freq='D')
        index.name = 'date'
        df = DataFrame(np.random.rand(10, 2), columns=list('AB'), index=index)
        df_col = df.reset_index()
        df_mult = df_col.copy()
        df_mult.index = pd.MultiIndex.from_arrays([range(10), df.index],
                                                  names=['index', 'date'])

        r = df.resample('2D')
        cases = [
            r,
            df_col.resample('2D', on='date'),
            df_mult.resample('2D', level='date'),
            df.groupby(pd.Grouper(freq='2D'))
        ]

        # passed lambda
        for t in cases:
            result = t.agg({'A': np.sum,
                            'B': lambda x: np.std(x, ddof=1)})
            rcustom = t['B'].apply(lambda x: np.std(x, ddof=1))
            expected = pd.concat([r['A'].sum(), rcustom], axis=1)
            assert_frame_equal(result, expected, check_like=True)

        # agg with renamers
        expected = pd.concat([t['A'].sum(),
                              t['B'].sum(),
                              t['A'].mean(),
                              t['B'].mean()],
                             axis=1)
        expected.columns = pd.MultiIndex.from_tuples([('result1', 'A'),
                                                      ('result1', 'B'),
                                                      ('result2', 'A'),
                                                      ('result2', 'B')])

        for t in cases:
            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                result = t[['A', 'B']].agg(OrderedDict([('result1', np.sum),
                                                        ('result2', np.mean)]))
            assert_frame_equal(result, expected, check_like=True)

        # agg with different hows
        expected = pd.concat([t['A'].sum(),
                              t['A'].std(),
                              t['B'].mean(),
                              t['B'].std()],
                             axis=1)
        expected.columns = pd.MultiIndex.from_tuples([('A', 'sum'),
                                                      ('A', 'std'),
                                                      ('B', 'mean'),
                                                      ('B', 'std')])
        for t in cases:
            result = t.agg(OrderedDict([('A', ['sum', 'std']),
                                        ('B', ['mean', 'std'])]))
            assert_frame_equal(result, expected, check_like=True)

        # equivalent of using a selection list / or not
        for t in cases:
            result = t[['A', 'B']].agg({'A': ['sum', 'std'],
                                        'B': ['mean', 'std']})
            assert_frame_equal(result, expected, check_like=True)

        # series like aggs
        for t in cases:
            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                result = t['A'].agg({'A': ['sum', 'std']})
            expected = pd.concat([t['A'].sum(),
                                  t['A'].std()],
                                 axis=1)
            expected.columns = pd.MultiIndex.from_tuples([('A', 'sum'),
                                                          ('A', 'std')])
            assert_frame_equal(result, expected, check_like=True)

            expected = pd.concat([t['A'].agg(['sum', 'std']),
                                  t['A'].agg(['mean', 'std'])],
                                 axis=1)
            expected.columns = pd.MultiIndex.from_tuples([('A', 'sum'),
                                                          ('A', 'std'),
                                                          ('B', 'mean'),
                                                          ('B', 'std')])
            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                result = t['A'].agg({'A': ['sum', 'std'],
                                     'B': ['mean', 'std']})
            assert_frame_equal(result, expected, check_like=True)

        # errors
        # invalid names in the agg specification
        for t in cases:
            def f():
                with tm.assert_produces_warning(FutureWarning,
                                                check_stacklevel=False):
                    t[['A']].agg({'A': ['sum', 'std'],
                                  'B': ['mean', 'std']})

            pytest.raises(KeyError, f)

    def test_agg_nested_dicts(self):

        np.random.seed(1234)
        index = date_range(datetime(2005, 1, 1),
                           datetime(2005, 1, 10), freq='D')
        index.name = 'date'
        df = DataFrame(np.random.rand(10, 2), columns=list('AB'), index=index)
        df_col = df.reset_index()
        df_mult = df_col.copy()
        df_mult.index = pd.MultiIndex.from_arrays([range(10), df.index],
                                                  names=['index', 'date'])
        r = df.resample('2D')
        cases = [
            r,
            df_col.resample('2D', on='date'),
            df_mult.resample('2D', level='date'),
            df.groupby(pd.Grouper(freq='2D'))
        ]

        for t in cases:
            def f():
                t.aggregate({'r1': {'A': ['mean', 'sum']},
                             'r2': {'B': ['mean', 'sum']}})
                pytest.raises(ValueError, f)

        for t in cases:
            expected = pd.concat([t['A'].mean(), t['A'].std(), t['B'].mean(),
                                  t['B'].std()], axis=1)
            expected.columns = pd.MultiIndex.from_tuples([('ra', 'mean'), (
                'ra', 'std'), ('rb', 'mean'), ('rb', 'std')])

            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                result = t[['A', 'B']].agg({'A': {'ra': ['mean', 'std']},
                                            'B': {'rb': ['mean', 'std']}})
            assert_frame_equal(result, expected, check_like=True)

            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                result = t.agg({'A': {'ra': ['mean', 'std']},
                                'B': {'rb': ['mean', 'std']}})
            assert_frame_equal(result, expected, check_like=True)

    def test_try_aggregate_non_existing_column(self):
        # GH 16766
        data = [
            {'dt': datetime(2017, 6, 1, 0), 'x': 1.0, 'y': 2.0},
            {'dt': datetime(2017, 6, 1, 1), 'x': 2.0, 'y': 2.0},
            {'dt': datetime(2017, 6, 1, 2), 'x': 3.0, 'y': 1.5}
        ]
        df = DataFrame(data).set_index('dt')

        # Error as we don't have 'z' column
        with pytest.raises(KeyError):
            df.resample('30T').agg({'x': ['mean'],
                                    'y': ['median'],
                                    'z': ['sum']})

    def test_selection_api_validation(self):
        # GH 13500
        index = date_range(datetime(2005, 1, 1),
                           datetime(2005, 1, 10), freq='D')

        rng = np.arange(len(index), dtype=np.int64)
        df = DataFrame({'date': index, 'a': rng},
                       index=pd.MultiIndex.from_arrays([rng, index],
                                                       names=['v', 'd']))
        df_exp = DataFrame({'a': rng}, index=index)

        # non DatetimeIndex
        with pytest.raises(TypeError):
            df.resample('2D', level='v')

        with pytest.raises(ValueError):
            df.resample('2D', on='date', level='d')

        with pytest.raises(TypeError):
            df.resample('2D', on=['a', 'date'])

        with pytest.raises(KeyError):
            df.resample('2D', level=['a', 'date'])

        # upsampling not allowed
        with pytest.raises(ValueError):
            df.resample('2D', level='d').asfreq()

        with pytest.raises(ValueError):
            df.resample('2D', on='date').asfreq()

        exp = df_exp.resample('2D').sum()
        exp.index.name = 'date'
        assert_frame_equal(exp, df.resample('2D', on='date').sum())

        exp.index.name = 'd'
        assert_frame_equal(exp, df.resample('2D', level='d').sum())


class Base(object):
    """
    base class for resampling testing, calling
    .create_series() generates a series of each index type
    """

    def create_index(self, *args, **kwargs):
        """ return the _index_factory created using the args, kwargs """
        factory = self._index_factory()
        return factory(*args, **kwargs)

    @pytest.fixture
    def _index_start(self):
        return datetime(2005, 1, 1)

    @pytest.fixture
    def _index_end(self):
        return datetime(2005, 1, 10)

    @pytest.fixture
    def _index_freq(self):
        return 'D'

    @pytest.fixture
    def index(self, _index_start, _index_end, _index_freq):
        return self.create_index(_index_start, _index_end, freq=_index_freq)

    @pytest.fixture
    def _series_name(self):
        raise AbstractMethodError(self)

    @pytest.fixture
    def _static_values(self, index):
        return np.arange(len(index))

    @pytest.fixture
    def series(self, index, _series_name, _static_values):
        return Series(_static_values, index=index, name=_series_name)

    @pytest.fixture
    def frame(self, index, _static_values):
        return DataFrame({'value': _static_values}, index=index)

    @pytest.fixture(params=[Series, DataFrame])
    def series_and_frame(self, request, index, _series_name, _static_values):
        if request.param == Series:
            return Series(_static_values, index=index, name=_series_name)
        if request.param == DataFrame:
            return DataFrame({'value': _static_values}, index=index)

    @pytest.mark.parametrize('freq', ['2D', '1H'])
    def test_asfreq(self, series_and_frame, freq):
        obj = series_and_frame

        result = obj.resample(freq).asfreq()
        new_index = self.create_index(obj.index[0], obj.index[-1], freq=freq)
        expected = obj.reindex(new_index)
        assert_almost_equal(result, expected)

    def test_asfreq_fill_value(self):
        # test for fill value during resampling, issue 3715

        s = self.create_series()

        result = s.resample('1H').asfreq()
        new_index = self.create_index(s.index[0], s.index[-1], freq='1H')
        expected = s.reindex(new_index)
        assert_series_equal(result, expected)

        frame = s.to_frame('value')
        frame.iloc[1] = None
        result = frame.resample('1H').asfreq(fill_value=4.0)
        new_index = self.create_index(frame.index[0],
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

    def test_resample_empty_series(self):
        # GH12771 & GH12868

        s = self.create_series()[:0]

        for freq in ['M', 'D', 'H']:
            # need to test for ohlc from GH13083
            methods = [method for method in resample_methods
                       if method != 'ohlc']
            for method in methods:
                result = getattr(s.resample(freq), method)()

                expected = s.copy()
                expected.index = s.index._shallow_copy(freq=freq)
                assert_index_equal(result.index, expected.index)
                assert result.index.freq == expected.index.freq
                assert_series_equal(result, expected, check_dtype=False)

    def test_resample_empty_dataframe(self):
        # GH13212
        index = self.create_series().index[:0]
        f = DataFrame(index=index)

        for freq in ['M', 'D', 'H']:
            # count retains dimensions too
            methods = downsample_methods + upsample_methods
            for method in methods:
                result = getattr(f.resample(freq), method)()
                if method != 'size':
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
    def test_resample_empty_dtypes(self, index, dtype):

        # Empty series were sometimes causing a segfault (for the functions
        # with Cython bounds-checking disabled) or an IndexError.  We just run
        # them to ensure they no longer do.  (GH #10228)
        for how in downsample_methods + upsample_methods:
            empty_series = Series([], index, dtype)
            try:
                getattr(empty_series.resample('d'), how)()
            except DataError:
                # Ignore these since some combinations are invalid
                # (ex: doing mean with dtype of np.object)
                pass

    def test_resample_loffset_arg_type(self):
        # GH 13218, 15002
        df = self.create_series().to_frame('value')
        expected_means = [df.values[i:i + 2].mean()
                          for i in range(0, len(df.values), 2)]
        expected_index = self.create_index(df.index[0],
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
