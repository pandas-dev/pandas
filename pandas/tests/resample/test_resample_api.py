# pylint: disable=E1101

from datetime import datetime

import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas import DataFrame, Series, Timestamp
from pandas.compat import range, OrderedDict
from pandas.core.base import SpecificationError
from pandas.core.dtypes.generic import ABCDataFrame, ABCSeries
from pandas.core.indexes.datetimes import date_range
from pandas.core.resample import DatetimeIndex, DatetimeIndexResampler
from pandas.util.testing import assert_frame_equal, assert_series_equal


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

    def test_api_changes_v018(self):

        # change from .resample(....., how=...)
        # to .resample(......).how()

        r = self.series.resample('H')
        assert isinstance(r, DatetimeIndexResampler)

        for how in ['sum', 'mean', 'prod', 'min', 'max', 'var', 'std']:
            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                result = self.series.resample('H', how=how)
                expected = getattr(self.series.resample('H'), how)()
                tm.assert_series_equal(result, expected)

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            result = self.series.resample('H', how='ohlc')
            expected = self.series.resample('H').ohlc()
            tm.assert_frame_equal(result, expected)

        # compat for pandas-like methods
        for how in ['sort_values', 'isna']:
            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                getattr(r, how)()

        # invalids as these can be setting operations
        r = self.series.resample('H')
        pytest.raises(ValueError, lambda: r.iloc[0])
        pytest.raises(ValueError, lambda: r.iat[0])
        pytest.raises(ValueError, lambda: r.loc[0])
        pytest.raises(ValueError, lambda: r.loc[
            Timestamp('2013-01-01 00:00:00', offset='H')])
        pytest.raises(ValueError, lambda: r.at[
            Timestamp('2013-01-01 00:00:00', offset='H')])

        def f():
            r[0] = 5

        pytest.raises(ValueError, f)

        # str/repr
        r = self.series.resample('H')
        with tm.assert_produces_warning(None):
            str(r)
        with tm.assert_produces_warning(None):
            repr(r)

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            tm.assert_numpy_array_equal(np.array(r), np.array(r.mean()))

        # masquerade as Series/DataFrame as needed for API compat
        assert isinstance(self.series.resample('H'), ABCSeries)
        assert not isinstance(self.frame.resample('H'), ABCSeries)
        assert not isinstance(self.series.resample('H'), ABCDataFrame)
        assert isinstance(self.frame.resample('H'), ABCDataFrame)

        # bin numeric ops
        for op in ['__add__', '__mul__', '__truediv__', '__div__', '__sub__']:

            if getattr(self.series, op, None) is None:
                continue
            r = self.series.resample('H')

            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                assert isinstance(getattr(r, op)(2), pd.Series)

        # unary numeric ops
        for op in ['__pos__', '__neg__', '__abs__', '__inv__']:

            if getattr(self.series, op, None) is None:
                continue
            r = self.series.resample('H')

            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                assert isinstance(getattr(r, op)(), pd.Series)

        # comparison ops
        for op in ['__lt__', '__le__', '__gt__', '__ge__', '__eq__', '__ne__']:
            r = self.series.resample('H')

            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                assert isinstance(getattr(r, op)(2), pd.Series)

        # IPython introspection shouldn't trigger warning GH 13618
        for op in ['_repr_json', '_repr_latex',
                   '_ipython_canary_method_should_not_exist_']:
            r = self.series.resample('H')
            with tm.assert_produces_warning(None):
                getattr(r, op, None)

        # getitem compat
        df = self.series.to_frame('foo')

        # same as prior versions for DataFrame
        pytest.raises(KeyError, lambda: df.resample('H')[0])

        # compat for Series
        # but we cannot be sure that we need a warning here
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            result = self.series.resample('H')[0]
            expected = self.series.resample('H').mean()[0]
            assert result == expected

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            result = self.series.resample('H')['2005-01-09 23:00:00']
            expected = self.series.resample('H').mean()['2005-01-09 23:00:00']
            assert result == expected

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
        df = pd.DataFrame({'key': ['A', 'B'] * 5,
                           'dates': pd.date_range('2016-01-01', periods=10),
                           'values': np.random.randn(10)})

        expected = df.set_index('dates').groupby('key').resample('D').mean()

        result = df.groupby('key').resample('D', on='dates').mean()
        assert_frame_equal(result, expected)

    def test_plot_api(self):
        tm._skip_if_no_mpl()

        # .resample(....).plot(...)
        # hitting warnings
        # GH 12448
        s = Series(np.random.randn(60),
                   index=date_range('2016-01-01', periods=60, freq='1min'))
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            result = s.resample('15min').plot()
            tm.assert_is_valid_plot_return_object(result)

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            result = s.resample('15min', how='sum').plot()
            tm.assert_is_valid_plot_return_object(result)

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
        with tm.assert_raises_regex(KeyError, '^[^A]+$'):
            # A should not be referenced as a bad column...
            # will have to rethink regex if you change message!
            g[['A', 'D']]

    def test_attribute_access(self):

        r = self.frame.resample('H')
        tm.assert_series_equal(r.A.sum(), r['A'].sum())

        # getting
        pytest.raises(AttributeError, lambda: r.F)

        # setting
        def f():
            r.F = 'bah'

        pytest.raises(ValueError, f)

    def test_api_compat_before_use(self):

        # make sure that we are setting the binner
        # on these attributes
        for attr in ['groups', 'ngroups', 'indices']:
            rng = pd.date_range('1/1/2012', periods=100, freq='S')
            ts = pd.Series(np.arange(len(rng)), index=rng)
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
        ts = pd.Series(np.arange(len(rng), dtype='int64'), index=rng)
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
        ts = pd.Series(np.arange(len(rng)), index=rng)
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
        ts = pd.Series(np.arange(len(rng), dtype='int64'), index=rng)
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
        df = pd.DataFrame(np.random.rand(10, 2),
                          columns=list('AB'),
                          index=index)
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
        df = pd.DataFrame(np.random.rand(10, 2),
                          columns=list('AB'),
                          index=index)
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

            pytest.raises(SpecificationError, f)

    def test_agg_nested_dicts(self):

        np.random.seed(1234)
        index = date_range(datetime(2005, 1, 1),
                           datetime(2005, 1, 10), freq='D')
        index.name = 'date'
        df = pd.DataFrame(np.random.rand(10, 2),
                          columns=list('AB'),
                          index=index)
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

    def test_selection_api_validation(self):
        # GH 13500
        index = date_range(datetime(2005, 1, 1),
                           datetime(2005, 1, 10), freq='D')
        df = pd.DataFrame({'date': index,
                           'a': np.arange(len(index), dtype=np.int64)},
                          index=pd.MultiIndex.from_arrays([
                              np.arange(len(index), dtype=np.int64),
                              index], names=['v', 'd']))
        df_exp = pd.DataFrame({'a': np.arange(len(index), dtype=np.int64)},
                              index=index)

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
