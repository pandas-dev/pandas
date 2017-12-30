# -*- coding: utf-8 -*-

"""
test .agg behavior / note that .apply is tested generally in test_groupby.py
"""

import pytest

import numpy as np
import pandas as pd

from pandas import concat, DataFrame, Index, MultiIndex, Series
from pandas.core.groupby import SpecificationError
from pandas.compat import OrderedDict
import pandas.util.testing as tm


class TestGroupByAggregate(object):

    def setup_method(self, method):
        self.ts = tm.makeTimeSeries()

        self.seriesd = tm.getSeriesData()
        self.tsd = tm.getTimeSeriesData()
        self.frame = DataFrame(self.seriesd)
        self.tsframe = DataFrame(self.tsd)

        self.df = DataFrame(
            {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
             'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
             'C': np.random.randn(8),
             'D': np.random.randn(8)})

        self.df_mixed_floats = DataFrame(
            {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
             'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
             'C': np.random.randn(8),
             'D': np.array(np.random.randn(8), dtype='float32')})

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        self.mframe = DataFrame(np.random.randn(10, 3), index=index,
                                columns=['A', 'B', 'C'])

        self.three_group = DataFrame(
            {'A': ['foo', 'foo', 'foo', 'foo', 'bar', 'bar', 'bar', 'bar',
                   'foo', 'foo', 'foo'],
             'B': ['one', 'one', 'one', 'two', 'one', 'one', 'one', 'two',
                   'two', 'two', 'one'],
             'C': ['dull', 'dull', 'shiny', 'dull', 'dull', 'shiny', 'shiny',
                   'dull', 'shiny', 'shiny', 'shiny'],
             'D': np.random.randn(11),
             'E': np.random.randn(11),
             'F': np.random.randn(11)})

    def test_agg_regression1(self):
        grouped = self.tsframe.groupby([lambda x: x.year, lambda x: x.month])
        result = grouped.agg(np.mean)
        expected = grouped.mean()
        tm.assert_frame_equal(result, expected)

    def test_agg_must_agg(self):
        grouped = self.df.groupby('A')['C']

        msg = "Must produce aggregated value"
        with tm.assert_raises_regex(Exception, msg):
            grouped.agg(lambda x: x.describe())
        with tm.assert_raises_regex(Exception, msg):
            grouped.agg(lambda x: x.index[:2])

    def test_agg_ser_multi_key(self):
        # TODO(wesm): unused
        ser = self.df.C  # noqa

        f = lambda x: x.sum()
        results = self.df.C.groupby([self.df.A, self.df.B]).aggregate(f)
        expected = self.df.groupby(['A', 'B']).sum()['C']
        tm.assert_series_equal(results, expected)

    def test_agg_apply_corner(self):
        # nothing to group, all NA
        grouped = self.ts.groupby(self.ts * np.nan)
        assert self.ts.dtype == np.float64

        # groupby float64 values results in Float64Index
        exp = Series([], dtype=np.float64,
                     index=pd.Index([], dtype=np.float64))
        tm.assert_series_equal(grouped.sum(), exp)
        tm.assert_series_equal(grouped.agg(np.sum), exp)
        tm.assert_series_equal(grouped.apply(np.sum), exp,
                               check_index_type=False)

        # DataFrame
        grouped = self.tsframe.groupby(self.tsframe['A'] * np.nan)
        exp_df = DataFrame(columns=self.tsframe.columns, dtype=float,
                           index=pd.Index([], dtype=np.float64))
        tm.assert_frame_equal(grouped.sum(), exp_df, check_names=False)
        tm.assert_frame_equal(grouped.agg(np.sum), exp_df, check_names=False)
        tm.assert_frame_equal(grouped.apply(np.sum), exp_df.iloc[:, :0],
                              check_names=False)

    def test_agg_grouping_is_list_tuple(self):
        from pandas.core.groupby import Grouping

        df = tm.makeTimeDataFrame()

        grouped = df.groupby(lambda x: x.year)
        grouper = grouped.grouper.groupings[0].grouper
        grouped.grouper.groupings[0] = Grouping(self.ts.index, list(grouper))

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        tm.assert_frame_equal(result, expected)

        grouped.grouper.groupings[0] = Grouping(self.ts.index, tuple(grouper))

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        tm.assert_frame_equal(result, expected)

    def test_agg_python_multiindex(self):
        grouped = self.mframe.groupby(['A', 'B'])

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('groupbyfunc', [
        lambda x: x.weekday(),
        [lambda x: x.month, lambda x: x.weekday()],
    ])
    def test_aggregate_str_func(self, groupbyfunc):
        grouped = self.tsframe.groupby(groupbyfunc)

        # single series
        result = grouped['A'].agg('std')
        expected = grouped['A'].std()
        tm.assert_series_equal(result, expected)

        # group frame by function name
        result = grouped.aggregate('var')
        expected = grouped.var()
        tm.assert_frame_equal(result, expected)

        # group frame by function dict
        result = grouped.agg(OrderedDict([['A', 'var'],
                                          ['B', 'std'],
                                          ['C', 'mean'],
                                          ['D', 'sem']]))
        expected = DataFrame(OrderedDict([['A', grouped['A'].var()],
                                          ['B', grouped['B'].std()],
                                          ['C', grouped['C'].mean()],
                                          ['D', grouped['D'].sem()]]))
        tm.assert_frame_equal(result, expected)

    def test_aggregate_item_by_item(self):
        df = self.df.copy()
        df['E'] = ['a'] * len(self.df)
        grouped = self.df.groupby('A')

        aggfun = lambda ser: ser.size
        result = grouped.agg(aggfun)
        foo = (self.df.A == 'foo').sum()
        bar = (self.df.A == 'bar').sum()
        K = len(result.columns)

        # GH5782
        # odd comparisons can result here, so cast to make easy
        exp = pd.Series(np.array([foo] * K), index=list('BCD'),
                        dtype=np.float64, name='foo')
        tm.assert_series_equal(result.xs('foo'), exp)

        exp = pd.Series(np.array([bar] * K), index=list('BCD'),
                        dtype=np.float64, name='bar')
        tm.assert_almost_equal(result.xs('bar'), exp)

        def aggfun(ser):
            return ser.size

        result = DataFrame().groupby(self.df.A).agg(aggfun)
        assert isinstance(result, DataFrame)
        assert len(result) == 0

    def test_wrap_agg_out(self):
        grouped = self.three_group.groupby(['A', 'B'])

        def func(ser):
            if ser.dtype == np.object:
                raise TypeError
            else:
                return ser.sum()

        result = grouped.aggregate(func)
        exp_grouped = self.three_group.loc[:, self.three_group.columns != 'C']
        expected = exp_grouped.groupby(['A', 'B']).aggregate(func)
        tm.assert_frame_equal(result, expected)

    def test_agg_multiple_functions_maintain_order(self):
        # GH #610
        funcs = [('mean', np.mean), ('max', np.max), ('min', np.min)]
        result = self.df.groupby('A')['C'].agg(funcs)
        exp_cols = Index(['mean', 'max', 'min'])

        tm.assert_index_equal(result.columns, exp_cols)

    def test_multiple_functions_tuples_and_non_tuples(self):
        # #1359
        funcs = [('foo', 'mean'), 'std']
        ex_funcs = [('foo', 'mean'), ('std', 'std')]

        result = self.df.groupby('A')['C'].agg(funcs)
        expected = self.df.groupby('A')['C'].agg(ex_funcs)
        tm.assert_frame_equal(result, expected)

        result = self.df.groupby('A').agg(funcs)
        expected = self.df.groupby('A').agg(ex_funcs)
        tm.assert_frame_equal(result, expected)

    def test_agg_multiple_functions_too_many_lambdas(self):
        grouped = self.df.groupby('A')
        funcs = ['mean', lambda x: x.mean(), lambda x: x.std()]

        msg = 'Function names must be unique, found multiple named <lambda>'
        with tm.assert_raises_regex(SpecificationError, msg):
            grouped.agg(funcs)

    def test_more_flexible_frame_multi_function(self):
        grouped = self.df.groupby('A')

        exmean = grouped.agg(OrderedDict([['C', np.mean], ['D', np.mean]]))
        exstd = grouped.agg(OrderedDict([['C', np.std], ['D', np.std]]))

        expected = concat([exmean, exstd], keys=['mean', 'std'], axis=1)
        expected = expected.swaplevel(0, 1, axis=1).sort_index(level=0, axis=1)

        d = OrderedDict([['C', [np.mean, np.std]], ['D', [np.mean, np.std]]])
        result = grouped.aggregate(d)

        tm.assert_frame_equal(result, expected)

        # be careful
        result = grouped.aggregate(OrderedDict([['C', np.mean],
                                                ['D', [np.mean, np.std]]]))
        expected = grouped.aggregate(OrderedDict([['C', np.mean],
                                                  ['D', [np.mean, np.std]]]))
        tm.assert_frame_equal(result, expected)

        def foo(x):
            return np.mean(x)

        def bar(x):
            return np.std(x, ddof=1)

        # this uses column selection & renaming
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            d = OrderedDict([['C', np.mean],
                             ['D', OrderedDict([['foo', np.mean],
                                                ['bar', np.std]])]])
            result = grouped.aggregate(d)

        d = OrderedDict([['C', [np.mean]], ['D', [foo, bar]]])
        expected = grouped.aggregate(d)

        tm.assert_frame_equal(result, expected)

    def test_multi_function_flexible_mix(self):
        # GH #1268
        grouped = self.df.groupby('A')

        # Expected
        d = OrderedDict([['C', OrderedDict([['foo', 'mean'], ['bar', 'std']])],
                         ['D', {'sum': 'sum'}]])
        # this uses column selection & renaming
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            expected = grouped.aggregate(d)

        # Test 1
        d = OrderedDict([['C', OrderedDict([['foo', 'mean'], ['bar', 'std']])],
                         ['D', 'sum']])
        # this uses column selection & renaming
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = grouped.aggregate(d)
        tm.assert_frame_equal(result, expected)

        # Test 2
        d = OrderedDict([['C', OrderedDict([['foo', 'mean'], ['bar', 'std']])],
                         ['D', ['sum']]])
        # this uses column selection & renaming
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = grouped.aggregate(d)
        tm.assert_frame_equal(result, expected)
