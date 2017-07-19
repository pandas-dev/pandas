# -*- coding: utf-8 -*-
from __future__ import print_function

import datetime

from string import ascii_lowercase

from numpy import nan

from pandas import (
    Index,
    MultiIndex,
    DataFrame,
    Series,
    NaT,
    date_range,
)
from pandas.util.testing import assert_frame_equal, assert_series_equal
from pandas.compat import (
    product as cart_product,
)
import numpy as np

import pandas.util.testing as tm
import pandas as pd

from .common import MixIn


class TestGroupBy(MixIn, tm.TestCase):
    """
    Grouping operations when keeping NaN values in the group
    with the dropna=False optione (GH 3729)

    Adapted from test_groupby.TestGroupBy
    """

    def setUp(self):
        MixIn.setUp(self)
        self.df_nan = DataFrame(
            {'A': [np.NaN, 1, np.NaN, 1, np.NaN, 1, np.NaN, np.NaN],
             'B': [1, 1, np.NaN, 3, np.NaN, np.NaN, 1, 3],
             'C': np.random.randn(8),
             'D': np.random.randn(8)})

    def test_basic_with_nan(self):
        """
        Adapted from TestGroupBy.test_basic
        """

        def checkit(dtype):
            data = Series(np.arange(9) // 3, index=np.arange(9), dtype=dtype)

            index = np.arange(9)
            np.random.shuffle(index)
            data = data.reindex(index)

            # same test as previously but with nan replacing 0 in groupby
            grouped = data.groupby(lambda x: x // 3 if x // 3 else np.nan,
                                   dropna=False)

            for k, v in grouped:
                self.assertEqual(len(v), 3)

            agged = grouped.aggregate(np.mean)
            self.assertEqual(agged[1], 1)

            assert_series_equal(agged, grouped.agg(np.mean))  # shorthand
            assert_series_equal(agged, grouped.mean())
            assert_series_equal(grouped.agg(np.sum), grouped.sum())

            data_nan = data.copy()
            data_nan[data_nan == 0] = np.nan
            value_grouped = data.groupby(data_nan, dropna=False)
            assert_series_equal(value_grouped.aggregate(np.mean), agged,
                                check_index_type=False)

            # complex agg
            agged = grouped.aggregate([np.mean, np.std])
            agged = grouped.aggregate({'one': np.mean, 'two': np.std})

            group_constants = {'nan': 10, '1.0': 20, '2.0': 30}
            agged = grouped.agg(
                lambda x: group_constants[str(x.name)] + x.mean())
            self.assertEqual(agged[1], 21)

            # corner cases
            self.assertRaises(Exception, grouped.aggregate, lambda x: x * 2)

        for dtype in ['int64', 'int32', 'float64', 'float32']:
            checkit(dtype)

    def test_first_last_nth_with_nan(self):
        """
        Adapted from TestGroupBy.test_first_last_nth_with
        """

        grouped = self.df_nan.groupby('A', dropna=False)
        first = grouped.first()

        expected = self.df_nan.ix[[1, 0], ['B', 'C', 'D']]
        expected.index = Index([1, np.NaN], name='A')

        assert_frame_equal(first, expected)

        nth = grouped.nth(0)
        assert_frame_equal(nth, expected)

        nth = grouped.nth(-1)
        expected = self.df_nan.ix[[5, 7], ['B', 'C', 'D']]
        expected.index = Index([1, np.NaN], name='A')

        assert_frame_equal(nth, expected)

        last = grouped.last()
        expected.iloc[0, 0] = 3.0  # there is a bug in first/last, cf GH 8427
        assert_frame_equal(last, expected)

        nth = grouped.nth(1)
        expected = self.df_nan.ix[[2, 3], ['B', 'C', 'D']].copy()
        expected.index = Index([np.NaN, 1], name='A')
        expected = expected.sort_index()
        assert_frame_equal(nth, expected)

    def test_nth_nan(self):
        """
        Adapted from TestGroupBy.test_nth
        """

        df = DataFrame([[np.NaN, np.nan], [np.NaN, 4], [5, 6]],
                       columns=['A', 'B'])
        g = df.groupby('A', dropna=False)

        assert_frame_equal(g.nth(0), df.iloc[[2, 0]].set_index('A'))
        assert_frame_equal(g.nth(1), df.iloc[[1]].set_index('A'))
        assert_frame_equal(g.nth(2), df.loc[[]].set_index('A'))
        assert_frame_equal(g.nth(-1), df.iloc[[2, 1]].set_index('A'))
        assert_frame_equal(g.nth(-2), df.iloc[[0]].set_index('A'))
        assert_frame_equal(g.nth(-3), df.loc[[]].set_index('A'))
        assert_series_equal(g.B.nth(0), df.set_index('A').B.iloc[[2, 0]])
        assert_series_equal(g.B.nth(1), df.set_index('A').B.iloc[[1]])
        assert_frame_equal(g[['B']].nth(0),
                           df.ix[[2, 0], ['A', 'B']].set_index('A'))

        # test dropna parameter for .nth()
        # this might be a border-line scenario: df is
        #      A    B
        # 0  NaN  NaN
        # 1  NaN  4.0
        # 2  5.0  6.0
        # line 1 is dropped, even though the NaN is in the grouping column
        exp = df.set_index('A')
        exp.iloc[1] = np.NaN
        expected = exp.iloc[[2, 1]]
        assert_frame_equal(g.nth(0, dropna='any'), expected)
        assert_frame_equal(g.nth(-1, dropna='any'), exp.iloc[[2, 1]])

        exp['B'] = np.nan
        assert_frame_equal(g.nth(7, dropna='any'), exp.iloc[[2, 1]])
        assert_frame_equal(g.nth(2, dropna='any'), exp.iloc[[2, 1]])

    def test_nth_multi_index_nan(self):
        """
        Adapted from TestGroupBy.test_nth_multi_index
        """

        df = self.three_group.replace('foo', np.NaN)
        grouped = df.groupby(['A', 'B'], dropna=False)
        result = grouped.nth(0)
        expected = grouped.first()
        assert_frame_equal(result, expected)

    def test_nth_multi_index_as_expected_nan(self):
        """
        Adapted from TestGroupBy.test_nth_multi_index_as_expected
        """

        three_group = DataFrame(
            {'A': ['foo', 'foo', 'foo', 'foo', np.NaN, np.NaN, np.NaN, np.NaN,
                   'foo', 'foo', 'foo'],
             'B': [np.NaN, np.NaN, np.NaN, 'two', np.NaN, np.NaN, np.NaN,
                   'two', 'two', 'two', np.NaN],
             'C': ['dull', 'dull', 'shiny', 'dull', 'dull', 'shiny', 'shiny',
                   'dull', 'shiny', 'shiny', 'shiny']})
        grouped = three_group.groupby(['A', 'B'], dropna=False)
        result = grouped.nth(0)
        expected = DataFrame(
            {'C': ['dull', 'dull', 'dull', 'dull']},
            index=MultiIndex(levels=[[nan, u'foo'], [nan, u'two']],
                             labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                             names=[u'A', u'B']))
        assert_frame_equal(result, expected)

    def test_groupby_series_dropna(self):
        """
        Basic test for grouping a series with the dropna=False option
        """

        s = Series([1, 2, 3, 4])
        result = s.groupby([1, np.NaN, np.NaN, 1], dropna=False).sum()

        # expected result
        expected = Series([5, 5], index=[1, np.NaN])

        assert_series_equal(result, expected)

    def test_groupby_frame_dropna(self):
        """
        Basic test for grouping a dataframe with the dropna=False option
        """

        data = [[1, 1], [np.NaN, 2], [1, 3], [np.NaN, 4]]
        df = DataFrame(data, columns=['a', 'b'])
        result = df.groupby('a', dropna=False).sum()

        # expected result
        expected = DataFrame([[4], [6]], index=[1, np.NaN], columns=['b'])
        expected.index.name = 'a'

        assert_frame_equal(result, expected)

    def test_groupby_multi_dropna(self):
        """
        Basic test for grouping a dataframe with a multiindex with the
        dropna=False option
        """

        data = [['a', 'b', 12, 12, 12],
                ['a', np.nan, 12.3, 233., 12],
                ['b', 'a', 123.23, 123, 1],
                ['a', 'b', 1, 1, 1.]]
        df = DataFrame(data, columns=['a', 'b', 'c', 'd', 'e'])
        result = df.groupby(['a', 'b'], dropna=False).sum()

        # expected result
        data = [[12.30, 233., 12.], [13., 13., 13.], [123.23, 123., 1.]]
        index = MultiIndex(levels=[['a', 'b'], [np.nan, 'a', 'b']],
                           labels=[[0, 0, 1], [0, 2, 1]],
                           names=['a', 'b'])
        expected = DataFrame(data, index=index, columns=['c', 'd', 'e'])

        assert_frame_equal(result, expected)

    def test_with_na(self):
        """
        Compare behavior with and without the dropna=False option
        """

        index = Index(np.arange(10))

        for dtype in ['float64', 'float32', 'int64', 'int32', 'int16', 'int8',
                      'datetime64[ns]']:
            values = Series(np.ones(10), index, dtype=dtype)
            labels = Series([nan, 'foo', 'bar', 'bar', nan, nan, 'bar',
                             'bar', nan, 'foo'], index=index)

            # this SHOULD be an int
            grouped = values.groupby(labels)
            agged = grouped.agg(len)
            expected = Series([4, 2], index=['bar', 'foo'])

            assert_series_equal(agged, expected, check_dtype=False)

            # self.assertTrue(issubclass(agged.dtype.type, np.integer))

            # explicity return a float from my function
            def f(x):
                return float(len(x))

            agged = grouped.agg(f)
            expected = Series([4, 2], index=['bar', 'foo'])

            assert_series_equal(agged, expected, check_dtype=False)
            # self.assertTrue(issubclass(agged.dtype.type, np.float))

            # GH 3729: keeping the NaNs

            # this SHOULD be an int
            grouped = values.groupby(labels, dropna=False)
            agged = grouped.agg(len)
            expected = Series([4, 4, 2], index=[nan, 'bar', 'foo'])

            assert_series_equal(agged, expected, check_dtype=False)
            # self.assertTrue(issubclass(agged.dtype.type, np.integer))

            # explicity return a float from my function
            def f(x):
                return float(len(x))

            agged = grouped.agg(f)
            expected = Series([4, 4, 2], index=[nan, 'bar', 'foo'])

            assert_series_equal(agged, expected, check_dtype=False)
            # self.assertTrue(issubclass(agged.dtype.type, np.float)

    def test_groupby_transform_with_nan_group(self):
        # GH 9941
        df = pd.DataFrame({'a': np.arange(10, dtype='int64'),
                           'b': [1, 1, 2, 3, np.nan, 4, 4, 5, 5, 5]})
        result = df.groupby(df.b)['a'].transform(max)
        expected = pd.Series([1., 1., 2., 3., np.nan, 6., 6., 9., 9., 9.],
                             name='a', dtype='float64')
        assert_series_equal(result, expected)

        result = df.groupby(df.b, dropna=False)['a'].transform(max)
        expected = pd.Series([1, 1, 2, 3, 4, 6, 6, 9, 9, 9],
                             name='a', dtype='int64')
        assert_series_equal(result, expected)

    def test__cython_agg_general_nan(self):
        """
        Adapted from TestGroupBy.test__cython_agg_general
        """

        ops = [('mean', np.mean),
               ('median', np.median),
               ('var', np.var),
               ('add', np.sum),
               ('prod', np.prod),
               ('min', np.min),
               ('max', np.max),
               ('first', lambda x: x.iloc[0]),
               ('last', lambda x: x.iloc[-1]), ]
        df = DataFrame(np.random.randn(1000))
        labels = np.random.randint(0, 50, size=1000).astype(float)
        labels[labels == 0] = np.nan

        for op, targop in ops:
            result = df.groupby(labels, dropna=False)._cython_agg_general(op)
            expected = df.groupby(labels, dropna=False).agg(targop)
            try:
                tm.assert_frame_equal(result, expected)
            except BaseException as exc:
                exc.args += ('operation: %s' % op, )
                raise

    def test_ops_general_nan(self):
        """
        Adapted from TestGroupBy.test_ops_general
        """

        ops = [('mean', np.mean),
               ('median', np.median),
               ('std', np.std),
               ('var', np.var),
               ('sum', np.sum),
               ('prod', np.prod),
               ('min', np.min),
               ('max', np.max),
               ('first', lambda x: x.iloc[0]),
               ('last', lambda x: x.iloc[-1]),
               ('count', np.size), ]
        try:
            from scipy.stats import sem
        except ImportError:
            pass
        else:
            ops.append(('sem', sem))
        df = DataFrame(np.random.randn(1000))
        labels = np.random.randint(0, 50, size=1000).astype(float)
        labels[labels == 0] = np.nan

        for op, targop in ops:
            result = getattr(df.groupby(labels, dropna=False), op)()
            result = result.astype(float)
            expected = df.groupby(labels, dropna=False).agg(targop)
            try:
                tm.assert_frame_equal(result, expected)
            except BaseException as exc:
                exc.args += ('operation: %s' % op, )
                raise

    def test_datetime_timedelta_nat(self):
        now = datetime.datetime.now()
        df = DataFrame([
            [now, now + datetime.timedelta(days=1)],
            [now, now + datetime.timedelta(days=2)],
            [NaT, now + datetime.timedelta(days=1)],
            [NaT, NaT],
        ], columns=['a', 'b'])
        expected = DataFrame([
            [now + datetime.timedelta(days=1)],
            [now + datetime.timedelta(days=1)],
        ], index=[now, NaT], columns=['b'])
        expected.index.name = 'a'

        # default behavior
        result = df.groupby('a').first()
        tm.assert_frame_equal(
            result.sort_index(),
            expected.iloc[:1],
        )

        # keeping NaT
        result = df.groupby('a', dropna=False).first()
        tm.assert_frame_equal(
            result.sort_index(),
            expected,
        )

        # now we test timedelta values
        df = df.applymap(lambda _: _ - now)  # workaround GH 8554
        expected = expected.applymap(lambda _: _ - now)
        expected.index = expected.index - now

        # default behavior
        result = df.groupby('a').first()
        tm.assert_frame_equal(
            result.sort_index(),
            expected.iloc[:1],
        )

        # keeping NaT
        result = df.groupby('a', dropna=False).first()
        tm.assert_frame_equal(
            result.sort_index(),
            expected,
        )

    def test_strings(self):
        df = DataFrame([
            ['a', 1],
            ['a', 2],
            [None, 1],
            [nan, 2],
            [nan, 3],
            [None, 4],
        ], columns=['x', 'y'])

        # default behavior: nans are dropped
        expected = Series([1], index=['a'], name='y')
        expected.index.name = 'x'

        result = df.groupby('x').y.first()
        tm.assert_series_equal(
            expected.sort_index(),
            result.sort_index(),
        )

        result = df.groupby(df.x).y.first()
        tm.assert_series_equal(
            expected.sort_index(),
            result.sort_index(),
        )

        result = df.groupby(Index(df.x)).y.first()
        tm.assert_series_equal(
            expected.sort_index(),
            result.sort_index(),
        )

        # not dropping the NaNs
        expected = Series([1, 1], index=['a', nan], name='y')
        expected.index.name = 'x'

        result = df.groupby('x', dropna=False).y.first()
        tm.assert_series_equal(
            expected.sort_index(),
            result.sort_index(),
        )

        result = df.groupby(df.x, dropna=False).y.first()
        tm.assert_series_equal(
            expected.sort_index(),
            result.sort_index(),
        )

        result = df.groupby(Index(df.x), dropna=False).y.first()
        tm.assert_series_equal(
            expected.sort_index(),
            result.sort_index(),
        )

    def test_series_groupby_nunique(self):

        def check_nunique(df, keys, as_index=True):
            choices = cart_product((False, True), repeat=3)
            for sort, nunique_dropna, groupby_dropna in choices:
                gr = df.groupby(
                    keys,
                    as_index=as_index,
                    sort=sort,
                    dropna=groupby_dropna,
                )
                left = gr['julie'].nunique(
                    dropna=nunique_dropna,
                )

                gr = df.groupby(
                    keys,
                    as_index=as_index,
                    sort=sort,
                    dropna=groupby_dropna,
                )
                right = gr['julie'].apply(
                    Series.nunique,
                    dropna=nunique_dropna,
                )
                if not as_index:
                    right = right.reset_index(drop=True)

                assert_series_equal(
                    left,
                    right,
                    check_names=False,
                )

        days = date_range('2015-08-23', periods=10)

        for n, m in cart_product(10 ** np.arange(2, 6), (10, 100, 1000)):
            frame = DataFrame({
                'jim': np.random.choice(
                    list(ascii_lowercase), n),
                'joe': np.random.choice(days, n),
                'julie': np.random.randint(0, m, n)
            })

            check_nunique(frame, ['jim'])
            check_nunique(frame, ['jim', 'joe'])

            frame.loc[1::17, 'jim'] = None
            frame.loc[3::37, 'joe'] = None
            frame.loc[7::19, 'julie'] = None
            frame.loc[8::19, 'julie'] = None
            frame.loc[9::19, 'julie'] = None

            check_nunique(frame, ['jim'])
            check_nunique(frame, ['jim', 'joe'])
            check_nunique(frame, ['jim'], as_index=False)
            check_nunique(frame, ['jim', 'joe'], as_index=False)
