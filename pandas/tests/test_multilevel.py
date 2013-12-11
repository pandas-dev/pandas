# pylint: disable-msg=W0612,E1101,W0141
import nose

from numpy.random import randn
import numpy as np

from pandas.core.index import Index, MultiIndex
from pandas import Panel, DataFrame, Series, notnull, isnull

from pandas.util.testing import (assert_almost_equal,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 assertRaisesRegexp)
import pandas.core.common as com
import pandas.util.testing as tm
from pandas.compat import (range, lrange, StringIO, lzip, u, cPickle,
                                product as cart_product, zip)
import pandas as pd

import pandas.index as _index


class TestMultiLevel(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        import warnings
        warnings.filterwarnings(action='ignore', category=FutureWarning)

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        self.frame = DataFrame(np.random.randn(10, 3), index=index,
                               columns=Index(['A', 'B', 'C'], name='exp'))

        self.single_level = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux']],
                                       labels=[[0, 1, 2, 3]],
                                       names=['first'])

        # create test series object
        arrays = [['bar', 'bar', 'baz', 'baz', 'qux', 'qux', 'foo', 'foo'],
                  ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        tuples = lzip(*arrays)
        index = MultiIndex.from_tuples(tuples)
        s = Series(randn(8), index=index)
        s[3] = np.NaN
        self.series = s

        tm.N = 100
        self.tdf = tm.makeTimeDataFrame()
        self.ymd = self.tdf.groupby([lambda x: x.year, lambda x: x.month,
                                     lambda x: x.day]).sum()

        # use Int64Index, to make sure things work
        self.ymd.index.set_levels([lev.astype('i8')
                                 for lev in self.ymd.index.levels],
                                 inplace=True)
        self.ymd.index.set_names(['year', 'month', 'day'],
                                 inplace=True)

    def test_append(self):
        a, b = self.frame[:5], self.frame[5:]

        result = a.append(b)
        tm.assert_frame_equal(result, self.frame)

        result = a['A'].append(b['A'])
        tm.assert_series_equal(result, self.frame['A'])

    def test_dataframe_constructor(self):
        multi = DataFrame(np.random.randn(4, 4),
                          index=[np.array(['a', 'a', 'b', 'b']),
                                 np.array(['x', 'y', 'x', 'y'])])
        tm.assert_isinstance(multi.index, MultiIndex)
        self.assert_(not isinstance(multi.columns, MultiIndex))

        multi = DataFrame(np.random.randn(4, 4),
                          columns=[['a', 'a', 'b', 'b'],
                                   ['x', 'y', 'x', 'y']])
        tm.assert_isinstance(multi.columns, MultiIndex)

    def test_series_constructor(self):
        multi = Series(1., index=[np.array(['a', 'a', 'b', 'b']),
                                  np.array(['x', 'y', 'x', 'y'])])
        tm.assert_isinstance(multi.index, MultiIndex)

        multi = Series(1., index=[['a', 'a', 'b', 'b'],
                                  ['x', 'y', 'x', 'y']])
        tm.assert_isinstance(multi.index, MultiIndex)

        multi = Series(lrange(4), index=[['a', 'a', 'b', 'b'],
                                        ['x', 'y', 'x', 'y']])
        tm.assert_isinstance(multi.index, MultiIndex)

    def test_reindex_level(self):
        # axis=0
        month_sums = self.ymd.sum(level='month')
        result = month_sums.reindex(self.ymd.index, level=1)
        expected = self.ymd.groupby(level='month').transform(np.sum)

        assert_frame_equal(result, expected)

        # Series
        result = month_sums['A'].reindex(self.ymd.index, level=1)
        expected = self.ymd['A'].groupby(level='month').transform(np.sum)
        assert_series_equal(result, expected)

        # axis=1
        month_sums = self.ymd.T.sum(axis=1, level='month')
        result = month_sums.reindex(columns=self.ymd.index, level=1)
        expected = self.ymd.groupby(level='month').transform(np.sum).T
        assert_frame_equal(result, expected)

    def test_binops_level(self):
        def _check_op(opname):
            op = getattr(DataFrame, opname)
            month_sums = self.ymd.sum(level='month')
            result = op(self.ymd, month_sums, level='month')

            broadcasted = self.ymd.groupby(level='month').transform(np.sum)
            expected = op(self.ymd, broadcasted)
            assert_frame_equal(result, expected)

            # Series
            op = getattr(Series, opname)
            result = op(self.ymd['A'], month_sums['A'], level='month')
            broadcasted = self.ymd['A'].groupby(
                level='month').transform(np.sum)
            expected = op(self.ymd['A'], broadcasted)
            assert_series_equal(result, expected)

        _check_op('sub')
        _check_op('add')
        _check_op('mul')
        _check_op('div')

    def test_pickle(self):

        def _test_roundtrip(frame):
            pickled = cPickle.dumps(frame)
            unpickled = cPickle.loads(pickled)
            assert_frame_equal(frame, unpickled)

        _test_roundtrip(self.frame)
        _test_roundtrip(self.frame.T)
        _test_roundtrip(self.ymd)
        _test_roundtrip(self.ymd.T)

    def test_reindex(self):
        reindexed = self.frame.ix[[('foo', 'one'), ('bar', 'one')]]
        expected = self.frame.ix[[0, 3]]
        assert_frame_equal(reindexed, expected)

    def test_reindex_preserve_levels(self):
        new_index = self.ymd.index[::10]
        chunk = self.ymd.reindex(new_index)
        self.assert_(chunk.index is new_index)

        chunk = self.ymd.ix[new_index]
        self.assert_(chunk.index is new_index)

        ymdT = self.ymd.T
        chunk = ymdT.reindex(columns=new_index)
        self.assert_(chunk.columns is new_index)

        chunk = ymdT.ix[:, new_index]
        self.assert_(chunk.columns is new_index)

    def test_sort_index_preserve_levels(self):
        result = self.frame.sort_index()
        self.assertEquals(result.index.names, self.frame.index.names)

    def test_repr_to_string(self):
        repr(self.frame)
        repr(self.ymd)
        repr(self.frame.T)
        repr(self.ymd.T)

        buf = StringIO()
        self.frame.to_string(buf=buf)
        self.ymd.to_string(buf=buf)
        self.frame.T.to_string(buf=buf)
        self.ymd.T.to_string(buf=buf)

    def test_repr_name_coincide(self):
        index = MultiIndex.from_tuples([('a', 0, 'foo'), ('b', 1, 'bar')],
                                       names=['a', 'b', 'c'])

        df = DataFrame({'value': [0, 1]}, index=index)

        lines = repr(df).split('\n')
        self.assert_(lines[2].startswith('a 0 foo'))

    def test_getitem_simple(self):
        df = self.frame.T

        col = df['foo', 'one']
        assert_almost_equal(col.values, df.values[:, 0])
        self.assertRaises(KeyError, df.__getitem__, ('foo', 'four'))
        self.assertRaises(KeyError, df.__getitem__, 'foobar')

    def test_series_getitem(self):
        s = self.ymd['A']

        result = s[2000, 3]
        result2 = s.ix[2000, 3]
        expected = s.reindex(s.index[42:65])
        expected.index = expected.index.droplevel(0).droplevel(0)
        assert_series_equal(result, expected)

        result = s[2000, 3, 10]
        expected = s[49]
        self.assertEquals(result, expected)

        # fancy
        result = s.ix[[(2000, 3, 10), (2000, 3, 13)]]
        expected = s.reindex(s.index[49:51])
        assert_series_equal(result, expected)

        # key error
        self.assertRaises(KeyError, s.__getitem__, (2000, 3, 4))

    def test_series_getitem_corner(self):
        s = self.ymd['A']

        # don't segfault, GH #495
        # out of bounds access
        self.assertRaises(IndexError, s.__getitem__, len(self.ymd))

        # generator
        result = s[(x > 0 for x in s)]
        expected = s[s > 0]
        assert_series_equal(result, expected)

    def test_series_setitem(self):
        s = self.ymd['A']

        s[2000, 3] = np.nan
        self.assert_(isnull(s.values[42:65]).all())
        self.assert_(notnull(s.values[:42]).all())
        self.assert_(notnull(s.values[65:]).all())

        s[2000, 3, 10] = np.nan
        self.assert_(isnull(s[49]))

    def test_series_slice_partial(self):
        pass

    def test_frame_getitem_setitem_boolean(self):
        df = self.frame.T.copy()
        values = df.values

        result = df[df > 0]
        expected = df.where(df > 0)
        assert_frame_equal(result, expected)

        df[df > 0] = 5
        values[values > 0] = 5
        assert_almost_equal(df.values, values)

        df[df == 5] = 0
        values[values == 5] = 0
        assert_almost_equal(df.values, values)

        # a df that needs alignment first
        df[df[:-1] < 0] = 2
        np.putmask(values[:-1], values[:-1] < 0, 2)
        assert_almost_equal(df.values, values)

        with assertRaisesRegexp(TypeError, 'boolean values only'):
            df[df * 0] = 2

    def test_frame_getitem_setitem_slice(self):
        # getitem
        result = self.frame.ix[:4]
        expected = self.frame[:4]
        assert_frame_equal(result, expected)

        # setitem
        cp = self.frame.copy()
        cp.ix[:4] = 0

        self.assert_((cp.values[:4] == 0).all())
        self.assert_((cp.values[4:] != 0).all())

    def test_frame_getitem_setitem_multislice(self):
        levels = [['t1', 't2'], ['a', 'b', 'c']]
        labels = [[0, 0, 0, 1, 1], [0, 1, 2, 0, 1]]
        midx = MultiIndex(labels=labels, levels=levels, names=[None, 'id'])
        df = DataFrame({'value': [1, 2, 3, 7, 8]}, index=midx)

        result = df.ix[:, 'value']
        assert_series_equal(df['value'], result)

        result = df.ix[1:3, 'value']
        assert_series_equal(df['value'][1:3], result)

        result = df.ix[:, :]
        assert_frame_equal(df, result)

        result = df
        df.ix[:, 'value'] = 10
        result['value'] = 10
        assert_frame_equal(df, result)

        df.ix[:, :] = 10
        assert_frame_equal(df, result)

    def test_frame_getitem_multicolumn_empty_level(self):
        f = DataFrame({'a': ['1', '2', '3'],
                       'b': ['2', '3', '4']})
        f.columns = [['level1 item1', 'level1 item2'],
                     ['', 'level2 item2'],
                     ['level3 item1', 'level3 item2']]

        result = f['level1 item1']
        expected = DataFrame([['1'], ['2'], ['3']], index=f.index,
                             columns=['level3 item1'])
        assert_frame_equal(result, expected)

    def test_frame_setitem_multi_column(self):
        df = DataFrame(randn(10, 4), columns=[['a', 'a', 'b', 'b'],
                                              [0, 1, 0, 1]])

        cp = df.copy()
        cp['a'] = cp['b']
        assert_frame_equal(cp['a'], cp['b'])

        # set with ndarray
        cp = df.copy()
        cp['a'] = cp['b'].values
        assert_frame_equal(cp['a'], cp['b'])

        #----------------------------------------
        # #1803
        columns = MultiIndex.from_tuples([('A', '1'), ('A', '2'), ('B', '1')])
        df = DataFrame(index=[1, 3, 5], columns=columns)

        # Works, but adds a column instead of updating the two existing ones
        df['A'] = 0.0  # Doesn't work
        self.assertTrue((df['A'].values == 0).all())

        # it broadcasts
        df['B', '1'] = [1, 2, 3]
        df['A'] = df['B', '1']
        assert_series_equal(df['A', '1'], df['B', '1'])
        assert_series_equal(df['A', '2'], df['B', '1'])

    def test_getitem_tuple_plus_slice(self):
        # GH #671
        df = DataFrame({'a': lrange(10),
                        'b': lrange(10),
                        'c': np.random.randn(10),
                        'd': np.random.randn(10)})

        idf = df.set_index(['a', 'b'])

        result = idf.ix[(0, 0), :]
        expected = idf.ix[0, 0]
        expected2 = idf.xs((0, 0))

        assert_series_equal(result, expected)
        assert_series_equal(result, expected2)

    def test_getitem_setitem_tuple_plus_columns(self):
        # GH #1013

        df = self.ymd[:5]

        result = df.ix[(2000, 1, 6), ['A', 'B', 'C']]
        expected = df.ix[2000, 1, 6][['A', 'B', 'C']]
        assert_series_equal(result, expected)

    def test_getitem_multilevel_index_tuple_unsorted(self):
        index_columns = list("abc")
        df = DataFrame([[0, 1, 0, "x"], [0, 0, 1, "y"]],
                       columns=index_columns + ["data"])
        df = df.set_index(index_columns)
        query_index = df.index[:1]
        rs = df.ix[query_index, "data"]
        xp = Series(['x'], index=MultiIndex.from_tuples([(0, 1, 0)]))
        assert_series_equal(rs, xp)

    def test_xs(self):
        xs = self.frame.xs(('bar', 'two'))
        xs2 = self.frame.ix[('bar', 'two')]

        assert_series_equal(xs, xs2)
        assert_almost_equal(xs.values, self.frame.values[4])

    def test_xs_partial(self):
        result = self.frame.xs('foo')
        result2 = self.frame.ix['foo']
        expected = self.frame.T['foo'].T
        assert_frame_equal(result, expected)
        assert_frame_equal(result, result2)

        result = self.ymd.xs((2000, 4))
        expected = self.ymd.ix[2000, 4]
        assert_frame_equal(result, expected)

        # ex from #1796
        index = MultiIndex(levels=[['foo', 'bar'], ['one', 'two'], [-1, 1]],
                           labels=[[0, 0, 0, 0, 1, 1, 1, 1],
                                   [0, 0, 1, 1, 0, 0, 1, 1],
                                   [0, 1, 0, 1, 0, 1, 0, 1]])
        df = DataFrame(np.random.randn(8, 4), index=index,
                       columns=list('abcd'))

        result = df.xs(['foo', 'one'])
        expected = df.ix['foo', 'one']
        assert_frame_equal(result, expected)

    def test_xs_level(self):
        result = self.frame.xs('two', level='second')
        expected = self.frame[self.frame.index.get_level_values(1) == 'two']
        expected.index = expected.index.droplevel(1)

        assert_frame_equal(result, expected)

        index = MultiIndex.from_tuples([('x', 'y', 'z'), ('a', 'b', 'c'),
                                        ('p', 'q', 'r')])
        df = DataFrame(np.random.randn(3, 5), index=index)
        result = df.xs('c', level=2)
        expected = df[1:2]
        expected.index = expected.index.droplevel(2)
        assert_frame_equal(result, expected)
        # can't produce a view of a multiindex with a level without copying
        with assertRaisesRegexp(ValueError, 'Cannot retrieve view'):
            self.frame.xs('two', level='second', copy=False)

    def test_xs_level_multiple(self):
        from pandas import read_table
        text = """                      A       B       C       D        E
one two three   four
a   b   10.0032 5    -0.5109 -2.3358 -0.4645  0.05076  0.3640
a   q   20      4     0.4473  1.4152  0.2834  1.00661  0.1744
x   q   30      3    -0.6662 -0.5243 -0.3580  0.89145  2.5838"""

        df = read_table(StringIO(text), sep='\s+')

        result = df.xs(('a', 4), level=['one', 'four'])
        expected = df.xs('a').xs(4, level='four')
        assert_frame_equal(result, expected)
        with assertRaisesRegexp(ValueError, 'Cannot retrieve view'):
            df.xs(('a', 4), level=['one', 'four'], copy=False)

        # GH2107
        dates = lrange(20111201, 20111205)
        ids = 'abcde'
        idx = MultiIndex.from_tuples([x for x in cart_product(dates, ids)])
        idx.names = ['date', 'secid']
        df = DataFrame(np.random.randn(len(idx), 3), idx, ['X', 'Y', 'Z'])
        rs = df.xs(20111201, level='date')
        xp = df.ix[20111201, :]
        assert_frame_equal(rs, xp)

    def test_xs_level0(self):
        from pandas import read_table
        text = """                      A       B       C       D        E
one two three   four
a   b   10.0032 5    -0.5109 -2.3358 -0.4645  0.05076  0.3640
a   q   20      4     0.4473  1.4152  0.2834  1.00661  0.1744
x   q   30      3    -0.6662 -0.5243 -0.3580  0.89145  2.5838"""

        df = read_table(StringIO(text), sep='\s+')

        result = df.xs('a', level=0)
        expected = df.xs('a')
        self.assertEqual(len(result), 2)
        assert_frame_equal(result, expected)

    def test_xs_level_series(self):
        s = self.frame['A']
        result = s[:, 'two']
        expected = self.frame.xs('two', level=1)['A']
        assert_series_equal(result, expected)

        s = self.ymd['A']
        result = s[2000, 5]
        expected = self.ymd.ix[2000, 5]['A']
        assert_series_equal(result, expected)

        # not implementing this for now

        self.assertRaises(TypeError, s.__getitem__, (2000, slice(3, 4)))

        # result = s[2000, 3:4]
        # lv =s.index.get_level_values(1)
        # expected = s[(lv == 3) | (lv == 4)]
        # expected.index = expected.index.droplevel(0)
        # assert_series_equal(result, expected)

        # can do this though

    def test_get_loc_single_level(self):
        s = Series(np.random.randn(len(self.single_level)),
                   index=self.single_level)
        for k in self.single_level.values:
            s[k]

    def test_getitem_toplevel(self):
        df = self.frame.T

        result = df['foo']
        expected = df.reindex(columns=df.columns[:3])
        expected.columns = expected.columns.droplevel(0)
        assert_frame_equal(result, expected)

        result = df['bar']
        result2 = df.ix[:, 'bar']

        expected = df.reindex(columns=df.columns[3:5])
        expected.columns = expected.columns.droplevel(0)
        assert_frame_equal(result, expected)
        assert_frame_equal(result, result2)

    def test_getitem_setitem_slice_integers(self):
        index = MultiIndex(levels=[[0, 1, 2], [0, 2]],
                           labels=[[0, 0, 1, 1, 2, 2],
                                   [0, 1, 0, 1, 0, 1]])

        frame = DataFrame(np.random.randn(len(index), 4), index=index,
                          columns=['a', 'b', 'c', 'd'])
        res = frame.ix[1:2]
        exp = frame.reindex(frame.index[2:])
        assert_frame_equal(res, exp)

        frame.ix[1:2] = 7
        self.assert_((frame.ix[1:2] == 7).values.all())

        series = Series(np.random.randn(len(index)), index=index)

        res = series.ix[1:2]
        exp = series.reindex(series.index[2:])
        assert_series_equal(res, exp)

        series.ix[1:2] = 7
        self.assert_((series.ix[1:2] == 7).values.all())

    def test_getitem_int(self):
        levels = [[0, 1], [0, 1, 2]]
        labels = [[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]]
        index = MultiIndex(levels=levels, labels=labels)

        frame = DataFrame(np.random.randn(6, 2), index=index)

        result = frame.ix[1]
        expected = frame[-3:]
        expected.index = expected.index.droplevel(0)
        assert_frame_equal(result, expected)

        # raises exception
        self.assertRaises(KeyError, frame.ix.__getitem__, 3)

        # however this will work
        result = self.frame.ix[2]
        expected = self.frame.xs(self.frame.index[2])
        assert_series_equal(result, expected)

    def test_getitem_partial(self):
        ymd = self.ymd.T
        result = ymd[2000, 2]

        expected = ymd.reindex(columns=ymd.columns[ymd.columns.labels[1] == 1])
        expected.columns = expected.columns.droplevel(0).droplevel(0)
        assert_frame_equal(result, expected)

    def test_getitem_slice_not_sorted(self):
        df = self.frame.sortlevel(1).T

        # buglet with int typechecking
        result = df.ix[:, :np.int32(3)]
        expected = df.reindex(columns=df.columns[:3])
        assert_frame_equal(result, expected)

    def test_setitem_change_dtype(self):
        dft = self.frame.T
        s = dft['foo', 'two']
        dft['foo', 'two'] = s > s.median()
        assert_series_equal(dft['foo', 'two'], s > s.median())
        tm.assert_isinstance(dft._data.blocks[1].items, MultiIndex)

        reindexed = dft.reindex(columns=[('foo', 'two')])
        assert_series_equal(reindexed['foo', 'two'], s > s.median())

    def test_frame_setitem_ix(self):
        self.frame.ix[('bar', 'two'), 'B'] = 5
        self.assertEquals(self.frame.ix[('bar', 'two'), 'B'], 5)

        # with integer labels
        df = self.frame.copy()
        df.columns = lrange(3)
        df.ix[('bar', 'two'), 1] = 7
        self.assertEquals(df.ix[('bar', 'two'), 1], 7)

    def test_fancy_slice_partial(self):
        result = self.frame.ix['bar':'baz']
        expected = self.frame[3:7]
        assert_frame_equal(result, expected)

        result = self.ymd.ix[(2000, 2):(2000, 4)]
        lev = self.ymd.index.labels[1]
        expected = self.ymd[(lev >= 1) & (lev <= 3)]
        assert_frame_equal(result, expected)

    def test_getitem_partial_column_select(self):
        idx = MultiIndex(labels=[[0, 0, 0], [0, 1, 1], [1, 0, 1]],
                         levels=[['a', 'b'], ['x', 'y'], ['p', 'q']])
        df = DataFrame(np.random.rand(3, 2), index=idx)

        result = df.ix[('a', 'y'), :]
        expected = df.ix[('a', 'y')]
        assert_frame_equal(result, expected)

        result = df.ix[('a', 'y'), [1, 0]]
        expected = df.ix[('a', 'y')][[1, 0]]
        assert_frame_equal(result, expected)

        self.assertRaises(KeyError, df.ix.__getitem__,
                          (('a', 'foo'), slice(None, None)))

    def test_sortlevel(self):
        df = self.frame.copy()
        df.index = np.arange(len(df))
        assertRaisesRegexp(TypeError, 'hierarchical index', df.sortlevel, 0)

        # axis=1

        # series
        a_sorted = self.frame['A'].sortlevel(0)
        with assertRaisesRegexp(TypeError, 'hierarchical index'):
            self.frame.reset_index()['A'].sortlevel()

        # preserve names
        self.assertEquals(a_sorted.index.names, self.frame.index.names)

        # inplace
        rs = self.frame.copy()
        rs.sortlevel(0, inplace=True)
        assert_frame_equal(rs, self.frame.sortlevel(0))

    def test_sortlevel_large_cardinality(self):

        # #2684 (int64)
        index = MultiIndex.from_arrays([np.arange(4000)]*3)
        df = DataFrame(np.random.randn(4000), index=index, dtype = np.int64)

        # it works!
        result = df.sortlevel(0)
        self.assertTrue(result.index.lexsort_depth == 3)

        # #2684 (int32)
        index = MultiIndex.from_arrays([np.arange(4000)]*3)
        df = DataFrame(np.random.randn(4000), index=index, dtype = np.int32)

        # it works!
        result = df.sortlevel(0)
        self.assert_((result.dtypes.values == df.dtypes.values).all() == True)
        self.assertTrue(result.index.lexsort_depth == 3)

    def test_delevel_infer_dtype(self):
        tuples = [tuple for tuple in cart_product(['foo', 'bar'],
                                                  [10, 20], [1.0, 1.1])]
        index = MultiIndex.from_tuples(tuples,
                                       names=['prm0', 'prm1', 'prm2'])
        df = DataFrame(np.random.randn(8, 3), columns=['A', 'B', 'C'],
                       index=index)
        deleveled = df.reset_index()
        self.assert_(com.is_integer_dtype(deleveled['prm1']))
        self.assert_(com.is_float_dtype(deleveled['prm2']))

    def test_reset_index_with_drop(self):
        deleveled = self.ymd.reset_index(drop=True)
        self.assertEquals(len(deleveled.columns), len(self.ymd.columns))

        deleveled = self.series.reset_index()
        tm.assert_isinstance(deleveled, DataFrame)
        self.assert_(
            len(deleveled.columns) == len(self.series.index.levels) + 1)

        deleveled = self.series.reset_index(drop=True)
        tm.assert_isinstance(deleveled, Series)

    def test_sortlevel_by_name(self):
        self.frame.index.names = ['first', 'second']
        result = self.frame.sortlevel(level='second')
        expected = self.frame.sortlevel(level=1)
        assert_frame_equal(result, expected)

    def test_sortlevel_mixed(self):
        sorted_before = self.frame.sortlevel(1)

        df = self.frame.copy()
        df['foo'] = 'bar'
        sorted_after = df.sortlevel(1)
        assert_frame_equal(sorted_before, sorted_after.drop(['foo'], axis=1))

        dft = self.frame.T
        sorted_before = dft.sortlevel(1, axis=1)
        dft['foo', 'three'] = 'bar'

        sorted_after = dft.sortlevel(1, axis=1)
        assert_frame_equal(sorted_before.drop([('foo', 'three')], axis=1),
                           sorted_after.drop([('foo', 'three')], axis=1))

    def test_count_level(self):
        def _check_counts(frame, axis=0):
            index = frame._get_axis(axis)
            for i in range(index.nlevels):
                result = frame.count(axis=axis, level=i)
                expected = frame.groupby(axis=axis, level=i).count(axis=axis)
                expected = expected.reindex_like(result).astype('i8')
                assert_frame_equal(result, expected)

        self.frame.ix[1, [1, 2]] = np.nan
        self.frame.ix[7, [0, 1]] = np.nan
        self.ymd.ix[1, [1, 2]] = np.nan
        self.ymd.ix[7, [0, 1]] = np.nan

        _check_counts(self.frame)
        _check_counts(self.ymd)
        _check_counts(self.frame.T, axis=1)
        _check_counts(self.ymd.T, axis=1)

        # can't call with level on regular DataFrame
        df = tm.makeTimeDataFrame()
        assertRaisesRegexp(TypeError, 'hierarchical', df.count, level=0)

        self.frame['D'] = 'foo'
        result = self.frame.count(level=0, numeric_only=True)
        assert_almost_equal(result.columns, ['A', 'B', 'C'])

    def test_count_level_series(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz'],
                                   ['one', 'two', 'three', 'four']],
                           labels=[[0, 0, 0, 2, 2],
                                   [2, 0, 1, 1, 2]])

        s = Series(np.random.randn(len(index)), index=index)

        result = s.count(level=0)
        expected = s.groupby(level=0).count()
        assert_series_equal(result.astype('f8'),
                            expected.reindex(result.index).fillna(0))

        result = s.count(level=1)
        expected = s.groupby(level=1).count()
        assert_series_equal(result.astype('f8'),
                            expected.reindex(result.index).fillna(0))

    def test_count_level_corner(self):
        s = self.frame['A'][:0]
        result = s.count(level=0)
        expected = Series(0, index=s.index.levels[0])
        assert_series_equal(result, expected)

        df = self.frame[:0]
        result = df.count(level=0)
        expected = DataFrame({}, index=s.index.levels[0],
                             columns=df.columns).fillna(0).astype(np.int64)
        assert_frame_equal(result, expected)

    def test_unstack(self):
        # just check that it works for now
        unstacked = self.ymd.unstack()
        unstacked2 = unstacked.unstack()

        # test that ints work
        unstacked = self.ymd.astype(int).unstack()

        # test that int32 work
        unstacked = self.ymd.astype(np.int32).unstack()

    def test_unstack_multiple_no_empty_columns(self):
        index = MultiIndex.from_tuples([(0, 'foo', 0), (0, 'bar', 0),
                                        (1, 'baz', 1), (1, 'qux', 1)])

        s = Series(np.random.randn(4), index=index)

        unstacked = s.unstack([1, 2])
        expected = unstacked.dropna(axis=1, how='all')
        assert_frame_equal(unstacked, expected)

    def test_stack(self):
        # regular roundtrip
        unstacked = self.ymd.unstack()
        restacked = unstacked.stack()
        assert_frame_equal(restacked, self.ymd)

        unlexsorted = self.ymd.sortlevel(2)

        unstacked = unlexsorted.unstack(2)
        restacked = unstacked.stack()
        assert_frame_equal(restacked.sortlevel(0), self.ymd)

        unlexsorted = unlexsorted[::-1]
        unstacked = unlexsorted.unstack(1)
        restacked = unstacked.stack().swaplevel(1, 2)
        assert_frame_equal(restacked.sortlevel(0), self.ymd)

        unlexsorted = unlexsorted.swaplevel(0, 1)
        unstacked = unlexsorted.unstack(0).swaplevel(0, 1, axis=1)
        restacked = unstacked.stack(0).swaplevel(1, 2)
        assert_frame_equal(restacked.sortlevel(0), self.ymd)

        # columns unsorted
        unstacked = self.ymd.unstack()
        unstacked = unstacked.sort(axis=1, ascending=False)
        restacked = unstacked.stack()
        assert_frame_equal(restacked, self.ymd)

        # more than 2 levels in the columns
        unstacked = self.ymd.unstack(1).unstack(1)

        result = unstacked.stack(1)
        expected = self.ymd.unstack()
        assert_frame_equal(result, expected)

        result = unstacked.stack(2)
        expected = self.ymd.unstack(1)
        assert_frame_equal(result, expected)

        result = unstacked.stack(0)
        expected = self.ymd.stack().unstack(1).unstack(1)
        assert_frame_equal(result, expected)

        # not all levels present in each echelon
        unstacked = self.ymd.unstack(2).ix[:, ::3]
        stacked = unstacked.stack().stack()
        ymd_stacked = self.ymd.stack()
        assert_series_equal(stacked, ymd_stacked.reindex(stacked.index))

        # stack with negative number
        result = self.ymd.unstack(0).stack(-2)
        expected = self.ymd.unstack(0).stack(0)

    def test_unstack_odd_failure(self):
        data = """day,time,smoker,sum,len
Fri,Dinner,No,8.25,3.
Fri,Dinner,Yes,27.03,9
Fri,Lunch,No,3.0,1
Fri,Lunch,Yes,13.68,6
Sat,Dinner,No,139.63,45
Sat,Dinner,Yes,120.77,42
Sun,Dinner,No,180.57,57
Sun,Dinner,Yes,66.82,19
Thur,Dinner,No,3.0,1
Thur,Lunch,No,117.32,44
Thur,Lunch,Yes,51.51,17"""

        df = pd.read_csv(StringIO(data)).set_index(['day', 'time', 'smoker'])

        # it works, #2100
        result = df.unstack(2)

        recons = result.stack()
        assert_frame_equal(recons, df)

    def test_stack_mixed_dtype(self):
        df = self.frame.T
        df['foo', 'four'] = 'foo'
        df = df.sortlevel(1, axis=1)

        stacked = df.stack()
        assert_series_equal(stacked['foo'], df['foo'].stack())
        self.assert_(stacked['bar'].dtype == np.float_)

    def test_unstack_bug(self):
        df = DataFrame({'state': ['naive', 'naive', 'naive',
                                  'activ', 'activ', 'activ'],
                        'exp': ['a', 'b', 'b', 'b', 'a', 'a'],
                        'barcode': [1, 2, 3, 4, 1, 3],
                        'v': ['hi', 'hi', 'bye', 'bye', 'bye', 'peace'],
                        'extra': np.arange(6.)})

        result = df.groupby(['state', 'exp', 'barcode', 'v']).apply(len)

        unstacked = result.unstack()
        restacked = unstacked.stack()
        assert_series_equal(restacked,
                            result.reindex(restacked.index).astype(float))

    def test_stack_unstack_preserve_names(self):
        unstacked = self.frame.unstack()
        self.assertEquals(unstacked.index.name, 'first')
        self.assertEquals(unstacked.columns.names, ['exp', 'second'])

        restacked = unstacked.stack()
        self.assertEquals(restacked.index.names, self.frame.index.names)

    def test_unstack_level_name(self):
        result = self.frame.unstack('second')
        expected = self.frame.unstack(level=1)
        assert_frame_equal(result, expected)

    def test_stack_level_name(self):
        unstacked = self.frame.unstack('second')
        result = unstacked.stack('exp')
        expected = self.frame.unstack().stack(0)
        assert_frame_equal(result, expected)

        result = self.frame.stack('exp')
        expected = self.frame.stack()
        assert_series_equal(result, expected)

    def test_stack_unstack_multiple(self):
        unstacked = self.ymd.unstack(['year', 'month'])
        expected = self.ymd.unstack('year').unstack('month')
        assert_frame_equal(unstacked, expected)
        self.assertEquals(unstacked.columns.names,
                          expected.columns.names)

        # series
        s = self.ymd['A']
        s_unstacked = s.unstack(['year', 'month'])
        assert_frame_equal(s_unstacked, expected['A'])

        restacked = unstacked.stack(['year', 'month'])
        restacked = restacked.swaplevel(0, 1).swaplevel(1, 2)
        restacked = restacked.sortlevel(0)

        assert_frame_equal(restacked, self.ymd)
        self.assertEquals(restacked.index.names, self.ymd.index.names)

        # GH #451
        unstacked = self.ymd.unstack([1, 2])
        expected = self.ymd.unstack(1).unstack(1).dropna(axis=1, how='all')
        assert_frame_equal(unstacked, expected)

        unstacked = self.ymd.unstack([2, 1])
        expected = self.ymd.unstack(2).unstack(1).dropna(axis=1, how='all')
        assert_frame_equal(unstacked, expected.ix[:, unstacked.columns])

    def test_stack_multiple_bug(self):
        """ bug when some uniques are not present in the data #3170"""
        id_col = ([1] * 3) + ([2] * 3)
        name = (['a'] * 3) + (['b'] * 3)
        date = pd.to_datetime(['2013-01-03', '2013-01-04', '2013-01-05'] * 2)
        var1 = np.random.randint(0, 100, 6)
        df = DataFrame(dict(ID=id_col, NAME=name, DATE=date, VAR1=var1))

        multi = df.set_index(['DATE', 'ID'])
        multi.columns.name = 'Params'
        unst = multi.unstack('ID')
        down = unst.resample('W-THU')

        rs = down.stack('ID')
        xp = unst.ix[:, ['VAR1']].resample('W-THU').stack('ID')
        xp.columns.name = 'Params'
        assert_frame_equal(rs, xp)

    def test_stack_dropna(self):
        # GH #3997
        df = pd.DataFrame({'A': ['a1', 'a2'],
                           'B': ['b1', 'b2'],
                           'C': [1, 1]})
        df = df.set_index(['A', 'B'])

        stacked = df.unstack().stack(dropna=False)
        self.assertTrue(len(stacked) > len(stacked.dropna()))

        stacked = df.unstack().stack(dropna=True)
        assert_frame_equal(stacked, stacked.dropna())

    def test_unstack_multiple_hierarchical(self):
        df = DataFrame(index=[[0, 0, 0, 0, 1, 1, 1, 1],
                              [0, 0, 1, 1, 0, 0, 1, 1],
                              [0, 1, 0, 1, 0, 1, 0, 1]],
                       columns=[[0, 0, 1, 1], [0, 1, 0, 1]])

        df.index.names = ['a', 'b', 'c']
        df.columns.names = ['d', 'e']

        # it works!
        df.unstack(['b', 'c'])

    def test_groupby_transform(self):
        s = self.frame['A']
        grouper = s.index.get_level_values(0)

        grouped = s.groupby(grouper)

        applied = grouped.apply(lambda x: x * 2)
        expected = grouped.transform(lambda x: x * 2)
        assert_series_equal(applied.reindex(expected.index), expected)

    def test_unstack_sparse_keyspace(self):
        # memory problems with naive impl #2278
        # Generate Long File & Test Pivot
        NUM_ROWS = 1000

        df = DataFrame({'A': np.random.randint(100, size=NUM_ROWS),
                        'B': np.random.randint(300, size=NUM_ROWS),
                        'C': np.random.randint(-7, 7, size=NUM_ROWS),
                        'D': np.random.randint(-19, 19, size=NUM_ROWS),
                        'E': np.random.randint(3000, size=NUM_ROWS),
                        'F': np.random.randn(NUM_ROWS)})

        idf = df.set_index(['A', 'B', 'C', 'D', 'E'])

        # it works! is sufficient
        idf.unstack('E')

    def test_unstack_unobserved_keys(self):
        # related to #2278 refactoring
        levels = [[0, 1], [0, 1, 2, 3]]
        labels = [[0, 0, 1, 1], [0, 2, 0, 2]]

        index = MultiIndex(levels, labels)

        df = DataFrame(np.random.randn(4, 2), index=index)

        result = df.unstack()
        self.assertEquals(len(result.columns), 4)

        recons = result.stack()
        assert_frame_equal(recons, df)

    def test_groupby_corner(self):
        midx = MultiIndex(levels=[['foo'], ['bar'], ['baz']],
                          labels=[[0], [0], [0]], names=['one', 'two', 'three'])
        df = DataFrame([np.random.rand(4)], columns=['a', 'b', 'c', 'd'],
                       index=midx)
        # should work
        df.groupby(level='three')

    def test_groupby_level_no_obs(self):
        # #1697
        midx = MultiIndex.from_tuples([('f1', 's1'), ('f1', 's2'),
                                       ('f2', 's1'), ('f2', 's2'),
                                       ('f3', 's1'), ('f3', 's2')])
        df = DataFrame(
            [[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]], columns=midx)
        df1 = df.select(lambda u: u[0] in ['f2', 'f3'], axis=1)

        grouped = df1.groupby(axis=1, level=0)
        result = grouped.sum()
        self.assert_((result.columns == ['f2', 'f3']).all())

    def test_join(self):
        a = self.frame.ix[:5, ['A']]
        b = self.frame.ix[2:, ['B', 'C']]

        joined = a.join(b, how='outer').reindex(self.frame.index)
        expected = self.frame.copy()
        expected.values[np.isnan(joined.values)] = np.nan

        self.assert_(not np.isnan(joined.values).all())

        assert_frame_equal(joined, expected, check_names=False)  # TODO what should join do with names ?

    def test_swaplevel(self):
        swapped = self.frame['A'].swaplevel(0, 1)
        swapped2 = self.frame['A'].swaplevel('first', 'second')
        self.assert_(not swapped.index.equals(self.frame.index))
        assert_series_equal(swapped, swapped2)

        back = swapped.swaplevel(0, 1)
        back2 = swapped.swaplevel('second', 'first')
        self.assert_(back.index.equals(self.frame.index))
        assert_series_equal(back, back2)

        ft = self.frame.T
        swapped = ft.swaplevel('first', 'second', axis=1)
        exp = self.frame.swaplevel('first', 'second').T
        assert_frame_equal(swapped, exp)

    def test_swaplevel_panel(self):
        panel = Panel({'ItemA': self.frame,
                       'ItemB': self.frame * 2})

        result = panel.swaplevel(0, 1, axis='major')
        expected = panel.copy()
        expected.major_axis = expected.major_axis.swaplevel(0, 1)
        tm.assert_panel_equal(result, expected)

    def test_reorder_levels(self):
        result = self.ymd.reorder_levels(['month', 'day', 'year'])
        expected = self.ymd.swaplevel(0, 1).swaplevel(1, 2)
        assert_frame_equal(result, expected)

        result = self.ymd['A'].reorder_levels(['month', 'day', 'year'])
        expected = self.ymd['A'].swaplevel(0, 1).swaplevel(1, 2)
        assert_series_equal(result, expected)

        result = self.ymd.T.reorder_levels(['month', 'day', 'year'], axis=1)
        expected = self.ymd.T.swaplevel(0, 1, axis=1).swaplevel(1, 2, axis=1)
        assert_frame_equal(result, expected)

        with assertRaisesRegexp(TypeError, 'hierarchical axis'):
            self.ymd.reorder_levels([1, 2], axis=1)

        with assertRaisesRegexp(IndexError, 'Too many levels'):
            self.ymd.index.reorder_levels([1, 2, 3])

    def test_insert_index(self):
        df = self.ymd[:5].T
        df[2000, 1, 10] = df[2000, 1, 7]
        tm.assert_isinstance(df.columns, MultiIndex)
        self.assert_((df[2000, 1, 10] == df[2000, 1, 7]).all())

    def test_alignment(self):
        x = Series(data=[1, 2, 3],
                   index=MultiIndex.from_tuples([("A", 1), ("A", 2), ("B", 3)]))

        y = Series(data=[4, 5, 6],
                   index=MultiIndex.from_tuples([("Z", 1), ("Z", 2), ("B", 3)]))

        res = x - y
        exp_index = x.index.union(y.index)
        exp = x.reindex(exp_index) - y.reindex(exp_index)
        assert_series_equal(res, exp)

        # hit non-monotonic code path
        res = x[::-1] - y[::-1]
        exp_index = x.index.union(y.index)
        exp = x.reindex(exp_index) - y.reindex(exp_index)
        assert_series_equal(res, exp)

    def test_is_lexsorted(self):
        levels = [[0, 1], [0, 1, 2]]

        index = MultiIndex(levels=levels,
                           labels=[[0, 0, 0, 1, 1, 1],
                                   [0, 1, 2, 0, 1, 2]])
        self.assert_(index.is_lexsorted())

        index = MultiIndex(levels=levels,
                           labels=[[0, 0, 0, 1, 1, 1],
                                   [0, 1, 2, 0, 2, 1]])
        self.assert_(not index.is_lexsorted())

        index = MultiIndex(levels=levels,
                           labels=[[0, 0, 1, 0, 1, 1],
                                   [0, 1, 0, 2, 2, 1]])
        self.assert_(not index.is_lexsorted())
        self.assert_(index.lexsort_depth == 0)

    def test_frame_getitem_view(self):
        df = self.frame.T.copy()

        # this works because we are modifying the underlying array
        # really a no-no
        df['foo'].values[:] = 0
        self.assert_((df['foo'].values == 0).all())

        # but not if it's mixed-type
        df['foo', 'four'] = 'foo'
        df = df.sortlevel(0, axis=1)

        # this will work, but will raise/warn as its chained assignment
        def f():
            df['foo']['one'] = 2
            return df
        self.assertRaises(com.SettingWithCopyError, f)

        try:
            df = f()
        except:
            pass
        self.assert_((df['foo', 'one'] == 0).all())

    def test_frame_getitem_not_sorted(self):
        df = self.frame.T
        df['foo', 'four'] = 'foo'

        arrays = [np.array(x) for x in zip(*df.columns._tuple_index)]

        result = df['foo']
        result2 = df.ix[:, 'foo']
        expected = df.reindex(columns=df.columns[arrays[0] == 'foo'])
        expected.columns = expected.columns.droplevel(0)
        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

        df = df.T
        result = df.xs('foo')
        result2 = df.ix['foo']
        expected = df.reindex(df.index[arrays[0] == 'foo'])
        expected.index = expected.index.droplevel(0)
        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

    def test_series_getitem_not_sorted(self):
        arrays = [['bar', 'bar', 'baz', 'baz', 'qux', 'qux', 'foo', 'foo'],
                 ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        tuples = lzip(*arrays)
        index = MultiIndex.from_tuples(tuples)
        s = Series(randn(8), index=index)

        arrays = [np.array(x) for x in zip(*index._tuple_index)]

        result = s['qux']
        result2 = s.ix['qux']
        expected = s[arrays[0] == 'qux']
        expected.index = expected.index.droplevel(0)
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

    def test_count(self):
        frame = self.frame.copy()
        frame.index.names = ['a', 'b']

        result = frame.count(level='b')
        expect = self.frame.count(level=1)
        assert_frame_equal(result, expect, check_names=False)

        result = frame.count(level='a')
        expect = self.frame.count(level=0)
        assert_frame_equal(result, expect, check_names=False)

        series = self.series.copy()
        series.index.names = ['a', 'b']

        result = series.count(level='b')
        expect = self.series.count(level=1)
        assert_series_equal(result, expect)

        result = series.count(level='a')
        expect = self.series.count(level=0)
        assert_series_equal(result, expect)

        self.assertRaises(KeyError, series.count, 'x')
        self.assertRaises(KeyError, frame.count, level='x')

    AGG_FUNCTIONS = ['sum', 'prod', 'min', 'max', 'median', 'mean', 'skew',
                     'mad', 'std', 'var']

    def test_series_group_min_max(self):
        for op, level, skipna in cart_product(self.AGG_FUNCTIONS,
                                              lrange(2),
                                              [False, True]):
            grouped = self.series.groupby(level=level)
            aggf = lambda x: getattr(x, op)(skipna=skipna)
            # skipna=True
            leftside = grouped.agg(aggf)
            rightside = getattr(self.series, op)(level=level, skipna=skipna)
            assert_series_equal(leftside, rightside)

    def test_frame_group_ops(self):
        self.frame.ix[1, [1, 2]] = np.nan
        self.frame.ix[7, [0, 1]] = np.nan

        for op, level, axis, skipna in cart_product(self.AGG_FUNCTIONS,
                                                    lrange(2), lrange(2),
                                                    [False, True]):
            if axis == 0:
                frame = self.frame
            else:
                frame = self.frame.T

            grouped = frame.groupby(level=level, axis=axis)

            pieces = []

            def aggf(x):
                pieces.append(x)
                return getattr(x, op)(skipna=skipna, axis=axis)
            leftside = grouped.agg(aggf)
            rightside = getattr(frame, op)(level=level, axis=axis,
                                           skipna=skipna)

            # for good measure, groupby detail
            level_index = frame._get_axis(axis).levels[level]

            self.assert_(leftside._get_axis(axis).equals(level_index))
            self.assert_(rightside._get_axis(axis).equals(level_index))

            assert_frame_equal(leftside, rightside)

    def test_stat_op_corner(self):
        obj = Series([10.0], index=MultiIndex.from_tuples([(2, 3)]))

        result = obj.sum(level=0)
        expected = Series([10.0], index=[2])
        assert_series_equal(result, expected)

    def test_frame_any_all_group(self):
        df = DataFrame(
            {'data': [False, False, True, False, True, False, True]},
            index=[
                ['one', 'one', 'two', 'one', 'two', 'two', 'two'],
                [0, 1, 0, 2, 1, 2, 3]])

        result = df.any(level=0)
        ex = DataFrame({'data': [False, True]}, index=['one', 'two'])
        assert_frame_equal(result, ex)

        result = df.all(level=0)
        ex = DataFrame({'data': [False, False]}, index=['one', 'two'])
        assert_frame_equal(result, ex)

    def test_std_var_pass_ddof(self):
        index = MultiIndex.from_arrays([np.arange(5).repeat(10),
                                        np.tile(np.arange(10), 5)])
        df = DataFrame(np.random.randn(len(index), 5), index=index)

        for meth in ['var', 'std']:
            ddof = 4
            alt = lambda x: getattr(x, meth)(ddof=ddof)

            result = getattr(df[0], meth)(level=0, ddof=ddof)
            expected = df[0].groupby(level=0).agg(alt)
            assert_series_equal(result, expected)

            result = getattr(df, meth)(level=0, ddof=ddof)
            expected = df.groupby(level=0).agg(alt)
            assert_frame_equal(result, expected)

    def test_frame_series_agg_multiple_levels(self):
        result = self.ymd.sum(level=['year', 'month'])
        expected = self.ymd.groupby(level=['year', 'month']).sum()
        assert_frame_equal(result, expected)

        result = self.ymd['A'].sum(level=['year', 'month'])
        expected = self.ymd['A'].groupby(level=['year', 'month']).sum()
        assert_series_equal(result, expected)

    def test_groupby_multilevel(self):
        result = self.ymd.groupby(level=[0, 1]).mean()

        k1 = self.ymd.index.get_level_values(0)
        k2 = self.ymd.index.get_level_values(1)

        expected = self.ymd.groupby([k1, k2]).mean()

        assert_frame_equal(result, expected, check_names=False)  # TODO groupby with level_values drops names
        self.assertEquals(result.index.names, self.ymd.index.names[:2])

        result2 = self.ymd.groupby(level=self.ymd.index.names[:2]).mean()
        assert_frame_equal(result, result2)

    def test_groupby_multilevel_with_transform(self):
        pass

    def test_multilevel_consolidate(self):
        index = MultiIndex.from_tuples([('foo', 'one'), ('foo', 'two'),
                                        ('bar', 'one'), ('bar', 'two')])
        df = DataFrame(np.random.randn(4, 4), index=index, columns=index)
        df['Totals', ''] = df.sum(1)
        df = df.consolidate()

    def test_ix_preserve_names(self):
        result = self.ymd.ix[2000]
        result2 = self.ymd['A'].ix[2000]
        self.assertEquals(result.index.names, self.ymd.index.names[1:])
        self.assertEquals(result2.index.names, self.ymd.index.names[1:])

        result = self.ymd.ix[2000, 2]
        result2 = self.ymd['A'].ix[2000, 2]
        self.assertEquals(result.index.name, self.ymd.index.names[2])
        self.assertEquals(result2.index.name, self.ymd.index.names[2])

    def test_partial_set(self):
        # GH #397
        df = self.ymd.copy()
        exp = self.ymd.copy()
        df.ix[2000, 4] = 0
        exp.ix[2000, 4].values[:] = 0
        assert_frame_equal(df, exp)

        df['A'].ix[2000, 4] = 1
        exp['A'].ix[2000, 4].values[:] = 1
        assert_frame_equal(df, exp)

        df.ix[2000] = 5
        exp.ix[2000].values[:] = 5
        assert_frame_equal(df, exp)

        # this works...for now
        df['A'].ix[14] = 5
        self.assertEquals(df['A'][14], 5)

    def test_unstack_preserve_types(self):
        # GH #403
        self.ymd['E'] = 'foo'
        self.ymd['F'] = 2

        unstacked = self.ymd.unstack('month')
        self.assert_(unstacked['A', 1].dtype == np.float64)
        self.assert_(unstacked['E', 1].dtype == np.object_)
        self.assert_(unstacked['F', 1].dtype == np.float64)

    def test_unstack_group_index_overflow(self):
        labels = np.tile(np.arange(500), 2)
        level = np.arange(500)

        index = MultiIndex(levels=[level] * 8 + [[0, 1]],
                           labels=[labels] * 8 + [np.arange(2).repeat(500)])

        s = Series(np.arange(1000), index=index)
        result = s.unstack()
        self.assertEqual(result.shape, (500, 2))

        # test roundtrip
        stacked = result.stack()
        assert_series_equal(s,
                            stacked.reindex(s.index))

        # put it at beginning
        index = MultiIndex(levels=[[0, 1]] + [level] * 8,
                           labels=[np.arange(2).repeat(500)] + [labels] * 8)

        s = Series(np.arange(1000), index=index)
        result = s.unstack(0)
        self.assertEqual(result.shape, (500, 2))

        # put it in middle
        index = MultiIndex(levels=[level] * 4 + [[0, 1]] + [level] * 4,
                           labels=([labels] * 4 + [np.arange(2).repeat(500)]
                                   + [labels] * 4))

        s = Series(np.arange(1000), index=index)
        result = s.unstack(4)
        self.assertEqual(result.shape, (500, 2))

    def test_getitem_lowerdim_corner(self):
        self.assertRaises(KeyError, self.frame.ix.__getitem__,
                          (('bar', 'three'), 'B'))


        # in theory should be inserting in a sorted space????
        self.frame.ix[('bar','three'),'B'] = 0
        self.assert_(self.frame.sortlevel().ix[('bar','three'),'B'] == 0)

    #----------------------------------------------------------------------
    # AMBIGUOUS CASES!

    def test_partial_ix_missing(self):
        raise nose.SkipTest("skipping for now")

        result = self.ymd.ix[2000, 0]
        expected = self.ymd.ix[2000]['A']
        assert_series_equal(result, expected)

        # need to put in some work here

        # self.ymd.ix[2000, 0] = 0
        # self.assert_((self.ymd.ix[2000]['A'] == 0).all())

        # Pretty sure the second (and maybe even the first) is already wrong.
        self.assertRaises(Exception, self.ymd.ix.__getitem__, (2000, 6))
        self.assertRaises(Exception, self.ymd.ix.__getitem__, (2000, 6), 0)

    #----------------------------------------------------------------------

    def test_to_html(self):
        self.ymd.columns.name = 'foo'
        self.ymd.to_html()
        self.ymd.T.to_html()

    def test_level_with_tuples(self):
        index = MultiIndex(levels=[[('foo', 'bar', 0), ('foo', 'baz', 0),
                                    ('foo', 'qux', 0)],
                                   [0, 1]],
                           labels=[[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]])

        series = Series(np.random.randn(6), index=index)
        frame = DataFrame(np.random.randn(6, 4), index=index)

        result = series[('foo', 'bar', 0)]
        result2 = series.ix[('foo', 'bar', 0)]
        expected = series[:2]
        expected.index = expected.index.droplevel(0)
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

        self.assertRaises(KeyError, series.__getitem__, (('foo', 'bar', 0), 2))

        result = frame.ix[('foo', 'bar', 0)]
        result2 = frame.xs(('foo', 'bar', 0))
        expected = frame[:2]
        expected.index = expected.index.droplevel(0)
        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

        index = MultiIndex(levels=[[('foo', 'bar'), ('foo', 'baz'),
                                    ('foo', 'qux')],
                                   [0, 1]],
                           labels=[[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]])

        series = Series(np.random.randn(6), index=index)
        frame = DataFrame(np.random.randn(6, 4), index=index)

        result = series[('foo', 'bar')]
        result2 = series.ix[('foo', 'bar')]
        expected = series[:2]
        expected.index = expected.index.droplevel(0)
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

        result = frame.ix[('foo', 'bar')]
        result2 = frame.xs(('foo', 'bar'))
        expected = frame[:2]
        expected.index = expected.index.droplevel(0)
        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

    def test_int_series_slicing(self):
        s = self.ymd['A']
        result = s[5:]
        expected = s.reindex(s.index[5:])
        assert_series_equal(result, expected)

        exp = self.ymd['A'].copy()
        s[5:] = 0
        exp.values[5:] = 0
        self.assert_(np.array_equal(s.values, exp.values))

        result = self.ymd[5:]
        expected = self.ymd.reindex(s.index[5:])
        assert_frame_equal(result, expected)

    def test_mixed_depth_get(self):
        arrays = [['a', 'top', 'top', 'routine1', 'routine1', 'routine2'],
                  ['', 'OD', 'OD', 'result1', 'result2', 'result1'],
                  ['', 'wx', 'wy', '', '', '']]

        tuples = sorted(zip(*arrays))
        index = MultiIndex.from_tuples(tuples)
        df = DataFrame(randn(4, 6), columns=index)

        result = df['a']
        expected = df['a', '', '']
        assert_series_equal(result, expected)
        self.assertEquals(result.name, 'a')

        result = df['routine1', 'result1']
        expected = df['routine1', 'result1', '']
        assert_series_equal(result, expected)
        self.assertEquals(result.name, ('routine1', 'result1'))

    def test_mixed_depth_insert(self):
        arrays = [['a', 'top', 'top', 'routine1', 'routine1', 'routine2'],
                  ['', 'OD', 'OD', 'result1', 'result2', 'result1'],
                  ['', 'wx', 'wy', '', '', '']]

        tuples = sorted(zip(*arrays))
        index = MultiIndex.from_tuples(tuples)
        df = DataFrame(randn(4, 6), columns=index)

        result = df.copy()
        expected = df.copy()
        result['b'] = [1, 2, 3, 4]
        expected['b', '', ''] = [1, 2, 3, 4]
        assert_frame_equal(result, expected)

    def test_mixed_depth_drop(self):
        arrays = [['a', 'top', 'top', 'routine1', 'routine1', 'routine2'],
                  ['', 'OD', 'OD', 'result1', 'result2', 'result1'],
                  ['', 'wx', 'wy', '', '', '']]

        tuples = sorted(zip(*arrays))
        index = MultiIndex.from_tuples(tuples)
        df = DataFrame(randn(4, 6), columns=index)

        result = df.drop('a', axis=1)
        expected = df.drop([('a', '', '')], axis=1)
        assert_frame_equal(expected, result)

        result = df.drop(['top'], axis=1)
        expected = df.drop([('top', 'OD', 'wx')], axis=1)
        expected = expected.drop([('top', 'OD', 'wy')], axis=1)
        assert_frame_equal(expected, result)

        result = df.drop(('top', 'OD', 'wx'), axis=1)
        expected = df.drop([('top', 'OD', 'wx')], axis=1)
        assert_frame_equal(expected, result)

        expected = df.drop([('top', 'OD', 'wy')], axis=1)
        expected = df.drop('top', axis=1)

        result = df.drop('result1', level=1, axis=1)
        expected = df.drop([('routine1', 'result1', ''),
                            ('routine2', 'result1', '')], axis=1)
        assert_frame_equal(expected, result)

    def test_drop_nonunique(self):
        df = DataFrame([["x-a", "x", "a", 1.5], ["x-a", "x", "a", 1.2],
                        ["z-c", "z", "c", 3.1], ["x-a", "x", "a", 4.1],
                        ["x-b", "x", "b", 5.1], ["x-b", "x", "b", 4.1],
                        ["x-b", "x", "b", 2.2],
                        ["y-a", "y", "a", 1.2], ["z-b", "z", "b", 2.1]],
                       columns=["var1", "var2", "var3", "var4"])

        grp_size = df.groupby("var1").size()
        drop_idx = grp_size.ix[grp_size == 1]

        idf = df.set_index(["var1", "var2", "var3"])

        # it works! #2101
        result = idf.drop(drop_idx.index, level=0).reset_index()
        expected = df[-df.var1.isin(drop_idx.index)]

        result.index = expected.index

        assert_frame_equal(result, expected)

    def test_mixed_depth_pop(self):
        arrays = [['a', 'top', 'top', 'routine1', 'routine1', 'routine2'],
                  ['', 'OD', 'OD', 'result1', 'result2', 'result1'],
                  ['', 'wx', 'wy', '', '', '']]

        tuples = sorted(zip(*arrays))
        index = MultiIndex.from_tuples(tuples)
        df = DataFrame(randn(4, 6), columns=index)

        df1 = df.copy()
        df2 = df.copy()
        result = df1.pop('a')
        expected = df2.pop(('a', '', ''))
        assert_series_equal(expected, result)
        assert_frame_equal(df1, df2)
        self.assertEquals(result.name, 'a')

        expected = df1['top']
        df1 = df1.drop(['top'], axis=1)
        result = df2.pop('top')
        assert_frame_equal(expected, result)
        assert_frame_equal(df1, df2)

    def test_reindex_level_partial_selection(self):
        result = self.frame.reindex(['foo', 'qux'], level=0)
        expected = self.frame.ix[[0, 1, 2, 7, 8, 9]]
        assert_frame_equal(result, expected)

        result = self.frame.T.reindex_axis(['foo', 'qux'], axis=1, level=0)
        assert_frame_equal(result, expected.T)

        result = self.frame.ix[['foo', 'qux']]
        assert_frame_equal(result, expected)

        result = self.frame['A'].ix[['foo', 'qux']]
        assert_series_equal(result, expected['A'])

        result = self.frame.T.ix[:, ['foo', 'qux']]
        assert_frame_equal(result, expected.T)

    def test_setitem_multiple_partial(self):
        expected = self.frame.copy()
        result = self.frame.copy()
        result.ix[['foo', 'bar']] = 0
        expected.ix['foo'] = 0
        expected.ix['bar'] = 0
        assert_frame_equal(result, expected)

        expected = self.frame.copy()
        result = self.frame.copy()
        result.ix['foo':'bar'] = 0
        expected.ix['foo'] = 0
        expected.ix['bar'] = 0
        assert_frame_equal(result, expected)

        expected = self.frame['A'].copy()
        result = self.frame['A'].copy()
        result.ix[['foo', 'bar']] = 0
        expected.ix['foo'] = 0
        expected.ix['bar'] = 0
        assert_series_equal(result, expected)

        expected = self.frame['A'].copy()
        result = self.frame['A'].copy()
        result.ix['foo':'bar'] = 0
        expected.ix['foo'] = 0
        expected.ix['bar'] = 0
        assert_series_equal(result, expected)

    def test_drop_level(self):
        result = self.frame.drop(['bar', 'qux'], level='first')
        expected = self.frame.ix[[0, 1, 2, 5, 6]]
        assert_frame_equal(result, expected)

        result = self.frame.drop(['two'], level='second')
        expected = self.frame.ix[[0, 2, 3, 6, 7, 9]]
        assert_frame_equal(result, expected)

        result = self.frame.T.drop(['bar', 'qux'], axis=1, level='first')
        expected = self.frame.ix[[0, 1, 2, 5, 6]].T
        assert_frame_equal(result, expected)

        result = self.frame.T.drop(['two'], axis=1, level='second')
        expected = self.frame.ix[[0, 2, 3, 6, 7, 9]].T
        assert_frame_equal(result, expected)

    def test_drop_preserve_names(self):
        index = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1],
                                        [1, 2, 3, 1, 2, 3]],
                                       names=['one', 'two'])

        df = DataFrame(np.random.randn(6, 3), index=index)

        result = df.drop([(0, 2)])
        self.assert_(result.index.names == ('one', 'two'))

    def test_unicode_repr_issues(self):
        levels = [Index([u('a/\u03c3'), u('b/\u03c3'), u('c/\u03c3')]),
                  Index([0, 1])]
        labels = [np.arange(3).repeat(2), np.tile(np.arange(2), 3)]
        index = MultiIndex(levels=levels, labels=labels)

        repr(index.levels)

        # NumPy bug
        # repr(index.get_level_values(1))

    def test_unicode_repr_level_names(self):
        index = MultiIndex.from_tuples([(0, 0), (1, 1)],
                                       names=[u('\u0394'), 'i1'])

        s = Series(lrange(2), index=index)
        df = DataFrame(np.random.randn(2, 4), index=index)
        repr(s)
        repr(df)

    def test_dataframe_insert_column_all_na(self):
        # GH #1534
        mix = MultiIndex.from_tuples(
            [('1a', '2a'), ('1a', '2b'), ('1a', '2c')])
        df = DataFrame([[1, 2], [3, 4], [5, 6]], index=mix)
        s = Series({(1, 1): 1, (1, 2): 2})
        df['new'] = s
        self.assert_(df['new'].isnull().all())

    def test_join_segfault(self):
        # 1532
        df1 = DataFrame({'a': [1, 1], 'b': [1, 2], 'x': [1, 2]})
        df2 = DataFrame({'a': [2, 2], 'b': [1, 2], 'y': [1, 2]})
        df1 = df1.set_index(['a', 'b'])
        df2 = df2.set_index(['a', 'b'])
        # it works!
        for how in ['left', 'right', 'outer']:
            df1.join(df2, how=how)

    def test_set_column_scalar_with_ix(self):
        subset = self.frame.index[[1, 4, 5]]

        self.frame.ix[subset] = 99
        self.assert_((self.frame.ix[subset].values == 99).all())

        col = self.frame['B']
        col[subset] = 97
        self.assert_((self.frame.ix[subset, 'B'] == 97).all())

    def test_frame_dict_constructor_empty_series(self):
        s1 = Series([1, 2, 3, 4], index=MultiIndex.from_tuples([(1, 2), (1, 3),
                                                              (2, 2), (2, 4)]))
        s2 = Series([1, 2, 3, 4],
                    index=MultiIndex.from_tuples([(1, 2), (1, 3), (3, 2), (3, 4)]))
        s3 = Series()

        # it works!
        df = DataFrame({'foo': s1, 'bar': s2, 'baz': s3})
        df = DataFrame.from_dict({'foo': s1, 'baz': s3, 'bar': s2})

    def test_indexing_ambiguity_bug_1678(self):
        columns = MultiIndex.from_tuples([('Ohio', 'Green'), ('Ohio', 'Red'),
                                          ('Colorado', 'Green')])
        index = MultiIndex.from_tuples(
            [('a', 1), ('a', 2), ('b', 1), ('b', 2)])

        frame = DataFrame(np.arange(12).reshape((4, 3)), index=index,
                          columns=columns)

        result = frame.ix[:, 1]
        exp = frame.icol(1)
        tm.assert_isinstance(result, Series)
        assert_series_equal(result, exp)

    def test_nonunique_assignment_1750(self):
        df = DataFrame([[1, 1, "x", "X"], [1, 1, "y", "Y"], [1, 2, "z", "Z"]],
                       columns=list("ABCD"))

        df = df.set_index(['A', 'B'])
        ix = MultiIndex.from_tuples([(1, 1)])

        df.ix[ix, "C"] = '_'

        self.assert_((df.xs((1, 1))['C'] == '_').all())

    def test_indexing_over_hashtable_size_cutoff(self):
        n = 10000

        old_cutoff = _index._SIZE_CUTOFF
        _index._SIZE_CUTOFF = 20000

        s = Series(np.arange(n),
                   MultiIndex.from_arrays((["a"] * n, np.arange(n))))

        # hai it works!
        self.assertEquals(s[("a", 5)], 5)
        self.assertEquals(s[("a", 6)], 6)
        self.assertEquals(s[("a", 7)], 7)

        _index._SIZE_CUTOFF = old_cutoff

    def test_multiindex_na_repr(self):
        # only an issue with long columns

        from numpy import nan
        df3 = DataFrame({
            'A' * 30: {('A', 'A0006000', 'nuit'): 'A0006000'},
            'B' * 30: {('A', 'A0006000', 'nuit'): nan},
            'C' * 30: {('A', 'A0006000', 'nuit'): nan},
            'D' * 30: {('A', 'A0006000', 'nuit'): nan},
            'E' * 30: {('A', 'A0006000', 'nuit'): 'A'},
            'F' * 30: {('A', 'A0006000', 'nuit'): nan},
        })

        idf = df3.set_index(['A' * 30, 'C' * 30])
        repr(idf)

    def test_assign_index_sequences(self):
        # #2200
        df = DataFrame({"a": [1, 2, 3],
                        "b": [4, 5, 6],
                        "c": [7, 8, 9]}).set_index(["a", "b"])
        l = list(df.index)
        l[0] = ("faz", "boo")
        df.index = l
        repr(df)

        # this travels an improper code path
        l[0] = ["faz", "boo"]
        df.index = l
        repr(df)

    def test_tuples_have_na(self):
        index = MultiIndex(levels=[[1, 0], [0, 1, 2, 3]],
                           labels=[[1, 1, 1, 1, -1, 0, 0, 0],
                                   [0, 1, 2, 3, 0, 1, 2, 3]])

        self.assertTrue(isnull(index[4][0]))
        self.assertTrue(isnull(index.values[4][0]))

    def test_duplicate_groupby_issues(self):
        idx_tp = [('600809', '20061231'), ('600809', '20070331'),
                  ('600809', '20070630'), ('600809', '20070331')]
        dt = ['demo','demo','demo','demo']

        idx = MultiIndex.from_tuples(idx_tp,names = ['STK_ID','RPT_Date'])
        s = Series(dt, index=idx)

        result = s.groupby(s.index).first()
        self.assertEquals(len(result), 3)

    def test_duplicate_mi(self):
        # GH 4516
        df = DataFrame([['foo','bar',1.0,1],['foo','bar',2.0,2],['bah','bam',3.0,3],
                        ['bah','bam',4.0,4],['foo','bar',5.0,5],['bah','bam',6.0,6]],
                       columns=list('ABCD'))
        df = df.set_index(['A','B'])
        df = df.sortlevel(0)
        expected = DataFrame([['foo','bar',1.0,1],['foo','bar',2.0,2],['foo','bar',5.0,5]],
                             columns=list('ABCD')).set_index(['A','B'])
        result = df.loc[('foo','bar')]
        assert_frame_equal(result,expected)

    def test_multiindex_set_index(self):
        # segfault in #3308
        d = {'t1': [2, 2.5, 3], 't2': [4, 5, 6]}
        df = DataFrame(d)
        tuples = [(0, 1), (0, 2), (1, 2)]
        df['tuples'] = tuples

        index = MultiIndex.from_tuples(df['tuples'])
        # it works!
        df.set_index(index)

if __name__ == '__main__':

    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
