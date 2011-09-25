# pylint: disable-msg=W0612,E1101
from cStringIO import StringIO
import operator
import unittest

from numpy import random, nan
from numpy.random import randn
import numpy as np

import pandas.core.datetools as datetools
from pandas.core.index import MultiIndex, NULL_INDEX
from pandas import Panel, DataFrame, Index, Series, notnull, isnull

from pandas.util.testing import (assert_almost_equal,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 randn)

import pandas.util.testing as tm

class TestMultiLevel(unittest.TestCase):

    def setUp(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]])
        self.frame = DataFrame(np.random.randn(10, 3), index=index,
                               columns=['A', 'B', 'C'])

        tm.N = 100
        self.tdf = tm.makeTimeDataFrame()
        self.ymd = self.tdf.groupby([lambda x: x.year, lambda x: x.month,
                                     lambda x: x.day]).sum()

    def test_pickle(self):
        import cPickle
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
        expected = s[42:65]
        expected.index = expected.index.droplevel(0).droplevel(0)
        assert_series_equal(result, expected)

        result = s[2000, 3, 10]
        expected = s[49]
        self.assertEquals(result, expected)

        # fancy
        result = s.ix[[(2000, 3, 10), (2000, 3, 13)]]
        expected = s[49:51]
        assert_series_equal(result, expected)

        # key error
        self.assertRaises(KeyError, s.__getitem__, (2000, 3, 4))

    def test_series_setitem(self):
        s = self.ymd['A']

        s[2000, 3] = np.nan
        self.assert_(isnull(s[42:65]).all())
        self.assert_(notnull(s[:42]).all())
        self.assert_(notnull(s[65:]).all())

        s[2000, 3, 10] = np.nan
        self.assert_(isnull(s[49]))

    def test_series_slice_partial(self):
        pass

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

    def test_fancy_2d(self):
        result = self.frame.ix['foo', 'B']
        expected = self.frame.xs('foo')['B']
        assert_series_equal(result, expected)

        ft = self.frame.T
        result = ft.ix['B', 'foo']
        expected = ft.xs('B')['foo']
        assert_series_equal(result, expected)

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

    def test_setitem_change_dtype(self):
        dft = self.frame.T
        s = dft['foo', 'two']
        dft['foo', 'two'] = s > s.median()
        assert_series_equal(dft['foo', 'two'], s > s.median())
        self.assert_(isinstance(dft._data.blocks[1].items, MultiIndex))

        reindexed = dft.reindex(columns=[('foo', 'two')])
        assert_series_equal(reindexed['foo', 'two'], s > s.median())

    def test_fancy_slice_partial(self):
        result = self.frame.ix['bar':'baz']
        expected = self.frame[3:7]
        assert_frame_equal(result, expected)

        result = self.ymd.ix[(2000,2):(2000,4)]
        lev = self.ymd.index.labels[1]
        expected = self.ymd[(lev >= 1) & (lev <= 3)]
        assert_frame_equal(result, expected)

    def test_sortlevel(self):
        df = self.frame.copy()
        df.index = np.arange(len(df))
        self.assertRaises(Exception, df.sortlevel, 0)

        # axis=1

        # series
        a_sorted = self.frame['A'].sortlevel(0)
        self.assertRaises(Exception,
                          self.frame.delevel()['A'].sortlevel)

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

        _check_counts(self.frame)
        _check_counts(self.ymd)
        _check_counts(self.frame.T, axis=1)
        _check_counts(self.ymd.T, axis=1)

        # can't call with level on regular DataFrame
        df = tm.makeTimeDataFrame()
        self.assertRaises(Exception, df.count, level=0)

    def test_count_level_corner(self):
        s = self.frame['A'][:0]
        result = s.count(level=0)
        expected = Series(0, index=s.index.levels[0])
        assert_series_equal(result, expected)

        df = self.frame[:0]
        result = df.count(level=0)
        expected = DataFrame({}, index=s.index.levels[0],
                             columns=df.columns).fillna(0).astype(int)
        assert_frame_equal(result, expected)

    def test_unstack(self):
        # just check that it works for now
        unstacked = self.ymd.unstack()
        unstacked2 = unstacked.unstack()

        # test that ints work
        unstacked = self.ymd.astype(int).unstack()

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

    def test_stack_mixed_dtype(self):
        df = self.frame.T
        df['foo', 'four'] = 'foo'
        df = df.sortlevel(1, axis=1)

        stacked = df.stack()
        assert_series_equal(stacked['foo'], df['foo'].stack())
        self.assert_(stacked['bar'].dtype == np.float_)

    def test_swaplevel(self):
        swapped = self.frame['A'].swaplevel(0, 1)
        self.assert_(not swapped.index.equals(self.frame.index))

        back = swapped.swaplevel(0, 1)
        self.assert_(back.index.equals(self.frame.index))

    def test_swaplevel_panel(self):
        panel = Panel({'ItemA' : self.frame,
                       'ItemB' : self.frame * 2})

        result = panel.swaplevel(0, 1, axis='major')
        expected = panel.copy()
        expected.major_axis = expected.major_axis.swaplevel(0, 1)

    def test_insert_index(self):
        df = self.ymd[:5].T
        df[2000, 1, 10] = df[2000, 1, 7]
        self.assert_(isinstance(df.columns, MultiIndex))
        self.assert_((df[2000, 1, 10] == df[2000, 1, 7]).all())

    def test_alignment(self):
        pass

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
        df = self.frame.T
        df['foo'].values[:] = 0
        self.assert_((df['foo'].values == 0).all())

        # but not if it's mixed-type
        df['foo', 'four'] = 'foo'
        df = df.sortlevel(0, axis=1)
        df['foo']['one'] = 2
        self.assert_((df['foo', 'one'] == 0).all())

    def test_frame_getitem_not_sorted(self):
        df = self.frame.T
        df['foo', 'four'] = 'foo'

        arrays = [np.array(x) for x in zip(*df.columns.get_tuple_index())]

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
        tuples = zip(*arrays)
        index = MultiIndex.from_tuples(tuples)
        s = Series(randn(8), index=index)

        arrays = [np.array(x) for x in zip(*index.get_tuple_index())]

        result = s['qux']
        result2 = s.ix['qux']
        expected = s[arrays[0] == 'qux']
        expected.index = expected.index.droplevel(0)
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

if __name__ == '__main__':

    # unittest.main()
    import nose
    # nose.runmodule(argv=[__file__,'-vvs','-x', '--pdb-failure'],
    #                exit=False)
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

