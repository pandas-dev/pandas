# pylint: disable-msg=W0612,E1101
from copy import deepcopy
from datetime import datetime, timedelta
from cStringIO import StringIO
import cPickle as pickle
import operator
import os
import unittest

from numpy import random, nan
from numpy.random import randn
import numpy as np

import pandas.core.datetools as datetools
from pandas.core.index import MultiIndex, NULL_INDEX
from pandas.core.api import (DataFrame, Index, Series, notnull, isnull)

from pandas.util.testing import (assert_almost_equal,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 randn)

import pandas.util.testing as tm

class TestSeriesMultiLevel(unittest.TestCase):
    pass

class TestDataFrameMultiLevel(unittest.TestCase):

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

    def test_repr_to_string(self):
        repr(self.frame)
        repr(self.ymd)
        repr(self.frame.T)
        repr(self.ymd.T)

        buf = StringIO()
        self.frame.toString(buf=buf)
        self.ymd.toString(buf=buf)
        self.frame.T.toString(buf=buf)
        self.ymd.T.toString(buf=buf)

    def test_getitem_simple(self):
        df = self.frame.T

        col = df['foo', 'one']
        assert_almost_equal(col.values, df.values[:, 0])
        self.assertRaises(KeyError, df.__getitem__, ('foo', 'four'))
        self.assertRaises(KeyError, df.__getitem__, 'foobar')

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
        assert_frame_equal(sorted_before, sorted_after.drop(['foo'], axis=1))

    def test_alignment(self):
        pass


if __name__ == '__main__':
    # unittest.main()
    import nose
    # nose.runmodule(argv=[__file__,'-vvs','-x', '--pdb-failure'],
    #                exit=False)
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

