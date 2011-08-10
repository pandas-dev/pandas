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

        self.tdf = tm.makeTimeDataFrame()
        self.ymd = self.tdf.groupby([lambda x: x.year, lambda x: x.month,
                                     lambda x: x.day]).sum()

    def test_getitem_simple(self):
        df = self.frame.T

        col = df['foo', 'one']
        assert_almost_equal(col.values, df.values[:, 0])
        self.assertRaises(KeyError, df.__getitem__, ('foo', 'four'))

    def test_xs(self):
        xs = self.frame.xs(('bar', 'two'))
        assert_almost_equal(xs.values, self.frame.values[4])

    def test_xs_partial(self):
        result = self.frame.xs('foo')
        expected = self.frame.T['foo'].T
        assert_frame_equal(result, expected)

    def test_getitem_toplevel(self):
        df = self.frame.T

        result = df['foo']
        expected = df.reindex(columns=df.columns[:3])
        assert_frame_equal(result, expected)

        result = df['bar']
        expected = df.reindex(columns=df.columns[3:5])
        assert_frame_equal(result, expected)

    def test_getitem_partial(self):
        ymd = self.ymd.T
        result = ymd[2000, 2]
        expected = ymd.reindex(columns=ymd.columns[ymd.columns.labels[1] == 1])
        assert_frame_equal(result, expected)

    def test_fancy_slice_partial(self):
        pass

    def test_fancy_select_toplevel(self):
        pass

    def test_alignment(self):
        pass


if __name__ == '__main__':
    # unittest.main()
    import nose
    # nose.runmodule(argv=[__file__,'-vvs','-x', '--pdb-failure'],
    #                exit=False)
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

