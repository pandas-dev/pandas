import unittest

import numpy as np
from pandas import Index
from pandas.util.testing import assert_almost_equal
import pandas.util.testing as common
import pandas._tseries as lib

class TestTseriesUtil(unittest.TestCase):

    def test_combineFunc(self):
        pass

    def test_reindex(self):
        pass

    def test_isnull(self):
        pass

    def test_groupby(self):
        pass

    def test_groupby_withnull(self):
        pass

    def test_getMergeVec(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler, mask = lib.getFillVec(old, new, old.indexMap,
                                          new.indexMap, None)

        expect_filler = [-1, 0, -1, -1, -1, 1, -1, -1, -1, -1, 2, -1]
        expect_mask = np.zeros(12, dtype=bool)
        expect_mask[[1, 5, 10]] = True

        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

        # corner case
        old = Index([1, 4])
        new = Index(range(5, 10))
        filler, mask = lib.getFillVec(old, new, old.indexMap,
                                          new.indexMap, None)

        expect_filler = [-1, -1, -1, -1, -1]
        expect_mask = np.zeros(5, dtype=bool)
        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

    def test_backfill(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler, mask = lib.getFillVec(old, new, old.indexMap,
                                          new.indexMap, 'BACKFILL')

        expect_filler = [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, -1]
        expect_mask = np.ones(12, dtype=bool)
        expect_mask[-1] = False

        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

        # corner case
        old = Index([1, 4])
        new = Index(range(5, 10))
        filler, mask = lib.getFillVec(old, new, old.indexMap,
                                          new.indexMap, 'BACKFILL')

        expect_filler = [-1, -1, -1, -1, -1]
        expect_mask = np.zeros(5, dtype=bool)
        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

    def test_pad(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler, mask = lib.getFillVec(old, new, old.indexMap,
                                          new.indexMap, 'PAD')

        expect_filler = [-1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2]
        expect_mask = np.ones(12, dtype=bool)
        expect_mask[0] = False

        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

        # corner case
        old = Index([5, 10])
        new = Index(range(5))
        filler, mask = lib.getFillVec(old, new, old.indexMap,
                                          new.indexMap, 'PAD')

        expect_filler = [-1, -1, -1, -1, -1]
        expect_mask = np.zeros(5, dtype=bool)
        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

def test_inner_join_indexer():
    a = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    b = np.array([0, 3, 5, 7, 9], dtype=np.int64)

    index, ares, bres = lib.inner_join_indexer(a, b)

    index_exp = np.array([3, 5], dtype=np.int64)
    assert_almost_equal(index, index_exp)

    aexp = np.array([2, 4])
    bexp = np.array([1, 2])
    assert_almost_equal(ares, aexp)
    assert_almost_equal(bres, bexp)

def test_outer_join_indexer():
    a = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    b = np.array([0, 3, 5, 7, 9], dtype=np.int64)

    index, ares, bres = lib.outer_join_indexer(a, b)

    index_exp = np.array([0, 1, 2, 3, 4, 5, 7, 9], dtype=np.int64)
    assert_almost_equal(index, index_exp)

    aexp = np.array([-1, 0, 1, 2, 3, 4, -1, -1], dtype=np.int32)
    bexp = np.array([0, -1, -1, 1, -1, 2, 3, 4])
    assert_almost_equal(ares, aexp)
    assert_almost_equal(bres, bexp)

class TestMoments(unittest.TestCase):
    pass

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

