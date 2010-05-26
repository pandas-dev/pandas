import unittest

import numpy as np
from pandas import Index
import pandas.util.testing as common
import pandas.lib.tseries as tseries

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

        filler, mask = tseries.getFillVec(old, new, old.indexMap,
                                          new.indexMap, None)

        expect_filler = [-1, 0, -1, -1, -1, 1, -1, -1, -1, -1, 2, -1]
        expect_mask = np.zeros(12, dtype=bool)
        expect_mask[[1, 5, 10]] = True

        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

        # corner case
        old = Index([1, 4])
        new = Index(range(5, 10))
        filler, mask = tseries.getFillVec(old, new, old.indexMap,
                                          new.indexMap, None)

        expect_filler = [-1, -1, -1, -1, -1]
        expect_mask = np.zeros(5, dtype=bool)
        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

    def test_backfill(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler, mask = tseries.getFillVec(old, new, old.indexMap,
                                          new.indexMap, 'BACKFILL')

        expect_filler = [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, -1]
        expect_mask = np.ones(12, dtype=bool)
        expect_mask[-1] = False

        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

        # corner case
        old = Index([1, 4])
        new = Index(range(5, 10))
        filler, mask = tseries.getFillVec(old, new, old.indexMap,
                                          new.indexMap, 'BACKFILL')

        expect_filler = [-1, -1, -1, -1, -1]
        expect_mask = np.zeros(5, dtype=bool)
        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

    def test_pad(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler, mask = tseries.getFillVec(old, new, old.indexMap,
                                          new.indexMap, 'PAD')

        expect_filler = [-1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2]
        expect_mask = np.ones(12, dtype=bool)
        expect_mask[0] = False

        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

        # corner case
        old = Index([5, 10])
        new = Index(range(5))
        filler, mask = tseries.getFillVec(old, new, old.indexMap,
                                          new.indexMap, 'PAD')

        expect_filler = [-1, -1, -1, -1, -1]
        expect_mask = np.zeros(5, dtype=bool)
        self.assert_(np.array_equal(filler, expect_filler))
        self.assert_(np.array_equal(mask, expect_mask))

class TestMoments(unittest.TestCase):
    pass
