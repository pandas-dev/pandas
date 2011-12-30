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

    def test_merge_indexer(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler = lib.merge_indexer_object(new, old.indexMap)

        expect_filler = [-1, 0, -1, -1, -1, 1, -1, -1, -1, -1, 2, -1]
        self.assert_(np.array_equal(filler, expect_filler))

        # corner case
        old = Index([1, 4])
        new = Index(range(5, 10))
        filler = lib.merge_indexer_object(new, old.indexMap)
        expect_filler = [-1, -1, -1, -1, -1]
        self.assert_(np.array_equal(filler, expect_filler))

    def test_backfill(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler = lib.backfill_object(old, new, old.indexMap, new.indexMap)

        expect_filler = [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, -1]
        self.assert_(np.array_equal(filler, expect_filler))

        # corner case
        old = Index([1, 4])
        new = Index(range(5, 10))
        filler = lib.backfill_object(old, new, old.indexMap, new.indexMap)

        expect_filler = [-1, -1, -1, -1, -1]
        self.assert_(np.array_equal(filler, expect_filler))

    def test_pad(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler = lib.pad_object(old, new, old.indexMap, new.indexMap)

        expect_filler = [-1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2]
        self.assert_(np.array_equal(filler, expect_filler))

        # corner case
        old = Index([5, 10])
        new = Index(range(5))
        filler = lib.pad_object(old, new, old.indexMap, new.indexMap)
        expect_filler = [-1, -1, -1, -1, -1]
        self.assert_(np.array_equal(filler, expect_filler))

def test_inner_join_indexer():
    a = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    b = np.array([0, 3, 5, 7, 9], dtype=np.int64)

    index, ares, bres = lib.inner_join_indexer_int64(a, b)

    index_exp = np.array([3, 5], dtype=np.int64)
    assert_almost_equal(index, index_exp)

    aexp = np.array([2, 4])
    bexp = np.array([1, 2])
    assert_almost_equal(ares, aexp)
    assert_almost_equal(bres, bexp)

def test_outer_join_indexer():
    a = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    b = np.array([0, 3, 5, 7, 9], dtype=np.int64)

    index, ares, bres = lib.outer_join_indexer_int64(a, b)

    index_exp = np.array([0, 1, 2, 3, 4, 5, 7, 9], dtype=np.int64)
    assert_almost_equal(index, index_exp)

    aexp = np.array([-1, 0, 1, 2, 3, 4, -1, -1], dtype=np.int32)
    bexp = np.array([0, -1, -1, 1, -1, 2, 3, 4])
    assert_almost_equal(ares, aexp)
    assert_almost_equal(bres, bexp)

def test_is_lexsorted():
    failure = [
        np.array([3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                  3, 3,
               3, 3,
           3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
           2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0]),
        np.array([30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16,
                  15, 14,
       13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0, 30, 29, 28,
       27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11,
       10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0, 30, 29, 28, 27, 26, 25,
       24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10,  9,  8,
        7,  6,  5,  4,  3,  2,  1,  0, 30, 29, 28, 27, 26, 25, 24, 23, 22,
       21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,
        4,  3,  2,  1,  0])]

    assert(not lib.is_lexsorted(failure))

# def test_get_group_index():
#     a = np.array([0, 1, 2, 0, 2, 1, 0, 0], dtype='i4')
#     b = np.array([1, 0, 3, 2, 0, 2, 3, 0], dtype='i4')
#     expected = np.array([1, 4, 11, 2, 8, 6, 3, 0], dtype='i4')

#     result = lib.get_group_index([a, b], (3, 4))

#     assert(np.array_equal(result, expected))

def test_groupsort_indexer():
    a = np.random.randint(0, 1000, 100).astype('i4')
    b = np.random.randint(0, 1000, 100).astype('i4')

    result = lib.groupsort_indexer(a, 1000)[0]

    # need to use a stable sort
    expected = np.argsort(a, kind='mergesort')
    assert(np.array_equal(result, expected))

    # compare with lexsort
    key = a * 1000 + b
    result = lib.groupsort_indexer(key, 1000000)[0]
    expected = np.lexsort((b, a))
    assert(np.array_equal(result, expected))


def test_duplicated_with_nas():
    keys = [0, 1, np.nan, 0, 2, np.nan]

    result = lib.duplicated(keys)
    expected = [False, False, False, True, False, True]
    assert(np.array_equal(result, expected))

    result = lib.duplicated(keys, take_last=True)
    expected = [True, False, True, False, False, False]
    assert(np.array_equal(result, expected))

    keys = [(0, 0), (0, np.nan), (np.nan, 0), (np.nan, np.nan)] * 2

    result = lib.duplicated(keys)
    falses = [False] * 4
    trues = [True] * 4
    expected = falses + trues
    assert(np.array_equal(result, expected))

    result = lib.duplicated(keys, take_last=True)
    expected = trues + falses
    assert(np.array_equal(result, expected))

def test_convert_objects():
    arr = np.array(['a', 'b', np.nan, np.nan, 'd', 'e', 'f'], dtype='O')
    result = lib.maybe_convert_objects(arr)
    assert(result.dtype == np.object_)

def test_convert_objects_ints():
    # test that we can detect many kinds of integers
    dtypes = ['i1', 'i2', 'i4', 'i8', 'u1', 'u2', 'u4', 'u8']

    for dtype_str in dtypes:
        arr = np.array(list(np.arange(20, dtype=dtype_str)), dtype='O')
        assert(arr[0].dtype == np.dtype(dtype_str))
        result = lib.maybe_convert_objects(arr)
        assert(issubclass(result.dtype.type, np.integer))

def test_rank():
    from scipy.stats import rankdata
    from numpy import nan
    def _check(arr):
        mask = -np.isfinite(arr)
        arr = arr.copy()
        result = lib.rank_1d_float64(arr)
        arr[mask] = np.inf
        exp = rankdata(arr)
        exp[mask] = np.nan
        assert_almost_equal(result, exp)

    _check(np.array([nan, nan, 5., 5., 5., nan, 1, 2, 3, nan]))
    _check(np.array([4., nan, 5., 5., 5., nan, 1, 2, 4., nan]))

def test_get_reverse_indexer():
    indexer = np.array([-1, -1, 1, 2, 0, -1, 3, 4], dtype='i4')
    result = lib.get_reverse_indexer(indexer, 5)
    expected = np.array([4, 2, 3, 6, 7], dtype='i4')
    assert(np.array_equal(result, expected))

class TestMoments(unittest.TestCase):
    pass

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

