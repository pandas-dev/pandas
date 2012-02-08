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

        filler = lib.merge_indexer_int64(new, old.indexMap)

        expect_filler = [-1, 0, -1, -1, -1, 1, -1, -1, -1, -1, 2, -1]
        self.assert_(np.array_equal(filler, expect_filler))

        # corner case
        old = Index([1, 4])
        new = Index(range(5, 10))
        filler = lib.merge_indexer_int64(new, old.indexMap)
        expect_filler = [-1, -1, -1, -1, -1]
        self.assert_(np.array_equal(filler, expect_filler))

    def test_backfill(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler = lib.backfill_int64(old, new, old.indexMap, new.indexMap)

        expect_filler = [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, -1]
        self.assert_(np.array_equal(filler, expect_filler))

        # corner case
        old = Index([1, 4])
        new = Index(range(5, 10))
        filler = lib.backfill_int64(old, new, old.indexMap, new.indexMap)

        expect_filler = [-1, -1, -1, -1, -1]
        self.assert_(np.array_equal(filler, expect_filler))

    def test_pad(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler = lib.pad_int64(old, new, old.indexMap, new.indexMap)

        expect_filler = [-1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2]
        self.assert_(np.array_equal(filler, expect_filler))

        # corner case
        old = Index([5, 10])
        new = Index(range(5))
        filler = lib.pad_int64(old, new, old.indexMap, new.indexMap)
        expect_filler = [-1, -1, -1, -1, -1]
        self.assert_(np.array_equal(filler, expect_filler))

def test_left_join_indexer():
    a = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    b = np.array([2, 2, 3, 4, 4], dtype=np.int64)

    result = lib.left_join_indexer_int64(b, a)
    expected = np.array([1, 1, 2, 3, 3], dtype='i4')
    assert(np.array_equal(result, expected))

def test_left_outer_join_bug():
    left = np.array([0, 1, 0, 1, 1, 2, 3, 1, 0, 2, 1, 2, 0, 1, 1, 2, 3, 2, 3,
                     2, 1, 1, 3, 0, 3, 2, 3, 0, 0, 2, 3, 2, 0, 3, 1, 3, 0, 1,
                     3, 0, 0, 1, 0, 3, 1, 0, 1, 0, 1, 1, 0, 2, 2, 2, 2, 2, 0,
                     3, 1, 2, 0, 0, 3, 1, 3, 2, 2, 0, 1, 3, 0, 2, 3, 2, 3, 3,
                     2, 3, 3, 1, 3, 2, 0, 0, 3, 1, 1, 1, 0, 2, 3, 3, 1, 2, 0,
                     3, 1, 2, 0, 2], dtype=np.int32)

    right = np.array([3, 1], dtype=np.int32)
    max_groups = 4

    lidx, ridx = lib.left_outer_join(left, right, max_groups, sort=False)

    exp_lidx = np.arange(len(left))
    exp_ridx = -np.ones(len(left))
    exp_ridx[left == 1] = 1
    exp_ridx[left == 3] = 0

    assert(np.array_equal(lidx, exp_lidx))
    assert(np.array_equal(ridx, exp_ridx))

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

def test_pad_backfill_object_segfault():
    from datetime import datetime
    old = np.array([], dtype='O')
    new = np.array([datetime(2010, 12, 31)], dtype='O')

    result = lib.pad_object(old, new, lib.map_indices_object(old),
                            lib.map_indices_object(new))
    expected = np.array([-1], dtype='i4')
    assert(np.array_equal(result, expected))

    result = lib.pad_object(new, old, lib.map_indices_object(new),
                            lib.map_indices_object(old))
    expected = np.array([], dtype='i4')
    assert(np.array_equal(result, expected))

    result = lib.backfill_object(old, new, lib.map_indices_object(old),
                                 lib.map_indices_object(new))
    expected = np.array([-1], dtype='i4')
    assert(np.array_equal(result, expected))

    result = lib.backfill_object(new, old, lib.map_indices_object(new),
                            lib.map_indices_object(old))
    expected = np.array([], dtype='i4')
    assert(np.array_equal(result, expected))

def test_arrmap():
    values = np.array(['foo', 'foo', 'bar', 'bar', 'baz', 'qux'], dtype='O')
    result = lib.arrmap_object(values, lambda x: x in ['foo', 'bar'])
    assert(result.dtype == np.bool_)

def test_series_grouper():
    from pandas import Series
    obj = Series(np.random.randn(10))
    dummy = obj[:0]

    labels = np.array([-1, -1, -1, 0, 0, 0, 1, 1, 1, 1], dtype='i4')

    grouper = lib.SeriesGrouper(obj, np.mean, labels, 2, dummy)
    result, counts = grouper.get_result()

    expected = np.array([obj[3:6].mean(), obj[6:].mean()])
    assert_almost_equal(result, expected)

    exp_counts = np.array([3, 4], dtype=np.int32)
    assert_almost_equal(counts, exp_counts)

class TestTypeInference(unittest.TestCase):

    def test_length_zero(self):
        result = lib.infer_dtype(np.array([], dtype='i4'))
        self.assertEqual(result, 'empty')

        result = lib.infer_dtype(np.array([], dtype='O'))
        self.assertEqual(result, 'empty')

    def test_integers(self):
        arr = np.array([1, 2, 3, np.int64(4), np.int32(5)], dtype='O')
        result = lib.infer_dtype(arr)
        self.assertEqual(result, 'integer')

        arr = np.array([1, 2, 3, np.int64(4), np.int32(5), 'foo'],
                       dtype='O')
        result = lib.infer_dtype(arr)
        self.assertEqual(result, 'mixed')

        arr = np.array([1, 2, 3, 4, 5], dtype='i4')
        result = lib.infer_dtype(arr)
        self.assertEqual(result, 'integer')

    def test_bools(self):
        arr = np.array([True, False, True, True, True], dtype='O')
        result = lib.infer_dtype(arr)
        self.assertEqual(result, 'boolean')

        arr = np.array([np.bool_(True), np.bool_(False)], dtype='O')
        result = lib.infer_dtype(arr)
        self.assertEqual(result, 'boolean')

        arr = np.array([True, False, True, 'foo'], dtype='O')
        result = lib.infer_dtype(arr)
        self.assertEqual(result, 'mixed')

        arr = np.array([True, False, True], dtype=bool)
        result = lib.infer_dtype(arr)
        self.assertEqual(result, 'boolean')

    def test_floats(self):
        arr = np.array([1., 2., 3., np.float64(4), np.float32(5)], dtype='O')
        result = lib.infer_dtype(arr)
        self.assertEqual(result, 'floating')

        arr = np.array([1, 2, 3, np.float64(4), np.float32(5), 'foo'],
                       dtype='O')
        result = lib.infer_dtype(arr)
        self.assertEqual(result, 'mixed')

        arr = np.array([1, 2, 3, 4, 5], dtype='f4')
        result = lib.infer_dtype(arr)
        self.assertEqual(result, 'floating')

        arr = np.array([1, 2, 3, 4, 5], dtype='f8')
        result = lib.infer_dtype(arr)
        self.assertEqual(result, 'floating')

    def test_string(self):
        pass

    def test_unicode(self):
        pass

    def test_datetime(self):
        pass

    def test_to_object_array_tuples(self):
        r = (5,6)
        values = [r]
        result = lib.to_object_array_tuples(values)

        try:
            # make sure record array works
            from collections import namedtuple
            record = namedtuple('record', 'x y')
            r = record(5, 6)
            values = [r]
            result = lib.to_object_array_tuples(values)
        except ImportError:
            pass

class TestMoments(unittest.TestCase):
    pass

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

