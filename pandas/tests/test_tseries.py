import unittest

from numpy import nan
import numpy as np
from pandas import Index, isnull, Timestamp
from pandas.util.testing import assert_almost_equal
import pandas.util.testing as common
import pandas.lib as lib
import pandas.algos as algos
from datetime import datetime


class TestTseriesUtil(unittest.TestCase):
    _multiprocess_can_split_ = True

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

    def test_backfill(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler = algos.backfill_int64(old, new)

        expect_filler = [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, -1]
        self.assert_(np.array_equal(filler, expect_filler))

        # corner case
        old = Index([1, 4])
        new = Index(range(5, 10))
        filler = algos.backfill_int64(old, new)

        expect_filler = [-1, -1, -1, -1, -1]
        self.assert_(np.array_equal(filler, expect_filler))

    def test_pad(self):
        old = Index([1, 5, 10])
        new = Index(range(12))

        filler = algos.pad_int64(old, new)

        expect_filler = [-1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2]
        self.assert_(np.array_equal(filler, expect_filler))

        # corner case
        old = Index([5, 10])
        new = Index(range(5))
        filler = algos.pad_int64(old, new)
        expect_filler = [-1, -1, -1, -1, -1]
        self.assert_(np.array_equal(filler, expect_filler))


def test_left_join_indexer_unique():
    a = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    b = np.array([2, 2, 3, 4, 4], dtype=np.int64)

    result = algos.left_join_indexer_unique_int64(b, a)
    expected = np.array([1, 1, 2, 3, 3], dtype=np.int64)
    assert(np.array_equal(result, expected))


def test_left_outer_join_bug():
    left = np.array([0, 1, 0, 1, 1, 2, 3, 1, 0, 2, 1, 2, 0, 1, 1, 2, 3, 2, 3,
                     2, 1, 1, 3, 0, 3, 2, 3, 0, 0, 2, 3, 2, 0, 3, 1, 3, 0, 1,
                     3, 0, 0, 1, 0, 3, 1, 0, 1, 0, 1, 1, 0, 2, 2, 2, 2, 2, 0,
                     3, 1, 2, 0, 0, 3, 1, 3, 2, 2, 0, 1, 3, 0, 2, 3, 2, 3, 3,
                     2, 3, 3, 1, 3, 2, 0, 0, 3, 1, 1, 1, 0, 2, 3, 3, 1, 2, 0,
                     3, 1, 2, 0, 2], dtype=np.int64)

    right = np.array([3, 1], dtype=np.int64)
    max_groups = 4

    lidx, ridx = algos.left_outer_join(left, right, max_groups, sort=False)

    exp_lidx = np.arange(len(left))
    exp_ridx = -np.ones(len(left))
    exp_ridx[left == 1] = 1
    exp_ridx[left == 3] = 0

    assert(np.array_equal(lidx, exp_lidx))
    assert(np.array_equal(ridx, exp_ridx))


def test_inner_join_indexer():
    a = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    b = np.array([0, 3, 5, 7, 9], dtype=np.int64)

    index, ares, bres = algos.inner_join_indexer_int64(a, b)

    index_exp = np.array([3, 5], dtype=np.int64)
    assert_almost_equal(index, index_exp)

    aexp = np.array([2, 4])
    bexp = np.array([1, 2])
    assert_almost_equal(ares, aexp)
    assert_almost_equal(bres, bexp)

    a = np.array([5], dtype=np.int64)
    b = np.array([5], dtype=np.int64)

    index, ares, bres = algos.inner_join_indexer_int64(a, b)
    assert_almost_equal(index, [5])
    assert_almost_equal(ares, [0])
    assert_almost_equal(bres, [0])


def test_outer_join_indexer():
    a = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    b = np.array([0, 3, 5, 7, 9], dtype=np.int64)

    index, ares, bres = algos.outer_join_indexer_int64(a, b)

    index_exp = np.array([0, 1, 2, 3, 4, 5, 7, 9], dtype=np.int64)
    assert_almost_equal(index, index_exp)

    aexp = np.array([-1, 0, 1, 2, 3, 4, -1, -1], dtype=np.int64)
    bexp = np.array([0, -1, -1, 1, -1, 2, 3, 4])
    assert_almost_equal(ares, aexp)
    assert_almost_equal(bres, bexp)

    a = np.array([5], dtype=np.int64)
    b = np.array([5], dtype=np.int64)

    index, ares, bres = algos.outer_join_indexer_int64(a, b)
    assert_almost_equal(index, [5])
    assert_almost_equal(ares, [0])
    assert_almost_equal(bres, [0])


def test_left_join_indexer():
    a = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    b = np.array([0, 3, 5, 7, 9], dtype=np.int64)

    index, ares, bres = algos.left_join_indexer_int64(a, b)

    assert_almost_equal(index, a)

    aexp = np.array([0, 1, 2, 3, 4], dtype=np.int64)
    bexp = np.array([-1, -1, 1, -1, 2], dtype=np.int64)
    assert_almost_equal(ares, aexp)
    assert_almost_equal(bres, bexp)

    a = np.array([5], dtype=np.int64)
    b = np.array([5], dtype=np.int64)

    index, ares, bres = algos.left_join_indexer_int64(a, b)
    assert_almost_equal(index, [5])
    assert_almost_equal(ares, [0])
    assert_almost_equal(bres, [0])


def test_left_join_indexer2():
    idx = Index([1, 1, 2, 5])
    idx2 = Index([1, 2, 5, 7, 9])

    res, lidx, ridx = algos.left_join_indexer_int64(idx2, idx)

    exp_res = np.array([1, 1, 2, 5, 7, 9], dtype=np.int64)
    assert_almost_equal(res, exp_res)

    exp_lidx = np.array([0, 0, 1, 2, 3, 4], dtype=np.int64)
    assert_almost_equal(lidx, exp_lidx)

    exp_ridx = np.array([0, 1, 2, 3, -1, -1], dtype=np.int64)
    assert_almost_equal(ridx, exp_ridx)


def test_outer_join_indexer2():
    idx = Index([1, 1, 2, 5])
    idx2 = Index([1, 2, 5, 7, 9])

    res, lidx, ridx = algos.outer_join_indexer_int64(idx2, idx)

    exp_res = np.array([1, 1, 2, 5, 7, 9], dtype=np.int64)
    assert_almost_equal(res, exp_res)

    exp_lidx = np.array([0, 0, 1, 2, 3, 4], dtype=np.int64)
    assert_almost_equal(lidx, exp_lidx)

    exp_ridx = np.array([0, 1, 2, 3, -1, -1], dtype=np.int64)
    assert_almost_equal(ridx, exp_ridx)


def test_inner_join_indexer2():
    idx = Index([1, 1, 2, 5])
    idx2 = Index([1, 2, 5, 7, 9])

    res, lidx, ridx = algos.inner_join_indexer_int64(idx2, idx)

    exp_res = np.array([1, 1, 2, 5], dtype=np.int64)
    assert_almost_equal(res, exp_res)

    exp_lidx = np.array([0, 0, 1, 2], dtype=np.int64)
    assert_almost_equal(lidx, exp_lidx)

    exp_ridx = np.array([0, 1, 2, 3], dtype=np.int64)
    assert_almost_equal(ridx, exp_ridx)


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
                  13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 30, 29, 28,
                  27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11,
                  10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 30, 29, 28, 27, 26, 25,
                  24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8,
                  7, 6, 5, 4, 3, 2, 1, 0, 30, 29, 28, 27, 26, 25, 24, 23, 22,
                  21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5,
                  4, 3, 2, 1, 0])]

    assert(not algos.is_lexsorted(failure))

# def test_get_group_index():
#     a = np.array([0, 1, 2, 0, 2, 1, 0, 0], dtype=np.int64)
#     b = np.array([1, 0, 3, 2, 0, 2, 3, 0], dtype=np.int64)
#     expected = np.array([1, 4, 11, 2, 8, 6, 3, 0], dtype=np.int64)

#     result = lib.get_group_index([a, b], (3, 4))

#     assert(np.array_equal(result, expected))


def test_groupsort_indexer():
    a = np.random.randint(0, 1000, 100).astype(np.int64)
    b = np.random.randint(0, 1000, 100).astype(np.int64)

    result = algos.groupsort_indexer(a, 1000)[0]

    # need to use a stable sort
    expected = np.argsort(a, kind='mergesort')
    assert(np.array_equal(result, expected))

    # compare with lexsort
    key = a * 1000 + b
    result = algos.groupsort_indexer(key, 1000000)[0]
    expected = np.lexsort((b, a))
    assert(np.array_equal(result, expected))


def test_ensure_platform_int():
    arr = np.arange(100)

    result = algos.ensure_platform_int(arr)
    assert(result is arr)


def test_duplicated_with_nas():
    keys = np.array([0, 1, nan, 0, 2, nan], dtype=object)

    result = lib.duplicated(keys)
    expected = [False, False, False, True, False, True]
    assert(np.array_equal(result, expected))

    result = lib.duplicated(keys, take_last=True)
    expected = [True, False, True, False, False, False]
    assert(np.array_equal(result, expected))

    keys = np.empty(8, dtype=object)
    for i, t in enumerate(zip([0, 0, nan, nan] * 2, [0, nan, 0, nan] * 2)):
        keys[i] = t

    result = lib.duplicated(keys)
    falses = [False] * 4
    trues = [True] * 4
    expected = falses + trues
    assert(np.array_equal(result, expected))

    result = lib.duplicated(keys, take_last=True)
    expected = trues + falses
    assert(np.array_equal(result, expected))


def test_maybe_booleans_to_slice():
    arr = np.array([0, 0, 1, 1, 1, 0, 1], dtype=np.uint8)
    result = lib.maybe_booleans_to_slice(arr)
    assert(result.dtype == np.bool_)

    result = lib.maybe_booleans_to_slice(arr[:0])
    assert(result == slice(0, 0))


def test_convert_objects():
    arr = np.array(['a', 'b', nan, nan, 'd', 'e', 'f'], dtype='O')
    result = lib.maybe_convert_objects(arr)
    assert(result.dtype == np.object_)


def test_convert_infs():
    arr = np.array(['inf', 'inf', 'inf'], dtype='O')
    result = lib.maybe_convert_numeric(arr, set(), False)
    assert(result.dtype == np.float64)

    arr = np.array(['-inf', '-inf', '-inf'], dtype='O')
    result = lib.maybe_convert_numeric(arr, set(), False)
    assert(result.dtype == np.float64)


def test_convert_objects_ints():
    # test that we can detect many kinds of integers
    dtypes = ['i1', 'i2', 'i4', 'i8', 'u1', 'u2', 'u4', 'u8']

    for dtype_str in dtypes:
        arr = np.array(list(np.arange(20, dtype=dtype_str)), dtype='O')
        assert(arr[0].dtype == np.dtype(dtype_str))
        result = lib.maybe_convert_objects(arr)
        assert(issubclass(result.dtype.type, np.integer))


def test_convert_objects_complex_number():
    for dtype in np.sctypes['complex']:
        arr = np.array(list(1j * np.arange(20, dtype=dtype)), dtype='O')
        assert(arr[0].dtype == np.dtype(dtype))
        result = lib.maybe_convert_objects(arr)
        assert(issubclass(result.dtype.type, np.complexfloating))


def test_rank():
    from pandas.compat.scipy import rankdata

    def _check(arr):
        mask = -np.isfinite(arr)
        arr = arr.copy()
        result = algos.rank_1d_float64(arr)
        arr[mask] = np.inf
        exp = rankdata(arr)
        exp[mask] = nan
        assert_almost_equal(result, exp)

    _check(np.array([nan, nan, 5., 5., 5., nan, 1, 2, 3, nan]))
    _check(np.array([4., nan, 5., 5., 5., nan, 1, 2, 4., nan]))


def test_get_reverse_indexer():
    indexer = np.array([-1, -1, 1, 2, 0, -1, 3, 4], dtype=np.int64)
    result = lib.get_reverse_indexer(indexer, 5)
    expected = np.array([4, 2, 3, 6, 7], dtype=np.int64)
    assert(np.array_equal(result, expected))


def test_pad_backfill_object_segfault():
    from datetime import datetime
    old = np.array([], dtype='O')
    new = np.array([datetime(2010, 12, 31)], dtype='O')

    result = algos.pad_object(old, new)
    expected = np.array([-1], dtype=np.int64)
    assert(np.array_equal(result, expected))

    result = algos.pad_object(new, old)
    expected = np.array([], dtype=np.int64)
    assert(np.array_equal(result, expected))

    result = algos.backfill_object(old, new)
    expected = np.array([-1], dtype=np.int64)
    assert(np.array_equal(result, expected))

    result = algos.backfill_object(new, old)
    expected = np.array([], dtype=np.int64)
    assert(np.array_equal(result, expected))


def test_arrmap():
    values = np.array(['foo', 'foo', 'bar', 'bar', 'baz', 'qux'], dtype='O')
    result = algos.arrmap_object(values, lambda x: x in ['foo', 'bar'])
    assert(result.dtype == np.bool_)


def test_series_grouper():
    from pandas import Series
    obj = Series(np.random.randn(10))
    dummy = obj[:0]

    labels = np.array([-1, -1, -1, 0, 0, 0, 1, 1, 1, 1], dtype=np.int64)

    grouper = lib.SeriesGrouper(obj, np.mean, labels, 2, dummy)
    result, counts = grouper.get_result()

    expected = np.array([obj[3:6].mean(), obj[6:].mean()])
    assert_almost_equal(result, expected)

    exp_counts = np.array([3, 4], dtype=np.int64)
    assert_almost_equal(counts, exp_counts)


def test_series_bin_grouper():
    from pandas import Series
    obj = Series(np.random.randn(10))
    dummy = obj[:0]

    bins = np.array([3, 6])

    grouper = lib.SeriesBinGrouper(obj, np.mean, bins, dummy)
    result, counts = grouper.get_result()

    expected = np.array([obj[:3].mean(), obj[3:6].mean(), obj[6:].mean()])
    assert_almost_equal(result, expected)

    exp_counts = np.array([3, 3, 4], dtype=np.int64)
    assert_almost_equal(counts, exp_counts)


class TestBinGroupers(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.obj = np.random.randn(10, 1)
        self.labels = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 2], dtype=np.int64)
        self.bins = np.array([3, 6], dtype=np.int64)

    def test_generate_bins(self):
        from pandas.core.groupby import generate_bins_generic
        values = np.array([1, 2, 3, 4, 5, 6], dtype=np.int64)
        binner = np.array([0, 3, 6, 9], dtype=np.int64)

        for func in [lib.generate_bins_dt64, generate_bins_generic]:
            bins = func(values, binner, closed='left')
            assert((bins == np.array([2, 5, 6])).all())

            bins = func(values, binner, closed='right')
            assert((bins == np.array([3, 6, 6])).all())

        for func in [lib.generate_bins_dt64, generate_bins_generic]:
            values = np.array([1, 2, 3, 4, 5, 6], dtype=np.int64)
            binner = np.array([0, 3, 6], dtype=np.int64)

            bins = func(values, binner, closed='right')
            assert((bins == np.array([3, 6])).all())

        self.assertRaises(ValueError, generate_bins_generic, values, [],
                          'right')
        self.assertRaises(ValueError, generate_bins_generic, values[:0],
                          binner, 'right')

        self.assertRaises(ValueError, generate_bins_generic,
                          values, [4], 'right')
        self.assertRaises(ValueError, generate_bins_generic,
                          values, [-3, -1], 'right')

    def test_group_bin_functions(self):

        dtypes = ['float32','float64']
        funcs  = ['add', 'mean', 'prod', 'min', 'max', 'var']

        np_funcs = {
            'add': np.sum,
            'mean': np.mean,
            'prod': np.prod,
            'min': np.min,
            'max': np.max,
            'var': lambda x: x.var(ddof=1) if len(x) >= 2 else np.nan
        }

        for fname in funcs:
            for d in dtypes:
                check_less_precise = False
                if d == 'float32':
                    check_less_precise = True
                args = [getattr(algos, 'group_%s_%s' % (fname,d)),
                        getattr(algos, 'group_%s_bin_%s' % (fname,d)),
                        np_funcs[fname],
                        d,
                        check_less_precise]
                self._check_versions(*args)

    def _check_versions(self, irr_func, bin_func, np_func, dtype, check_less_precise):
        obj = self.obj.astype(dtype)

        cts = np.zeros(3, dtype=np.int64)
        exp = np.zeros((3, 1), dtype)
        irr_func(exp, cts, obj, self.labels)

        # bin-based version
        bins = np.array([3, 6], dtype=np.int64)
        out = np.zeros((3, 1), dtype)
        counts = np.zeros(len(out), dtype=np.int64)
        bin_func(out, counts, obj, bins)

        assert_almost_equal(out, exp, check_less_precise=check_less_precise)

        bins = np.array([3, 9, 10], dtype=np.int64)
        out = np.zeros((3, 1), dtype)
        counts = np.zeros(len(out), dtype=np.int64)
        bin_func(out, counts, obj, bins)
        exp = np.array([np_func(obj[:3]), np_func(obj[3:9]),
                        np_func(obj[9:])],
                       dtype=dtype)
        assert_almost_equal(out.squeeze(), exp, check_less_precise=check_less_precise)

        # duplicate bins
        bins = np.array([3, 6, 10, 10], dtype=np.int64)
        out = np.zeros((4, 1), dtype)
        counts = np.zeros(len(out), dtype=np.int64)
        bin_func(out, counts, obj, bins)
        exp = np.array([np_func(obj[:3]), np_func(obj[3:6]),
                        np_func(obj[6:10]), np.nan],
                       dtype=dtype)
        assert_almost_equal(out.squeeze(), exp, check_less_precise=check_less_precise)


def test_group_ohlc():

    def _check(dtype):
        obj = np.array(np.random.randn(20),dtype=dtype)

        bins = np.array([6, 12], dtype=np.int64)
        out = np.zeros((3, 4), dtype)
        counts = np.zeros(len(out), dtype=np.int64)
        
        func = getattr(algos,'group_ohlc_%s' % dtype)
        func(out, counts, obj[:, None], bins)

        def _ohlc(group):
            if isnull(group).all():
                return np.repeat(nan, 4)
            return [group[0], group.max(), group.min(), group[-1]]

        expected = np.array([_ohlc(obj[:6]), _ohlc(obj[6:12]),
                             _ohlc(obj[12:])])

        assert_almost_equal(out, expected)
        assert_almost_equal(counts, [6, 6, 8])

        obj[:6] = nan
        func(out, counts, obj[:, None], bins)
        expected[0] = nan
        assert_almost_equal(out, expected)

    _check('float32')
    _check('float64')

def test_try_parse_dates():
    from dateutil.parser import parse

    arr = np.array(['5/1/2000', '6/1/2000', '7/1/2000'], dtype=object)

    result = lib.try_parse_dates(arr, dayfirst=True)
    expected = [parse(d, dayfirst=True) for d in arr]
    assert(np.array_equal(result, expected))


class TestTypeInference(unittest.TestCase):
    _multiprocess_can_split_ = True

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
        self.assertEqual(result, 'mixed-integer')

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
        self.assertEqual(result, 'mixed-integer')

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
        import datetime
        dates = [datetime.datetime(2012, 1, x) for x in range(1, 20)]
        index = Index(dates)
        self.assert_(index.inferred_type == 'datetime64')

    def test_date(self):
        import datetime
        dates = [datetime.date(2012, 1, x) for x in range(1, 20)]
        index = Index(dates)
        self.assert_(index.inferred_type == 'date')

    def test_to_object_array_tuples(self):
        r = (5, 6)
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


class TestReducer(unittest.TestCase):

    def test_int_index(self):
        from pandas.core.series import Series

        arr = np.random.randn(100, 4)

        result = lib.reduce(arr, np.sum, labels=Index(np.arange(4)))
        expected = arr.sum(0)
        assert_almost_equal(result, expected)

        result = lib.reduce(arr, np.sum, axis=1, labels=Index(np.arange(100)))
        expected = arr.sum(1)
        assert_almost_equal(result, expected)

        dummy = Series(0., index=np.arange(100))
        result = lib.reduce(
            arr, np.sum, dummy=dummy, labels=Index(np.arange(4)))
        expected = arr.sum(0)
        assert_almost_equal(result, expected)

        dummy = Series(0., index=np.arange(4))
        result = lib.reduce(arr, np.sum, axis=1,
                            dummy=dummy, labels=Index(np.arange(100)))
        expected = arr.sum(1)
        assert_almost_equal(result, expected)


class TestTsUtil(unittest.TestCase):
    def test_min_valid(self):
        # Ensure that Timestamp.min is a valid Timestamp
        Timestamp(Timestamp.min)

    def test_max_valid(self):
        # Ensure that Timestamp.max is a valid Timestamp
        Timestamp(Timestamp.max)

    def test_to_datetime_bijective(self):
        # Ensure that converting to datetime and back only loses precision
        # by going from nanoseconds to microseconds.
        self.assertEqual(Timestamp(Timestamp.max.to_pydatetime()).value/1000, Timestamp.max.value/1000)
        self.assertEqual(Timestamp(Timestamp.min.to_pydatetime()).value/1000, Timestamp.min.value/1000)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
