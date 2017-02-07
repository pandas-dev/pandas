# -*- coding: utf-8 -*-
from pandas.compat import range

import numpy as np
from numpy.random import RandomState
from numpy import nan
from datetime import datetime
from itertools import permutations
from pandas import Series, Categorical, CategoricalIndex, Index
import pandas as pd

from pandas import compat
import pandas.algos as _algos
from pandas.compat import lrange
import pandas.core.algorithms as algos
import pandas.util.testing as tm
import pandas.hashtable as hashtable
from pandas.compat.numpy import np_array_datetime64_compat
from pandas.util.testing import assert_almost_equal


class TestMatch(tm.TestCase):

    def test_ints(self):
        values = np.array([0, 2, 1])
        to_match = np.array([0, 1, 2, 2, 0, 1, 3, 0])

        result = algos.match(to_match, values)
        expected = np.array([0, 2, 1, 1, 0, 2, -1, 0], dtype=np.int64)
        self.assert_numpy_array_equal(result, expected)

        result = Series(algos.match(to_match, values, np.nan))
        expected = Series(np.array([0, 2, 1, 1, 0, 2, np.nan, 0]))
        tm.assert_series_equal(result, expected)

        s = pd.Series(np.arange(5), dtype=np.float32)
        result = algos.match(s, [2, 4])
        expected = np.array([-1, -1, 0, -1, 1], dtype=np.int64)
        self.assert_numpy_array_equal(result, expected)

        result = Series(algos.match(s, [2, 4], np.nan))
        expected = Series(np.array([np.nan, np.nan, 0, np.nan, 1]))
        tm.assert_series_equal(result, expected)

    def test_strings(self):
        values = ['foo', 'bar', 'baz']
        to_match = ['bar', 'foo', 'qux', 'foo', 'bar', 'baz', 'qux']

        result = algos.match(to_match, values)
        expected = np.array([1, 0, -1, 0, 1, 2, -1], dtype=np.int64)
        self.assert_numpy_array_equal(result, expected)

        result = Series(algos.match(to_match, values, np.nan))
        expected = Series(np.array([1, 0, np.nan, 0, 1, 2, np.nan]))
        tm.assert_series_equal(result, expected)


class TestSafeSort(tm.TestCase):

    def test_basic_sort(self):
        values = [3, 1, 2, 0, 4]
        result = algos.safe_sort(values)
        expected = np.array([0, 1, 2, 3, 4])
        tm.assert_numpy_array_equal(result, expected)

        values = list("baaacb")
        result = algos.safe_sort(values)
        expected = np.array(list("aaabbc"))
        tm.assert_numpy_array_equal(result, expected)

        values = []
        result = algos.safe_sort(values)
        expected = np.array([])
        tm.assert_numpy_array_equal(result, expected)

    def test_labels(self):
        values = [3, 1, 2, 0, 4]
        expected = np.array([0, 1, 2, 3, 4])

        labels = [0, 1, 1, 2, 3, 0, -1, 4]
        result, result_labels = algos.safe_sort(values, labels)
        expected_labels = np.array([3, 1, 1, 2, 0, 3, -1, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(result_labels, expected_labels)

        # na_sentinel
        labels = [0, 1, 1, 2, 3, 0, 99, 4]
        result, result_labels = algos.safe_sort(values, labels,
                                                na_sentinel=99)
        expected_labels = np.array([3, 1, 1, 2, 0, 3, 99, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(result_labels, expected_labels)

        # out of bound indices
        labels = [0, 101, 102, 2, 3, 0, 99, 4]
        result, result_labels = algos.safe_sort(values, labels)
        expected_labels = np.array([3, -1, -1, 2, 0, 3, -1, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(result_labels, expected_labels)

        labels = []
        result, result_labels = algos.safe_sort(values, labels)
        expected_labels = np.array([], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(result_labels, expected_labels)

    def test_mixed_integer(self):
        values = np.array(['b', 1, 0, 'a', 0, 'b'], dtype=object)
        result = algos.safe_sort(values)
        expected = np.array([0, 0, 1, 'a', 'b', 'b'], dtype=object)
        tm.assert_numpy_array_equal(result, expected)

        values = np.array(['b', 1, 0, 'a'], dtype=object)
        labels = [0, 1, 2, 3, 0, -1, 1]
        result, result_labels = algos.safe_sort(values, labels)
        expected = np.array([0, 1, 'a', 'b'], dtype=object)
        expected_labels = np.array([3, 1, 0, 2, 3, -1, 1], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(result_labels, expected_labels)

    def test_unsortable(self):
        # GH 13714
        arr = np.array([1, 2, datetime.now(), 0, 3], dtype=object)
        if compat.PY2 and not pd._np_version_under1p10:
            # RuntimeWarning: tp_compare didn't return -1 or -2 for exception
            with tm.assert_produces_warning(RuntimeWarning):
                tm.assertRaises(TypeError, algos.safe_sort, arr)
        else:
            tm.assertRaises(TypeError, algos.safe_sort, arr)

    def test_exceptions(self):
        with tm.assertRaisesRegexp(TypeError,
                                   "Only list-like objects are allowed"):
            algos.safe_sort(values=1)

        with tm.assertRaisesRegexp(TypeError,
                                   "Only list-like objects or None"):
            algos.safe_sort(values=[0, 1, 2], labels=1)

        with tm.assertRaisesRegexp(ValueError, "values should be unique"):
            algos.safe_sort(values=[0, 1, 2, 1], labels=[0, 1])


class TestFactorize(tm.TestCase):

    def test_basic(self):

        labels, uniques = algos.factorize(['a', 'b', 'b', 'a', 'a', 'c', 'c',
                                           'c'])
        self.assert_numpy_array_equal(
            uniques, np.array(['a', 'b', 'c'], dtype=object))

        labels, uniques = algos.factorize(['a', 'b', 'b', 'a',
                                           'a', 'c', 'c', 'c'], sort=True)
        exp = np.array([0, 1, 1, 0, 0, 2, 2, 2], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        exp = np.array(['a', 'b', 'c'], dtype=object)
        self.assert_numpy_array_equal(uniques, exp)

        labels, uniques = algos.factorize(list(reversed(range(5))))
        exp = np.array([0, 1, 2, 3, 4], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        exp = np.array([4, 3, 2, 1, 0], dtype=np.int64)
        self.assert_numpy_array_equal(uniques, exp)

        labels, uniques = algos.factorize(list(reversed(range(5))), sort=True)

        exp = np.array([4, 3, 2, 1, 0], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        exp = np.array([0, 1, 2, 3, 4], dtype=np.int64)
        self.assert_numpy_array_equal(uniques, exp)

        labels, uniques = algos.factorize(list(reversed(np.arange(5.))))
        exp = np.array([0, 1, 2, 3, 4], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        exp = np.array([4., 3., 2., 1., 0.], dtype=np.float64)
        self.assert_numpy_array_equal(uniques, exp)

        labels, uniques = algos.factorize(list(reversed(np.arange(5.))),
                                          sort=True)
        exp = np.array([4, 3, 2, 1, 0], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        exp = np.array([0., 1., 2., 3., 4.], dtype=np.float64)
        self.assert_numpy_array_equal(uniques, exp)

    def test_mixed(self):

        # doc example reshaping.rst
        x = Series(['A', 'A', np.nan, 'B', 3.14, np.inf])
        labels, uniques = algos.factorize(x)

        exp = np.array([0, 0, -1, 1, 2, 3], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        exp = pd.Index(['A', 'B', 3.14, np.inf])
        tm.assert_index_equal(uniques, exp)

        labels, uniques = algos.factorize(x, sort=True)
        exp = np.array([2, 2, -1, 3, 0, 1], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        exp = pd.Index([3.14, np.inf, 'A', 'B'])
        tm.assert_index_equal(uniques, exp)

    def test_datelike(self):

        # M8
        v1 = pd.Timestamp('20130101 09:00:00.00004')
        v2 = pd.Timestamp('20130101')
        x = Series([v1, v1, v1, v2, v2, v1])
        labels, uniques = algos.factorize(x)

        exp = np.array([0, 0, 0, 1, 1, 0], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        exp = pd.DatetimeIndex([v1, v2])
        self.assert_index_equal(uniques, exp)

        labels, uniques = algos.factorize(x, sort=True)
        exp = np.array([1, 1, 1, 0, 0, 1], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        exp = pd.DatetimeIndex([v2, v1])
        self.assert_index_equal(uniques, exp)

        # period
        v1 = pd.Period('201302', freq='M')
        v2 = pd.Period('201303', freq='M')
        x = Series([v1, v1, v1, v2, v2, v1])

        # periods are not 'sorted' as they are converted back into an index
        labels, uniques = algos.factorize(x)
        exp = np.array([0, 0, 0, 1, 1, 0], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        self.assert_index_equal(uniques, pd.PeriodIndex([v1, v2]))

        labels, uniques = algos.factorize(x, sort=True)
        exp = np.array([0, 0, 0, 1, 1, 0], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        self.assert_index_equal(uniques, pd.PeriodIndex([v1, v2]))

        # GH 5986
        v1 = pd.to_timedelta('1 day 1 min')
        v2 = pd.to_timedelta('1 day')
        x = Series([v1, v2, v1, v1, v2, v2, v1])
        labels, uniques = algos.factorize(x)
        exp = np.array([0, 1, 0, 0, 1, 1, 0], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        self.assert_index_equal(uniques, pd.to_timedelta([v1, v2]))

        labels, uniques = algos.factorize(x, sort=True)
        exp = np.array([1, 0, 1, 1, 0, 0, 1], dtype=np.intp)
        self.assert_numpy_array_equal(labels, exp)
        self.assert_index_equal(uniques, pd.to_timedelta([v2, v1]))

    def test_factorize_nan(self):
        # nan should map to na_sentinel, not reverse_indexer[na_sentinel]
        # rizer.factorize should not raise an exception if na_sentinel indexes
        # outside of reverse_indexer
        key = np.array([1, 2, 1, np.nan], dtype='O')
        rizer = hashtable.Factorizer(len(key))
        for na_sentinel in (-1, 20):
            ids = rizer.factorize(key, sort=True, na_sentinel=na_sentinel)
            expected = np.array([0, 1, 0, na_sentinel], dtype='int32')
            self.assertEqual(len(set(key)), len(set(expected)))
            self.assertTrue(np.array_equal(
                pd.isnull(key), expected == na_sentinel))

        # nan still maps to na_sentinel when sort=False
        key = np.array([0, np.nan, 1], dtype='O')
        na_sentinel = -1

        # TODO(wesm): unused?
        ids = rizer.factorize(key, sort=False, na_sentinel=na_sentinel)  # noqa

        expected = np.array([2, -1, 0], dtype='int32')
        self.assertEqual(len(set(key)), len(set(expected)))
        self.assertTrue(
            np.array_equal(pd.isnull(key), expected == na_sentinel))

    def test_complex_sorting(self):
        # gh 12666 - check no segfault
        # Test not valid numpy versions older than 1.11
        if pd._np_version_under1p11:
            self.skipTest("Test valid only for numpy 1.11+")

        x17 = np.array([complex(i) for i in range(17)], dtype=object)

        self.assertRaises(TypeError, algos.factorize, x17[::-1], sort=True)

    def test_uint64_factorize(self):
        data = np.array([2**63, 1, 2**63], dtype=np.uint64)
        exp_labels = np.array([0, 1, 0], dtype=np.intp)
        exp_uniques = np.array([2**63, 1], dtype=np.uint64)

        labels, uniques = algos.factorize(data)
        tm.assert_numpy_array_equal(labels, exp_labels)
        tm.assert_numpy_array_equal(uniques, exp_uniques)

        data = np.array([2**63, -1, 2**63], dtype=object)
        exp_labels = np.array([0, 1, 0], dtype=np.intp)
        exp_uniques = np.array([2**63, -1], dtype=object)

        labels, uniques = algos.factorize(data)
        tm.assert_numpy_array_equal(labels, exp_labels)
        tm.assert_numpy_array_equal(uniques, exp_uniques)


class TestUnique(tm.TestCase):

    def test_ints(self):
        arr = np.random.randint(0, 100, size=50)

        result = algos.unique(arr)
        tm.assertIsInstance(result, np.ndarray)

    def test_objects(self):
        arr = np.random.randint(0, 100, size=50).astype('O')

        result = algos.unique(arr)
        tm.assertIsInstance(result, np.ndarray)

    def test_object_refcount_bug(self):
        lst = ['A', 'B', 'C', 'D', 'E']
        for i in range(1000):
            len(algos.unique(lst))

    def test_on_index_object(self):

        mindex = pd.MultiIndex.from_arrays([np.arange(5).repeat(5), np.tile(
            np.arange(5), 5)])
        expected = mindex.values
        expected.sort()

        mindex = mindex.repeat(2)

        result = pd.unique(mindex)
        result.sort()

        tm.assert_almost_equal(result, expected)

    def test_datetime64_dtype_array_returned(self):
        # GH 9431
        expected = np_array_datetime64_compat(
            ['2015-01-03T00:00:00.000000000+0000',
             '2015-01-01T00:00:00.000000000+0000'],
            dtype='M8[ns]')

        dt_index = pd.to_datetime(['2015-01-03T00:00:00.000000000+0000',
                                   '2015-01-01T00:00:00.000000000+0000',
                                   '2015-01-01T00:00:00.000000000+0000'])
        result = algos.unique(dt_index)
        tm.assert_numpy_array_equal(result, expected)
        self.assertEqual(result.dtype, expected.dtype)

        s = pd.Series(dt_index)
        result = algos.unique(s)
        tm.assert_numpy_array_equal(result, expected)
        self.assertEqual(result.dtype, expected.dtype)

        arr = s.values
        result = algos.unique(arr)
        tm.assert_numpy_array_equal(result, expected)
        self.assertEqual(result.dtype, expected.dtype)

    def test_timedelta64_dtype_array_returned(self):
        # GH 9431
        expected = np.array([31200, 45678, 10000], dtype='m8[ns]')

        td_index = pd.to_timedelta([31200, 45678, 31200, 10000, 45678])
        result = algos.unique(td_index)
        tm.assert_numpy_array_equal(result, expected)
        self.assertEqual(result.dtype, expected.dtype)

        s = pd.Series(td_index)
        result = algos.unique(s)
        tm.assert_numpy_array_equal(result, expected)
        self.assertEqual(result.dtype, expected.dtype)

        arr = s.values
        result = algos.unique(arr)
        tm.assert_numpy_array_equal(result, expected)
        self.assertEqual(result.dtype, expected.dtype)

    def test_uint64_overflow(self):
        s = pd.Series([1, 2, 2**63, 2**63], dtype=np.uint64)
        exp = np.array([1, 2, 2**63], dtype=np.uint64)
        tm.assert_numpy_array_equal(algos.unique(s), exp)


class TestIsin(tm.TestCase):

    def test_invalid(self):

        self.assertRaises(TypeError, lambda: algos.isin(1, 1))
        self.assertRaises(TypeError, lambda: algos.isin(1, [1]))
        self.assertRaises(TypeError, lambda: algos.isin([1], 1))

    def test_basic(self):

        result = algos.isin([1, 2], [1])
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(np.array([1, 2]), [1])
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(pd.Series([1, 2]), [1])
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(pd.Series([1, 2]), pd.Series([1]))
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(pd.Series([1, 2]), set([1]))
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(['a', 'b'], ['a'])
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(pd.Series(['a', 'b']), pd.Series(['a']))
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(pd.Series(['a', 'b']), set(['a']))
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(['a', 'b'], [1])
        expected = np.array([False, False])
        tm.assert_numpy_array_equal(result, expected)

        arr = pd.date_range('20130101', periods=3).values
        result = algos.isin(arr, [arr[0]])
        expected = np.array([True, False, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(arr, arr[0:2])
        expected = np.array([True, True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(arr, set(arr[0:2]))
        expected = np.array([True, True, False])
        tm.assert_numpy_array_equal(result, expected)

        arr = pd.timedelta_range('1 day', periods=3).values
        result = algos.isin(arr, [arr[0]])
        expected = np.array([True, False, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(arr, arr[0:2])
        expected = np.array([True, True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.isin(arr, set(arr[0:2]))
        expected = np.array([True, True, False])
        tm.assert_numpy_array_equal(result, expected)

    def test_large(self):

        s = pd.date_range('20000101', periods=2000000, freq='s').values
        result = algos.isin(s, s[0:2])
        expected = np.zeros(len(s), dtype=bool)
        expected[0] = True
        expected[1] = True
        tm.assert_numpy_array_equal(result, expected)


class TestValueCounts(tm.TestCase):

    def test_value_counts(self):
        np.random.seed(1234)
        from pandas.tools.tile import cut

        arr = np.random.randn(4)
        factor = cut(arr, 4)

        tm.assertIsInstance(factor, Categorical)
        result = algos.value_counts(factor)
        cats = ['(-1.194, -0.535]', '(-0.535, 0.121]', '(0.121, 0.777]',
                '(0.777, 1.433]']
        expected_index = CategoricalIndex(cats, cats, ordered=True)
        expected = Series([1, 1, 1, 1], index=expected_index)
        tm.assert_series_equal(result.sort_index(), expected.sort_index())

    def test_value_counts_bins(self):
        s = [1, 2, 3, 4]
        result = algos.value_counts(s, bins=1)
        self.assertEqual(result.tolist(), [4])
        self.assertEqual(result.index[0], 0.997)

        result = algos.value_counts(s, bins=2, sort=False)
        self.assertEqual(result.tolist(), [2, 2])
        self.assertEqual(result.index[0], 0.997)
        self.assertEqual(result.index[1], 2.5)

    def test_value_counts_dtypes(self):
        result = algos.value_counts([1, 1.])
        self.assertEqual(len(result), 1)

        result = algos.value_counts([1, 1.], bins=1)
        self.assertEqual(len(result), 1)

        result = algos.value_counts(Series([1, 1., '1']))  # object
        self.assertEqual(len(result), 2)

        self.assertRaises(TypeError, lambda s: algos.value_counts(s, bins=1),
                          ['1', 1])

    def test_value_counts_nat(self):
        td = Series([np.timedelta64(10000), pd.NaT], dtype='timedelta64[ns]')
        dt = pd.to_datetime(['NaT', '2014-01-01'])

        for s in [td, dt]:
            vc = algos.value_counts(s)
            vc_with_na = algos.value_counts(s, dropna=False)
            self.assertEqual(len(vc), 1)
            self.assertEqual(len(vc_with_na), 2)

        exp_dt = pd.Series({pd.Timestamp('2014-01-01 00:00:00'): 1})
        tm.assert_series_equal(algos.value_counts(dt), exp_dt)
        # TODO same for (timedelta)

    def test_value_counts_datetime_outofbounds(self):
        # GH 13663
        s = pd.Series([datetime(3000, 1, 1), datetime(5000, 1, 1),
                       datetime(5000, 1, 1), datetime(6000, 1, 1),
                       datetime(3000, 1, 1), datetime(3000, 1, 1)])
        res = s.value_counts()

        exp_index = pd.Index([datetime(3000, 1, 1), datetime(5000, 1, 1),
                              datetime(6000, 1, 1)], dtype=object)
        exp = pd.Series([3, 2, 1], index=exp_index)
        tm.assert_series_equal(res, exp)

        # GH 12424
        res = pd.to_datetime(pd.Series(['2362-01-01', np.nan]),
                             errors='ignore')
        exp = pd.Series(['2362-01-01', np.nan], dtype=object)
        tm.assert_series_equal(res, exp)

    def test_categorical(self):
        s = Series(pd.Categorical(list('aaabbc')))
        result = s.value_counts()
        expected = pd.Series([3, 2, 1],
                             index=pd.CategoricalIndex(['a', 'b', 'c']))
        tm.assert_series_equal(result, expected, check_index_type=True)

        # preserve order?
        s = s.cat.as_ordered()
        result = s.value_counts()
        expected.index = expected.index.as_ordered()
        tm.assert_series_equal(result, expected, check_index_type=True)

    def test_categorical_nans(self):
        s = Series(pd.Categorical(list('aaaaabbbcc')))  # 4,3,2,1 (nan)
        s.iloc[1] = np.nan
        result = s.value_counts()
        expected = pd.Series([4, 3, 2], index=pd.CategoricalIndex(
            ['a', 'b', 'c'], categories=['a', 'b', 'c']))
        tm.assert_series_equal(result, expected, check_index_type=True)
        result = s.value_counts(dropna=False)
        expected = pd.Series([
            4, 3, 2, 1
        ], index=pd.CategoricalIndex(['a', 'b', 'c', np.nan]))
        tm.assert_series_equal(result, expected, check_index_type=True)

        # out of order
        s = Series(pd.Categorical(
            list('aaaaabbbcc'), ordered=True, categories=['b', 'a', 'c']))
        s.iloc[1] = np.nan
        result = s.value_counts()
        expected = pd.Series([4, 3, 2], index=pd.CategoricalIndex(
            ['a', 'b', 'c'], categories=['b', 'a', 'c'], ordered=True))
        tm.assert_series_equal(result, expected, check_index_type=True)

        result = s.value_counts(dropna=False)
        expected = pd.Series([4, 3, 2, 1], index=pd.CategoricalIndex(
            ['a', 'b', 'c', np.nan], categories=['b', 'a', 'c'], ordered=True))
        tm.assert_series_equal(result, expected, check_index_type=True)

    def test_categorical_zeroes(self):
        # keep the `d` category with 0
        s = Series(pd.Categorical(
            list('bbbaac'), categories=list('abcd'), ordered=True))
        result = s.value_counts()
        expected = Series([3, 2, 1, 0], index=pd.Categorical(
            ['b', 'a', 'c', 'd'], categories=list('abcd'), ordered=True))
        tm.assert_series_equal(result, expected, check_index_type=True)

    def test_dropna(self):
        # https://github.com/pandas-dev/pandas/issues/9443#issuecomment-73719328

        tm.assert_series_equal(
            pd.Series([True, True, False]).value_counts(dropna=True),
            pd.Series([2, 1], index=[True, False]))
        tm.assert_series_equal(
            pd.Series([True, True, False]).value_counts(dropna=False),
            pd.Series([2, 1], index=[True, False]))

        tm.assert_series_equal(
            pd.Series([True, True, False, None]).value_counts(dropna=True),
            pd.Series([2, 1], index=[True, False]))
        tm.assert_series_equal(
            pd.Series([True, True, False, None]).value_counts(dropna=False),
            pd.Series([2, 1, 1], index=[True, False, np.nan]))
        tm.assert_series_equal(
            pd.Series([10.3, 5., 5.]).value_counts(dropna=True),
            pd.Series([2, 1], index=[5., 10.3]))
        tm.assert_series_equal(
            pd.Series([10.3, 5., 5.]).value_counts(dropna=False),
            pd.Series([2, 1], index=[5., 10.3]))

        tm.assert_series_equal(
            pd.Series([10.3, 5., 5., None]).value_counts(dropna=True),
            pd.Series([2, 1], index=[5., 10.3]))

        # 32-bit linux has a different ordering
        if not compat.is_platform_32bit():
            tm.assert_series_equal(
                pd.Series([10.3, 5., 5., None]).value_counts(dropna=False),
                pd.Series([2, 1, 1], index=[5., 10.3, np.nan]))

    def test_value_counts_normalized(self):
        # GH12558
        s = Series([1, 2, np.nan, np.nan, np.nan])
        dtypes = (np.float64, np.object, 'M8[ns]')
        for t in dtypes:
            s_typed = s.astype(t)
            result = s_typed.value_counts(normalize=True, dropna=False)
            expected = Series([0.6, 0.2, 0.2],
                              index=Series([np.nan, 2.0, 1.0], dtype=t))
            tm.assert_series_equal(result, expected)

            result = s_typed.value_counts(normalize=True, dropna=True)
            expected = Series([0.5, 0.5],
                              index=Series([2.0, 1.0], dtype=t))
            tm.assert_series_equal(result, expected)

    def test_value_counts_uint64(self):
        arr = np.array([2**63], dtype=np.uint64)
        expected = Series([1], index=[2**63])
        result = algos.value_counts(arr)

        tm.assert_series_equal(result, expected)

        arr = np.array([-1, 2**63], dtype=object)
        expected = Series([1, 1], index=[-1, 2**63])
        result = algos.value_counts(arr)

        tm.assert_series_equal(result, expected)


class TestDuplicated(tm.TestCase):

    def test_duplicated_with_nas(self):
        keys = np.array([0, 1, np.nan, 0, 2, np.nan], dtype=object)

        result = algos.duplicated(keys)
        expected = np.array([False, False, False, True, False, True])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.duplicated(keys, keep='first')
        expected = np.array([False, False, False, True, False, True])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.duplicated(keys, keep='last')
        expected = np.array([True, False, True, False, False, False])
        tm.assert_numpy_array_equal(result, expected)

        result = algos.duplicated(keys, keep=False)
        expected = np.array([True, False, True, True, False, True])
        tm.assert_numpy_array_equal(result, expected)

        keys = np.empty(8, dtype=object)
        for i, t in enumerate(zip([0, 0, np.nan, np.nan] * 2,
                                  [0, np.nan, 0, np.nan] * 2)):
            keys[i] = t

        result = algos.duplicated(keys)
        falses = [False] * 4
        trues = [True] * 4
        expected = np.array(falses + trues)
        tm.assert_numpy_array_equal(result, expected)

        result = algos.duplicated(keys, keep='last')
        expected = np.array(trues + falses)
        tm.assert_numpy_array_equal(result, expected)

        result = algos.duplicated(keys, keep=False)
        expected = np.array(trues + trues)
        tm.assert_numpy_array_equal(result, expected)

    def test_numeric_object_likes(self):
        cases = [np.array([1, 2, 1, 5, 3,
                           2, 4, 1, 5, 6]),
                 np.array([1.1, 2.2, 1.1, np.nan, 3.3,
                           2.2, 4.4, 1.1, np.nan, 6.6]),
                 np.array([1 + 1j, 2 + 2j, 1 + 1j, 5 + 5j, 3 + 3j,
                           2 + 2j, 4 + 4j, 1 + 1j, 5 + 5j, 6 + 6j]),
                 np.array(['a', 'b', 'a', 'e', 'c',
                           'b', 'd', 'a', 'e', 'f'], dtype=object),
                 np.array([1, 2**63, 1, 3**5, 10,
                           2**63, 39, 1, 3**5, 7], dtype=np.uint64)]

        exp_first = np.array([False, False, True, False, False,
                              True, False, True, True, False])
        exp_last = np.array([True, True, True, True, False,
                             False, False, False, False, False])
        exp_false = exp_first | exp_last

        for case in cases:
            res_first = algos.duplicated(case, keep='first')
            tm.assert_numpy_array_equal(res_first, exp_first)

            res_last = algos.duplicated(case, keep='last')
            tm.assert_numpy_array_equal(res_last, exp_last)

            res_false = algos.duplicated(case, keep=False)
            tm.assert_numpy_array_equal(res_false, exp_false)

            # index
            for idx in [pd.Index(case), pd.Index(case, dtype='category')]:
                res_first = idx.duplicated(keep='first')
                tm.assert_numpy_array_equal(res_first, exp_first)

                res_last = idx.duplicated(keep='last')
                tm.assert_numpy_array_equal(res_last, exp_last)

                res_false = idx.duplicated(keep=False)
                tm.assert_numpy_array_equal(res_false, exp_false)

            # series
            for s in [pd.Series(case), pd.Series(case, dtype='category')]:
                res_first = s.duplicated(keep='first')
                tm.assert_series_equal(res_first, pd.Series(exp_first))

                res_last = s.duplicated(keep='last')
                tm.assert_series_equal(res_last, pd.Series(exp_last))

                res_false = s.duplicated(keep=False)
                tm.assert_series_equal(res_false, pd.Series(exp_false))

    def test_datetime_likes(self):

        dt = ['2011-01-01', '2011-01-02', '2011-01-01', 'NaT', '2011-01-03',
              '2011-01-02', '2011-01-04', '2011-01-01', 'NaT', '2011-01-06']
        td = ['1 days', '2 days', '1 days', 'NaT', '3 days',
              '2 days', '4 days', '1 days', 'NaT', '6 days']

        cases = [np.array([pd.Timestamp(d) for d in dt]),
                 np.array([pd.Timestamp(d, tz='US/Eastern') for d in dt]),
                 np.array([pd.Period(d, freq='D') for d in dt]),
                 np.array([np.datetime64(d) for d in dt]),
                 np.array([pd.Timedelta(d) for d in td])]

        exp_first = np.array([False, False, True, False, False,
                              True, False, True, True, False])
        exp_last = np.array([True, True, True, True, False,
                             False, False, False, False, False])
        exp_false = exp_first | exp_last

        for case in cases:
            res_first = algos.duplicated(case, keep='first')
            tm.assert_numpy_array_equal(res_first, exp_first)

            res_last = algos.duplicated(case, keep='last')
            tm.assert_numpy_array_equal(res_last, exp_last)

            res_false = algos.duplicated(case, keep=False)
            tm.assert_numpy_array_equal(res_false, exp_false)

            # index
            for idx in [pd.Index(case), pd.Index(case, dtype='category'),
                        pd.Index(case, dtype=object)]:
                res_first = idx.duplicated(keep='first')
                tm.assert_numpy_array_equal(res_first, exp_first)

                res_last = idx.duplicated(keep='last')
                tm.assert_numpy_array_equal(res_last, exp_last)

                res_false = idx.duplicated(keep=False)
                tm.assert_numpy_array_equal(res_false, exp_false)

            # series
            for s in [pd.Series(case), pd.Series(case, dtype='category'),
                      pd.Series(case, dtype=object)]:
                res_first = s.duplicated(keep='first')
                tm.assert_series_equal(res_first, pd.Series(exp_first))

                res_last = s.duplicated(keep='last')
                tm.assert_series_equal(res_last, pd.Series(exp_last))

                res_false = s.duplicated(keep=False)
                tm.assert_series_equal(res_false, pd.Series(exp_false))

    def test_unique_index(self):
        cases = [pd.Index([1, 2, 3]), pd.RangeIndex(0, 3)]
        for case in cases:
            self.assertTrue(case.is_unique)
            tm.assert_numpy_array_equal(case.duplicated(),
                                        np.array([False, False, False]))


class GroupVarTestMixin(object):

    def test_group_var_generic_1d(self):
        prng = RandomState(1234)

        out = (np.nan * np.ones((5, 1))).astype(self.dtype)
        counts = np.zeros(5, dtype='int64')
        values = 10 * prng.rand(15, 1).astype(self.dtype)
        labels = np.tile(np.arange(5), (3, )).astype('int64')

        expected_out = (np.squeeze(values)
                        .reshape((5, 3), order='F')
                        .std(axis=1, ddof=1) ** 2)[:, np.newaxis]
        expected_counts = counts + 3

        self.algo(out, counts, values, labels)
        self.assertTrue(np.allclose(out, expected_out, self.rtol))
        tm.assert_numpy_array_equal(counts, expected_counts)

    def test_group_var_generic_1d_flat_labels(self):
        prng = RandomState(1234)

        out = (np.nan * np.ones((1, 1))).astype(self.dtype)
        counts = np.zeros(1, dtype='int64')
        values = 10 * prng.rand(5, 1).astype(self.dtype)
        labels = np.zeros(5, dtype='int64')

        expected_out = np.array([[values.std(ddof=1) ** 2]])
        expected_counts = counts + 5

        self.algo(out, counts, values, labels)

        self.assertTrue(np.allclose(out, expected_out, self.rtol))
        tm.assert_numpy_array_equal(counts, expected_counts)

    def test_group_var_generic_2d_all_finite(self):
        prng = RandomState(1234)

        out = (np.nan * np.ones((5, 2))).astype(self.dtype)
        counts = np.zeros(5, dtype='int64')
        values = 10 * prng.rand(10, 2).astype(self.dtype)
        labels = np.tile(np.arange(5), (2, )).astype('int64')

        expected_out = np.std(values.reshape(2, 5, 2), ddof=1, axis=0) ** 2
        expected_counts = counts + 2

        self.algo(out, counts, values, labels)
        self.assertTrue(np.allclose(out, expected_out, self.rtol))
        tm.assert_numpy_array_equal(counts, expected_counts)

    def test_group_var_generic_2d_some_nan(self):
        prng = RandomState(1234)

        out = (np.nan * np.ones((5, 2))).astype(self.dtype)
        counts = np.zeros(5, dtype='int64')
        values = 10 * prng.rand(10, 2).astype(self.dtype)
        values[:, 1] = np.nan
        labels = np.tile(np.arange(5), (2, )).astype('int64')

        expected_out = np.vstack([values[:, 0]
                                  .reshape(5, 2, order='F')
                                  .std(ddof=1, axis=1) ** 2,
                                  np.nan * np.ones(5)]).T.astype(self.dtype)
        expected_counts = counts + 2

        self.algo(out, counts, values, labels)
        tm.assert_almost_equal(out, expected_out, check_less_precise=6)
        tm.assert_numpy_array_equal(counts, expected_counts)

    def test_group_var_constant(self):
        # Regression test from GH 10448.

        out = np.array([[np.nan]], dtype=self.dtype)
        counts = np.array([0], dtype='int64')
        values = 0.832845131556193 * np.ones((3, 1), dtype=self.dtype)
        labels = np.zeros(3, dtype='int64')

        self.algo(out, counts, values, labels)

        self.assertEqual(counts[0], 3)
        self.assertTrue(out[0, 0] >= 0)
        tm.assert_almost_equal(out[0, 0], 0.0)


class TestGroupVarFloat64(tm.TestCase, GroupVarTestMixin):
    __test__ = True

    algo = algos.algos.group_var_float64
    dtype = np.float64
    rtol = 1e-5

    def test_group_var_large_inputs(self):

        prng = RandomState(1234)

        out = np.array([[np.nan]], dtype=self.dtype)
        counts = np.array([0], dtype='int64')
        values = (prng.rand(10 ** 6) + 10 ** 12).astype(self.dtype)
        values.shape = (10 ** 6, 1)
        labels = np.zeros(10 ** 6, dtype='int64')

        self.algo(out, counts, values, labels)

        self.assertEqual(counts[0], 10 ** 6)
        tm.assert_almost_equal(out[0, 0], 1.0 / 12, check_less_precise=True)


class TestGroupVarFloat32(tm.TestCase, GroupVarTestMixin):
    __test__ = True

    algo = algos.algos.group_var_float32
    dtype = np.float32
    rtol = 1e-2


class TestHashTable(tm.TestCase):

    def test_lookup_nan(self):
        xs = np.array([2.718, 3.14, np.nan, -7, 5, 2, 3])
        m = hashtable.Float64HashTable()
        m.map_locations(xs)
        self.assert_numpy_array_equal(m.lookup(xs),
                                      np.arange(len(xs), dtype=np.int64))

    def test_lookup_overflow(self):
        xs = np.array([1, 2, 2**63], dtype=np.uint64)
        m = hashtable.UInt64HashTable()
        m.map_locations(xs)
        self.assert_numpy_array_equal(m.lookup(xs),
                                      np.arange(len(xs), dtype=np.int64))

    def test_get_unique(self):
        s = pd.Series([1, 2, 2**63, 2**63], dtype=np.uint64)
        exp = np.array([1, 2, 2**63], dtype=np.uint64)
        self.assert_numpy_array_equal(s.unique(), exp)

    def test_vector_resize(self):
        # Test for memory errors after internal vector
        # reallocations (pull request #7157)

        def _test_vector_resize(htable, uniques, dtype, nvals):
            vals = np.array(np.random.randn(1000), dtype=dtype)
            # get_labels appends to the vector
            htable.get_labels(vals[:nvals], uniques, 0, -1)
            # to_array resizes the vector
            uniques.to_array()
            htable.get_labels(vals, uniques, 0, -1)

        test_cases = [
            (hashtable.PyObjectHashTable, hashtable.ObjectVector, 'object'),
            (hashtable.StringHashTable, hashtable.ObjectVector, 'object'),
            (hashtable.Float64HashTable, hashtable.Float64Vector, 'float64'),
            (hashtable.Int64HashTable, hashtable.Int64Vector, 'int64'),
            (hashtable.UInt64HashTable, hashtable.UInt64Vector, 'uint64')]

        for (tbl, vect, dtype) in test_cases:
            # resizing to empty is a special case
            _test_vector_resize(tbl(), vect(), dtype, 0)
            _test_vector_resize(tbl(), vect(), dtype, 10)


def test_quantile():
    s = Series(np.random.randn(100))

    result = algos.quantile(s, [0, .25, .5, .75, 1.])
    expected = algos.quantile(s.values, [0, .25, .5, .75, 1.])
    tm.assert_almost_equal(result, expected)


def test_unique_label_indices():
    from pandas.hashtable import unique_label_indices

    a = np.random.randint(1, 1 << 10, 1 << 15).astype('i8')

    left = unique_label_indices(a)
    right = np.unique(a, return_index=True)[1]

    tm.assert_numpy_array_equal(left, right,
                                check_dtype=False)

    a[np.random.choice(len(a), 10)] = -1
    left = unique_label_indices(a)
    right = np.unique(a, return_index=True)[1][1:]
    tm.assert_numpy_array_equal(left, right,
                                check_dtype=False)


class TestRank(tm.TestCase):

    def test_scipy_compat(self):
        tm._skip_if_no_scipy()
        from scipy.stats import rankdata

        def _check(arr):
            mask = ~np.isfinite(arr)
            arr = arr.copy()
            result = _algos.rank_1d_float64(arr)
            arr[mask] = np.inf
            exp = rankdata(arr)
            exp[mask] = nan
            assert_almost_equal(result, exp)

        _check(np.array([nan, nan, 5., 5., 5., nan, 1, 2, 3, nan]))
        _check(np.array([4., nan, 5., 5., 5., nan, 1, 2, 4., nan]))

    def test_basic(self):
        exp = np.array([1, 2], dtype=np.float64)

        for dtype in np.typecodes['AllInteger']:
            s = Series([1, 100], dtype=dtype)
            tm.assert_numpy_array_equal(algos.rank(s), exp)

    def test_uint64_overflow(self):
        exp = np.array([1, 2], dtype=np.float64)

        for dtype in [np.float64, np.uint64]:
            s = Series([1, 2**63], dtype=dtype)
            tm.assert_numpy_array_equal(algos.rank(s), exp)

    def test_too_many_ndims(self):
        arr = np.array([[[1, 2, 3], [4, 5, 6], [7, 8, 9]]])
        msg = "Array with ndim > 2 are not supported"

        with tm.assertRaisesRegexp(TypeError, msg):
            algos.rank(arr)


def test_pad_backfill_object_segfault():

    old = np.array([], dtype='O')
    new = np.array([datetime(2010, 12, 31)], dtype='O')

    result = _algos.pad_object(old, new)
    expected = np.array([-1], dtype=np.int64)
    assert (np.array_equal(result, expected))

    result = _algos.pad_object(new, old)
    expected = np.array([], dtype=np.int64)
    assert (np.array_equal(result, expected))

    result = _algos.backfill_object(old, new)
    expected = np.array([-1], dtype=np.int64)
    assert (np.array_equal(result, expected))

    result = _algos.backfill_object(new, old)
    expected = np.array([], dtype=np.int64)
    assert (np.array_equal(result, expected))


def test_arrmap():
    values = np.array(['foo', 'foo', 'bar', 'bar', 'baz', 'qux'], dtype='O')
    result = _algos.arrmap_object(values, lambda x: x in ['foo', 'bar'])
    assert (result.dtype == np.bool_)


class TestTseriesUtil(tm.TestCase):

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
        new = Index(lrange(12))

        filler = _algos.backfill_int64(old.values, new.values)

        expect_filler = np.array([0, 0, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, -1], dtype=np.int64)
        self.assert_numpy_array_equal(filler, expect_filler)

        # corner case
        old = Index([1, 4])
        new = Index(lrange(5, 10))
        filler = _algos.backfill_int64(old.values, new.values)

        expect_filler = np.array([-1, -1, -1, -1, -1], dtype=np.int64)
        self.assert_numpy_array_equal(filler, expect_filler)

    def test_pad(self):
        old = Index([1, 5, 10])
        new = Index(lrange(12))

        filler = _algos.pad_int64(old.values, new.values)

        expect_filler = np.array([-1, 0, 0, 0, 0, 1,
                                  1, 1, 1, 1, 2, 2], dtype=np.int64)
        self.assert_numpy_array_equal(filler, expect_filler)

        # corner case
        old = Index([5, 10])
        new = Index(lrange(5))
        filler = _algos.pad_int64(old.values, new.values)
        expect_filler = np.array([-1, -1, -1, -1, -1], dtype=np.int64)
        self.assert_numpy_array_equal(filler, expect_filler)


def test_is_lexsorted():
    failure = [
        np.array([3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                  3, 3,
                  3, 3,
                  3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2,
                  2, 2, 2, 2, 2, 2, 2,
                  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                  1, 1, 1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1, 1, 1,
                  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0]),
        np.array([30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16,
                  15, 14,
                  13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 30, 29, 28,
                  27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13,
                  12, 11,
                  10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 30, 29, 28, 27, 26, 25,
                  24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10,
                  9, 8,
                  7, 6, 5, 4, 3, 2, 1, 0, 30, 29, 28, 27, 26, 25, 24, 23, 22,
                  21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7,
                  6, 5,
                  4, 3, 2, 1, 0])]

    assert (not _algos.is_lexsorted(failure))

# def test_get_group_index():
#     a = np.array([0, 1, 2, 0, 2, 1, 0, 0], dtype=np.int64)
#     b = np.array([1, 0, 3, 2, 0, 2, 3, 0], dtype=np.int64)
#     expected = np.array([1, 4, 11, 2, 8, 6, 3, 0], dtype=np.int64)

#     result = lib.get_group_index([a, b], (3, 4))

#     assert(np.array_equal(result, expected))


def test_groupsort_indexer():
    a = np.random.randint(0, 1000, 100).astype(np.int64)
    b = np.random.randint(0, 1000, 100).astype(np.int64)

    result = _algos.groupsort_indexer(a, 1000)[0]

    # need to use a stable sort
    expected = np.argsort(a, kind='mergesort')
    assert (np.array_equal(result, expected))

    # compare with lexsort
    key = a * 1000 + b
    result = _algos.groupsort_indexer(key, 1000000)[0]
    expected = np.lexsort((b, a))
    assert (np.array_equal(result, expected))


def test_infinity_sort():
    # GH 13445
    # numpy's argsort can be unhappy if something is less than
    # itself.  Instead, let's give our infinities a self-consistent
    # ordering, but outside the float extended real line.

    Inf = _algos.Infinity()
    NegInf = _algos.NegInfinity()

    ref_nums = [NegInf, float("-inf"), -1e100, 0, 1e100, float("inf"), Inf]

    assert all(Inf >= x for x in ref_nums)
    assert all(Inf > x or x is Inf for x in ref_nums)
    assert Inf >= Inf and Inf == Inf
    assert not Inf < Inf and not Inf > Inf

    assert all(NegInf <= x for x in ref_nums)
    assert all(NegInf < x or x is NegInf for x in ref_nums)
    assert NegInf <= NegInf and NegInf == NegInf
    assert not NegInf < NegInf and not NegInf > NegInf

    for perm in permutations(ref_nums):
        assert sorted(perm) == ref_nums

    # smoke tests
    np.array([_algos.Infinity()] * 32).argsort()
    np.array([_algos.NegInfinity()] * 32).argsort()


def test_ensure_platform_int():
    arr = np.arange(100, dtype=np.intp)

    result = _algos.ensure_platform_int(arr)
    assert (result is arr)


def test_int64_add_overflow():
    # see gh-14068
    msg = "Overflow in int64 addition"
    m = np.iinfo(np.int64).max
    n = np.iinfo(np.int64).min

    with tm.assertRaisesRegexp(OverflowError, msg):
        algos.checked_add_with_arr(np.array([m, m]), m)
    with tm.assertRaisesRegexp(OverflowError, msg):
        algos.checked_add_with_arr(np.array([m, m]), np.array([m, m]))
    with tm.assertRaisesRegexp(OverflowError, msg):
        algos.checked_add_with_arr(np.array([n, n]), n)
    with tm.assertRaisesRegexp(OverflowError, msg):
        algos.checked_add_with_arr(np.array([n, n]), np.array([n, n]))
    with tm.assertRaisesRegexp(OverflowError, msg):
        algos.checked_add_with_arr(np.array([m, n]), np.array([n, n]))
    with tm.assertRaisesRegexp(OverflowError, msg):
        algos.checked_add_with_arr(np.array([m, m]), np.array([m, m]),
                                   arr_mask=np.array([False, True]))
    with tm.assertRaisesRegexp(OverflowError, msg):
        algos.checked_add_with_arr(np.array([m, m]), np.array([m, m]),
                                   b_mask=np.array([False, True]))
    with tm.assertRaisesRegexp(OverflowError, msg):
        algos.checked_add_with_arr(np.array([m, m]), np.array([m, m]),
                                   arr_mask=np.array([False, True]),
                                   b_mask=np.array([False, True]))
    with tm.assertRaisesRegexp(OverflowError, msg):
        with tm.assert_produces_warning(RuntimeWarning):
            algos.checked_add_with_arr(np.array([m, m]),
                                       np.array([np.nan, m]))

    # Check that the nan boolean arrays override whether or not
    # the addition overflows. We don't check the result but just
    # the fact that an OverflowError is not raised.
    with tm.assertRaises(AssertionError):
        with tm.assertRaisesRegexp(OverflowError, msg):
            algos.checked_add_with_arr(np.array([m, m]), np.array([m, m]),
                                       arr_mask=np.array([True, True]))
    with tm.assertRaises(AssertionError):
        with tm.assertRaisesRegexp(OverflowError, msg):
            algos.checked_add_with_arr(np.array([m, m]), np.array([m, m]),
                                       b_mask=np.array([True, True]))
    with tm.assertRaises(AssertionError):
        with tm.assertRaisesRegexp(OverflowError, msg):
            algos.checked_add_with_arr(np.array([m, m]), np.array([m, m]),
                                       arr_mask=np.array([True, False]),
                                       b_mask=np.array([False, True]))


class TestMode(tm.TestCase):

    def test_no_mode(self):
        exp = Series([], dtype=np.float64)
        tm.assert_series_equal(algos.mode([]), exp)

        exp = Series([], dtype=np.int)
        tm.assert_series_equal(algos.mode([1]), exp)

        exp = Series([], dtype=np.object)
        tm.assert_series_equal(algos.mode(['a', 'b', 'c']), exp)

    def test_number_mode(self):
        exp_single = [1]
        data_single = [1] * 5 + [2] * 3

        exp_multi = [1, 3]
        data_multi = [1] * 5 + [2] * 3 + [3] * 5

        for dt in np.typecodes['AllInteger'] + np.typecodes['Float']:
            s = Series(data_single, dtype=dt)
            exp = Series(exp_single, dtype=dt)
            tm.assert_series_equal(algos.mode(s), exp)

            s = Series(data_multi, dtype=dt)
            exp = Series(exp_multi, dtype=dt)
            tm.assert_series_equal(algos.mode(s), exp)

    def test_strobj_mode(self):
        exp = ['b']
        data = ['a'] * 2 + ['b'] * 3

        s = Series(data, dtype='c')
        exp = Series(exp, dtype='c')
        tm.assert_series_equal(algos.mode(s), exp)

        exp = ['bar']
        data = ['foo'] * 2 + ['bar'] * 3

        for dt in [str, object]:
            s = Series(data, dtype=dt)
            exp = Series(exp, dtype=dt)
            tm.assert_series_equal(algos.mode(s), exp)

    def test_datelike_mode(self):
        exp = Series([], dtype="M8[ns]")
        s = Series(['2011-01-03', '2013-01-02',
                    '1900-05-03'], dtype='M8[ns]')
        tm.assert_series_equal(algos.mode(s), exp)

        exp = Series(['2011-01-03', '2013-01-02'], dtype='M8[ns]')
        s = Series(['2011-01-03', '2013-01-02', '1900-05-03',
                    '2011-01-03', '2013-01-02'], dtype='M8[ns]')
        tm.assert_series_equal(algos.mode(s), exp)

    def test_timedelta_mode(self):
        exp = Series([], dtype='timedelta64[ns]')
        s = Series(['1 days', '-1 days', '0 days'],
                   dtype='timedelta64[ns]')
        tm.assert_series_equal(algos.mode(s), exp)

        exp = Series(['2 min', '1 day'], dtype='timedelta64[ns]')
        s = Series(['1 day', '1 day', '-1 day', '-1 day 2 min',
                    '2 min', '2 min'], dtype='timedelta64[ns]')
        tm.assert_series_equal(algos.mode(s), exp)

    def test_mixed_dtype(self):
        exp = Series(['foo'])
        s = Series([1, 'foo', 'foo'])
        tm.assert_series_equal(algos.mode(s), exp)

    def test_uint64_overflow(self):
        exp = Series([2**63], dtype=np.uint64)
        s = Series([1, 2**63, 2**63], dtype=np.uint64)
        tm.assert_series_equal(algos.mode(s), exp)

        exp = Series([], dtype=np.uint64)
        s = Series([1, 2**63], dtype=np.uint64)
        tm.assert_series_equal(algos.mode(s), exp)

    def test_categorical(self):
        c = Categorical([1, 2])
        exp = Series([], dtype=np.int64)
        tm.assert_series_equal(algos.mode(c), exp)

        c = Categorical([1, 'a', 'a'])
        exp = Series(['a'], dtype=object)
        tm.assert_series_equal(algos.mode(c), exp)

        c = Categorical([1, 1, 2, 3, 3])
        exp = Series([1, 3], dtype=np.int64)
        tm.assert_series_equal(algos.mode(c), exp)

    def test_index(self):
        idx = Index([1, 2, 3])
        exp = Series([], dtype=np.int64)
        tm.assert_series_equal(algos.mode(idx), exp)

        idx = Index([1, 'a', 'a'])
        exp = Series(['a'], dtype=object)
        tm.assert_series_equal(algos.mode(idx), exp)

        idx = Index([1, 1, 2, 3, 3])
        exp = Series([1, 3], dtype=np.int64)
        tm.assert_series_equal(algos.mode(idx), exp)

        exp = Series(['2 min', '1 day'], dtype='timedelta64[ns]')
        idx = Index(['1 day', '1 day', '-1 day', '-1 day 2 min',
                     '2 min', '2 min'], dtype='timedelta64[ns]')
        tm.assert_series_equal(algos.mode(idx), exp)
