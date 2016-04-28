# -*- coding: utf-8 -*-
from pandas.compat import range

import numpy as np
from numpy.random import RandomState

from pandas.core.api import Series, Categorical, CategoricalIndex
import pandas as pd

from pandas import compat
import pandas.core.algorithms as algos
import pandas.util.testing as tm
import pandas.hashtable as hashtable
from pandas.compat.numpy import np_array_datetime64_compat


class TestMatch(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_ints(self):
        values = np.array([0, 2, 1])
        to_match = np.array([0, 1, 2, 2, 0, 1, 3, 0])

        result = algos.match(to_match, values)
        expected = np.array([0, 2, 1, 1, 0, 2, -1, 0])
        self.assert_numpy_array_equal(result, expected)

        result = Series(algos.match(to_match, values, np.nan))
        expected = Series(np.array([0, 2, 1, 1, 0, 2, np.nan, 0]))
        tm.assert_series_equal(result, expected)

        s = pd.Series(np.arange(5), dtype=np.float32)
        result = algos.match(s, [2, 4])
        expected = np.array([-1, -1, 0, -1, 1])
        self.assert_numpy_array_equal(result, expected)

        result = Series(algos.match(s, [2, 4], np.nan))
        expected = Series(np.array([np.nan, np.nan, 0, np.nan, 1]))
        tm.assert_series_equal(result, expected)

    def test_strings(self):
        values = ['foo', 'bar', 'baz']
        to_match = ['bar', 'foo', 'qux', 'foo', 'bar', 'baz', 'qux']

        result = algos.match(to_match, values)
        expected = np.array([1, 0, -1, 0, 1, 2, -1])
        self.assert_numpy_array_equal(result, expected)

        result = Series(algos.match(to_match, values, np.nan))
        expected = Series(np.array([1, 0, np.nan, 0, 1, 2, np.nan]))
        tm.assert_series_equal(result, expected)


class TestFactorize(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_basic(self):

        labels, uniques = algos.factorize(['a', 'b', 'b', 'a', 'a', 'c', 'c',
                                           'c'])
        self.assert_numpy_array_equal(
            uniques, np.array(['a', 'b', 'c'], dtype=object))

        labels, uniques = algos.factorize(['a', 'b', 'b', 'a',
                                           'a', 'c', 'c', 'c'], sort=True)
        self.assert_numpy_array_equal(labels, np.array(
            [0, 1, 1, 0, 0, 2, 2, 2], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, np.array(
            ['a', 'b', 'c'], dtype=object))

        labels, uniques = algos.factorize(list(reversed(range(5))))
        self.assert_numpy_array_equal(labels, np.array(
            [0, 1, 2, 3, 4], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, np.array(
            [4, 3, 2, 1, 0], dtype=np.int64))

        labels, uniques = algos.factorize(list(reversed(range(5))), sort=True)
        self.assert_numpy_array_equal(labels, np.array(
            [4, 3, 2, 1, 0], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, np.array(
            [0, 1, 2, 3, 4], dtype=np.int64))

        labels, uniques = algos.factorize(list(reversed(np.arange(5.))))
        self.assert_numpy_array_equal(labels, np.array(
            [0., 1., 2., 3., 4.], dtype=np.float64))
        self.assert_numpy_array_equal(uniques, np.array(
            [4, 3, 2, 1, 0], dtype=np.int64))

        labels, uniques = algos.factorize(
            list(reversed(np.arange(5.))), sort=True)
        self.assert_numpy_array_equal(labels, np.array(
            [4, 3, 2, 1, 0], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, np.array(
            [0., 1., 2., 3., 4.], dtype=np.float64))

    def test_mixed(self):

        # doc example reshaping.rst
        x = Series(['A', 'A', np.nan, 'B', 3.14, np.inf])
        labels, uniques = algos.factorize(x)

        self.assert_numpy_array_equal(labels, np.array(
            [0, 0, -1, 1, 2, 3], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, np.array(
            ['A', 'B', 3.14, np.inf], dtype=object))

        labels, uniques = algos.factorize(x, sort=True)
        self.assert_numpy_array_equal(labels, np.array(
            [2, 2, -1, 3, 0, 1], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, np.array(
            [3.14, np.inf, 'A', 'B'], dtype=object))

    def test_datelike(self):

        # M8
        v1 = pd.Timestamp('20130101 09:00:00.00004')
        v2 = pd.Timestamp('20130101')
        x = Series([v1, v1, v1, v2, v2, v1])
        labels, uniques = algos.factorize(x)
        self.assert_numpy_array_equal(labels, np.array(
            [0, 0, 0, 1, 1, 0], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, np.array(
            [v1.value, v2.value], dtype='M8[ns]'))

        labels, uniques = algos.factorize(x, sort=True)
        self.assert_numpy_array_equal(labels, np.array(
            [1, 1, 1, 0, 0, 1], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, np.array(
            [v2.value, v1.value], dtype='M8[ns]'))

        # period
        v1 = pd.Period('201302', freq='M')
        v2 = pd.Period('201303', freq='M')
        x = Series([v1, v1, v1, v2, v2, v1])

        # periods are not 'sorted' as they are converted back into an index
        labels, uniques = algos.factorize(x)
        self.assert_numpy_array_equal(labels, np.array(
            [0, 0, 0, 1, 1, 0], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, pd.PeriodIndex([v1, v2]))

        labels, uniques = algos.factorize(x, sort=True)
        self.assert_numpy_array_equal(labels, np.array(
            [0, 0, 0, 1, 1, 0], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, pd.PeriodIndex([v1, v2]))

        # GH 5986
        v1 = pd.to_timedelta('1 day 1 min')
        v2 = pd.to_timedelta('1 day')
        x = Series([v1, v2, v1, v1, v2, v2, v1])
        labels, uniques = algos.factorize(x)
        self.assert_numpy_array_equal(labels, np.array(
            [0, 1, 0, 0, 1, 1, 0], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, pd.to_timedelta([v1, v2]))

        labels, uniques = algos.factorize(x, sort=True)
        self.assert_numpy_array_equal(labels, np.array(
            [1, 0, 1, 1, 0, 0, 1], dtype=np.int64))
        self.assert_numpy_array_equal(uniques, pd.to_timedelta([v2, v1]))

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
            (hashtable.Float64HashTable, hashtable.Float64Vector, 'float64'),
            (hashtable.Int64HashTable, hashtable.Int64Vector, 'int64')]

        for (tbl, vect, dtype) in test_cases:
            # resizing to empty is a special case
            _test_vector_resize(tbl(), vect(), dtype, 0)
            _test_vector_resize(tbl(), vect(), dtype, 10)


class TestIndexer(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_outer_join_indexer(self):
        typemap = [('int32', algos.algos.outer_join_indexer_int32),
                   ('int64', algos.algos.outer_join_indexer_int64),
                   ('float32', algos.algos.outer_join_indexer_float32),
                   ('float64', algos.algos.outer_join_indexer_float64),
                   ('object', algos.algos.outer_join_indexer_object)]

        for dtype, indexer in typemap:
            left = np.arange(3, dtype=dtype)
            right = np.arange(2, 5, dtype=dtype)
            empty = np.array([], dtype=dtype)

            result, lindexer, rindexer = indexer(left, right)
            tm.assertIsInstance(result, np.ndarray)
            tm.assertIsInstance(lindexer, np.ndarray)
            tm.assertIsInstance(rindexer, np.ndarray)
            tm.assert_numpy_array_equal(result, np.arange(5, dtype=dtype))
            tm.assert_numpy_array_equal(lindexer, np.array([0, 1, 2, -1, -1]))
            tm.assert_numpy_array_equal(rindexer, np.array([-1, -1, 0, 1, 2]))

            result, lindexer, rindexer = indexer(empty, right)
            tm.assert_numpy_array_equal(result, right)
            tm.assert_numpy_array_equal(lindexer, np.array([-1, -1, -1]))
            tm.assert_numpy_array_equal(rindexer, np.array([0, 1, 2]))

            result, lindexer, rindexer = indexer(left, empty)
            tm.assert_numpy_array_equal(result, left)
            tm.assert_numpy_array_equal(lindexer, np.array([0, 1, 2]))
            tm.assert_numpy_array_equal(rindexer, np.array([-1, -1, -1]))


class TestUnique(tm.TestCase):
    _multiprocess_can_split_ = True

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


class TestIsin(tm.TestCase):
    _multiprocess_can_split_ = True

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
    _multiprocess_can_split_ = True

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
        # https://github.com/pydata/pandas/issues/9443#issuecomment-73719328

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
        np.testing.assert_allclose(out, expected_out, self.rtol)
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

        np.testing.assert_allclose(out, expected_out, self.rtol)
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
        np.testing.assert_allclose(out, expected_out, self.rtol)
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
                                  np.nan * np.ones(5)]).T
        expected_counts = counts + 2

        self.algo(out, counts, values, labels)
        np.testing.assert_allclose(out, expected_out, self.rtol)
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
    _multiprocess_can_split_ = True

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
    _multiprocess_can_split_ = True

    algo = algos.algos.group_var_float32
    dtype = np.float32
    rtol = 1e-2


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

    tm.assert_numpy_array_equal(left, right)

    a[np.random.choice(len(a), 10)] = -1
    left = unique_label_indices(a)
    right = np.unique(a, return_index=True)[1][1:]
    tm.assert_numpy_array_equal(left, right)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
