#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas as pd
import unittest
import warnings
import nose
import numpy as np
import sys
from pandas import Series, DataFrame
import pandas.util.testing as tm
from pandas.util.testing import (
    assert_almost_equal, assertRaisesRegexp, raise_with_traceback,
    assert_index_equal, assert_series_equal, assert_frame_equal,
    assert_numpy_array_equal, assert_isinstance, RNGContext
)
from pandas.compat import is_platform_windows

# let's get meta.

class TestAssertAlmostEqual(tm.TestCase):
    _multiprocess_can_split_ = True

    def _assert_almost_equal_both(self, a, b, **kwargs):
        assert_almost_equal(a, b, **kwargs)
        assert_almost_equal(b, a, **kwargs)

    def _assert_not_almost_equal_both(self, a, b, **kwargs):
        self.assertRaises(AssertionError, assert_almost_equal, a, b, **kwargs)
        self.assertRaises(AssertionError, assert_almost_equal, b, a, **kwargs)

    def test_assert_almost_equal_numbers(self):
        self._assert_almost_equal_both(1.1, 1.1)
        self._assert_almost_equal_both(1.1, 1.100001)
        self._assert_almost_equal_both(np.int16(1), 1.000001)
        self._assert_almost_equal_both(np.float64(1.1), 1.1)
        self._assert_almost_equal_both(np.uint32(5), 5)

        self._assert_not_almost_equal_both(1.1, 1)
        self._assert_not_almost_equal_both(1.1, True)
        self._assert_not_almost_equal_both(1, 2)
        self._assert_not_almost_equal_both(1.0001, np.int16(1))

    def test_assert_almost_equal_numbers_with_zeros(self):
        self._assert_almost_equal_both(0, 0)
        self._assert_almost_equal_both(0.000001, 0)

        self._assert_not_almost_equal_both(0.001, 0)
        self._assert_not_almost_equal_both(1, 0)

    def test_assert_almost_equal_numbers_with_mixed(self):
        self._assert_not_almost_equal_both(1, 'abc')
        self._assert_not_almost_equal_both(1, [1,])
        self._assert_not_almost_equal_both(1, object())

    def test_assert_almost_equal_edge_case_ndarrays(self):
        self._assert_almost_equal_both(np.array([], dtype='M8[ns]'),
                                       np.array([], dtype='float64'))
        self._assert_almost_equal_both(np.array([], dtype=str),
                                       np.array([], dtype='int64'))

    def test_assert_almost_equal_dicts(self):
        self._assert_almost_equal_both({'a': 1, 'b': 2}, {'a': 1, 'b': 2})

        self._assert_not_almost_equal_both({'a': 1, 'b': 2}, {'a': 1, 'b': 3})
        self._assert_not_almost_equal_both(
            {'a': 1, 'b': 2}, {'a': 1, 'b': 2, 'c': 3}
        )
        self._assert_not_almost_equal_both({'a': 1}, 1)
        self._assert_not_almost_equal_both({'a': 1}, 'abc')
        self._assert_not_almost_equal_both({'a': 1}, [1,])

    def test_assert_almost_equal_dict_like_object(self):
        class DictLikeObj(object):
            def keys(self):
                return ('a',)

            def __getitem__(self, item):
                if item == 'a':
                    return 1

        self._assert_almost_equal_both({'a': 1}, DictLikeObj())

        self._assert_not_almost_equal_both({'a': 2}, DictLikeObj())

    def test_assert_almost_equal_strings(self):
        self._assert_almost_equal_both('abc', 'abc')

        self._assert_not_almost_equal_both('abc', 'abcd')
        self._assert_not_almost_equal_both('abc', 'abd')
        self._assert_not_almost_equal_both('abc', 1)
        self._assert_not_almost_equal_both('abc', [1,])

    def test_assert_almost_equal_iterables(self):
        self._assert_almost_equal_both([1, 2, 3], [1, 2, 3])
        self._assert_almost_equal_both(np.array([1, 2, 3]), [1, 2, 3])

        # Can't compare generators
        self._assert_not_almost_equal_both(iter([1, 2, 3]), [1, 2, 3])

        self._assert_not_almost_equal_both([1, 2, 3], [1, 2, 4])
        self._assert_not_almost_equal_both([1, 2, 3], [1, 2, 3, 4])
        self._assert_not_almost_equal_both([1, 2, 3], 1)

    def test_assert_almost_equal_null(self):
        self._assert_almost_equal_both(None, None)
        self._assert_almost_equal_both(None, np.NaN)

        self._assert_not_almost_equal_both(None, 0)
        self._assert_not_almost_equal_both(np.NaN, 0)

    def test_assert_almost_equal_inf(self):
        self._assert_almost_equal_both(np.inf, np.inf)
        self._assert_almost_equal_both(np.inf, float("inf"))

        self._assert_not_almost_equal_both(np.inf, 0)


class TestUtilTesting(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_raise_with_traceback(self):
        with assertRaisesRegexp(LookupError, "error_text"):
            try:
                raise ValueError("THIS IS AN ERROR")
            except ValueError as e:
                e = LookupError("error_text")
                raise_with_traceback(e)
        with assertRaisesRegexp(LookupError, "error_text"):
            try:
                raise ValueError("This is another error")
            except ValueError:
                e = LookupError("error_text")
                _, _, traceback = sys.exc_info()
                raise_with_traceback(e, traceback)


class TestAssertNumpyArrayEqual(tm.TestCase):

    def test_numpy_array_equal_message(self):

        if is_platform_windows():
            raise nose.SkipTest("windows has incomparable line-endings and uses L on the shape")

        expected = """numpy array are different

numpy array shapes are different
\\[left\\]:  \\(2,\\)
\\[right\\]: \\(3,\\)"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_numpy_array_equal(np.array([1, 2]), np.array([3, 4, 5]))

        with assertRaisesRegexp(AssertionError, expected):
            assert_almost_equal(np.array([1, 2]), np.array([3, 4, 5]))

        # scalar comparison
        expected = """: 1 != 2"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_numpy_array_equal(1, 2)
        expected = """expected 2\\.00000 but got 1\\.00000, with decimal 5"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_almost_equal(1, 2)

        # array / scalar array comparison
        expected = """(numpy array|Iterable) are different

First object is iterable, second isn't
\\[left\\]:  \\[1\\]
\\[right\\]: 1"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_numpy_array_equal(np.array([1]), 1)
        with assertRaisesRegexp(AssertionError, expected):
            assert_almost_equal(np.array([1]), 1)

        # scalar / array comparison
        expected = """(numpy array|Iterable) are different

Second object is iterable, first isn't
\\[left\\]:  1
\\[right\\]: \\[1\\]"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_numpy_array_equal(1, np.array([1]))
        with assertRaisesRegexp(AssertionError, expected):
            assert_almost_equal(1, np.array([1]))

        expected = """numpy array are different

numpy array values are different \\(66\\.66667 %\\)
\\[left\\]:  \\[nan, 2\\.0, 3\\.0\\]
\\[right\\]: \\[1\\.0, nan, 3\\.0\\]"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_numpy_array_equal(np.array([np.nan, 2, 3]), np.array([1, np.nan, 3]))
        with assertRaisesRegexp(AssertionError, expected):
            assert_almost_equal(np.array([np.nan, 2, 3]), np.array([1, np.nan, 3]))

        expected = """numpy array are different

numpy array values are different \\(50\\.0 %\\)
\\[left\\]:  \\[1, 2\\]
\\[right\\]: \\[1, 3\\]"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_numpy_array_equal(np.array([1, 2]), np.array([1, 3]))
        with assertRaisesRegexp(AssertionError, expected):
            assert_almost_equal(np.array([1, 2]), np.array([1, 3]))


        expected = """numpy array are different

numpy array values are different \\(50\\.0 %\\)
\\[left\\]:  \\[1\\.1, 2\\.000001\\]
\\[right\\]: \\[1\\.1, 2.0\\]"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_numpy_array_equal(np.array([1.1, 2.000001]), np.array([1.1, 2.0]))

        # must pass
        assert_almost_equal(np.array([1.1, 2.000001]), np.array([1.1, 2.0]))

        expected = """numpy array are different

numpy array values are different \\(16\\.66667 %\\)
\\[left\\]:  \\[\\[1, 2\\], \\[3, 4\\], \\[5, 6\\]\\]
\\[right\\]: \\[\\[1, 3\\], \\[3, 4\\], \\[5, 6\\]\\]"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_numpy_array_equal(np.array([[1, 2], [3, 4], [5, 6]]),
                                     np.array([[1, 3], [3, 4], [5, 6]]))
        with assertRaisesRegexp(AssertionError, expected):
            assert_almost_equal(np.array([[1, 2], [3, 4], [5, 6]]),
                                np.array([[1, 3], [3, 4], [5, 6]]))

        expected = """numpy array are different

numpy array values are different \\(25\\.0 %\\)
\\[left\\]:  \\[\\[1, 2\\], \\[3, 4\\]\\]
\\[right\\]: \\[\\[1, 3\\], \\[3, 4\\]\\]"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_numpy_array_equal(np.array([[1, 2], [3, 4]]),
                                     np.array([[1, 3], [3, 4]]))
        with assertRaisesRegexp(AssertionError, expected):
            assert_almost_equal(np.array([[1, 2], [3, 4]]),
                                np.array([[1, 3], [3, 4]]))

        # allow to overwrite message
        expected = """Index are different

Index shapes are different
\\[left\\]:  \\(2,\\)
\\[right\\]: \\(3,\\)"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_numpy_array_equal(np.array([1, 2]), np.array([3, 4, 5]),
                                     obj='Index')
        with assertRaisesRegexp(AssertionError, expected):
            assert_almost_equal(np.array([1, 2]), np.array([3, 4, 5]),
                                obj='Index')

    def test_assert_almost_equal_iterable_message(self):

        expected = """Iterable are different

Iterable length are different
\\[left\\]:  2
\\[right\\]: 3"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_almost_equal([1, 2], [3, 4, 5])

        expected = """Iterable are different

Iterable values are different \\(50\\.0 %\\)
\\[left\\]:  \\[1, 2\\]
\\[right\\]: \\[1, 3\\]"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_almost_equal([1, 2], [1, 3])


class TestAssertIndexEqual(unittest.TestCase):
    _multiprocess_can_split_ = True

    def test_index_equal_message(self):

        expected = """Index are different

Index levels are different
\\[left\\]:  1, Int64Index\\(\\[1, 2, 3\\], dtype='int64'\\)
\\[right\\]: 2, MultiIndex\\(levels=\\[\\[u?'A', u?'B'\\], \\[1, 2, 3, 4\\]\\],
           labels=\\[\\[0, 0, 1, 1\\], \\[0, 1, 2, 3\\]\\]\\)"""
        idx1 = pd.Index([1, 2, 3])
        idx2 = pd.MultiIndex.from_tuples([('A', 1), ('A', 2), ('B', 3), ('B', 4)])
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2)
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2, exact=False)


        expected = """MultiIndex level \\[1\\] are different

MultiIndex level \\[1\\] values are different \\(25\\.0 %\\)
\\[left\\]:  Int64Index\\(\\[2, 2, 3, 4\\], dtype='int64'\\)
\\[right\\]: Int64Index\\(\\[1, 2, 3, 4\\], dtype='int64'\\)"""
        idx1 = pd.MultiIndex.from_tuples([('A', 2), ('A', 2), ('B', 3), ('B', 4)])
        idx2 = pd.MultiIndex.from_tuples([('A', 1), ('A', 2), ('B', 3), ('B', 4)])
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2)
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2, check_exact=False)

        expected = """Index are different

Index length are different
\\[left\\]:  3, Int64Index\\(\\[1, 2, 3\\], dtype='int64'\\)
\\[right\\]: 4, Int64Index\\(\\[1, 2, 3, 4\\], dtype='int64'\\)"""
        idx1 = pd.Index([1, 2, 3])
        idx2 = pd.Index([1, 2, 3, 4])
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2)
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2, check_exact=False)

        expected = """Index are different

Index classes are different
\\[left\\]:  Int64Index\\(\\[1, 2, 3\\], dtype='int64'\\)
\\[right\\]: Float64Index\\(\\[1\\.0, 2\\.0, 3\\.0\\], dtype='float64'\\)"""
        idx1 = pd.Index([1, 2, 3])
        idx2 = pd.Index([1, 2, 3.0])
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2, exact=True)
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2, exact=True, check_exact=False)

        expected = """Index are different

Index values are different \\(33\\.33333 %\\)
\\[left\\]:  Float64Index\\(\\[1.0, 2.0, 3.0], dtype='float64'\\)
\\[right\\]: Float64Index\\(\\[1.0, 2.0, 3.0000000001\\], dtype='float64'\\)"""
        idx1 = pd.Index([1, 2, 3.])
        idx2 = pd.Index([1, 2, 3.0000000001])
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2)

        # must success
        assert_index_equal(idx1, idx2, check_exact=False)

        expected = """Index are different

Index values are different \\(33\\.33333 %\\)
\\[left\\]:  Float64Index\\(\\[1.0, 2.0, 3.0], dtype='float64'\\)
\\[right\\]: Float64Index\\(\\[1.0, 2.0, 3.0001\\], dtype='float64'\\)"""
        idx1 = pd.Index([1, 2, 3.])
        idx2 = pd.Index([1, 2, 3.0001])
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2)
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2, check_exact=False)
        # must success
        assert_index_equal(idx1, idx2, check_exact=False, check_less_precise=True)

        expected = """Index are different

Index values are different \\(33\\.33333 %\\)
\\[left\\]:  Int64Index\\(\\[1, 2, 3\\], dtype='int64'\\)
\\[right\\]: Int64Index\\(\\[1, 2, 4\\], dtype='int64'\\)"""
        idx1 = pd.Index([1, 2, 3])
        idx2 = pd.Index([1, 2, 4])
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2)
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2, check_less_precise=True)

        expected = """MultiIndex level \\[1\\] are different

MultiIndex level \\[1\\] values are different \\(25\\.0 %\\)
\\[left\\]:  Int64Index\\(\\[2, 2, 3, 4\\], dtype='int64'\\)
\\[right\\]: Int64Index\\(\\[1, 2, 3, 4\\], dtype='int64'\\)"""
        idx1 = pd.MultiIndex.from_tuples([('A', 2), ('A', 2), ('B', 3), ('B', 4)])
        idx2 = pd.MultiIndex.from_tuples([('A', 1), ('A', 2), ('B', 3), ('B', 4)])
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2)
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2, check_exact=False)

    def test_index_equal_metadata_message(self):

        expected = """Index are different

Attribute "names" are different
\\[left\\]:  \\[None\\]
\\[right\\]: \\[u?'x'\\]"""
        idx1 = pd.Index([1, 2, 3])
        idx2 = pd.Index([1, 2, 3], name='x')
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2)

        # same name, should pass
        assert_index_equal(pd.Index([1, 2, 3], name=np.nan),
                              pd.Index([1, 2, 3], name=np.nan))
        assert_index_equal(pd.Index([1, 2, 3], name=pd.NaT),
                              pd.Index([1, 2, 3], name=pd.NaT))


        expected = """Index are different

Attribute "names" are different
\\[left\\]:  \\[nan\\]
\\[right\\]: \\[NaT\\]"""
        idx1 = pd.Index([1, 2, 3], name=np.nan)
        idx2 = pd.Index([1, 2, 3], name=pd.NaT)
        with assertRaisesRegexp(AssertionError, expected):
            assert_index_equal(idx1, idx2)


class TestAssertSeriesEqual(tm.TestCase):
    _multiprocess_can_split_ = True

    def _assert_equal(self, x, y, **kwargs):
        assert_series_equal(x,y,**kwargs)
        assert_series_equal(y,x,**kwargs)

    def _assert_not_equal(self, a, b, **kwargs):
        self.assertRaises(AssertionError, assert_series_equal, a, b, **kwargs)
        self.assertRaises(AssertionError, assert_series_equal, b, a, **kwargs)

    def test_equal(self):
        self._assert_equal(Series(range(3)),Series(range(3)))
        self._assert_equal(Series(list('abc')),Series(list('abc')))

    def test_not_equal(self):
        self._assert_not_equal(Series(range(3)),Series(range(3))+1)
        self._assert_not_equal(Series(list('abc')),Series(list('xyz')))
        self._assert_not_equal(Series(range(3)),Series(range(4)))
        self._assert_not_equal(Series(range(3)),Series(range(3),dtype='float64'))
        self._assert_not_equal(Series(range(3)),Series(range(3),index=[1,2,4]))

        # ATM meta data is not checked in assert_series_equal
        # self._assert_not_equal(Series(range(3)),Series(range(3),name='foo'),check_names=True)

    def test_less_precise(self):
        s1 =  Series([0.12345],dtype='float64')
        s2 =  Series([0.12346],dtype='float64')

        self.assertRaises(AssertionError, assert_series_equal, s1, s2)
        self._assert_equal(s1,s2,check_less_precise=True)

        s1 =  Series([0.12345],dtype='float32')
        s2 =  Series([0.12346],dtype='float32')

        self.assertRaises(AssertionError, assert_series_equal, s1, s2)
        self._assert_equal(s1,s2,check_less_precise=True)

        # even less than less precise
        s1 =  Series([0.1235],dtype='float32')
        s2 =  Series([0.1236],dtype='float32')

        self.assertRaises(AssertionError, assert_series_equal, s1, s2)
        self.assertRaises(AssertionError, assert_series_equal, s1, s2, True)

    def test_index_dtype(self):
        df1 = DataFrame.from_records(
                {'a':[1,2],'c':['l1','l2']}, index=['a'])
        df2 = DataFrame.from_records(
                {'a':[1.0,2.0],'c':['l1','l2']}, index=['a'])
        self._assert_not_equal(df1.c, df2.c, check_index_type=True)

    def test_multiindex_dtype(self):
        df1 = DataFrame.from_records(
                {'a':[1,2],'b':[2.1,1.5],'c':['l1','l2']}, index=['a','b'])
        df2 = DataFrame.from_records(
                {'a':[1.0,2.0],'b':[2.1,1.5],'c':['l1','l2']}, index=['a','b'])
        self._assert_not_equal(df1.c, df2.c, check_index_type=True)

    def test_series_equal_message(self):

        expected = """Series are different

Series length are different
\\[left\\]:  3, Int64Index\\(\\[0, 1, 2\\], dtype='int64'\\)
\\[right\\]: 4, Int64Index\\(\\[0, 1, 2, 3\\], dtype='int64'\\)"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_series_equal(pd.Series([1, 2, 3]), pd.Series([1, 2, 3, 4]))


        expected = """Series are different

Series values are different \\(33\\.33333 %\\)
\\[left\\]:  \\[1, 2, 3\\]
\\[right\\]: \\[1, 2, 4\\]"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_series_equal(pd.Series([1, 2, 3]), pd.Series([1, 2, 4]))
        with assertRaisesRegexp(AssertionError, expected):
            assert_series_equal(pd.Series([1, 2, 3]), pd.Series([1, 2, 4]),
                                   check_less_precise=True)


class TestAssertFrameEqual(tm.TestCase):
    _multiprocess_can_split_ = True

    def _assert_equal(self, x, y, **kwargs):
        assert_frame_equal(x,y,**kwargs)
        assert_frame_equal(y,x,**kwargs)

    def _assert_not_equal(self, a, b, **kwargs):
        self.assertRaises(AssertionError, assert_frame_equal, a, b, **kwargs)
        self.assertRaises(AssertionError, assert_frame_equal, b, a, **kwargs)

    def test_index_dtype(self):
        df1 = DataFrame.from_records(
                {'a':[1,2],'c':['l1','l2']}, index=['a'])
        df2 = DataFrame.from_records(
                {'a':[1.0,2.0],'c':['l1','l2']}, index=['a'])
        self._assert_not_equal(df1, df2, check_index_type=True)

    def test_multiindex_dtype(self):
        df1 = DataFrame.from_records(
                {'a':[1,2],'b':[2.1,1.5],'c':['l1','l2']}, index=['a','b'])
        df2 = DataFrame.from_records(
                {'a':[1.0,2.0],'b':[2.1,1.5],'c':['l1','l2']}, index=['a','b'])
        self._assert_not_equal(df1, df2, check_index_type=True)

    def test_empty_dtypes(self):
        df1=pd.DataFrame(columns=["col1","col2"])
        df1["col1"] = df1["col1"].astype('int64')
        df2=pd.DataFrame(columns=["col1","col2"])
        self._assert_equal(df1, df2, check_dtype=False)
        self._assert_not_equal(df1, df2, check_dtype=True)

    def test_frame_equal_message(self):

        expected = """DataFrame are different

DataFrame shape \\(number of rows\\) are different
\\[left\\]:  3, Int64Index\\(\\[0, 1, 2\\], dtype='int64'\\)
\\[right\\]: 4, Int64Index\\(\\[0, 1, 2, 3\\], dtype='int64'\\)"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_frame_equal(pd.DataFrame({'A':[1, 2, 3]}),
                                  pd.DataFrame({'A':[1, 2, 3, 4]}))


        expected = """DataFrame are different

DataFrame shape \\(number of columns\\) are different
\\[left\\]:  2, Index\\(\\[u?'A', u?'B'\\], dtype='object'\\)
\\[right\\]: 1, Index\\(\\[u?'A'\\], dtype='object'\\)"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_frame_equal(pd.DataFrame({'A':[1, 2, 3], 'B':[4, 5, 6]}),
                                  pd.DataFrame({'A':[1, 2, 3]}))


        expected = """DataFrame\\.index are different

DataFrame\\.index values are different \\(33\\.33333 %\\)
\\[left\\]:  Index\\(\\[u?'a', u?'b', u?'c'\\], dtype='object'\\)
\\[right\\]: Index\\(\\[u?'a', u?'b', u?'d'\\], dtype='object'\\)"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_frame_equal(pd.DataFrame({'A':[1, 2, 3], 'B':[4, 5, 6]},
                                               index=['a', 'b', 'c']),
                                  pd.DataFrame({'A':[1, 2, 3], 'B':[4, 5, 6]},
                                               index=['a', 'b', 'd']))

        expected = """DataFrame\\.columns are different

DataFrame\\.columns values are different \\(50\\.0 %\\)
\\[left\\]:  Index\\(\\[u?'A', u?'B'\\], dtype='object'\\)
\\[right\\]: Index\\(\\[u?'A', u?'b'\\], dtype='object'\\)"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_frame_equal(pd.DataFrame({'A':[1, 2, 3], 'B':[4, 5, 6]},
                                               index=['a', 'b', 'c']),
                                  pd.DataFrame({'A':[1, 2, 3], 'b':[4, 5, 6]},
                                               index=['a', 'b', 'c']))


        expected = """DataFrame\\.iloc\\[:, 1\\] are different

DataFrame\\.iloc\\[:, 1\\] values are different \\(33\\.33333 %\\)
\\[left\\]:  \\[4, 5, 6\\]
\\[right\\]: \\[4, 5, 7\\]"""
        with assertRaisesRegexp(AssertionError, expected):
            assert_frame_equal(pd.DataFrame({'A':[1, 2, 3], 'B':[4, 5, 6]}),
                               pd.DataFrame({'A':[1, 2, 3], 'B':[4, 5, 7]}))

        with assertRaisesRegexp(AssertionError, expected):
            assert_frame_equal(pd.DataFrame({'A':[1, 2, 3], 'B':[4, 5, 6]}),
                                  pd.DataFrame({'A':[1, 2, 3], 'B':[4, 5, 7]}),
                                  by_blocks=True)


class TestRNGContext(unittest.TestCase):

    def test_RNGContext(self):
        expected0 = 1.764052345967664
        expected1 = 1.6243453636632417

        with RNGContext(0):
            with RNGContext(1):
                self.assertEqual(np.random.randn(), expected1)
            self.assertEqual(np.random.randn(), expected0)



class TestDeprecatedTests(tm.TestCase):

    def test_warning(self):

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.assertEquals(1, 1)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.assertNotEquals(1, 2)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.assert_(True)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.assertAlmostEquals(1.0, 1.0000000001)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.assertNotAlmostEquals(1, 2)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            tm.assert_isinstance(Series([1, 2]), Series, msg='xxx')


class TestLocale(tm.TestCase):

    def test_locale(self):
        if sys.platform == 'win32':
            raise nose.SkipTest("skipping on win platforms as locale not available")

        #GH9744
        locales = tm.get_locales()
        self.assertTrue(len(locales) >= 1)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
