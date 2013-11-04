#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas as pd
import unittest
import warnings
import nose
import numpy as np
import sys
from pandas import Series
from pandas.util.testing import (
    assert_almost_equal, assertRaisesRegexp, raise_with_traceback, assert_series_equal
)

# let's get meta.

class TestAssertAlmostEqual(unittest.TestCase):
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

class TestUtilTesting(unittest.TestCase):
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

class TestAssertSeriesEqual(unittest.TestCase):
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

