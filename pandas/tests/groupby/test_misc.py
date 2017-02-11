""" misc non-groupby routines, as they are defined in core/groupby.py """

import pytest
import numpy as np
from numpy import nan
from pandas.util import testing as tm
from pandas.core.groupby import _nargsort, _lexsort_indexer


class TestSorting(tm.TestCase):

    def test_lexsort_indexer(self):
        keys = [[nan] * 5 + list(range(100)) + [nan] * 5]
        # orders=True, na_position='last'
        result = _lexsort_indexer(keys, orders=True, na_position='last')
        exp = list(range(5, 105)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

        # orders=True, na_position='first'
        result = _lexsort_indexer(keys, orders=True, na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(5, 105))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

        # orders=False, na_position='last'
        result = _lexsort_indexer(keys, orders=False, na_position='last')
        exp = list(range(104, 4, -1)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

        # orders=False, na_position='first'
        result = _lexsort_indexer(keys, orders=False, na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(104, 4, -1))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

    def test_nargsort(self):
        # np.argsort(items) places NaNs last
        items = [nan] * 5 + list(range(100)) + [nan] * 5
        # np.argsort(items2) may not place NaNs first
        items2 = np.array(items, dtype='O')

        try:
            # GH 2785; due to a regression in NumPy1.6.2
            np.argsort(np.array([[1, 2], [1, 3], [1, 2]], dtype='i'))
            np.argsort(items2, kind='mergesort')
        except TypeError:
            pytest.skip('requested sort not available for type')

        # mergesort is the most difficult to get right because we want it to be
        # stable.

        # According to numpy/core/tests/test_multiarray, """The number of
        # sorted items must be greater than ~50 to check the actual algorithm
        # because quick and merge sort fall over to insertion sort for small
        # arrays."""

        # mergesort, ascending=True, na_position='last'
        result = _nargsort(items, kind='mergesort', ascending=True,
                           na_position='last')
        exp = list(range(5, 105)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=True, na_position='first'
        result = _nargsort(items, kind='mergesort', ascending=True,
                           na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(5, 105))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='last'
        result = _nargsort(items, kind='mergesort', ascending=False,
                           na_position='last')
        exp = list(range(104, 4, -1)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='first'
        result = _nargsort(items, kind='mergesort', ascending=False,
                           na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(104, 4, -1))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=True, na_position='last'
        result = _nargsort(items2, kind='mergesort', ascending=True,
                           na_position='last')
        exp = list(range(5, 105)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=True, na_position='first'
        result = _nargsort(items2, kind='mergesort', ascending=True,
                           na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(5, 105))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='last'
        result = _nargsort(items2, kind='mergesort', ascending=False,
                           na_position='last')
        exp = list(range(104, 4, -1)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='first'
        result = _nargsort(items2, kind='mergesort', ascending=False,
                           na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(104, 4, -1))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)
