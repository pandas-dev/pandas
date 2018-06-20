# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest
from pandas import (CategoricalIndex, DatetimeIndex, Float64Index, Index,
                    Int64Index, IntervalIndex, MultiIndex, PeriodIndex,
                    RangeIndex, Series, TimedeltaIndex, UInt64Index, compat,
                    isna)
from pandas._libs.tslib import iNaT
from pandas.compat import PY3
from pandas.core.indexes.base import InvalidIndexError
from pandas.core.indexes.datetimelike import DatetimeIndexOpsMixin


def test_sortlevel(_index):
    import random

    tuples = list(_index)
    random.shuffle(tuples)

    index = MultiIndex.from_tuples(tuples)

    sorted_idx, _ = index.sortlevel(0)
    expected = MultiIndex.from_tuples(sorted(tuples))
    assert sorted_idx.equals(expected)

    sorted_idx, _ = index.sortlevel(0, ascending=False)
    assert sorted_idx.equals(expected[::-1])

    sorted_idx, _ = index.sortlevel(1)
    by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
    expected = MultiIndex.from_tuples(by1)
    assert sorted_idx.equals(expected)

    sorted_idx, _ = index.sortlevel(1, ascending=False)
    assert sorted_idx.equals(expected[::-1])


def test_sortlevel_not_sort_remaining():
    mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
    sorted_idx, _ = mi.sortlevel('A', sort_remaining=False)
    assert sorted_idx.equals(mi)


def test_sortlevel_deterministic():
    tuples = [('bar', 'one'), ('foo', 'two'), ('qux', 'two'),
              ('foo', 'one'), ('baz', 'two'), ('qux', 'one')]

    index = MultiIndex.from_tuples(tuples)

    sorted_idx, _ = index.sortlevel(0)
    expected = MultiIndex.from_tuples(sorted(tuples))
    assert sorted_idx.equals(expected)

    sorted_idx, _ = index.sortlevel(0, ascending=False)
    assert sorted_idx.equals(expected[::-1])

    sorted_idx, _ = index.sortlevel(1)
    by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
    expected = MultiIndex.from_tuples(by1)
    assert sorted_idx.equals(expected)

    sorted_idx, _ = index.sortlevel(1, ascending=False)
    assert sorted_idx.equals(expected[::-1])


def test_sort(indices):
    pytest.raises(TypeError, indices.sort)


def test_numpy_argsort(named_index):
    for k, ind in named_index.items():
        result = np.argsort(ind)
        expected = ind.argsort()
        tm.assert_numpy_array_equal(result, expected)

        # these are the only two types that perform
        # pandas compatibility input validation - the
        # rest already perform separate (or no) such
        # validation via their 'values' attribute as
        # defined in pandas.core.indexes/base.py - they
        # cannot be changed at the moment due to
        # backwards compatibility concerns
        if isinstance(type(ind), (CategoricalIndex, RangeIndex)):
            msg = "the 'axis' parameter is not supported"
            tm.assert_raises_regex(ValueError, msg,
                                   np.argsort, ind, axis=1)

            msg = "the 'kind' parameter is not supported"
            tm.assert_raises_regex(ValueError, msg, np.argsort,
                                   ind, kind='mergesort')

            msg = "the 'order' parameter is not supported"
            tm.assert_raises_regex(ValueError, msg, np.argsort,
                                   ind, order=('a', 'b'))
