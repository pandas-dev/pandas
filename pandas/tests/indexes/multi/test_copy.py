# -*- coding: utf-8 -*-

from copy import copy, deepcopy

import pandas.util.testing as tm
from pandas import (CategoricalIndex, IntervalIndex, MultiIndex, PeriodIndex,
                    RangeIndex, Series, compat)


def assert_multiindex_copied(copy, original):
    # Levels should be (at least, shallow copied)
    tm.assert_copy(copy.levels, original.levels)
    tm.assert_almost_equal(copy.labels, original.labels)

    # Labels doesn't matter which way copied
    tm.assert_almost_equal(copy.labels, original.labels)
    assert copy.labels is not original.labels

    # Names doesn't matter which way copied
    assert copy.names == original.names
    assert copy.names is not original.names

    # Sort order should be copied
    assert copy.sortorder == original.sortorder


def test_copy(idx):
    i_copy = idx.copy()

    assert_multiindex_copied(i_copy, idx)


def test_shallow_copy(idx):
    i_copy = idx._shallow_copy()

    assert_multiindex_copied(i_copy, idx)


def test_view(idx):
    i_view = idx.view()
    assert_multiindex_copied(i_view, idx)


def test_copy_name(idx):
    # gh-12309: Check that the "name" argument
    # passed at initialization is honored.

    # TODO: Remove or refactor MultiIndex not tested.
    for name, index in compat.iteritems({'idx': idx}):
        if isinstance(index, MultiIndex):
            continue

        first = index.__class__(index, copy=True, name='mario')
        second = first.__class__(first, copy=False)

        # Even though "copy=False", we want a new object.
        assert first is not second

        # Not using tm.assert_index_equal() since names differ.
        assert index.equals(first)

        assert first.name == 'mario'
        assert second.name == 'mario'

        s1 = Series(2, index=first)
        s2 = Series(3, index=second[:-1])

        if not isinstance(index, CategoricalIndex):
            # See gh-13365
            s3 = s1 * s2
            assert s3.index.name == 'mario'


def test_ensure_copied_data(idx):
    # Check the "copy" argument of each Index.__new__ is honoured
    # GH12309
    # TODO: REMOVE THIS TEST.  MultiIndex is tested seperately as noted below.

    for name, index in compat.iteritems({'idx': idx}):
        init_kwargs = {}
        if isinstance(index, PeriodIndex):
            # Needs "freq" specification:
            init_kwargs['freq'] = index.freq
        elif isinstance(index, (RangeIndex, MultiIndex, CategoricalIndex)):
            # RangeIndex cannot be initialized from data
            # MultiIndex and CategoricalIndex are tested separately
            continue

        index_type = index.__class__
        result = index_type(index.values, copy=True, **init_kwargs)
        tm.assert_index_equal(index, result)
        tm.assert_numpy_array_equal(index.values, result.values,
                                    check_same='copy')

        if isinstance(index, PeriodIndex):
            # .values an object array of Period, thus copied
            result = index_type(ordinal=index.asi8, copy=False,
                                **init_kwargs)
            tm.assert_numpy_array_equal(index._ndarray_values,
                                        result._ndarray_values,
                                        check_same='same')
        elif isinstance(index, IntervalIndex):
            # checked in test_interval.py
            pass
        else:
            result = index_type(index.values, copy=False, **init_kwargs)
            tm.assert_numpy_array_equal(index.values, result.values,
                                        check_same='same')
            tm.assert_numpy_array_equal(index._ndarray_values,
                                        result._ndarray_values,
                                        check_same='same')


def test_copy_and_deepcopy(indices):

    if isinstance(indices, MultiIndex):
        return
    for func in (copy, deepcopy):
        idx_copy = func(indices)
        assert idx_copy is not indices
        assert idx_copy.equals(indices)

    new_copy = indices.copy(deep=True, name="banana")
    assert new_copy.name == "banana"
