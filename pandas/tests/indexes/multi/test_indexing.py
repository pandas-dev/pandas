# -*- coding: utf-8 -*-


from datetime import timedelta

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest
from pandas import Index, MultiIndex
from pandas.compat import lrange
from pandas.core.indexes.base import InvalidIndexError
from pandas.util.testing import assert_almost_equal
from pandas import (CategoricalIndex, DatetimeIndex, Float64Index, Index,
                    Int64Index, IntervalIndex, MultiIndex, PeriodIndex,
                    RangeIndex, Series, TimedeltaIndex, UInt64Index, compat,
                    isna)


def test_slice_locs_partial(_index):
    sorted_idx, _ = _index.sortlevel(0)

    result = sorted_idx.slice_locs(('foo', 'two'), ('qux', 'one'))
    assert result == (1, 5)

    result = sorted_idx.slice_locs(None, ('qux', 'one'))
    assert result == (0, 5)

    result = sorted_idx.slice_locs(('foo', 'two'), None)
    assert result == (1, len(sorted_idx))

    result = sorted_idx.slice_locs('bar', 'baz')
    assert result == (2, 4)


def test_slice_locs():
    df = tm.makeTimeDataFrame()
    stacked = df.stack()
    idx = stacked.index

    slob = slice(*idx.slice_locs(df.index[5], df.index[15]))
    sliced = stacked[slob]
    expected = df[5:16].stack()
    tm.assert_almost_equal(sliced.values, expected.values)

    slob = slice(*idx.slice_locs(df.index[5] + timedelta(seconds=30),
                                 df.index[15] - timedelta(seconds=30)))
    sliced = stacked[slob]
    expected = df[6:15].stack()
    tm.assert_almost_equal(sliced.values, expected.values)


def test_slice_locs_with_type_mismatch():
    df = tm.makeTimeDataFrame()
    stacked = df.stack()
    idx = stacked.index
    tm.assert_raises_regex(TypeError, '^Level type mismatch',
                           idx.slice_locs, (1, 3))
    tm.assert_raises_regex(TypeError, '^Level type mismatch',
                           idx.slice_locs,
                           df.index[5] + timedelta(
                               seconds=30), (5, 2))
    df = tm.makeCustomDataframe(5, 5)
    stacked = df.stack()
    idx = stacked.index
    with tm.assert_raises_regex(TypeError, '^Level type mismatch'):
        idx.slice_locs(timedelta(seconds=30))
    # TODO: Try creating a UnicodeDecodeError in exception message
    with tm.assert_raises_regex(TypeError, '^Level type mismatch'):
        idx.slice_locs(df.index[1], (16, "a"))


def test_slice_locs_not_sorted():
    index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
        lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
            [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])

    tm.assert_raises_regex(KeyError, "[Kk]ey length.*greater than "
                           "MultiIndex lexsort depth",
                           index.slice_locs, (1, 0, 1), (2, 1, 0))

    # works
    sorted_index, _ = index.sortlevel(0)
    # should there be a test case here???
    sorted_index.slice_locs((1, 0, 1), (2, 1, 0))


def test_slice_locs_not_contained():
    # some searchsorted action

    index = MultiIndex(levels=[[0, 2, 4, 6], [0, 2, 4]],
                       labels=[[0, 0, 0, 1, 1, 2, 3, 3, 3],
                               [0, 1, 2, 1, 2, 2, 0, 1, 2]], sortorder=0)

    result = index.slice_locs((1, 0), (5, 2))
    assert result == (3, 6)

    result = index.slice_locs(1, 5)
    assert result == (3, 6)

    result = index.slice_locs((2, 2), (5, 2))
    assert result == (3, 6)

    result = index.slice_locs(2, 5)
    assert result == (3, 6)

    result = index.slice_locs((1, 0), (6, 3))
    assert result == (3, 8)

    result = index.slice_locs(-1, 10)
    assert result == (0, len(index))



def test_to_series(_index):
    # assert that we are creating a copy of the index

    idx = _index
    s = idx.to_series()
    assert s.values is not idx.values
    assert s.index is not idx
    assert s.name == idx.name


def test_to_series_with_arguments(_index):
    # GH18699

    # index kwarg
    idx = _index
    s = idx.to_series(index=idx)

    assert s.values is not idx.values
    assert s.index is idx
    assert s.name == idx.name

    # name kwarg
    idx = _index
    s = idx.to_series(name='__test')

    assert s.values is not idx.values
    assert s.index is not idx
    assert s.name != idx.name


def test_shift(_index):

    # GH8083 test the base class for shift
    idx = _index
    pytest.raises(NotImplementedError, idx.shift, 1)
    pytest.raises(NotImplementedError, idx.shift, 1, 2)


def test_insert_base(named_index):

    for name, idx in compat.iteritems(named_index):
        result = idx[1:4]

        if not len(idx):
            continue

        # test 0th element
        assert idx[0:4].equals(result.insert(0, idx[0]))


def test_delete_base(named_index):

    for name, idx in compat.iteritems(named_index):

        if not len(idx):
            continue

        if isinstance(idx, RangeIndex):
            # tested in class
            continue

        expected = idx[1:]
        result = idx.delete(0)
        assert result.equals(expected)
        assert result.name == expected.name

        expected = idx[:-1]
        result = idx.delete(-1)
        assert result.equals(expected)
        assert result.name == expected.name

        with pytest.raises((IndexError, ValueError)):
            # either depending on numpy version
            result = idx.delete(len(idx))


def test_fillna(named_index):
    # GH 11343
    for name, index in named_index.items():
        if len(index) == 0:
            pass
        elif isinstance(index, MultiIndex):
            idx = index.copy()
            msg = "isna is not defined for MultiIndex"
            with tm.assert_raises_regex(NotImplementedError, msg):
                idx.fillna(idx[0])
        else:
            idx = index.copy()
            result = idx.fillna(idx[0])
            tm.assert_index_equal(result, idx)
            assert result is not idx

            msg = "'value' must be a scalar, passed: "
            with tm.assert_raises_regex(TypeError, msg):
                idx.fillna([idx[0]])

            idx = index.copy()
            values = idx.values

            if isinstance(index, DatetimeIndexOpsMixin):
                values[1] = iNaT
            elif isinstance(index, (Int64Index, UInt64Index)):
                continue
            else:
                values[1] = np.nan

            if isinstance(index, PeriodIndex):
                idx = index.__class__(values, freq=index.freq)
            else:
                idx = index.__class__(values)

            expected = np.array([False] * len(idx), dtype=bool)
            expected[1] = True
            tm.assert_numpy_array_equal(idx._isnan, expected)
            assert idx.hasnans


def test_putmask_with_wrong_mask(_index):
    # GH18368
    index = _index

    with pytest.raises(ValueError):
        index.putmask(np.ones(len(index) + 1, np.bool), 1)

    with pytest.raises(ValueError):
        index.putmask(np.ones(len(index) - 1, np.bool), 1)

    with pytest.raises(ValueError):
        index.putmask('foo', 1)
