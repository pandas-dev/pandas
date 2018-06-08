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


def test_get_loc(_index):
    assert _index.get_loc(('foo', 'two')) == 1
    assert _index.get_loc(('baz', 'two')) == 3
    pytest.raises(KeyError, _index.get_loc, ('bar', 'two'))
    pytest.raises(KeyError, _index.get_loc, 'quux')

    pytest.raises(NotImplementedError, _index.get_loc, 'foo',
                  method='nearest')

    # 3 levels
    index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
        lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
            [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])
    pytest.raises(KeyError, index.get_loc, (1, 1))
    assert index.get_loc((2, 0)) == slice(3, 5)


def test_get_loc_duplicates():
    index = Index([2, 2, 2, 2])
    result = index.get_loc(2)
    expected = slice(0, 4)
    assert result == expected
    # pytest.raises(Exception, index.get_loc, 2)

    index = Index(['c', 'a', 'a', 'b', 'b'])
    rs = index.get_loc('c')
    xp = 0
    assert rs == xp


def test_get_loc_level():
    index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
        lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
            [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])

    loc, new_index = index.get_loc_level((0, 1))
    expected = slice(1, 2)
    exp_index = index[expected].droplevel(0).droplevel(0)
    assert loc == expected
    assert new_index.equals(exp_index)

    loc, new_index = index.get_loc_level((0, 1, 0))
    expected = 1
    assert loc == expected
    assert new_index is None

    pytest.raises(KeyError, index.get_loc_level, (2, 2))

    index = MultiIndex(levels=[[2000], lrange(4)], labels=[np.array(
        [0, 0, 0, 0]), np.array([0, 1, 2, 3])])
    result, new_index = index.get_loc_level((2000, slice(None, None)))
    expected = slice(None, None)
    assert result == expected
    assert new_index.equals(index.droplevel(0))


@pytest.mark.parametrize('level', [0, 1])
@pytest.mark.parametrize('null_val', [np.nan, pd.NaT, None])
def test_get_loc_nan(level, null_val):
    # GH 18485 : NaN in MultiIndex
    levels = [['a', 'b'], ['c', 'd']]
    key = ['b', 'd']
    levels[level] = np.array([0, null_val], dtype=type(null_val))
    key[level] = null_val
    idx = MultiIndex.from_product(levels)
    assert idx.get_loc(tuple(key)) == 3


def test_get_loc_missing_nan():
    # GH 8569
    idx = MultiIndex.from_arrays([[1.0, 2.0], [3.0, 4.0]])
    assert isinstance(idx.get_loc(1), slice)
    pytest.raises(KeyError, idx.get_loc, 3)
    pytest.raises(KeyError, idx.get_loc, np.nan)
    pytest.raises(KeyError, idx.get_loc, [np.nan])


@pytest.mark.parametrize('dtype1', [int, float, bool, str])
@pytest.mark.parametrize('dtype2', [int, float, bool, str])
def test_get_loc_multiple_dtypes(dtype1, dtype2):
    # GH 18520
    levels = [np.array([0, 1]).astype(dtype1),
              np.array([0, 1]).astype(dtype2)]
    idx = pd.MultiIndex.from_product(levels)
    assert idx.get_loc(idx[2]) == 2


@pytest.mark.parametrize('level', [0, 1])
@pytest.mark.parametrize('dtypes', [[int, float], [float, int]])
def test_get_loc_implicit_cast(level, dtypes):
    # GH 18818, GH 15994 : as flat index, cast int to float and vice-versa
    levels = [['a', 'b'], ['c', 'd']]
    key = ['b', 'd']
    lev_dtype, key_dtype = dtypes
    levels[level] = np.array([0, 1], dtype=lev_dtype)
    key[level] = key_dtype(1)
    idx = MultiIndex.from_product(levels)
    assert idx.get_loc(tuple(key)) == 3


def test_get_loc_cast_bool():
    # GH 19086 : int is casted to bool, but not vice-versa
    levels = [[False, True], np.arange(2, dtype='int64')]
    idx = MultiIndex.from_product(levels)

    assert idx.get_loc((0, 1)) == 1
    assert idx.get_loc((1, 0)) == 2

    pytest.raises(KeyError, idx.get_loc, (False, True))
    pytest.raises(KeyError, idx.get_loc, (True, False))


def test_get_indexer():
    major_axis = Index(lrange(4))
    minor_axis = Index(lrange(2))

    major_labels = np.array([0, 0, 1, 2, 2, 3, 3], dtype=np.intp)
    minor_labels = np.array([0, 1, 0, 0, 1, 0, 1], dtype=np.intp)

    index = MultiIndex(levels=[major_axis, minor_axis],
                       labels=[major_labels, minor_labels])
    idx1 = index[:5]
    idx2 = index[[1, 3, 5]]

    r1 = idx1.get_indexer(idx2)
    assert_almost_equal(r1, np.array([1, 3, -1], dtype=np.intp))

    r1 = idx2.get_indexer(idx1, method='pad')
    e1 = np.array([-1, 0, 0, 1, 1], dtype=np.intp)
    assert_almost_equal(r1, e1)

    r2 = idx2.get_indexer(idx1[::-1], method='pad')
    assert_almost_equal(r2, e1[::-1])

    rffill1 = idx2.get_indexer(idx1, method='ffill')
    assert_almost_equal(r1, rffill1)

    r1 = idx2.get_indexer(idx1, method='backfill')
    e1 = np.array([0, 0, 1, 1, 2], dtype=np.intp)
    assert_almost_equal(r1, e1)

    r2 = idx2.get_indexer(idx1[::-1], method='backfill')
    assert_almost_equal(r2, e1[::-1])

    rbfill1 = idx2.get_indexer(idx1, method='bfill')
    assert_almost_equal(r1, rbfill1)

    # pass non-MultiIndex
    r1 = idx1.get_indexer(idx2.values)
    rexp1 = idx1.get_indexer(idx2)
    assert_almost_equal(r1, rexp1)

    r1 = idx1.get_indexer([1, 2, 3])
    assert (r1 == [-1, -1, -1]).all()

    # create index with duplicates
    idx1 = Index(lrange(10) + lrange(10))
    idx2 = Index(lrange(20))

    msg = "Reindexing only valid with uniquely valued Index objects"
    with tm.assert_raises_regex(InvalidIndexError, msg):
        idx1.get_indexer(idx2)


def test_get_indexer_nearest():
    midx = MultiIndex.from_tuples([('a', 1), ('b', 2)])
    with pytest.raises(NotImplementedError):
        midx.get_indexer(['a'], method='nearest')
    with pytest.raises(NotImplementedError):
        midx.get_indexer(['a'], method='pad', tolerance=2)
