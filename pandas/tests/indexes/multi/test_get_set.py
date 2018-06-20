# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest
from pandas import CategoricalIndex, Index, MultiIndex
from pandas.compat import range
from pandas.compat import PY3, PYPY, lrange, lzip, range, u
from pandas.util.testing import assert_almost_equal
from pandas.core.indexes.base import InvalidIndexError
from pandas import (CategoricalIndex, DatetimeIndex, Float64Index, Index,
                    Int64Index, IntervalIndex, MultiIndex, PeriodIndex,
                    RangeIndex, Series, TimedeltaIndex, UInt64Index, compat,
                    isna)


def test_get_loc(idx):
    assert idx.get_loc(('foo', 'two')) == 1
    assert idx.get_loc(('baz', 'two')) == 3
    pytest.raises(KeyError, idx.get_loc, ('bar', 'two'))
    pytest.raises(KeyError, idx.get_loc, 'quux')

    pytest.raises(NotImplementedError, idx.get_loc, 'foo',
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


def test_set_name_methods(idx, index_names):
    # so long as these are synonyms, we don't need to test set_names
    assert idx.rename == idx.set_names
    new_names = [name + "SUFFIX" for name in index_names]
    ind = idx.set_names(new_names)
    assert idx.names == index_names
    assert ind.names == new_names
    with tm.assert_raises_regex(ValueError, "^Length"):
        ind.set_names(new_names + new_names)
    new_names2 = [name + "SUFFIX2" for name in new_names]
    res = ind.set_names(new_names2, inplace=True)
    assert res is None
    assert ind.names == new_names2

    # set names for specific level (# GH7792)
    ind = idx.set_names(new_names[0], level=0)
    assert idx.names == index_names
    assert ind.names == [new_names[0], index_names[1]]

    res = ind.set_names(new_names2[0], level=0, inplace=True)
    assert res is None
    assert ind.names == [new_names2[0], index_names[1]]

    # set names for multiple levels
    ind = idx.set_names(new_names, level=[0, 1])
    assert idx.names == index_names
    assert ind.names == new_names

    res = ind.set_names(new_names2, level=[0, 1], inplace=True)
    assert res is None
    assert ind.names == new_names2


def test_set_levels_labels_directly(idx):
    # setting levels/labels directly raises AttributeError

    levels = idx.levels
    new_levels = [[lev + 'a' for lev in level] for level in levels]

    labels = idx.labels
    major_labels, minor_labels = labels
    major_labels = [(x + 1) % 3 for x in major_labels]
    minor_labels = [(x + 1) % 1 for x in minor_labels]
    new_labels = [major_labels, minor_labels]

    with pytest.raises(AttributeError):
        idx.levels = new_levels

    with pytest.raises(AttributeError):
        idx.labels = new_labels


def test_set_levels(idx):
    # side note - you probably wouldn't want to use levels and labels
    # directly like this - but it is possible.
    levels = idx.levels
    new_levels = [[lev + 'a' for lev in level] for level in levels]

    def assert_matching(actual, expected, check_dtype=False):
        # avoid specifying internal representation
        # as much as possible
        assert len(actual) == len(expected)
        for act, exp in zip(actual, expected):
            act = np.asarray(act)
            exp = np.asarray(exp)
            tm.assert_numpy_array_equal(act, exp, check_dtype=check_dtype)

    # level changing [w/o mutation]
    ind2 = idx.set_levels(new_levels)
    assert_matching(ind2.levels, new_levels)
    assert_matching(idx.levels, levels)

    # level changing [w/ mutation]
    ind2 = idx.copy()
    inplace_return = ind2.set_levels(new_levels, inplace=True)
    assert inplace_return is None
    assert_matching(ind2.levels, new_levels)

    # level changing specific level [w/o mutation]
    ind2 = idx.set_levels(new_levels[0], level=0)
    assert_matching(ind2.levels, [new_levels[0], levels[1]])
    assert_matching(idx.levels, levels)

    ind2 = idx.set_levels(new_levels[1], level=1)
    assert_matching(ind2.levels, [levels[0], new_levels[1]])
    assert_matching(idx.levels, levels)

    # level changing multiple levels [w/o mutation]
    ind2 = idx.set_levels(new_levels, level=[0, 1])
    assert_matching(ind2.levels, new_levels)
    assert_matching(idx.levels, levels)

    # level changing specific level [w/ mutation]
    ind2 = idx.copy()
    inplace_return = ind2.set_levels(new_levels[0], level=0, inplace=True)
    assert inplace_return is None
    assert_matching(ind2.levels, [new_levels[0], levels[1]])
    assert_matching(idx.levels, levels)

    ind2 = idx.copy()
    inplace_return = ind2.set_levels(new_levels[1], level=1, inplace=True)
    assert inplace_return is None
    assert_matching(ind2.levels, [levels[0], new_levels[1]])
    assert_matching(idx.levels, levels)

    # level changing multiple levels [w/ mutation]
    ind2 = idx.copy()
    inplace_return = ind2.set_levels(new_levels, level=[0, 1],
                                     inplace=True)
    assert inplace_return is None
    assert_matching(ind2.levels, new_levels)
    assert_matching(idx.levels, levels)

    # illegal level changing should not change levels
    # GH 13754
    original_index = idx.copy()
    for inplace in [True, False]:
        with tm.assert_raises_regex(ValueError, "^On"):
            idx.set_levels(['c'], level=0, inplace=inplace)
        assert_matching(idx.levels, original_index.levels,
                        check_dtype=True)

        with tm.assert_raises_regex(ValueError, "^On"):
            idx.set_labels([0, 1, 2, 3, 4, 5], level=0,
                              inplace=inplace)
        assert_matching(idx.labels, original_index.labels,
                        check_dtype=True)

        with tm.assert_raises_regex(TypeError, "^Levels"):
            idx.set_levels('c', level=0, inplace=inplace)
        assert_matching(idx.levels, original_index.levels,
                        check_dtype=True)

        with tm.assert_raises_regex(TypeError, "^Labels"):
            idx.set_labels(1, level=0, inplace=inplace)
        assert_matching(idx.labels, original_index.labels,
                        check_dtype=True)


def test_set_labels(idx):
    # side note - you probably wouldn't want to use levels and labels
    # directly like this - but it is possible.
    labels = idx.labels
    major_labels, minor_labels = labels
    major_labels = [(x + 1) % 3 for x in major_labels]
    minor_labels = [(x + 1) % 1 for x in minor_labels]
    new_labels = [major_labels, minor_labels]

    def assert_matching(actual, expected):
        # avoid specifying internal representation
        # as much as possible
        assert len(actual) == len(expected)
        for act, exp in zip(actual, expected):
            act = np.asarray(act)
            exp = np.asarray(exp, dtype=np.int8)
            tm.assert_numpy_array_equal(act, exp)

    # label changing [w/o mutation]
    ind2 = idx.set_labels(new_labels)
    assert_matching(ind2.labels, new_labels)
    assert_matching(idx.labels, labels)

    # label changing [w/ mutation]
    ind2 = idx.copy()
    inplace_return = ind2.set_labels(new_labels, inplace=True)
    assert inplace_return is None
    assert_matching(ind2.labels, new_labels)

    # label changing specific level [w/o mutation]
    ind2 = idx.set_labels(new_labels[0], level=0)
    assert_matching(ind2.labels, [new_labels[0], labels[1]])
    assert_matching(idx.labels, labels)

    ind2 = idx.set_labels(new_labels[1], level=1)
    assert_matching(ind2.labels, [labels[0], new_labels[1]])
    assert_matching(idx.labels, labels)

    # label changing multiple levels [w/o mutation]
    ind2 = idx.set_labels(new_labels, level=[0, 1])
    assert_matching(ind2.labels, new_labels)
    assert_matching(idx.labels, labels)

    # label changing specific level [w/ mutation]
    ind2 = idx.copy()
    inplace_return = ind2.set_labels(new_labels[0], level=0, inplace=True)
    assert inplace_return is None
    assert_matching(ind2.labels, [new_labels[0], labels[1]])
    assert_matching(idx.labels, labels)

    ind2 = idx.copy()
    inplace_return = ind2.set_labels(new_labels[1], level=1, inplace=True)
    assert inplace_return is None
    assert_matching(ind2.labels, [labels[0], new_labels[1]])
    assert_matching(idx.labels, labels)

    # label changing multiple levels [w/ mutation]
    ind2 = idx.copy()
    inplace_return = ind2.set_labels(new_labels, level=[0, 1],
                                     inplace=True)
    assert inplace_return is None
    assert_matching(ind2.labels, new_labels)
    assert_matching(idx.labels, labels)

    # label changing for levels of different magnitude of categories
    ind = pd.MultiIndex.from_tuples([(0, i) for i in range(130)])
    new_labels = range(129, -1, -1)
    expected = pd.MultiIndex.from_tuples(
        [(0, i) for i in new_labels])

    # [w/o mutation]
    result = ind.set_labels(labels=new_labels, level=1)
    assert result.equals(expected)

    # [w/ mutation]
    result = ind.copy()
    result.set_labels(labels=new_labels, level=1, inplace=True)
    assert result.equals(expected)


def test_set_levels_labels_names_bad_input(idx):
    levels, labels = idx.levels, idx.labels
    names = idx.names

    with tm.assert_raises_regex(ValueError, 'Length of levels'):
        idx.set_levels([levels[0]])

    with tm.assert_raises_regex(ValueError, 'Length of labels'):
        idx.set_labels([labels[0]])

    with tm.assert_raises_regex(ValueError, 'Length of names'):
        idx.set_names([names[0]])

    # shouldn't scalar data error, instead should demand list-like
    with tm.assert_raises_regex(TypeError, 'list of lists-like'):
        idx.set_levels(levels[0])

    # shouldn't scalar data error, instead should demand list-like
    with tm.assert_raises_regex(TypeError, 'list of lists-like'):
        idx.set_labels(labels[0])

    # shouldn't scalar data error, instead should demand list-like
    with tm.assert_raises_regex(TypeError, 'list-like'):
        idx.set_names(names[0])

    # should have equal lengths
    with tm.assert_raises_regex(TypeError, 'list of lists-like'):
        idx.set_levels(levels[0], level=[0, 1])

    with tm.assert_raises_regex(TypeError, 'list-like'):
        idx.set_levels(levels, level=0)

    # should have equal lengths
    with tm.assert_raises_regex(TypeError, 'list of lists-like'):
        idx.set_labels(labels[0], level=[0, 1])

    with tm.assert_raises_regex(TypeError, 'list-like'):
        idx.set_labels(labels, level=0)

    # should have equal lengths
    with tm.assert_raises_regex(ValueError, 'Length of names'):
        idx.set_names(names[0], level=[0, 1])

    with tm.assert_raises_regex(TypeError, 'string'):
        idx.set_names(names, level=0)


@pytest.mark.parametrize('inplace', [True, False])
def test_set_names_with_nlevel_1(inplace):
    # GH 21149
    # Ensure that .set_names for MultiIndex with
    # nlevels == 1 does not raise any errors
    expected = pd.MultiIndex(levels=[[0, 1]],
                             labels=[[0, 1]],
                             names=['first'])
    m = pd.MultiIndex.from_product([[0, 1]])
    result = m.set_names('first', level=0, inplace=inplace)

    if inplace:
        result = m

    tm.assert_index_equal(result, expected)


def test_set_levels_categorical():
    # GH13854
    index = MultiIndex.from_arrays([list("xyzx"), [0, 1, 2, 3]])
    for ordered in [False, True]:
        cidx = CategoricalIndex(list("bac"), ordered=ordered)
        result = index.set_levels(cidx, 0)
        expected = MultiIndex(levels=[cidx, [0, 1, 2, 3]],
                              labels=index.labels)
        tm.assert_index_equal(result, expected)

        result_lvl = result.get_level_values(0)
        expected_lvl = CategoricalIndex(list("bacb"),
                                        categories=cidx.categories,
                                        ordered=cidx.ordered)
        tm.assert_index_equal(result_lvl, expected_lvl)


def test_set_value_keeps_names():
    # motivating example from #3742
    lev1 = ['hans', 'hans', 'hans', 'grethe', 'grethe', 'grethe']
    lev2 = ['1', '2', '3'] * 2
    idx = pd.MultiIndex.from_arrays([lev1, lev2], names=['Name', 'Number'])
    df = pd.DataFrame(
        np.random.randn(6, 4),
        columns=['one', 'two', 'three', 'four'],
        index=idx)
    df = df.sort_index()
    assert df._is_copy is None
    assert df.index.names == ('Name', 'Number')
    df.at[('grethe', '4'), 'one'] = 99.34
    assert df._is_copy is None
    assert df.index.names == ('Name', 'Number')


def test_get_level_number_integer(idx):
    idx.names = [1, 0]
    assert idx._get_level_number(1) == 0
    assert idx._get_level_number(0) == 1
    pytest.raises(IndexError, idx._get_level_number, 2)
    tm.assert_raises_regex(KeyError, 'Level fourth not found',
                           idx._get_level_number, 'fourth')


def test_get_level_values(idx):
    result = idx.get_level_values(0)
    expected = Index(['foo', 'foo', 'bar', 'baz', 'qux', 'qux'],
                     name='first')
    tm.assert_index_equal(result, expected)
    assert result.name == 'first'

    result = idx.get_level_values('first')
    expected = idx.get_level_values(0)
    tm.assert_index_equal(result, expected)

    # GH 10460
    index = MultiIndex(
        levels=[CategoricalIndex(['A', 'B']),
                CategoricalIndex([1, 2, 3])],
        labels=[np.array([0, 0, 0, 1, 1, 1]),
                np.array([0, 1, 2, 0, 1, 2])])

    exp = CategoricalIndex(['A', 'A', 'A', 'B', 'B', 'B'])
    tm.assert_index_equal(index.get_level_values(0), exp)
    exp = CategoricalIndex([1, 2, 3, 1, 2, 3])
    tm.assert_index_equal(index.get_level_values(1), exp)


def test_getitem(idx):
    # scalar
    assert idx[2] == ('bar', 'one')

    # slice
    result = idx[2:5]
    expected = idx[[2, 3, 4]]
    assert result.equals(expected)

    # boolean
    result = idx[[True, False, True, False, True, True]]
    result2 = idx[np.array([True, False, True, False, True, True])]
    expected = idx[[0, 2, 4, 5]]
    assert result.equals(expected)
    assert result2.equals(expected)


def test_getitem_group_select(idx):
    sorted_idx, _ = idx.sortlevel(0)
    assert sorted_idx.get_loc('baz') == slice(3, 4)
    assert sorted_idx.get_loc('foo') == slice(0, 2)


def test_get_value_duplicates():
    index = MultiIndex(levels=[['D', 'B', 'C'],
                               [0, 26, 27, 37, 57, 67, 75, 82]],
                       labels=[[0, 0, 0, 1, 2, 2, 2, 2, 2, 2],
                               [1, 3, 4, 6, 0, 2, 2, 3, 5, 7]],
                       names=['tag', 'day'])

    assert index.get_loc('D') == slice(0, 3)
    with pytest.raises(KeyError):
        index._engine.get_value(np.array([]), 'D')


def test_get_level_values_all_na():
    # GH 17924 when level entirely consists of nan
    arrays = [[np.nan, np.nan, np.nan], ['a', np.nan, 1]]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(0)
    expected = pd.Index([np.nan, np.nan, np.nan], dtype=np.float64)
    tm.assert_index_equal(result, expected)

    result = index.get_level_values(1)
    expected = pd.Index(['a', np.nan, 1], dtype=object)
    tm.assert_index_equal(result, expected)


def test_get_level_values_int_with_na():
    # GH 17924
    arrays = [['a', 'b', 'b'], [1, np.nan, 2]]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(1)
    expected = Index([1, np.nan, 2])
    tm.assert_index_equal(result, expected)

    arrays = [['a', 'b', 'b'], [np.nan, np.nan, 2]]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(1)
    expected = Index([np.nan, np.nan, 2])
    tm.assert_index_equal(result, expected)


def test_get_level_values_na():
    arrays = [[np.nan, np.nan, np.nan], ['a', np.nan, 1]]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(0)
    expected = pd.Index([np.nan, np.nan, np.nan])
    tm.assert_index_equal(result, expected)

    result = index.get_level_values(1)
    expected = pd.Index(['a', np.nan, 1])
    tm.assert_index_equal(result, expected)

    arrays = [['a', 'b', 'b'], pd.DatetimeIndex([0, 1, pd.NaT])]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(1)
    expected = pd.DatetimeIndex([0, 1, pd.NaT])
    tm.assert_index_equal(result, expected)

    arrays = [[], []]
    index = pd.MultiIndex.from_arrays(arrays)
    result = index.get_level_values(0)
    expected = pd.Index([], dtype=object)
    tm.assert_index_equal(result, expected)


def test_get_unique_index(idx):
    idx = idx[[0, 1, 0, 1, 1, 0, 0]]
    expected = idx._shallow_copy(idx[[0, 1]])

    for dropna in [False, True]:
        result = idx._get_unique_index(dropna=dropna)
        assert result.unique
        tm.assert_index_equal(result, expected)


def test_get_indexer_consistency(named_index):
    # See GH 16819
    for name, index in named_index.items():
        if isinstance(index, IntervalIndex):
            continue

        if index.is_unique or isinstance(index, CategoricalIndex):
            indexer = index.get_indexer(index[0:2])
            assert isinstance(indexer, np.ndarray)
            assert indexer.dtype == np.intp
        else:
            e = "Reindexing only valid with uniquely valued Index objects"
            with tm.assert_raises_regex(InvalidIndexError, e):
                indexer = index.get_indexer(index[0:2])

        indexer, _ = index.get_indexer_non_unique(index[0:2])
        assert isinstance(indexer, np.ndarray)
        assert indexer.dtype == np.intp
