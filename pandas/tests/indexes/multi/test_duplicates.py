# -*- coding: utf-8 -*-

from itertools import product
import pytest

import numpy as np

from pandas.compat import range, u
from pandas import MultiIndex, DatetimeIndex
from pandas._libs import hashtable
import pandas.util.testing as tm


@pytest.mark.parametrize('names', [None, ['first', 'second']])
def test_unique(names):
    mi = MultiIndex.from_arrays([[1, 2, 1, 2], [1, 1, 1, 2]], names=names)

    res = mi.unique()
    exp = MultiIndex.from_arrays([[1, 2, 2], [1, 1, 2]], names=mi.names)
    tm.assert_index_equal(res, exp)

    mi = MultiIndex.from_arrays([list('aaaa'), list('abab')],
                                names=names)
    res = mi.unique()
    exp = MultiIndex.from_arrays([list('aa'), list('ab')], names=mi.names)
    tm.assert_index_equal(res, exp)

    mi = MultiIndex.from_arrays([list('aaaa'), list('aaaa')], names=names)
    res = mi.unique()
    exp = MultiIndex.from_arrays([['a'], ['a']], names=mi.names)
    tm.assert_index_equal(res, exp)

    # GH #20568 - empty MI
    mi = MultiIndex.from_arrays([[], []], names=names)
    res = mi.unique()
    tm.assert_index_equal(mi, res)


def test_unique_datetimelike():
    idx1 = DatetimeIndex(['2015-01-01', '2015-01-01', '2015-01-01',
                          '2015-01-01', 'NaT', 'NaT'])
    idx2 = DatetimeIndex(['2015-01-01', '2015-01-01', '2015-01-02',
                          '2015-01-02', 'NaT', '2015-01-01'],
                         tz='Asia/Tokyo')
    result = MultiIndex.from_arrays([idx1, idx2]).unique()

    eidx1 = DatetimeIndex(['2015-01-01', '2015-01-01', 'NaT', 'NaT'])
    eidx2 = DatetimeIndex(['2015-01-01', '2015-01-02',
                           'NaT', '2015-01-01'],
                          tz='Asia/Tokyo')
    exp = MultiIndex.from_arrays([eidx1, eidx2])
    tm.assert_index_equal(result, exp)


@pytest.mark.parametrize('level', [0, 'first', 1, 'second'])
def test_unique_level(idx, level):
    # GH #17896 - with level= argument
    result = idx.unique(level=level)
    expected = idx.get_level_values(level).unique()
    tm.assert_index_equal(result, expected)

    # With already unique level
    mi = MultiIndex.from_arrays([[1, 3, 2, 4], [1, 3, 2, 5]],
                                names=['first', 'second'])
    result = mi.unique(level=level)
    expected = mi.get_level_values(level)
    tm.assert_index_equal(result, expected)

    # With empty MI
    mi = MultiIndex.from_arrays([[], []], names=['first', 'second'])
    result = mi.unique(level=level)
    expected = mi.get_level_values(level)


@pytest.mark.parametrize('dropna', [True, False])
def test_get_unique_index(idx, dropna):
    mi = idx[[0, 1, 0, 1, 1, 0, 0]]
    expected = mi._shallow_copy(mi[[0, 1]])

    result = mi._get_unique_index(dropna=dropna)
    assert result.unique
    tm.assert_index_equal(result, expected)


def test_duplicate_multiindex_labels():
    # GH 17464
    # Make sure that a MultiIndex with duplicate levels throws a ValueError
    with pytest.raises(ValueError):
        mi = MultiIndex([['A'] * 10, range(10)], [[0] * 10, range(10)])

    # And that using set_levels with duplicate levels fails
    mi = MultiIndex.from_arrays([['A', 'A', 'B', 'B', 'B'],
                                 [1, 2, 1, 2, 3]])
    with pytest.raises(ValueError):
        mi.set_levels([['A', 'B', 'A', 'A', 'B'], [2, 1, 3, -2, 5]],
                      inplace=True)


@pytest.mark.parametrize('names', [['a', 'b', 'a'], [1, 1, 2],
                                   [1, 'a', 1]])
def test_duplicate_level_names(names):
    # GH18872, GH19029
    mi = MultiIndex.from_product([[0, 1]] * 3, names=names)
    assert mi.names == names

    # With .rename()
    mi = MultiIndex.from_product([[0, 1]] * 3)
    mi = mi.rename(names)
    assert mi.names == names

    # With .rename(., level=)
    mi.rename(names[1], level=1, inplace=True)
    mi = mi.rename([names[0], names[2]], level=[0, 2])
    assert mi.names == names


def test_duplicate_meta_data():
    # GH 10115
    mi = MultiIndex(
        levels=[[0, 1], [0, 1, 2]],
        labels=[[0, 0, 0, 0, 1, 1, 1],
                [0, 1, 2, 0, 0, 1, 2]])

    for idx in [mi,
                mi.set_names([None, None]),
                mi.set_names([None, 'Num']),
                mi.set_names(['Upper', 'Num']), ]:
        assert idx.has_duplicates
        assert idx.drop_duplicates().names == idx.names


def test_has_duplicates(idx, idx_dup):
    # see fixtures
    assert idx.is_unique
    assert not idx.has_duplicates
    assert not idx_dup.is_unique
    assert idx_dup.has_duplicates

    mi = MultiIndex(levels=[[0, 1], [0, 1, 2]],
                    labels=[[0, 0, 0, 0, 1, 1, 1],
                            [0, 1, 2, 0, 0, 1, 2]])
    assert not mi.is_unique
    assert mi.has_duplicates


def test_has_duplicates_from_tuples():
    # GH 9075
    t = [(u('x'), u('out'), u('z'), 5, u('y'), u('in'), u('z'), 169),
         (u('x'), u('out'), u('z'), 7, u('y'), u('in'), u('z'), 119),
         (u('x'), u('out'), u('z'), 9, u('y'), u('in'), u('z'), 135),
         (u('x'), u('out'), u('z'), 13, u('y'), u('in'), u('z'), 145),
         (u('x'), u('out'), u('z'), 14, u('y'), u('in'), u('z'), 158),
         (u('x'), u('out'), u('z'), 16, u('y'), u('in'), u('z'), 122),
         (u('x'), u('out'), u('z'), 17, u('y'), u('in'), u('z'), 160),
         (u('x'), u('out'), u('z'), 18, u('y'), u('in'), u('z'), 180),
         (u('x'), u('out'), u('z'), 20, u('y'), u('in'), u('z'), 143),
         (u('x'), u('out'), u('z'), 21, u('y'), u('in'), u('z'), 128),
         (u('x'), u('out'), u('z'), 22, u('y'), u('in'), u('z'), 129),
         (u('x'), u('out'), u('z'), 25, u('y'), u('in'), u('z'), 111),
         (u('x'), u('out'), u('z'), 28, u('y'), u('in'), u('z'), 114),
         (u('x'), u('out'), u('z'), 29, u('y'), u('in'), u('z'), 121),
         (u('x'), u('out'), u('z'), 31, u('y'), u('in'), u('z'), 126),
         (u('x'), u('out'), u('z'), 32, u('y'), u('in'), u('z'), 155),
         (u('x'), u('out'), u('z'), 33, u('y'), u('in'), u('z'), 123),
         (u('x'), u('out'), u('z'), 12, u('y'), u('in'), u('z'), 144)]

    mi = MultiIndex.from_tuples(t)
    assert not mi.has_duplicates


def test_has_duplicates_overflow():
    # handle int64 overflow if possible
    def check(nlevels, with_nulls):
        labels = np.tile(np.arange(500), 2)
        level = np.arange(500)

        if with_nulls:  # inject some null values
            labels[500] = -1  # common nan value
            labels = [labels.copy() for i in range(nlevels)]
            for i in range(nlevels):
                labels[i][500 + i - nlevels // 2] = -1

            labels += [np.array([-1, 1]).repeat(500)]
        else:
            labels = [labels] * nlevels + [np.arange(2).repeat(500)]

        levels = [level] * nlevels + [[0, 1]]

        # no dups
        mi = MultiIndex(levels=levels, labels=labels)
        assert not mi.has_duplicates

        # with a dup
        if with_nulls:
            def f(a):
                return np.insert(a, 1000, a[0])
            labels = list(map(f, labels))
            mi = MultiIndex(levels=levels, labels=labels)
        else:
            values = mi.values.tolist()
            mi = MultiIndex.from_tuples(values + [values[0]])

        assert mi.has_duplicates

    # no overflow
    check(4, False)
    check(4, True)

    # overflow possible
    check(8, False)
    check(8, True)


@pytest.mark.parametrize('keep, expected', [
    ('first', np.array([False, False, False, True, True, False])),
    ('last', np.array([False, True, True, False, False, False])),
    (False, np.array([False, True, True, True, True, False]))
])
def test_duplicated(idx_dup, keep, expected):
    result = idx_dup.duplicated(keep=keep)
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize('keep', ['first', 'last', False])
def test_duplicated_large(keep):
    # GH 9125
    n, k = 200, 5000
    levels = [np.arange(n), tm.makeStringIndex(n), 1000 + np.arange(n)]
    labels = [np.random.choice(n, k * n) for lev in levels]
    mi = MultiIndex(levels=levels, labels=labels)

    result = mi.duplicated(keep=keep)
    expected = hashtable.duplicated_object(mi.values, keep=keep)
    tm.assert_numpy_array_equal(result, expected)


def test_get_duplicates():
    # GH5873
    for a in [101, 102]:
        mi = MultiIndex.from_arrays([[101, a], [3.5, np.nan]])
        assert not mi.has_duplicates

        with tm.assert_produces_warning(FutureWarning):
            # Deprecated - see GH20239
            assert mi.get_duplicates().equals(MultiIndex.from_arrays([[], []]))

        tm.assert_numpy_array_equal(mi.duplicated(),
                                    np.zeros(2, dtype='bool'))

    for n in range(1, 6):  # 1st level shape
        for m in range(1, 5):  # 2nd level shape
            # all possible unique combinations, including nan
            lab = product(range(-1, n), range(-1, m))
            mi = MultiIndex(levels=[list('abcde')[:n], list('WXYZ')[:m]],
                            labels=np.random.permutation(list(lab)).T)
            assert len(mi) == (n + 1) * (m + 1)
            assert not mi.has_duplicates

            with tm.assert_produces_warning(FutureWarning):
                # Deprecated - see GH20239
                assert mi.get_duplicates().equals(MultiIndex.from_arrays(
                    [[], []]))

            tm.assert_numpy_array_equal(mi.duplicated(),
                                        np.zeros(len(mi), dtype='bool'))
