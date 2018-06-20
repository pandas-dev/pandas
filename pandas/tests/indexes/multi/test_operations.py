# -*- coding: utf-8 -*-

import warnings
from itertools import product

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest
from pandas import (DataFrame, DatetimeIndex, Float64Index, Index, Int64Index,
                    MultiIndex, PeriodIndex, TimedeltaIndex, UInt64Index,
                    compat, date_range, period_range)
from pandas.compat import lrange, range, u
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.core.indexes.datetimelike import DatetimeIndexOpsMixin
from pandas.util.testing import assert_copy


def check_level_names(index, names):
    assert [level.name for level in index.levels] == list(names)


def test_insert(idx):
    # key contained in all levels
    new_index = idx.insert(0, ('bar', 'two'))
    assert new_index.equal_levels(idx)
    assert new_index[0] == ('bar', 'two')

    # key not contained in all levels
    new_index = idx.insert(0, ('abc', 'three'))

    exp0 = Index(list(idx.levels[0]) + ['abc'], name='first')
    tm.assert_index_equal(new_index.levels[0], exp0)

    exp1 = Index(list(idx.levels[1]) + ['three'], name='second')
    tm.assert_index_equal(new_index.levels[1], exp1)
    assert new_index[0] == ('abc', 'three')

    # key wrong length
    msg = "Item must have length equal to number of levels"
    with tm.assert_raises_regex(ValueError, msg):
        idx.insert(0, ('foo2',))

    left = pd.DataFrame([['a', 'b', 0], ['b', 'd', 1]],
                        columns=['1st', '2nd', '3rd'])
    left.set_index(['1st', '2nd'], inplace=True)
    ts = left['3rd'].copy(deep=True)

    left.loc[('b', 'x'), '3rd'] = 2
    left.loc[('b', 'a'), '3rd'] = -1
    left.loc[('b', 'b'), '3rd'] = 3
    left.loc[('a', 'x'), '3rd'] = 4
    left.loc[('a', 'w'), '3rd'] = 5
    left.loc[('a', 'a'), '3rd'] = 6

    ts.loc[('b', 'x')] = 2
    ts.loc['b', 'a'] = -1
    ts.loc[('b', 'b')] = 3
    ts.loc['a', 'x'] = 4
    ts.loc[('a', 'w')] = 5
    ts.loc['a', 'a'] = 6

    right = pd.DataFrame([['a', 'b', 0], ['b', 'd', 1], ['b', 'x', 2],
                          ['b', 'a', -1], ['b', 'b', 3], ['a', 'x', 4],
                          ['a', 'w', 5], ['a', 'a', 6]],
                         columns=['1st', '2nd', '3rd'])
    right.set_index(['1st', '2nd'], inplace=True)
    # FIXME data types changes to float because
    # of intermediate nan insertion;
    tm.assert_frame_equal(left, right, check_dtype=False)
    tm.assert_series_equal(ts, right['3rd'])

    # GH9250
    idx = [('test1', i) for i in range(5)] + \
        [('test2', i) for i in range(6)] + \
        [('test', 17), ('test', 18)]

    left = pd.Series(np.linspace(0, 10, 11),
                     pd.MultiIndex.from_tuples(idx[:-2]))

    left.loc[('test', 17)] = 11
    left.loc[('test', 18)] = 12

    right = pd.Series(np.linspace(0, 12, 13),
                      pd.MultiIndex.from_tuples(idx))

    tm.assert_series_equal(left, right)


def test_bounds(idx):
    idx._bounds


def test_append(idx):
    result = idx[:3].append(idx[3:])
    assert result.equals(idx)

    foos = [idx[:1], idx[1:3], idx[3:]]
    result = foos[0].append(foos[1:])
    assert result.equals(idx)

    # empty
    result = idx.append([])
    assert result.equals(idx)


def test_groupby(idx):
    groups = idx.groupby(np.array([1, 1, 1, 2, 2, 2]))
    labels = idx.get_values().tolist()
    exp = {1: labels[:3], 2: labels[3:]}
    tm.assert_dict_equal(groups, exp)

    # GH5620
    groups = idx.groupby(idx)
    exp = {key: [key] for key in idx}
    tm.assert_dict_equal(groups, exp)


def test_truncate():
    major_axis = Index(lrange(4))
    minor_axis = Index(lrange(2))

    major_labels = np.array([0, 0, 1, 2, 3, 3])
    minor_labels = np.array([0, 1, 0, 1, 0, 1])

    index = MultiIndex(levels=[major_axis, minor_axis],
                       labels=[major_labels, minor_labels])

    result = index.truncate(before=1)
    assert 'foo' not in result.levels[0]
    assert 1 in result.levels[0]

    result = index.truncate(after=1)
    assert 2 not in result.levels[0]
    assert 1 in result.levels[0]

    result = index.truncate(before=1, after=2)
    assert len(result.levels[0]) == 2

    # after < before
    pytest.raises(ValueError, index.truncate, 3, 1)


def test_where():
    i = MultiIndex.from_tuples([('A', 1), ('A', 2)])

    def f():
        i.where(True)

    pytest.raises(NotImplementedError, f)


def test_where_array_like():
    i = MultiIndex.from_tuples([('A', 1), ('A', 2)])
    klasses = [list, tuple, np.array, pd.Series]
    cond = [False, True]

    for klass in klasses:
        def f():
            return i.where(klass(cond))
        pytest.raises(NotImplementedError, f)


def test_reorder_levels(idx):
    # this blows up
    tm.assert_raises_regex(IndexError, '^Too many levels',
                           idx.reorder_levels, [2, 1, 0])


def test_astype(idx):
    expected = idx.copy()
    actual = idx.astype('O')
    assert_copy(actual.levels, expected.levels)
    assert_copy(actual.labels, expected.labels)
    check_level_names(actual, expected.names)

    with tm.assert_raises_regex(TypeError, "^Setting.*dtype.*object"):
        idx.astype(np.dtype(int))


@pytest.mark.parametrize('ordered', [True, False])
def test_astype_category(idx, ordered):
    # GH 18630
    msg = '> 1 ndim Categorical are not supported at this time'
    with tm.assert_raises_regex(NotImplementedError, msg):
        idx.astype(CategoricalDtype(ordered=ordered))

    if ordered is False:
        # dtype='category' defaults to ordered=False, so only test once
        with tm.assert_raises_regex(NotImplementedError, msg):
            idx.astype('category')


@pytest.mark.parametrize('first_type,second_type', [
    ('int64', 'int64'),
    ('datetime64[D]', 'str')])
def test_remove_unused_levels_large(first_type, second_type):
    # GH16556

    # because tests should be deterministic (and this test in particular
    # checks that levels are removed, which is not the case for every
    # random input):
    rng = np.random.RandomState(4)  # seed is arbitrary value that works

    size = 1 << 16
    df = DataFrame(dict(
        first=rng.randint(0, 1 << 13, size).astype(first_type),
        second=rng.randint(0, 1 << 10, size).astype(second_type),
        third=rng.rand(size)))
    df = df.groupby(['first', 'second']).sum()
    df = df[df.third < 0.1]

    result = df.index.remove_unused_levels()
    assert len(result.levels[0]) < len(df.index.levels[0])
    assert len(result.levels[1]) < len(df.index.levels[1])
    assert result.equals(df.index)

    expected = df.reset_index().set_index(['first', 'second']).index
    tm.assert_index_equal(result, expected)


def test_repeat():
    reps = 2
    numbers = [1, 2, 3]
    names = np.array(['foo', 'bar'])

    m = MultiIndex.from_product([
        numbers, names], names=names)
    expected = MultiIndex.from_product([
        numbers, names.repeat(reps)], names=names)
    tm.assert_index_equal(m.repeat(reps), expected)

    with tm.assert_produces_warning(FutureWarning):
        result = m.repeat(n=reps)
        tm.assert_index_equal(result, expected)


def test_numpy_repeat():
    reps = 2
    numbers = [1, 2, 3]
    names = np.array(['foo', 'bar'])

    m = MultiIndex.from_product([
        numbers, names], names=names)
    expected = MultiIndex.from_product([
        numbers, names.repeat(reps)], names=names)
    tm.assert_index_equal(np.repeat(m, reps), expected)

    msg = "the 'axis' parameter is not supported"
    tm.assert_raises_regex(
        ValueError, msg, np.repeat, m, reps, axis=1)


def test_append_mixed_dtypes():
    # GH 13660
    dti = date_range('2011-01-01', freq='M', periods=3, )
    dti_tz = date_range('2011-01-01', freq='M', periods=3, tz='US/Eastern')
    pi = period_range('2011-01', freq='M', periods=3)

    mi = MultiIndex.from_arrays([[1, 2, 3],
                                 [1.1, np.nan, 3.3],
                                 ['a', 'b', 'c'],
                                 dti, dti_tz, pi])
    assert mi.nlevels == 6

    res = mi.append(mi)
    exp = MultiIndex.from_arrays([[1, 2, 3, 1, 2, 3],
                                  [1.1, np.nan, 3.3, 1.1, np.nan, 3.3],
                                  ['a', 'b', 'c', 'a', 'b', 'c'],
                                  dti.append(dti),
                                  dti_tz.append(dti_tz),
                                  pi.append(pi)])
    tm.assert_index_equal(res, exp)

    other = MultiIndex.from_arrays([['x', 'y', 'z'], ['x', 'y', 'z'],
                                    ['x', 'y', 'z'], ['x', 'y', 'z'],
                                    ['x', 'y', 'z'], ['x', 'y', 'z']])

    res = mi.append(other)
    exp = MultiIndex.from_arrays([[1, 2, 3, 'x', 'y', 'z'],
                                  [1.1, np.nan, 3.3, 'x', 'y', 'z'],
                                  ['a', 'b', 'c', 'x', 'y', 'z'],
                                  dti.append(pd.Index(['x', 'y', 'z'])),
                                  dti_tz.append(pd.Index(['x', 'y', 'z'])),
                                  pi.append(pd.Index(['x', 'y', 'z']))])
    tm.assert_index_equal(res, exp)


def test_take(named_index):
    indexer = [4, 3, 0, 2]
    for k, ind in named_index.items():

        # separate
        if k in ['boolIndex', 'tuples', 'empty']:
            continue

        result = ind.take(indexer)
        expected = ind[indexer]
        assert result.equals(expected)

        if not isinstance(ind,
                          (DatetimeIndex, PeriodIndex, TimedeltaIndex)):
            # GH 10791
            with pytest.raises(AttributeError):
                ind.freq


def test_take_invalid_kwargs(idx):
    idx = idx
    indices = [1, 2]

    msg = r"take\(\) got an unexpected keyword argument 'foo'"
    tm.assert_raises_regex(TypeError, msg, idx.take,
                           indices, foo=2)

    msg = "the 'out' parameter is not supported"
    tm.assert_raises_regex(ValueError, msg, idx.take,
                           indices, out=indices)

    msg = "the 'mode' parameter is not supported"
    tm.assert_raises_regex(ValueError, msg, idx.take,
                           indices, mode='clip')


def test_take_fill_value():
    # GH 12631
    vals = [['A', 'B'],
            [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02')]]
    idx = pd.MultiIndex.from_product(vals, names=['str', 'dt'])

    result = idx.take(np.array([1, 0, -1]))
    exp_vals = [('A', pd.Timestamp('2011-01-02')),
                ('A', pd.Timestamp('2011-01-01')),
                ('B', pd.Timestamp('2011-01-02'))]
    expected = pd.MultiIndex.from_tuples(exp_vals, names=['str', 'dt'])
    tm.assert_index_equal(result, expected)

    # fill_value
    result = idx.take(np.array([1, 0, -1]), fill_value=True)
    exp_vals = [('A', pd.Timestamp('2011-01-02')),
                ('A', pd.Timestamp('2011-01-01')),
                (np.nan, pd.NaT)]
    expected = pd.MultiIndex.from_tuples(exp_vals, names=['str', 'dt'])
    tm.assert_index_equal(result, expected)

    # allow_fill=False
    result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                      fill_value=True)
    exp_vals = [('A', pd.Timestamp('2011-01-02')),
                ('A', pd.Timestamp('2011-01-01')),
                ('B', pd.Timestamp('2011-01-02'))]
    expected = pd.MultiIndex.from_tuples(exp_vals, names=['str', 'dt'])
    tm.assert_index_equal(result, expected)

    msg = ('When allow_fill=True and fill_value is not None, '
           'all indices must be >= -1')
    with tm.assert_raises_regex(ValueError, msg):
        idx.take(np.array([1, 0, -2]), fill_value=True)
    with tm.assert_raises_regex(ValueError, msg):
        idx.take(np.array([1, 0, -5]), fill_value=True)

    with pytest.raises(IndexError):
        idx.take(np.array([1, -5]))


def test_iter(idx):
    result = list(idx)
    expected = [('foo', 'one'), ('foo', 'two'), ('bar', 'one'),
                ('baz', 'two'), ('qux', 'one'), ('qux', 'two')]
    assert result == expected


def test_sub(idx):

    first = idx

    # - now raises (previously was set op difference)
    with pytest.raises(TypeError):
        first - idx[-3:]
    with pytest.raises(TypeError):
        idx[-3:] - first
    with pytest.raises(TypeError):
        idx[-3:] - first.tolist()
    with pytest.raises(TypeError):
        first.tolist() - idx[-3:]


def test_nlevels(idx):
    assert idx.nlevels == 2


def test_argsort(idx):
    result = idx.argsort()
    expected = idx.values.argsort()
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize('level0', [['a', 'd', 'b'],
                                    ['a', 'd', 'b', 'unused']])
@pytest.mark.parametrize('level1', [['w', 'x', 'y', 'z'],
                                    ['w', 'x', 'y', 'z', 'unused']])
def test_remove_unused_nan(level0, level1):
    # GH 18417
    mi = pd.MultiIndex(levels=[level0, level1],
                       labels=[[0, 2, -1, 1, -1], [0, 1, 2, 3, 2]])

    result = mi.remove_unused_levels()
    tm.assert_index_equal(result, mi)
    for level in 0, 1:
        assert('unused' not in result.levels[level])


@pytest.mark.parametrize('names', [None, ['first', 'second']])
def test_unique(names):
    mi = pd.MultiIndex.from_arrays([[1, 2, 1, 2], [1, 1, 1, 2]],
                                   names=names)

    res = mi.unique()
    exp = pd.MultiIndex.from_arrays([[1, 2, 2], [1, 1, 2]], names=mi.names)
    tm.assert_index_equal(res, exp)

    mi = pd.MultiIndex.from_arrays([list('aaaa'), list('abab')],
                                   names=names)
    res = mi.unique()
    exp = pd.MultiIndex.from_arrays([list('aa'), list('ab')],
                                    names=mi.names)
    tm.assert_index_equal(res, exp)

    mi = pd.MultiIndex.from_arrays([list('aaaa'), list('aaaa')],
                                   names=names)
    res = mi.unique()
    exp = pd.MultiIndex.from_arrays([['a'], ['a']], names=mi.names)
    tm.assert_index_equal(res, exp)

    # GH #20568 - empty MI
    mi = pd.MultiIndex.from_arrays([[], []], names=names)
    res = mi.unique()
    tm.assert_index_equal(mi, res)


def test_unique_datetimelike():
    idx1 = pd.DatetimeIndex(['2015-01-01', '2015-01-01', '2015-01-01',
                             '2015-01-01', 'NaT', 'NaT'])
    idx2 = pd.DatetimeIndex(['2015-01-01', '2015-01-01', '2015-01-02',
                             '2015-01-02', 'NaT', '2015-01-01'],
                            tz='Asia/Tokyo')
    result = pd.MultiIndex.from_arrays([idx1, idx2]).unique()

    eidx1 = pd.DatetimeIndex(['2015-01-01', '2015-01-01', 'NaT', 'NaT'])
    eidx2 = pd.DatetimeIndex(['2015-01-01', '2015-01-02',
                              'NaT', '2015-01-01'],
                             tz='Asia/Tokyo')
    exp = pd.MultiIndex.from_arrays([eidx1, eidx2])
    tm.assert_index_equal(result, exp)


@pytest.mark.parametrize('level', [0, 'first', 1, 'second'])
def test_unique_level(idx, level):
    # GH #17896 - with level= argument
    result = idx.unique(level=level)
    expected = idx.get_level_values(level).unique()
    tm.assert_index_equal(result, expected)

    # With already unique level
    mi = pd.MultiIndex.from_arrays([[1, 3, 2, 4], [1, 3, 2, 5]],
                                   names=['first', 'second'])
    result = mi.unique(level=level)
    expected = mi.get_level_values(level)
    tm.assert_index_equal(result, expected)

    # With empty MI
    mi = pd.MultiIndex.from_arrays([[], []], names=['first', 'second'])
    result = mi.unique(level=level)
    expected = mi.get_level_values(level)


def test_multiindex_compare():
        # GH 21149
        # Ensure comparison operations for MultiIndex with nlevels == 1
        # behave consistently with those for MultiIndex with nlevels > 1

    midx = pd.MultiIndex.from_product([[0, 1]])

    # Equality self-test: MultiIndex object vs self
    expected = pd.Series([True, True])
    result = pd.Series(midx == midx)
    tm.assert_series_equal(result, expected)

    # Greater than comparison: MultiIndex object vs self
    expected = pd.Series([False, False])
    result = pd.Series(midx > midx)
    tm.assert_series_equal(result, expected)


def test_duplicate_multiindex_labels():
    # GH 17464
    # Make sure that a MultiIndex with duplicate levels throws a ValueError
    with pytest.raises(ValueError):
        ind = pd.MultiIndex([['A'] * 10, range(10)], [[0] * 10, range(10)])

    # And that using set_levels with duplicate levels fails
    ind = MultiIndex.from_arrays([['A', 'A', 'B', 'B', 'B'],
                                  [1, 2, 1, 2, 3]])
    with pytest.raises(ValueError):
        ind.set_levels([['A', 'B', 'A', 'A', 'B'], [2, 1, 3, -2, 5]],
                       inplace=True)


@pytest.mark.parametrize('names', [['a', 'b', 'a'], ['1', '1', '2'],
                                   ['1', 'a', '1']])
def test_duplicate_level_names(names):
    # GH18872
    pytest.raises(ValueError, pd.MultiIndex.from_product,
                  [[0, 1]] * 3, names=names)

    # With .rename()
    mi = pd.MultiIndex.from_product([[0, 1]] * 3)
    tm.assert_raises_regex(ValueError, "Duplicated level name:",
                           mi.rename, names)

    # With .rename(., level=)
    mi.rename(names[0], level=1, inplace=True)
    tm.assert_raises_regex(ValueError, "Duplicated level name:",
                           mi.rename, names[:2], level=[0, 2])


def test_duplicate_meta_data():
    # GH 10115
    index = MultiIndex(
        levels=[[0, 1], [0, 1, 2]],
        labels=[[0, 0, 0, 0, 1, 1, 1],
                [0, 1, 2, 0, 0, 1, 2]])

    for idx in [index,
                index.set_names([None, None]),
                index.set_names([None, 'Num']),
                index.set_names(['Upper', 'Num']), ]:
        assert idx.has_duplicates
        assert idx.drop_duplicates().names == idx.names


def test_duplicates(idx):
    assert not idx.has_duplicates
    assert idx.append(idx).has_duplicates

    index = MultiIndex(levels=[[0, 1], [0, 1, 2]], labels=[
        [0, 0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 0, 1, 2]])
    assert index.has_duplicates

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

    index = pd.MultiIndex.from_tuples(t)
    assert not index.has_duplicates

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
        index = MultiIndex(levels=levels, labels=labels)
        assert not index.has_duplicates

        # with a dup
        if with_nulls:
            def f(a):
                return np.insert(a, 1000, a[0])
            labels = list(map(f, labels))
            index = MultiIndex(levels=levels, labels=labels)
        else:
            values = index.values.tolist()
            index = MultiIndex.from_tuples(values + [values[0]])

        assert index.has_duplicates

    # no overflow
    check(4, False)
    check(4, True)

    # overflow possible
    check(8, False)
    check(8, True)

    # GH 9125
    n, k = 200, 5000
    levels = [np.arange(n), tm.makeStringIndex(n), 1000 + np.arange(n)]
    labels = [np.random.choice(n, k * n) for lev in levels]
    mi = MultiIndex(levels=levels, labels=labels)

    for keep in ['first', 'last', False]:
        left = mi.duplicated(keep=keep)
        right = pd._libs.hashtable.duplicated_object(mi.values, keep=keep)
        tm.assert_numpy_array_equal(left, right)

    # GH5873
    for a in [101, 102]:
        mi = MultiIndex.from_arrays([[101, a], [3.5, np.nan]])
        assert not mi.has_duplicates

        with warnings.catch_warnings(record=True):
            # Deprecated - see GH20239
            assert mi.get_duplicates().equals(MultiIndex.from_arrays(
                [[], []]))

        tm.assert_numpy_array_equal(mi.duplicated(), np.zeros(
            2, dtype='bool'))

    for n in range(1, 6):  # 1st level shape
        for m in range(1, 5):  # 2nd level shape
            # all possible unique combinations, including nan
            lab = product(range(-1, n), range(-1, m))
            mi = MultiIndex(levels=[list('abcde')[:n], list('WXYZ')[:m]],
                            labels=np.random.permutation(list(lab)).T)
            assert len(mi) == (n + 1) * (m + 1)
            assert not mi.has_duplicates

            with warnings.catch_warnings(record=True):
                # Deprecated - see GH20239
                assert mi.get_duplicates().equals(MultiIndex.from_arrays(
                    [[], []]))

            tm.assert_numpy_array_equal(mi.duplicated(), np.zeros(
                len(mi), dtype='bool'))


def test_map(idx):
    # callable
    index = idx

    # we don't infer UInt64
    if isinstance(index, pd.UInt64Index):
        expected = index.astype('int64')
    else:
        expected = index

    result = index.map(lambda x: x)
    tm.assert_index_equal(result, expected)


@pytest.mark.parametrize(
    "mapper",
    [
        lambda values, idx: {i: e for e, i in zip(values, idx)},
        lambda values, idx: pd.Series(values, idx)])
def test_map_dictlike(idx, mapper):

    if isinstance(idx, (pd.CategoricalIndex, pd.IntervalIndex)):
        pytest.skip("skipping tests for {}".format(type(idx)))

    identity = mapper(idx.values, idx)

    # we don't infer to UInt64 for a dict
    if isinstance(idx, pd.UInt64Index) and isinstance(identity, dict):
        expected = idx.astype('int64')
    else:
        expected = idx

    result = idx.map(identity)
    tm.assert_index_equal(result, expected)

    # empty mappable
    expected = pd.Index([np.nan] * len(idx))
    result = idx.map(mapper(expected, idx))
    tm.assert_index_equal(result, expected)


def test_numpy_ufuncs(named_index):
    # test ufuncs of numpy 1.9.2. see:
    # http://docs.scipy.org/doc/numpy/reference/ufuncs.html

    # some functions are skipped because it may return different result
    # for unicode input depending on numpy version

    for name, idx in compat.iteritems(named_index):
        for func in [np.exp, np.exp2, np.expm1, np.log, np.log2, np.log10,
                     np.log1p, np.sqrt, np.sin, np.cos, np.tan, np.arcsin,
                     np.arccos, np.arctan, np.sinh, np.cosh, np.tanh,
                     np.arcsinh, np.arccosh, np.arctanh, np.deg2rad,
                     np.rad2deg]:
            if isinstance(idx, DatetimeIndexOpsMixin):
                # raise TypeError or ValueError (PeriodIndex)
                # PeriodIndex behavior should be changed in future version
                with pytest.raises(Exception):
                    with np.errstate(all='ignore'):
                        func(idx)
            elif isinstance(idx, (Float64Index, Int64Index, UInt64Index)):
                # coerces to float (e.g. np.sin)
                with np.errstate(all='ignore'):
                    result = func(idx)
                    exp = Index(func(idx.values), name=idx.name)

                tm.assert_index_equal(result, exp)
                assert isinstance(result, pd.Float64Index)
            else:
                # raise AttributeError or TypeError
                if len(idx) == 0:
                    continue
                else:
                    with pytest.raises(Exception):
                        with np.errstate(all='ignore'):
                            func(idx)

        for func in [np.isfinite, np.isinf, np.isnan, np.signbit]:
            if isinstance(idx, DatetimeIndexOpsMixin):
                # raise TypeError or ValueError (PeriodIndex)
                with pytest.raises(Exception):
                    func(idx)
            elif isinstance(idx, (Float64Index, Int64Index, UInt64Index)):
                # Results in bool array
                result = func(idx)
                assert isinstance(result, np.ndarray)
                assert not isinstance(result, Index)
            else:
                if len(idx) == 0:
                    continue
                else:
                    with pytest.raises(Exception):
                        func(idx)
