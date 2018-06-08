# -*- coding: utf-8 -*-

import warnings
from itertools import product

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest
from pandas import DataFrame, Index, MultiIndex, date_range, period_range
from pandas.compat import PYPY, lrange, lzip, range, u
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.util.testing import assert_copy


def check_level_names(index, names):
    assert [level.name for level in index.levels] == list(names)


def test_difference(_index):

    first = _index
    result = first.difference(_index[-3:])
    expected = MultiIndex.from_tuples(sorted(_index[:-3].values),
                                      sortorder=0,
                                      names=_index.names)

    assert isinstance(result, MultiIndex)
    assert result.equals(expected)
    assert result.names == _index.names

    # empty difference: reflexive
    result = _index.difference(_index)
    expected = _index[:0]
    assert result.equals(expected)
    assert result.names == _index.names

    # empty difference: superset
    result = _index[-3:].difference(_index)
    expected = _index[:0]
    assert result.equals(expected)
    assert result.names == _index.names

    # empty difference: degenerate
    result = _index[:0].difference(_index)
    expected = _index[:0]
    assert result.equals(expected)
    assert result.names == _index.names

    # names not the same
    chunklet = _index[-3:]
    chunklet.names = ['foo', 'baz']
    result = first.difference(chunklet)
    assert result.names == (None, None)

    # empty, but non-equal
    result = _index.difference(_index.sortlevel(1)[0])
    assert len(result) == 0

    # raise Exception called with non-MultiIndex
    result = first.difference(first.values)
    assert result.equals(first[:0])

    # name from empty array
    result = first.difference([])
    assert first.equals(result)
    assert first.names == result.names

    # name from non-empty array
    result = first.difference([('foo', 'one')])
    expected = pd.MultiIndex.from_tuples([('bar', 'one'), ('baz', 'two'), (
        'foo', 'two'), ('qux', 'one'), ('qux', 'two')])
    expected.names = first.names
    assert first.names == result.names
    tm.assert_raises_regex(TypeError, "other must be a MultiIndex "
                           "or a list of tuples",
                           first.difference, [1, 2, 3, 4, 5])


def test_union(_index):
    piece1 = _index[:5][::-1]
    piece2 = _index[3:]

    the_union = piece1 | piece2

    tups = sorted(_index.values)
    expected = MultiIndex.from_tuples(tups)

    assert the_union.equals(expected)

    # corner case, pass self or empty thing:
    the_union = _index.union(_index)
    assert the_union is _index

    the_union = _index.union(_index[:0])
    assert the_union is _index

    # won't work in python 3
    # tuples = _index.values
    # result = _index[:4] | tuples[4:]
    # assert result.equals(tuples)

    # not valid for python 3
    # def test_union_with_regular_index(self):
    #     other = Index(['A', 'B', 'C'])

    #     result = other.union(_index)
    #     assert ('foo', 'one') in result
    #     assert 'B' in result

    #     result2 = _index.union(other)
    #     assert result.equals(result2)


def test_intersection(_index):
    piece1 = _index[:5][::-1]
    piece2 = _index[3:]

    the_int = piece1 & piece2
    tups = sorted(_index[3:5].values)
    expected = MultiIndex.from_tuples(tups)
    assert the_int.equals(expected)

    # corner case, pass self
    the_int = _index.intersection(_index)
    assert the_int is _index

    # empty intersection: disjoint
    empty = _index[:2] & _index[2:]
    expected = _index[:0]
    assert empty.equals(expected)

    # can't do in python 3
    # tuples = _index.values
    # result = _index & tuples
    # assert result.equals(tuples)


def test_insert(_index):
    # key contained in all levels
    new_index = _index.insert(0, ('bar', 'two'))
    assert new_index.equal_levels(_index)
    assert new_index[0] == ('bar', 'two')

    # key not contained in all levels
    new_index = _index.insert(0, ('abc', 'three'))

    exp0 = Index(list(_index.levels[0]) + ['abc'], name='first')
    tm.assert_index_equal(new_index.levels[0], exp0)

    exp1 = Index(list(_index.levels[1]) + ['three'], name='second')
    tm.assert_index_equal(new_index.levels[1], exp1)
    assert new_index[0] == ('abc', 'three')

    # key wrong length
    msg = "Item must have length equal to number of levels"
    with tm.assert_raises_regex(ValueError, msg):
        _index.insert(0, ('foo2',))

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


def test_is_all_dates(_index):
    assert not _index.is_all_dates


def test_is_numeric(_index):
    # MultiIndex is never numeric
    assert not _index.is_numeric()


def test_bounds(_index):
    _index._bounds


def test_equals_multi(_index):
    assert _index.equals(_index)
    assert not _index.equals(_index.values)
    assert _index.equals(Index(_index.values))

    assert _index.equal_levels(_index)
    assert not _index.equals(_index[:-1])
    assert not _index.equals(_index[-1])

    # different number of levels
    index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
        lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
            [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])

    index2 = MultiIndex(levels=index.levels[:-1], labels=index.labels[:-1])
    assert not index.equals(index2)
    assert not index.equal_levels(index2)

    # levels are different
    major_axis = Index(lrange(4))
    minor_axis = Index(lrange(2))

    major_labels = np.array([0, 0, 1, 2, 2, 3])
    minor_labels = np.array([0, 1, 0, 0, 1, 0])

    index = MultiIndex(levels=[major_axis, minor_axis],
                       labels=[major_labels, minor_labels])
    assert not _index.equals(index)
    assert not _index.equal_levels(index)

    # some of the labels are different
    major_axis = Index(['foo', 'bar', 'baz', 'qux'])
    minor_axis = Index(['one', 'two'])

    major_labels = np.array([0, 0, 2, 2, 3, 3])
    minor_labels = np.array([0, 1, 0, 1, 0, 1])

    index = MultiIndex(levels=[major_axis, minor_axis],
                       labels=[major_labels, minor_labels])
    assert not _index.equals(index)


def test_identical(_index):
    mi = _index.copy()
    mi2 = _index.copy()
    assert mi.identical(mi2)

    mi = mi.set_names(['new1', 'new2'])
    assert mi.equals(mi2)
    assert not mi.identical(mi2)

    mi2 = mi2.set_names(['new1', 'new2'])
    assert mi.identical(mi2)

    mi3 = Index(mi.tolist(), names=mi.names)
    mi4 = Index(mi.tolist(), names=mi.names, tupleize_cols=False)
    assert mi.identical(mi3)
    assert not mi.identical(mi4)
    assert mi.equals(mi4)


def test_append(_index):
    result = _index[:3].append(_index[3:])
    assert result.equals(_index)

    foos = [_index[:1], _index[1:3], _index[3:]]
    result = foos[0].append(foos[1:])
    assert result.equals(_index)

    # empty
    result = _index.append([])
    assert result.equals(_index)


def test_groupby(_index):
    groups = _index.groupby(np.array([1, 1, 1, 2, 2, 2]))
    labels = _index.get_values().tolist()
    exp = {1: labels[:3], 2: labels[3:]}
    tm.assert_dict_equal(groups, exp)

    # GH5620
    groups = _index.groupby(_index)
    exp = {key: [key] for key in _index}
    tm.assert_dict_equal(groups, exp)


def test_equals_operator(_index):
    # GH9785
    assert (_index == _index).all()


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


def test_reorder_levels(_index):
    # this blows up
    tm.assert_raises_regex(IndexError, '^Too many levels',
                           _index.reorder_levels, [2, 1, 0])


def test_astype(_index):
    expected = _index.copy()
    actual = _index.astype('O')
    assert_copy(actual.levels, expected.levels)
    assert_copy(actual.labels, expected.labels)
    check_level_names(actual, expected.names)

    with tm.assert_raises_regex(TypeError, "^Setting.*dtype.*object"):
        _index.astype(np.dtype(int))


@pytest.mark.parametrize('ordered', [True, False])
def test_astype_category(_index, ordered):
    # GH 18630
    msg = '> 1 ndim Categorical are not supported at this time'
    with tm.assert_raises_regex(NotImplementedError, msg):
        _index.astype(CategoricalDtype(ordered=ordered))

    if ordered is False:
        # dtype='category' defaults to ordered=False, so only test once
        with tm.assert_raises_regex(NotImplementedError, msg):
            _index.astype('category')


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


def test_is_():
    mi = MultiIndex.from_tuples(lzip(range(10), range(10)))
    assert mi.is_(mi)
    assert mi.is_(mi.view())
    assert mi.is_(mi.view().view().view().view())
    mi2 = mi.view()
    # names are metadata, they don't change id
    mi2.names = ["A", "B"]
    assert mi2.is_(mi)
    assert mi.is_(mi2)

    assert mi.is_(mi.set_names(["C", "D"]))
    mi2 = mi.view()
    mi2.set_names(["E", "F"], inplace=True)
    assert mi.is_(mi2)
    # levels are inherent properties, they change identity
    mi3 = mi2.set_levels([lrange(10), lrange(10)])
    assert not mi3.is_(mi2)
    # shouldn't change
    assert mi2.is_(mi)
    mi4 = mi3.view()

    # GH 17464 - Remove duplicate MultiIndex levels
    mi4.set_levels([lrange(10), lrange(10)], inplace=True)
    assert not mi4.is_(mi3)
    mi5 = mi.view()
    mi5.set_levels(mi5.levels, inplace=True)
    assert not mi5.is_(mi)


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


def test_iter(_index):
    result = list(_index)
    expected = [('foo', 'one'), ('foo', 'two'), ('bar', 'one'),
                ('baz', 'two'), ('qux', 'one'), ('qux', 'two')]
    assert result == expected


def test_sub(_index):

    first = _index

    # - now raises (previously was set op difference)
    with pytest.raises(TypeError):
        first - _index[-3:]
    with pytest.raises(TypeError):
        _index[-3:] - first
    with pytest.raises(TypeError):
        _index[-3:] - first.tolist()
    with pytest.raises(TypeError):
        first.tolist() - _index[-3:]


def test_nlevels(_index):
    assert _index.nlevels == 2


def test_argsort(_index):
    result = _index.argsort()
    expected = _index.values.argsort()
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
def test_unique_level(_index, level):
    # GH #17896 - with level= argument
    result = _index.unique(level=level)
    expected = _index.get_level_values(level).unique()
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


@pytest.mark.skipif(not PYPY, reason="tuples cmp recursively on PyPy")
def test_isin_nan_pypy():
    idx = MultiIndex.from_arrays([['foo', 'bar'], [1.0, np.nan]])
    tm.assert_numpy_array_equal(idx.isin([('bar', np.nan)]),
                                np.array([False, True]))
    tm.assert_numpy_array_equal(idx.isin([('bar', float('nan'))]),
                                np.array([False, True]))


def test_isin():
    values = [('foo', 2), ('bar', 3), ('quux', 4)]

    idx = MultiIndex.from_arrays([['qux', 'baz', 'foo', 'bar'], np.arange(
        4)])
    result = idx.isin(values)
    expected = np.array([False, False, True, True])
    tm.assert_numpy_array_equal(result, expected)

    # empty, return dtype bool
    idx = MultiIndex.from_arrays([[], []])
    result = idx.isin(values)
    assert len(result) == 0
    assert result.dtype == np.bool_


@pytest.mark.skipif(PYPY, reason="tuples cmp recursively on PyPy")
def test_isin_nan_not_pypy():
    idx = MultiIndex.from_arrays([['foo', 'bar'], [1.0, np.nan]])
    tm.assert_numpy_array_equal(idx.isin([('bar', np.nan)]),
                                np.array([False, False]))
    tm.assert_numpy_array_equal(idx.isin([('bar', float('nan'))]),
                                np.array([False, False]))


def test_isin_level_kwarg():
    idx = MultiIndex.from_arrays([['qux', 'baz', 'foo', 'bar'], np.arange(
        4)])

    vals_0 = ['foo', 'bar', 'quux']
    vals_1 = [2, 3, 10]

    expected = np.array([False, False, True, True])
    tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level=0))
    tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level=-2))

    tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level=1))
    tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level=-1))

    pytest.raises(IndexError, idx.isin, vals_0, level=5)
    pytest.raises(IndexError, idx.isin, vals_0, level=-5)

    pytest.raises(KeyError, idx.isin, vals_0, level=1.0)
    pytest.raises(KeyError, idx.isin, vals_1, level=-1.0)
    pytest.raises(KeyError, idx.isin, vals_1, level='A')

    idx.names = ['A', 'B']
    tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level='A'))
    tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level='B'))

    pytest.raises(KeyError, idx.isin, vals_1, level='C')


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


def test_duplicates(_index):
    assert not _index.has_duplicates
    assert _index.append(_index).has_duplicates

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
