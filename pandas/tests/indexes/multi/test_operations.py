# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest
from pandas import (DatetimeIndex, Float64Index, Index, Int64Index, MultiIndex,
                    PeriodIndex, TimedeltaIndex, UInt64Index, date_range,
                    period_range)
from pandas.compat import lrange, range
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


def test_take(idx):
    indexer = [4, 3, 0, 2]
    result = idx.take(indexer)
    expected = idx[indexer]
    assert result.equals(expected)

    if not isinstance(idx,
                      (DatetimeIndex, PeriodIndex, TimedeltaIndex)):
        # GH 10791
        with pytest.raises(AttributeError):
            idx.freq


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


def test_argsort(idx):
    result = idx.argsort()
    expected = idx.values.argsort()
    tm.assert_numpy_array_equal(result, expected)


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


def test_numpy_ufuncs(idx):
    # test ufuncs of numpy 1.9.2. see:
    # http://docs.scipy.org/doc/numpy/reference/ufuncs.html

    # some functions are skipped because it may return different result
    # for unicode input depending on numpy version

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
