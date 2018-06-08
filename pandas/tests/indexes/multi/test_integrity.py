# -*- coding: utf-8 -*-

import re


import pytest

import numpy as np

import pandas as pd

from pandas import DataFrame, MultiIndex, date_range
from pandas.compat import long, lrange, range
from pandas.errors import PerformanceWarning, UnsortedIndexError
from pandas.core.dtypes.cast import construct_1d_object_array_from_listlike

import pandas.util.testing as tm


def test_labels_dtypes():

    # GH 8456
    i = MultiIndex.from_tuples([('A', 1), ('A', 2)])
    assert i.labels[0].dtype == 'int8'
    assert i.labels[1].dtype == 'int8'

    i = MultiIndex.from_product([['a'], range(40)])
    assert i.labels[1].dtype == 'int8'
    i = MultiIndex.from_product([['a'], range(400)])
    assert i.labels[1].dtype == 'int16'
    i = MultiIndex.from_product([['a'], range(40000)])
    assert i.labels[1].dtype == 'int32'

    i = pd.MultiIndex.from_product([['a'], range(1000)])
    assert (i.labels[0] >= 0).all()
    assert (i.labels[1] >= 0).all()


def test_values_boxed():
    tuples = [(1, pd.Timestamp('2000-01-01')), (2, pd.NaT),
              (3, pd.Timestamp('2000-01-03')),
              (1, pd.Timestamp('2000-01-04')),
              (2, pd.Timestamp('2000-01-02')),
              (3, pd.Timestamp('2000-01-03'))]
    result = pd.MultiIndex.from_tuples(tuples)
    expected = construct_1d_object_array_from_listlike(tuples)
    tm.assert_numpy_array_equal(result.values, expected)
    # Check that code branches for boxed values produce identical results
    tm.assert_numpy_array_equal(result.values[:4], result[:4].values)


def test_values_multiindex_datetimeindex():
    # Test to ensure we hit the boxing / nobox part of MI.values
    ints = np.arange(10 ** 18, 10 ** 18 + 5)
    naive = pd.DatetimeIndex(ints)
    aware = pd.DatetimeIndex(ints, tz='US/Central')

    idx = pd.MultiIndex.from_arrays([naive, aware])
    result = idx.values

    outer = pd.DatetimeIndex([x[0] for x in result])
    tm.assert_index_equal(outer, naive)

    inner = pd.DatetimeIndex([x[1] for x in result])
    tm.assert_index_equal(inner, aware)

    # n_lev > n_lab
    result = idx[:2].values

    outer = pd.DatetimeIndex([x[0] for x in result])
    tm.assert_index_equal(outer, naive[:2])

    inner = pd.DatetimeIndex([x[1] for x in result])
    tm.assert_index_equal(inner, aware[:2])


def test_values_multiindex_periodindex():
    # Test to ensure we hit the boxing / nobox part of MI.values
    ints = np.arange(2007, 2012)
    pidx = pd.PeriodIndex(ints, freq='D')

    idx = pd.MultiIndex.from_arrays([ints, pidx])
    result = idx.values

    outer = pd.Int64Index([x[0] for x in result])
    tm.assert_index_equal(outer, pd.Int64Index(ints))

    inner = pd.PeriodIndex([x[1] for x in result])
    tm.assert_index_equal(inner, pidx)

    # n_lev > n_lab
    result = idx[:2].values

    outer = pd.Int64Index([x[0] for x in result])
    tm.assert_index_equal(outer, pd.Int64Index(ints[:2]))

    inner = pd.PeriodIndex([x[1] for x in result])
    tm.assert_index_equal(inner, pidx[:2])


def test_consistency():
    # need to construct an overflow
    major_axis = lrange(70000)
    minor_axis = lrange(10)

    major_labels = np.arange(70000)
    minor_labels = np.repeat(lrange(10), 7000)

    # the fact that is works means it's consistent
    index = MultiIndex(levels=[major_axis, minor_axis],
                       labels=[major_labels, minor_labels])

    # inconsistent
    major_labels = np.array([0, 0, 1, 1, 1, 2, 2, 3, 3])
    minor_labels = np.array([0, 1, 0, 1, 1, 0, 1, 0, 1])
    index = MultiIndex(levels=[major_axis, minor_axis],
                       labels=[major_labels, minor_labels])

    assert not index.is_unique


def test_hash_collisions():
    # non-smoke test that we don't get hash collisions

    index = MultiIndex.from_product([np.arange(1000), np.arange(1000)],
                                    names=['one', 'two'])
    result = index.get_indexer(index.values)
    tm.assert_numpy_array_equal(result, np.arange(
        len(index), dtype='intp'))

    for i in [0, 1, len(index) - 2, len(index) - 1]:
        result = index.get_loc(index[i])
        assert result == i


def test_equals_missing_values():
    # make sure take is not using -1
    i = pd.MultiIndex.from_tuples([(0, pd.NaT),
                                   (0, pd.Timestamp('20130101'))])
    result = i[0:1].equals(i[0])
    assert not result
    result = i[1:2].equals(i[1])
    assert not result


def test_dims():
    pass


def take_invalid_kwargs():
    vals = [['A', 'B'],
            [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02')]]
    idx = pd.MultiIndex.from_product(vals, names=['str', 'dt'])
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


def test_isna_behavior(_index):
    # should not segfault GH5123
    # NOTE: if MI representation changes, may make sense to allow
    # isna(MI)
    with pytest.raises(NotImplementedError):
        pd.isna(_index)


def test_large_multiindex_error():
    # GH12527
    df_below_1000000 = pd.DataFrame(
        1, index=pd.MultiIndex.from_product([[1, 2], range(499999)]),
        columns=['dest'])
    with pytest.raises(KeyError):
        df_below_1000000.loc[(-1, 0), 'dest']
    with pytest.raises(KeyError):
        df_below_1000000.loc[(3, 0), 'dest']
    df_above_1000000 = pd.DataFrame(
        1, index=pd.MultiIndex.from_product([[1, 2], range(500001)]),
        columns=['dest'])
    with pytest.raises(KeyError):
        df_above_1000000.loc[(-1, 0), 'dest']
    with pytest.raises(KeyError):
        df_above_1000000.loc[(3, 0), 'dest']


def test_nan_stays_float():

    # GH 7031
    idx0 = pd.MultiIndex(levels=[["A", "B"], []],
                         labels=[[1, 0], [-1, -1]],
                         names=[0, 1])
    idx1 = pd.MultiIndex(levels=[["C"], ["D"]],
                         labels=[[0], [0]],
                         names=[0, 1])
    idxm = idx0.join(idx1, how='outer')
    assert pd.isna(idx0.get_level_values(1)).all()
    # the following failed in 0.14.1
    assert pd.isna(idxm.get_level_values(1)[:-1]).all()

    df0 = pd.DataFrame([[1, 2]], index=idx0)
    df1 = pd.DataFrame([[3, 4]], index=idx1)
    dfm = df0 - df1
    assert pd.isna(df0.index.get_level_values(1)).all()
    # the following failed in 0.14.1
    assert pd.isna(dfm.index.get_level_values(1)[:-1]).all()


def test_million_record_attribute_error():
    # GH 18165
    r = list(range(1000000))
    df = pd.DataFrame({'a': r, 'b': r},
                      index=pd.MultiIndex.from_tuples([(x, x) for x in r]))

    with tm.assert_raises_regex(AttributeError,
                                "'Series' object has no attribute 'foo'"):
        df['a'].foo()


def test_can_hold_identifiers(_index):
    idx = _index
    key = idx[0]
    assert idx._can_hold_identifiers_and_holds_name(key) is True


def test_metadata_immutable(_index):
    levels, labels = _index.levels, _index.labels
    # shouldn't be able to set at either the top level or base level
    mutable_regex = re.compile('does not support mutable operations')
    with tm.assert_raises_regex(TypeError, mutable_regex):
        levels[0] = levels[0]
    with tm.assert_raises_regex(TypeError, mutable_regex):
        levels[0][0] = levels[0][0]
    # ditto for labels
    with tm.assert_raises_regex(TypeError, mutable_regex):
        labels[0] = labels[0]
    with tm.assert_raises_regex(TypeError, mutable_regex):
        labels[0][0] = labels[0][0]
    # and for names
    names = _index.names
    with tm.assert_raises_regex(TypeError, mutable_regex):
        names[0] = names[0]


def test_boolean_context_compat2():

    # boolean context compat
    # GH7897
    i1 = MultiIndex.from_tuples([('A', 1), ('A', 2)])
    i2 = MultiIndex.from_tuples([('A', 1), ('A', 3)])
    common = i1.intersection(i2)

    def f():
        if common:
            pass

    tm.assert_raises_regex(ValueError, 'The truth value of a', f)


def test_inplace_mutation_resets_values():
    levels = [['a', 'b', 'c'], [4]]
    levels2 = [[1, 2, 3], ['a']]
    labels = [[0, 1, 0, 2, 2, 0], [0, 0, 0, 0, 0, 0]]

    mi1 = MultiIndex(levels=levels, labels=labels)
    mi2 = MultiIndex(levels=levels2, labels=labels)
    vals = mi1.values.copy()
    vals2 = mi2.values.copy()

    assert mi1._tuples is not None

    # Make sure level setting works
    new_vals = mi1.set_levels(levels2).values
    tm.assert_almost_equal(vals2, new_vals)

    # Non-inplace doesn't kill _tuples [implementation detail]
    tm.assert_almost_equal(mi1._tuples, vals)

    # ...and values is still same too
    tm.assert_almost_equal(mi1.values, vals)

    # Inplace should kill _tuples
    mi1.set_levels(levels2, inplace=True)
    tm.assert_almost_equal(mi1.values, vals2)

    # Make sure label setting works too
    labels2 = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
    exp_values = np.empty((6,), dtype=object)
    exp_values[:] = [(long(1), 'a')] * 6

    # Must be 1d array of tuples
    assert exp_values.shape == (6,)
    new_values = mi2.set_labels(labels2).values

    # Not inplace shouldn't change
    tm.assert_almost_equal(mi2._tuples, vals2)

    # Should have correct values
    tm.assert_almost_equal(exp_values, new_values)

    # ...and again setting inplace should kill _tuples, etc
    mi2.set_labels(labels2, inplace=True)
    tm.assert_almost_equal(mi2.values, new_values)


def test_level_setting_resets_attributes():
    ind = pd.MultiIndex.from_arrays([
        ['A', 'A', 'B', 'B', 'B'], [1, 2, 1, 2, 3]
    ])
    assert ind.is_monotonic
    ind.set_levels([['A', 'B'], [1, 3, 2]], inplace=True)
    # if this fails, probably didn't reset the cache correctly.
    assert not ind.is_monotonic


def test_partial_string_timestamp_multiindex():
    # GH10331
    dr = pd.date_range('2016-01-01', '2016-01-03', freq='12H')
    abc = ['a', 'b', 'c']
    ix = pd.MultiIndex.from_product([dr, abc])
    df = pd.DataFrame({'c1': range(0, 15)}, index=ix)
    idx = pd.IndexSlice

    #                        c1
    # 2016-01-01 00:00:00 a   0
    #                     b   1
    #                     c   2
    # 2016-01-01 12:00:00 a   3
    #                     b   4
    #                     c   5
    # 2016-01-02 00:00:00 a   6
    #                     b   7
    #                     c   8
    # 2016-01-02 12:00:00 a   9
    #                     b  10
    #                     c  11
    # 2016-01-03 00:00:00 a  12
    #                     b  13
    #                     c  14

    # partial string matching on a single index
    for df_swap in (df.swaplevel(),
                    df.swaplevel(0),
                    df.swaplevel(0, 1)):
        df_swap = df_swap.sort_index()
        just_a = df_swap.loc['a']
        result = just_a.loc['2016-01-01']
        expected = df.loc[idx[:, 'a'], :].iloc[0:2]
        expected.index = expected.index.droplevel(1)
        tm.assert_frame_equal(result, expected)

    # indexing with IndexSlice
    result = df.loc[idx['2016-01-01':'2016-02-01', :], :]
    expected = df
    tm.assert_frame_equal(result, expected)

    # match on secondary index
    result = df_swap.loc[idx[:, '2016-01-01':'2016-01-01'], :]
    expected = df_swap.iloc[[0, 1, 5, 6, 10, 11]]
    tm.assert_frame_equal(result, expected)

    # Even though this syntax works on a single index, this is somewhat
    # ambiguous and we don't want to extend this behavior forward to work
    # in multi-indexes. This would amount to selecting a scalar from a
    # column.
    with pytest.raises(KeyError):
        df['2016-01-01']

    # partial string match on year only
    result = df.loc['2016']
    expected = df
    tm.assert_frame_equal(result, expected)

    # partial string match on date
    result = df.loc['2016-01-01']
    expected = df.iloc[0:6]
    tm.assert_frame_equal(result, expected)

    # partial string match on date and hour, from middle
    result = df.loc['2016-01-02 12']
    expected = df.iloc[9:12]
    tm.assert_frame_equal(result, expected)

    # partial string match on secondary index
    result = df_swap.loc[idx[:, '2016-01-02'], :]
    expected = df_swap.iloc[[2, 3, 7, 8, 12, 13]]
    tm.assert_frame_equal(result, expected)

    # tuple selector with partial string match on date
    result = df.loc[('2016-01-01', 'a'), :]
    expected = df.iloc[[0, 3]]
    tm.assert_frame_equal(result, expected)

    # Slicing date on first level should break (of course)
    with pytest.raises(KeyError):
        df_swap.loc['2016-01-01']

    # GH12685 (partial string with daily resolution or below)
    dr = date_range('2013-01-01', periods=100, freq='D')
    ix = MultiIndex.from_product([dr, ['a', 'b']])
    df = DataFrame(np.random.randn(200, 1), columns=['A'], index=ix)

    result = df.loc[idx['2013-03':'2013-03', :], :]
    expected = df.iloc[118:180]
    tm.assert_frame_equal(result, expected)


def test_rangeindex_fallback_coercion_bug():
    # GH 12893
    foo = pd.DataFrame(np.arange(100).reshape((10, 10)))
    bar = pd.DataFrame(np.arange(100).reshape((10, 10)))
    df = pd.concat({'foo': foo.stack(), 'bar': bar.stack()}, axis=1)
    df.index.names = ['fizz', 'buzz']

    str(df)
    expected = pd.DataFrame({'bar': np.arange(100),
                             'foo': np.arange(100)},
                            index=pd.MultiIndex.from_product(
                                [range(10), range(10)],
                                names=['fizz', 'buzz']))
    tm.assert_frame_equal(df, expected, check_like=True)

    result = df.index.get_level_values('fizz')
    expected = pd.Int64Index(np.arange(10), name='fizz').repeat(10)
    tm.assert_index_equal(result, expected)

    result = df.index.get_level_values('buzz')
    expected = pd.Int64Index(np.tile(np.arange(10), 10), name='buzz')
    tm.assert_index_equal(result, expected)


def test_unsortedindex():
    # GH 11897
    mi = pd.MultiIndex.from_tuples([('z', 'a'), ('x', 'a'), ('y', 'b'),
                                    ('x', 'b'), ('y', 'a'), ('z', 'b')],
                                   names=['one', 'two'])
    df = pd.DataFrame([[i, 10 * i] for i in lrange(6)], index=mi,
                      columns=['one', 'two'])

    # GH 16734: not sorted, but no real slicing
    result = df.loc(axis=0)['z', 'a']
    expected = df.iloc[0]
    tm.assert_series_equal(result, expected)

    with pytest.raises(UnsortedIndexError):
        df.loc(axis=0)['z', slice('a')]
    df.sort_index(inplace=True)
    assert len(df.loc(axis=0)['z', :]) == 2

    with pytest.raises(KeyError):
        df.loc(axis=0)['q', :]


def test_unsortedindex_doc_examples():
    # http://pandas.pydata.org/pandas-docs/stable/advanced.html#sorting-a-multiindex  # noqa
    dfm = DataFrame({'jim': [0, 0, 1, 1],
                     'joe': ['x', 'x', 'z', 'y'],
                     'jolie': np.random.rand(4)})

    dfm = dfm.set_index(['jim', 'joe'])
    with tm.assert_produces_warning(PerformanceWarning):
        dfm.loc[(1, 'z')]

    with pytest.raises(UnsortedIndexError):
        dfm.loc[(0, 'y'):(1, 'z')]

    assert not dfm.index.is_lexsorted()
    assert dfm.index.lexsort_depth == 1

    # sort it
    dfm = dfm.sort_index()
    dfm.loc[(1, 'z')]
    dfm.loc[(0, 'y'):(1, 'z')]

    assert dfm.index.is_lexsorted()
    assert dfm.index.lexsort_depth == 2
