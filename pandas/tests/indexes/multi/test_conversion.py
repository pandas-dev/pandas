# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import pandas.util.testing as tm
from pandas import DataFrame, MultiIndex, date_range
from pandas.compat import range


def test_tolist(idx):
    result = idx.tolist()
    exp = list(idx.values)
    assert result == exp


def test_to_frame():
    tuples = [(1, 'one'), (1, 'two'), (2, 'one'), (2, 'two')]

    index = MultiIndex.from_tuples(tuples)
    result = index.to_frame(index=False)
    expected = DataFrame(tuples)
    tm.assert_frame_equal(result, expected)

    result = index.to_frame()
    expected.index = index
    tm.assert_frame_equal(result, expected)

    tuples = [(1, 'one'), (1, 'two'), (2, 'one'), (2, 'two')]
    index = MultiIndex.from_tuples(tuples, names=['first', 'second'])
    result = index.to_frame(index=False)
    expected = DataFrame(tuples)
    expected.columns = ['first', 'second']
    tm.assert_frame_equal(result, expected)

    result = index.to_frame()
    expected.index = index
    tm.assert_frame_equal(result, expected)

    index = MultiIndex.from_product([range(5),
                                     pd.date_range('20130101', periods=3)])
    result = index.to_frame(index=False)
    expected = DataFrame(
        {0: np.repeat(np.arange(5, dtype='int64'), 3),
            1: np.tile(pd.date_range('20130101', periods=3), 5)})
    tm.assert_frame_equal(result, expected)

    index = MultiIndex.from_product([range(5),
                                     pd.date_range('20130101', periods=3)])
    result = index.to_frame()
    expected.index = index
    tm.assert_frame_equal(result, expected)


def test_to_hierarchical():
    index = MultiIndex.from_tuples([(1, 'one'), (1, 'two'), (2, 'one'), (
        2, 'two')])
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        result = index.to_hierarchical(3)
    expected = MultiIndex(levels=[[1, 2], ['one', 'two']],
                          labels=[[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                                  [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]])
    tm.assert_index_equal(result, expected)
    assert result.names == index.names

    # K > 1
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        result = index.to_hierarchical(3, 2)
    expected = MultiIndex(levels=[[1, 2], ['one', 'two']],
                          labels=[[0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
                                  [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]])
    tm.assert_index_equal(result, expected)
    assert result.names == index.names

    # non-sorted
    index = MultiIndex.from_tuples([(2, 'c'), (1, 'b'),
                                    (2, 'a'), (2, 'b')],
                                   names=['N1', 'N2'])

    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        result = index.to_hierarchical(2)
    expected = MultiIndex.from_tuples([(2, 'c'), (2, 'c'), (1, 'b'),
                                       (1, 'b'),
                                       (2, 'a'), (2, 'a'),
                                       (2, 'b'), (2, 'b')],
                                      names=['N1', 'N2'])
    tm.assert_index_equal(result, expected)
    assert result.names == index.names


def test_roundtrip_pickle_with_tz():

    # GH 8367
    # round-trip of timezone
    index = MultiIndex.from_product(
        [[1, 2], ['a', 'b'], date_range('20130101', periods=3,
                                        tz='US/Eastern')
         ], names=['one', 'two', 'three'])
    unpickled = tm.round_trip_pickle(index)
    assert index.equal_levels(unpickled)


def test_pickle(indices):

    unpickled = tm.round_trip_pickle(indices)
    assert indices.equals(unpickled)
    original_name, indices.name = indices.name, 'foo'
    unpickled = tm.round_trip_pickle(indices)
    assert indices.equals(unpickled)
    indices.name = original_name


def test_to_series(idx):
    # assert that we are creating a copy of the index

    s = idx.to_series()
    assert s.values is not idx.values
    assert s.index is not idx
    assert s.name == idx.name


def test_to_series_with_arguments(idx):
    # GH18699

    # index kwarg
    s = idx.to_series(index=idx)

    assert s.values is not idx.values
    assert s.index is idx
    assert s.name == idx.name

    # name kwarg
    idx = idx
    s = idx.to_series(name='__test')

    assert s.values is not idx.values
    assert s.index is not idx
    assert s.name != idx.name
