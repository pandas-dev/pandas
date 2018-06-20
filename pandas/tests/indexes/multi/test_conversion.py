# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest
from pandas import DataFrame, MultiIndex, date_range
from pandas.compat import PY3, range
from pandas.util.testing import assert_almost_equal


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
    result = index.to_hierarchical(3)
    expected = MultiIndex(levels=[[1, 2], ['one', 'two']],
                          labels=[[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                                  [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]])
    tm.assert_index_equal(result, expected)
    assert result.names == index.names

    # K > 1
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

    result = index.to_hierarchical(2)
    expected = MultiIndex.from_tuples([(2, 'c'), (2, 'c'), (1, 'b'),
                                       (1, 'b'),
                                       (2, 'a'), (2, 'a'),
                                       (2, 'b'), (2, 'b')],
                                      names=['N1', 'N2'])
    tm.assert_index_equal(result, expected)
    assert result.names == index.names


def test_legacy_pickle():
    if PY3:
        pytest.skip("testing for legacy pickles not "
                    "support on py3")

    path = tm.get_data_path('multiindex_v1.pickle')
    obj = pd.read_pickle(path)

    obj2 = MultiIndex.from_tuples(obj.values)
    assert obj.equals(obj2)

    res = obj.get_indexer(obj)
    exp = np.arange(len(obj), dtype=np.intp)
    assert_almost_equal(res, exp)

    res = obj.get_indexer(obj2[::-1])
    exp = obj.get_indexer(obj[::-1])
    exp2 = obj2.get_indexer(obj2[::-1])
    assert_almost_equal(res, exp)
    assert_almost_equal(exp, exp2)


def test_legacy_v2_unpickle():

    # 0.7.3 -> 0.8.0 format manage
    path = tm.get_data_path('mindex_073.pickle')
    obj = pd.read_pickle(path)

    obj2 = MultiIndex.from_tuples(obj.values)
    assert obj.equals(obj2)

    res = obj.get_indexer(obj)
    exp = np.arange(len(obj), dtype=np.intp)
    assert_almost_equal(res, exp)

    res = obj.get_indexer(obj2[::-1])
    exp = obj.get_indexer(obj[::-1])
    exp2 = obj2.get_indexer(obj2[::-1])
    assert_almost_equal(res, exp)
    assert_almost_equal(exp, exp2)


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
