# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import pandas.util.testing as tm
from pandas import Index, MultiIndex


def check_level_names(index, names):
    assert [level.name for level in index.levels] == list(names)


def test_reindex(_index):
    result, indexer = _index.reindex(list(_index[:4]))
    assert isinstance(result, MultiIndex)
    check_level_names(result, _index[:4].names)

    result, indexer = _index.reindex(list(_index))
    assert isinstance(result, MultiIndex)
    assert indexer is None
    check_level_names(result, _index.names)


def test_reindex_level(_index):
    idx = Index(['one'])

    target, indexer = _index.reindex(idx, level='second')
    target2, indexer2 = idx.reindex(_index, level='second')

    exp_index = _index.join(idx, level='second', how='right')
    exp_index2 = _index.join(idx, level='second', how='left')

    assert target.equals(exp_index)
    exp_indexer = np.array([0, 2, 4])
    tm.assert_numpy_array_equal(indexer, exp_indexer, check_dtype=False)

    assert target2.equals(exp_index2)
    exp_indexer2 = np.array([0, -1, 0, -1, 0, -1])
    tm.assert_numpy_array_equal(indexer2, exp_indexer2, check_dtype=False)

    tm.assert_raises_regex(TypeError, "Fill method not supported",
                           _index.reindex, _index,
                           method='pad', level='second')

    tm.assert_raises_regex(TypeError, "Fill method not supported",
                           idx.reindex, idx, method='bfill',
                           level='first')


def test_reindex_preserves_names_when_target_is_list_or_ndarray(_index):
    # GH6552
    idx = _index.copy()
    target = idx.copy()
    idx.names = target.names = [None, None]

    other_dtype = pd.MultiIndex.from_product([[1, 2], [3, 4]])

    # list & ndarray cases
    assert idx.reindex([])[0].names == [None, None]
    assert idx.reindex(np.array([]))[0].names == [None, None]
    assert idx.reindex(target.tolist())[0].names == [None, None]
    assert idx.reindex(target.values)[0].names == [None, None]
    assert idx.reindex(other_dtype.tolist())[0].names == [None, None]
    assert idx.reindex(other_dtype.values)[0].names == [None, None]

    idx.names = ['foo', 'bar']
    assert idx.reindex([])[0].names == ['foo', 'bar']
    assert idx.reindex(np.array([]))[0].names == ['foo', 'bar']
    assert idx.reindex(target.tolist())[0].names == ['foo', 'bar']
    assert idx.reindex(target.values)[0].names == ['foo', 'bar']
    assert idx.reindex(other_dtype.tolist())[0].names == ['foo', 'bar']
    assert idx.reindex(other_dtype.values)[0].names == ['foo', 'bar']


def test_reindex_lvl_preserves_names_when_target_is_list_or_array():
    # GH7774
    idx = pd.MultiIndex.from_product([[0, 1], ['a', 'b']],
                                     names=['foo', 'bar'])
    assert idx.reindex([], level=0)[0].names == ['foo', 'bar']
    assert idx.reindex([], level=1)[0].names == ['foo', 'bar']


def test_reindex_lvl_preserves_type_if_target_is_empty_list_or_array():
    # GH7774
    idx = pd.MultiIndex.from_product([[0, 1], ['a', 'b']])
    assert idx.reindex([], level=0)[0].levels[0].dtype.type == np.int64
    assert idx.reindex([], level=1)[0].levels[1].dtype.type == np.object_
