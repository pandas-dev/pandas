# -*- coding: utf-8 -*-

import numpy as np

import pandas as pd

from pandas import (Index, MultiIndex)

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base


class TestReIndex(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def setup_method(self, method):
        major_axis = Index(['foo', 'bar', 'baz', 'qux'])
        minor_axis = Index(['one', 'two'])

        major_labels = np.array([0, 0, 1, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])
        self.index_names = ['first', 'second']
        self.indices = dict(index=MultiIndex(levels=[major_axis, minor_axis],
                                             labels=[major_labels, minor_labels
                                                     ], names=self.index_names,
                                             verify_integrity=False))
        self.setup_indices()

    def create_index(self):
        return self.index

    def test_reindex(self):
        result, indexer = self.index.reindex(list(self.index[:4]))
        assert isinstance(result, MultiIndex)
        self.check_level_names(result, self.index[:4].names)

        result, indexer = self.index.reindex(list(self.index))
        assert isinstance(result, MultiIndex)
        assert indexer is None
        self.check_level_names(result, self.index.names)

    def test_reindex_level(self):
        idx = Index(['one'])

        target, indexer = self.index.reindex(idx, level='second')
        target2, indexer2 = idx.reindex(self.index, level='second')

        exp_index = self.index.join(idx, level='second', how='right')
        exp_index2 = self.index.join(idx, level='second', how='left')

        assert target.equals(exp_index)
        exp_indexer = np.array([0, 2, 4])
        tm.assert_numpy_array_equal(indexer, exp_indexer, check_dtype=False)

        assert target2.equals(exp_index2)
        exp_indexer2 = np.array([0, -1, 0, -1, 0, -1])
        tm.assert_numpy_array_equal(indexer2, exp_indexer2, check_dtype=False)

        tm.assert_raises_regex(TypeError, "Fill method not supported",
                               self.index.reindex, self.index,
                               method='pad', level='second')

        tm.assert_raises_regex(TypeError, "Fill method not supported",
                               idx.reindex, idx, method='bfill',
                               level='first')

    def test_reindex_preserves_names_when_target_is_list_or_ndarray(self):
        # GH6552
        idx = self.index.copy()
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

    def test_reindex_lvl_preserves_names_when_target_is_list_or_array(self):
        # GH7774
        idx = pd.MultiIndex.from_product([[0, 1], ['a', 'b']],
                                         names=['foo', 'bar'])
        assert idx.reindex([], level=0)[0].names == ['foo', 'bar']
        assert idx.reindex([], level=1)[0].names == ['foo', 'bar']

    def test_reindex_lvl_preserves_type_if_target_is_empty_list_or_array(self):
        # GH7774
        idx = pd.MultiIndex.from_product([[0, 1], ['a', 'b']])
        assert idx.reindex([], level=0)[0].levels[0].dtype.type == np.int64
        assert idx.reindex([], level=1)[0].levels[1].dtype.type == np.object_
