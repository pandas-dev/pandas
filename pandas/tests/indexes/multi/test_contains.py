# -*- coding: utf-8 -*-

import numpy as np

import pandas as pd

from pandas import (Index, MultiIndex)

from pandas.tests.indexes.common import Base


class TestContains(Base):
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

    def test_contains(self):
        assert ('foo', 'two') in self.index
        assert ('bar', 'two') not in self.index
        assert None not in self.index

    def test_contains_top_level(self):
        midx = MultiIndex.from_product([['A', 'B'], [1, 2]])
        assert 'A' in midx
        assert 'A' not in midx._engine

    def test_contains_with_nat(self):
        # MI with a NaT
        mi = MultiIndex(levels=[['C'],
                                pd.date_range('2012-01-01', periods=5)],
                        labels=[[0, 0, 0, 0, 0, 0], [-1, 0, 1, 2, 3, 4]],
                        names=[None, 'B'])
        assert ('C', pd.Timestamp('2012-01-01')) in mi
        for val in mi.values:
            assert val in mi
