# -*- coding: utf-8 -*-

import warnings

import numpy as np

import pandas as pd

from pandas import (Index, MultiIndex)

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base


class TestFormat(Base):
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

    def test_format(self):
        self.index.format()
        self.index[:0].format()

    def test_format_integer_names(self):
        index = MultiIndex(levels=[[0, 1], [0, 1]],
                           labels=[[0, 0, 1, 1], [0, 1, 0, 1]], names=[0, 1])
        index.format(names=True)

    def test_format_sparse_display(self):
        index = MultiIndex(levels=[[0, 1], [0, 1], [0, 1], [0]],
                           labels=[[0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 0, 1],
                                   [0, 1, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0]])

        result = index.format()
        assert result[3] == '1  0  0  0'

    def test_format_sparse_config(self):
        warn_filters = warnings.filters
        warnings.filterwarnings('ignore', category=FutureWarning,
                                module=".*format")
        # GH1538
        pd.set_option('display.multi_sparse', False)

        result = self.index.format()
        assert result[1] == 'foo  two'

        tm.reset_display_options()

        warnings.filters = warn_filters
