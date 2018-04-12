# -*- coding: utf-8 -*-

import warnings


import pandas as pd

from pandas import MultiIndex

import pandas.util.testing as tm

from .common import Base


class TestFormat(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

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
