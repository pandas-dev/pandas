# -*- coding: utf-8 -*-

import re
import warnings

from datetime import timedelta
from itertools import product

import pytest

import numpy as np

import pandas as pd

from pandas import (CategoricalIndex, DataFrame, Index, MultiIndex,
                    compat, date_range, period_range)
from pandas.compat import PY3, long, lrange, lzip, range, u, PYPY
from pandas.errors import PerformanceWarning, UnsortedIndexError
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.core.indexes.base import InvalidIndexError
from pandas.core.dtypes.cast import construct_1d_object_array_from_listlike
from pandas._libs.tslib import Timestamp

import pandas.util.testing as tm

from pandas.util.testing import assert_almost_equal, assert_copy

from .common import Base


class TestContains(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

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
