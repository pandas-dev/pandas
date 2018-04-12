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


class TestTuples(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']
    
    def test_from_tuples(self):
        tm.assert_raises_regex(TypeError, 'Cannot infer number of levels '
                                          'from empty list',
                               MultiIndex.from_tuples, [])

        expected = MultiIndex(levels=[[1, 3], [2, 4]],
                              labels=[[0, 1], [0, 1]],
                              names=['a', 'b'])

        # input tuples
        result = MultiIndex.from_tuples(((1, 2), (3, 4)), names=['a', 'b'])
        tm.assert_index_equal(result, expected)

    def test_from_tuples_iterator(self):
        # GH 18434
        # input iterator for tuples
        expected = MultiIndex(levels=[[1, 3], [2, 4]],
                              labels=[[0, 1], [0, 1]],
                              names=['a', 'b'])

        result = MultiIndex.from_tuples(zip([1, 3], [2, 4]), names=['a', 'b'])
        tm.assert_index_equal(result, expected)

        # input non-iterables
        with tm.assert_raises_regex(
                TypeError, 'Input must be a list / sequence of tuple-likes.'):
            MultiIndex.from_tuples(0)

    def test_from_tuples_empty(self):
        # GH 16777
        result = MultiIndex.from_tuples([], names=['a', 'b'])
        expected = MultiIndex.from_arrays(arrays=[[], []],
                                          names=['a', 'b'])
        tm.assert_index_equal(result, expected)
