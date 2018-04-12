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


class TestSortLevel(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']
    

    def test_sortlevel(self):
        import random

        tuples = list(self.index)
        random.shuffle(tuples)

        index = MultiIndex.from_tuples(tuples)

        sorted_idx, _ = index.sortlevel(0)
        expected = MultiIndex.from_tuples(sorted(tuples))
        assert sorted_idx.equals(expected)

        sorted_idx, _ = index.sortlevel(0, ascending=False)
        assert sorted_idx.equals(expected[::-1])

        sorted_idx, _ = index.sortlevel(1)
        by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
        expected = MultiIndex.from_tuples(by1)
        assert sorted_idx.equals(expected)

        sorted_idx, _ = index.sortlevel(1, ascending=False)
        assert sorted_idx.equals(expected[::-1])

    def test_sortlevel_not_sort_remaining(self):
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        sorted_idx, _ = mi.sortlevel('A', sort_remaining=False)
        assert sorted_idx.equals(mi)

    def test_sortlevel_deterministic(self):
        tuples = [('bar', 'one'), ('foo', 'two'), ('qux', 'two'),
                  ('foo', 'one'), ('baz', 'two'), ('qux', 'one')]

        index = MultiIndex.from_tuples(tuples)

        sorted_idx, _ = index.sortlevel(0)
        expected = MultiIndex.from_tuples(sorted(tuples))
        assert sorted_idx.equals(expected)

        sorted_idx, _ = index.sortlevel(0, ascending=False)
        assert sorted_idx.equals(expected[::-1])

        sorted_idx, _ = index.sortlevel(1)
        by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
        expected = MultiIndex.from_tuples(by1)
        assert sorted_idx.equals(expected)

        sorted_idx, _ = index.sortlevel(1, ascending=False)
        assert sorted_idx.equals(expected[::-1])
