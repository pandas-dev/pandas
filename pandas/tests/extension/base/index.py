import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas.core.indexes.extension import ExtensionIndex

from .base import BaseExtensionTests


class BaseIndexTests(BaseExtensionTests):
    """Tests for ExtensionIndex."""

    def test_constructor(self, data):
        result = ExtensionIndex(data, name='test')
        assert result.name == 'test'
        self.assert_extension_array_equal(data, result._values)

    def test_series_constructor(self, data):
        result = pd.Series(range(len(data)), index=data)
        assert isinstance(result.index, ExtensionIndex)

    def test_asarray(self, data):
        idx = ExtensionIndex(data)
        tm.assert_numpy_array_equal(np.array(idx), np.array(data))

    def test_repr(self, data):
        idx = ExtensionIndex(data, name='test')
        repr(idx)
        s = pd.Series(range(len(data)), index=data)
        repr(s)

    def test_indexing_scalar(self, data):
        s = pd.Series(range(len(data)), index=data)
        label = data[1]
        assert s[label] == 1
        assert s.iloc[1] == 1
        assert s.loc[label] == 1

    def test_indexing_list(self, data):
        s = pd.Series(range(len(data)), index=data)
        labels = [data[1], data[3]]
        exp = pd.Series([1, 3], index=data[[1, 3]])
        self.assert_series_equal(s[labels], exp)
        self.assert_series_equal(s.loc[labels], exp)
        self.assert_series_equal(s.iloc[[1, 3]], exp)

    def test_contains(self, data_missing, data_for_sorting, na_value):
        idx = ExtensionIndex(data_missing)
        assert data_missing[0] in idx
        assert data_missing[1] in idx
        assert na_value in idx
        assert '__random' not in idx
        idx = ExtensionIndex(data_for_sorting)
        assert na_value not in idx

    def test_na(self, data_missing):
        idx = ExtensionIndex(data_missing)
        result = idx.isna()
        expected = np.array([True, False], dtype=bool)
        tm.assert_numpy_array_equal(result, expected)
        result = idx.notna()
        tm.assert_numpy_array_equal(result, ~expected)
        assert idx.hasnans #is True

    def test_monotonic(self, data_for_sorting):
        data = data_for_sorting
        idx = ExtensionIndex(data)
        assert idx.is_monotonic_increasing is False
        assert idx.is_monotonic_decreasing is False

        idx = ExtensionIndex(data[[2, 0, 1]])
        assert idx.is_monotonic_increasing is True
        assert idx.is_monotonic_decreasing is False

        idx = ExtensionIndex(data[[1, 0, 2]])
        assert idx.is_monotonic_increasing is False
        assert idx.is_monotonic_decreasing is True

    def test_is_unique(self, data_for_sorting, data_for_grouping):
        idx = ExtensionIndex(data_for_sorting)
        assert idx.is_unique is True

        idx = ExtensionIndex(data_for_grouping)
        assert idx.is_unique is False

    def test_take(self, data):
        idx = ExtensionIndex(data)
        expected = ExtensionIndex(data.take([0, 2, 3]))
        result = idx.take([0, 2, 3])
        tm.assert_index_equal(result, expected)

    def test_getitem(self, data):
        idx = ExtensionIndex(data)
        assert idx[0] == data[0]
        tm.assert_index_equal(idx[[0, 1]], ExtensionIndex(data[[0, 1]]))
