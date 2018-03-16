import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm

from .base import BaseExtensionTests


class BaseMethodsTests(BaseExtensionTests):
    """Various Series and DataFrame methods."""

    @pytest.mark.parametrize('dropna', [True, False])
    def test_value_counts(self, all_data, dropna):
        all_data = all_data[:10]
        if dropna:
            other = np.array(all_data[~all_data.isna()])
        else:
            other = all_data

        result = pd.Series(all_data).value_counts(dropna=dropna).sort_index()
        expected = pd.Series(other).value_counts(dropna=dropna).sort_index()

        self.assert_series_equal(result, expected)

    def test_count(self, data_missing):
        df = pd.DataFrame({"A": data_missing})
        result = df.count(axis='columns')
        expected = pd.Series([0, 1])
        self.assert_series_equal(result, expected)

    def test_apply_simple_series(self, data):
        result = pd.Series(data).apply(id)
        assert isinstance(result, pd.Series)

    @pytest.mark.parametrize('box', [pd.Series, lambda x: x])
    @pytest.mark.parametrize('method', [lambda x: x.unique(), pd.unique])
    def test_unique(self, data, box, method):
        duplicated = box(data._constructor_from_sequence([data[0], data[0]]))

        result = method(duplicated)

        assert len(result) == 1
        assert isinstance(result, type(data))
        assert result[0] == duplicated[0]

    @pytest.mark.parametrize('na_sentinel', [-1, -2])
    def test_factorize(self, data_for_grouping, na_sentinel):
        labels, uniques = pd.factorize(data_for_grouping,
                                       na_sentinel=na_sentinel)
        expected_labels = np.array([0, 0, na_sentinel,
                                   na_sentinel, 1, 1, 0, 2],
                                   dtype='int64')
        expected_uniques = data_for_grouping.take([0, 4, 7])

        tm.assert_numpy_array_equal(labels, expected_labels)
        self.assert_extension_array_equal(uniques, expected_uniques)

    @pytest.mark.parametrize('na_sentinel', [-1, -2])
    def test_factorize_equivalence(self, data_for_grouping, na_sentinel):
        l1, u1 = pd.factorize(data_for_grouping, na_sentinel=na_sentinel)
        l2, u2 = data_for_grouping.factorize(na_sentinel=na_sentinel)

        tm.assert_numpy_array_equal(l1, l2)
        self.assert_extension_array_equal(u1, u2)
