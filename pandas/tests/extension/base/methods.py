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
        expected = pd.Series(other).value_counts(
            dropna=dropna).sort_index()

        self.assert_series_equal(result, expected)

    def test_count(self, data_missing):
        df = pd.DataFrame({"A": data_missing})
        result = df.count(axis='columns')
        expected = pd.Series([0, 1])
        self.assert_series_equal(result, expected)

    def test_apply_simple_series(self, data):
        result = pd.Series(data).apply(id)
        assert isinstance(result, pd.Series)

    def test_argsort(self, data_for_sorting):
        result = pd.Series(data_for_sorting).argsort()
        expected = pd.Series(np.array([2, 0, 1], dtype=np.int64))
        self.assert_series_equal(result, expected)

    def test_argsort_missing(self, data_missing_for_sorting):
        result = pd.Series(data_missing_for_sorting).argsort()
        expected = pd.Series(np.array([1, -1, 0], dtype=np.int64))
        self.assert_series_equal(result, expected)

    @pytest.mark.parametrize('ascending', [True, False])
    def test_sort_values(self, data_for_sorting, ascending):
        ser = pd.Series(data_for_sorting)
        result = ser.sort_values(ascending=ascending)
        expected = ser.iloc[[2, 0, 1]]
        if not ascending:
            expected = expected[::-1]

        self.assert_series_equal(result, expected)

    @pytest.mark.parametrize('ascending', [True, False])
    def test_sort_values_missing(self, data_missing_for_sorting, ascending):
        ser = pd.Series(data_missing_for_sorting)
        result = ser.sort_values(ascending=ascending)
        if ascending:
            expected = ser.iloc[[2, 0, 1]]
        else:
            expected = ser.iloc[[0, 2, 1]]
        self.assert_series_equal(result, expected)

    @pytest.mark.parametrize('ascending', [True, False])
    def test_sort_values_frame(self, data_for_sorting, ascending):
        df = pd.DataFrame({"A": [1, 2, 1],
                           "B": data_for_sorting})
        result = df.sort_values(['A', 'B'])
        expected = pd.DataFrame({"A": [1, 1, 2],
                                 'B': data_for_sorting.take([2, 0, 1])},
                                index=[2, 0, 1])
        self.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('box', [pd.Series, lambda x: x])
    @pytest.mark.parametrize('method', [lambda x: x.unique(), pd.unique])
    def test_unique(self, data, box, method):
        duplicated = box(data._from_sequence([data[0], data[0]]))

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
                                   dtype=np.intp)
        expected_uniques = data_for_grouping.take([0, 4, 7])

        tm.assert_numpy_array_equal(labels, expected_labels)
        self.assert_extension_array_equal(uniques, expected_uniques)

    @pytest.mark.parametrize('na_sentinel', [-1, -2])
    def test_factorize_equivalence(self, data_for_grouping, na_sentinel):
        l1, u1 = pd.factorize(data_for_grouping, na_sentinel=na_sentinel)
        l2, u2 = data_for_grouping.factorize(na_sentinel=na_sentinel)

        tm.assert_numpy_array_equal(l1, l2)
        self.assert_extension_array_equal(u1, u2)

    def test_fillna_copy_frame(self, data_missing):
        arr = data_missing.take([1, 1])
        df = pd.DataFrame({"A": arr})

        filled_val = df.iloc[0, 0]
        result = df.fillna(filled_val)
        assert df.values.base is not result.values.base

        if isinstance(arr, pd.SparseArray):
            assert df.A._values.to_dense() is arr.to_dense()
        else:
            assert df.A._values is arr

    def test_fillna_copy_series(self, data_missing):
        arr = data_missing.take([1, 1])
        ser = pd.Series(arr)

        filled_val = ser[0]
        result = ser.fillna(filled_val)
        assert ser._values is not result._values

        if isinstance(arr, pd.SparseArray):
            assert ser._values.to_dense() is arr.to_dense()
        else:
            assert ser._values is arr

    def test_fillna_length_mismatch(self, data_missing):
        with (tm.assert_raises_regex(ValueError,
              "Length of 'value' does not match.")):
            data_missing.fillna(data_missing.take([1]))

    def test_combine_le(self, data_repeated):
        # GH 20825
        # Test that combine works when doing a <= (le) comparison
        orig_data1, orig_data2 = data_repeated(2)
        s1 = pd.Series(orig_data1)
        s2 = pd.Series(orig_data2)
        result = s1.combine(s2, lambda x1, x2: x1 <= x2)
        expected = pd.Series([a <= b for (a, b) in
                              zip(list(orig_data1), list(orig_data2))])
        self.assert_series_equal(result, expected)

        val = s1.iloc[0]
        result = s1.combine(val, lambda x1, x2: x1 <= x2)
        expected = pd.Series([a <= val for a in list(orig_data1)])
        self.assert_series_equal(result, expected)

    def test_combine_add(self, data_repeated):
        # GH 20825
        orig_data1, orig_data2 = data_repeated(2)
        s1 = pd.Series(orig_data1)
        s2 = pd.Series(orig_data2)
        result = s1.combine(s2, lambda x1, x2: x1 + x2)
        with np.errstate(over='ignore'):
            expected = pd.Series(
                orig_data1._from_sequence([a + b for (a, b) in
                                           zip(list(orig_data1),
                                               list(orig_data2))]))
        self.assert_series_equal(result, expected)

        val = s1.iloc[0]
        result = s1.combine(val, lambda x1, x2: x1 + x2)
        expected = pd.Series(
            orig_data1._from_sequence([a + val for a in list(orig_data1)]))
        self.assert_series_equal(result, expected)

    @pytest.mark.parametrize('frame', [True, False])
    @pytest.mark.parametrize('periods, indices', [
        (-2, [2, 3, 4, -1, -1]),
        (0, [0, 1, 2, 3, 4]),
        (2, [-1, -1, 0, 1, 2]),
    ])
    def test_container_shift(self, data, frame, periods, indices):
        # https://github.com/pandas-dev/pandas/issues/22386
        subset = data[:5]
        data = pd.Series(subset, name='A')
        expected = pd.Series(subset.take(indices, allow_fill=True), name='A')

        if frame:
            result = data.to_frame(name='A').assign(B=1).shift(periods)
            expected = pd.concat([
                expected,
                pd.Series([1] * 5, name='B').shift(periods)
            ], axis=1)
            compare = self.assert_frame_equal
        else:
            result = data.shift(periods)
            compare = self.assert_series_equal

        compare(result, expected)

    @pytest.mark.parametrize("as_frame", [True, False])
    def test_hash_pandas_object_works(self, data, as_frame):
        # https://github.com/pandas-dev/pandas/issues/23066
        data = pd.Series(data)
        if as_frame:
            data = data.to_frame()
        a = pd.util.hash_pandas_object(data)
        b = pd.util.hash_pandas_object(data)
        self.assert_equal(a, b)
