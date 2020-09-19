import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Index, Series
import pandas.testing as tm

class TestBasicComplexSupport:
    @pytest.mark.parametrize(
        "array,expected",
        [
            (
                [1 + 1j, 0, 1, 1j, 1 + 2j],
                Series([1, 1, 1, 1, 1], index=[1 + 2j, 1 + 1j, 1j, 1, 0]),
            ),
            (
                [1 + 2j, 0, 1j, 1, 1j, 1 + 1j],
                # index is sorted by value counts in descending order by default
                Series([2, 1, 1, 1, 1], index=[1j, 1 + 2j, 1 + 1j, 1, 0]),
            ),
        ],
    )
    def test_value_counts(self, array, expected):
        result = pd.value_counts(array)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "array,expected",
        [
            (
                [1 + 1j, 0, 1, 1j, 1 + 2j, 1 + 2j],
                np.array([(1 + 1j), 0j, (1 + 0j), 1j, (1 + 2j)]),
            )
        ],
    )
    def test_unique(self, array, expected):
        result = pd.unique(array)
        np.testing.assert_array_equal(result, expected)

    @pytest.mark.parametrize(
        "array,expected",
        [
            (
                [0, 1j, 1j, 1, 1 + 1j, 1 + 2j, 1 + 1j],
                Series([False, False, True, False, False, False, True], dtype=bool),
            )
        ],
    )
    def test_duplicated(self, array, expected):
        result = Series(array, dtype=np.complex64).duplicated()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "array,expected",
        [
            (
                [0, 1j, 1j, 1, 1 + 1j, 1 + 2j, 1 + 1j],
                Series([False, True, True, False, True, True, True], dtype=bool),
            )
        ],
    )
    def test_isin(self, array, expected):
        result = Series(array).isin([1j, 1 + 1j, 1 + 2j])
        tm.assert_series_equal(result, expected)

    def test_factorize(self):
        array = [1, 2, 2 + 1j]
        labels, uniques = pd.factorize(array)

        expected_labels = np.array([0, 1, 2], dtype=np.intp)
        np.testing.assert_array_equal(labels, expected_labels)

        expected_uniques = np.array([(1 + 0j), (2 + 0j), (2 + 1j)], dtype=np.complex64)
        np.testing.assert_array_equal(uniques, expected_uniques)

    @pytest.mark.parametrize(
        "frame,expected",
        [
            (
                DataFrame([{"a": 1, "b": 1 + 1j}, {"a": 1, "b": 1 + 2j}]),
                DataFrame(
                    np.array([1, 1], dtype=np.int64),
                    index=Index([(1 + 1j), (1 + 2j)], dtype="object", name="b"),
                    columns=Index(["a"], dtype="object"),
                ),
            )
        ],
    )
    def test_groupby(self, frame, expected):
        result = frame.groupby("b", sort=False).count()
        tm.assert_frame_equal(result, expected)

        # sorting of the index should fail since complex numbers are unordered
        with pytest.raises(
            TypeError,
            match="'<' not supported between instances of 'complex' and 'complex'",
        ):
            frame.groupby("b", sort=True).count()

    @pytest.mark.parametrize(
        "array,expected",
        [
            ([0, 1j, 1, 1, 1 + 1j, 1 + 2j], Series([1], dtype=np.complex128)),
            ([1 + 1j, 2j, 1 + 1j], Series([1 + 1j], dtype=np.complex128)),
        ],
    )
    def test_unimode(self, array, expected):
        result = Series(array).mode()
        tm.assert_series_equal(result, expected)

    # mode tries to sort multimodal series.
    # A warning will be raised since complex numbers
    # are not ordered.
    @pytest.mark.parametrize(
        "array,expected",
        [
            (
                # no modes
                [0, 1j, 1, 1 + 1j, 1 + 2j],
                Series([0, 1, 1j, 1 + 1j, 1 + 2j], dtype=np.complex128),
            ),
            ([1 + 1j, 2j, 1 + 1j, 2j, 3], Series([1 + 1j, 2j], dtype=np.complex128)),
        ],
    )
    def test_multimode(self, array, expected):
        with pytest.warns(UserWarning):
            result = Series(array).mode()
        tm.assert_series_equal(result, expected)
