import numpy as np
import pytest

import pandas as pd
from pandas import (
    Index,
    Series,
    DataFrame,
)
import pandas.core.algorithms as algos
import pandas.util.testing as tm


class TestComplexSupportBasic:
    @pytest.mark.parametrize("array,expected", [
        (
            [1 + 1j, 0, 1, 1j, 1 + 2j],
            Series([1, 1, 1, 1, 1], index=[1 + 2j, 1 + 1j, 1j, 1, 0])
        ),
        (
            [1 + 2j, 0, 1j, 1, 1j, 1 + 1j],
            # index is sorted by value counts in descending order by default
            Series([2, 1, 1, 1, 1], index=[1j, 1 + 2j, 1 + 1j, 1, 0])
        )
    ])
    def test_value_counts(self, array, expected):
        result = algos.value_counts(array)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("array,expected", [
        (
            [1 + 1j, 0, 1, 1j, 1 + 2j, 1 + 2j],
            np.array([(1 + 1j), 0j, (1 + 0j), 1j, (1 + 2j)], dtype=object)
        ),
    ])
    def test_unique(self, array, expected):
        result = algos.unique(array)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize("array,expected", [
        (
            [0, 1j, 1j, 1, 1 + 1j, 1 + 2j, 1 + 1j],
            Series([False, False, True, False, False, False, True], dtype=bool)
        ),
    ])
    def test_duplicated(self, array, expected):
        result = Series(array).duplicated()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("array,expected", [
        (
            [0, 1j, 1j, 1, 1 + 1j, 1 + 2j, 1 + 1j],
            Series([False, True, True, False, True, True, True], dtype=bool)
        ),
    ])
    def test_isin(self, array, expected):
        result = Series(array).isin([1j, 1 + 1j, 1 + 2j])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("array,expected", [
        (
            [1, 2, 2 + 1j],
            (np.array([0, 1, 2]), np.array([(1 + 0j), (2 + 0j), (2 + 1j)],
                                           dtype=object))
        ),
    ])
    def test_factorize(self, array, expected):
        result = pd.factorize(array)
        assert len(result) == 2

        for i, r in enumerate(result):
            tm.assert_numpy_array_equal(r, expected[i])

    @pytest.mark.parametrize("frame,expected", [
        (
            DataFrame([dict(a=1, b=1 + 1j), dict(a=1, b=1 + 2j)]),
            DataFrame(
                np.array([1, 1]),
                index=Index([(1 + 1j), (1 + 2j)], dtype='object', name='b'),
                columns=Index(['a'], dtype='object'))
        ),
    ])
    def test_groupby(self, frame, expected):
        result = frame.groupby("b", sort=False).count()
        tm.assert_frame_equal(result, expected)
