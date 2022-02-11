import numpy as np
import pytest

from pandas import (
    Int64Dtype,
    IntervalIndex,
    Series,
    period_range,
)
import pandas._testing as tm


class TestValues:
    @pytest.mark.parametrize(
        "data",
        [
            period_range("2000", periods=4),
            IntervalIndex.from_breaks([1, 2, 3, 4]),
        ],
    )
    def test_values_object_extension_dtypes(self, data):
        # https://github.com/pandas-dev/pandas/issues/23995
        result = Series(data).values
        expected = np.array(data.astype(object))
        tm.assert_numpy_array_equal(result, expected)

    def test_values(self, datetime_series):
        tm.assert_almost_equal(
            datetime_series.values, datetime_series, check_dtype=False
        )

    def test_values_multiplying_large_series_by_NA(self):

        Int64NA = Series([np.nan], dtype=Int64Dtype()).iloc[0]

        result = Int64NA * Series(np.zeros(10001))
        expected = Series(np.zeros(10001) * Int64NA)

        tm.assert_series_equal(result, expected)
