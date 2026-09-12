import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestValues:
    @pytest.mark.parametrize(
        "data",
        [
            pd.period_range("2000", periods=4),
            pd.IntervalIndex.from_breaks([1, 2, 3, 4]),
        ],
    )
    def test_values_object_extension_dtypes(self, data):
        # https://github.com/pandas-dev/pandas/issues/23995
        msg = "Series.values returning an object-dtype ndarray"
        with tm.assert_produces_warning(pd.errors.Pandas4Warning, match=msg):
            result = pd.Series(data).values
        expected = np.array(data.astype(object))
        tm.assert_numpy_array_equal(result, expected)

    def test_values(self, datetime_series):
        tm.assert_almost_equal(
            datetime_series.values, list(datetime_series), check_dtype=False
        )
