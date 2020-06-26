import numpy as np
import pytest

from pandas.core.dtypes.cast import safe_convert_to_ndarray

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "values, expected",
    [
        (pd.Series([1, 2, 3], dtype=int), np.array([1, 2, 3], dtype=int)),
        (
            # Nullable integer type cast to float to handle missing values
            pd.Series([1, np.NaN, 3], dtype="Int64"),
            np.array([1, np.NaN, 3], dtype=float),
        ),
        (
            # Nullable boolean type cast to float to handle missing values
            pd.Series([True, np.NaN, False], dtype="boolean"),
            np.array([1.0, np.NaN, 0.0], dtype=float),
        ),
        (
            # Normal datetime cast not changed
            pd.to_datetime([2001, None, 2003], format="%Y"),
            np.array(["2001", "NaT", "2003"], dtype="datetime64").astype(
                "datetime64[ns]"
            ),
        ),
        (
            # Extended datetime should be downcast to normal datetime
            pd.to_datetime([2001, None, 2003], format="%Y", utc=True),
            np.array(["2001", "NaT", "2003"], dtype="datetime64").astype(
                "datetime64[ns]"
            ),
        ),
        (
            # Downcast to naive datetime should result in local dates, not UTC
            pd.to_datetime([2001, None, 2003], format="%Y").tz_localize(
                tz="US/Eastern"
            ),
            np.array(["2001", "NaT", "2003"], dtype="datetime64").astype(
                "datetime64[ns]"
            ),
        ),
    ],
)
def test_safe_convert_to_ndarray(values, expected):
    result = safe_convert_to_ndarray(values)
    tm.assert_numpy_array_equal(result, expected)
