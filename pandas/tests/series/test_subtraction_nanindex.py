import numpy as np

import pandas as pd
import pandas.testing as tm


def test_series_subtraction_with_nan_and_levels():
    ix1 = pd.MultiIndex.from_arrays(
        [
            [np.nan, 81, 81, 82, 82],
            [np.nan] * 5,
            pd.to_datetime(
                [np.nan, "2018-06-01", "2018-07-01", "2018-07-01", "2018-08-01"]
            ),
        ],
        names=["foo", "bar", "date"],
    )

    s1 = pd.Series([np.nan, 25.058969, 22.519751, 20.847981, 21.625236], index=ix1)

    ix2 = pd.Index([81, 82, 83, 84, 85, 86, 87], name="foo")
    s2 = pd.Series(
        [28.2800, 25.2500, 22.2200, 16.7660, 14.0087, 14.9480, 29.2900], index=ix2
    )

    expected = pd.Series(
        [np.nan, -3.221031, -5.760249, -4.402019, -3.624764], index=ix1, dtype="float64"
    )

    result = s1 - s2

    result = result.astype("float64")

    tm.assert_series_equal(result, expected)
