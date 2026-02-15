import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "data",
    [
        ["a", None],  # Original case
        [None, None],  # All nulls
        [],  # Empty list
        ["", None],  # Empty string and null
    ],
)
def test_string_array_construction_consistency_with_series(data):
    # GH#57702: Ensure pd.array and pd.Series(data).values are consistent
    # for StringDtype(storage="python")
    dtype = pd.StringDtype(storage="python")

    result_array = pd.array(data, dtype=dtype)
    result_series_array = pd.Series(data, dtype=dtype).values

    tm.assert_extension_array_equal(result_array, result_series_array)
