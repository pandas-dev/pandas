import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestUnaryOps:
    def test_invert(self):
        a = pd.array([True, False, None], dtype="boolean")
        expected = pd.array([False, True, None], dtype="boolean")
        tm.assert_extension_array_equal(~a, expected)

        expected = pd.Series(expected, index=["a", "b", "c"], name="name")
        result = ~pd.Series(a, index=["a", "b", "c"], name="name")
        tm.assert_series_equal(result, expected)

        df = pd.DataFrame({"A": a, "B": [True, False, False]}, index=["a", "b", "c"])
        result = ~df
        expected = pd.DataFrame(
            {"A": expected, "B": [False, True, True]}, index=["a", "b", "c"]
        )
        tm.assert_frame_equal(result, expected)


def test_repr():
    df = pd.DataFrame({"A": pd.array([True, False, None], dtype="boolean")})
    expected = "       A\n0   True\n1  False\n2   <NA>"
    assert repr(df) == expected

    expected = "0     True\n1    False\n2     <NA>\nName: A, dtype: boolean"
    assert repr(df.A) == expected

    expected = "<BooleanArray>\n[True, False, <NA>]\nLength: 3, dtype: boolean"
    assert repr(df.A.array) == expected


@pytest.mark.parametrize("na", [None, np.nan, pd.NA])
def test_setitem_missing_values(na):
    arr = pd.array([True, False, None], dtype="boolean")
    expected = pd.array([True, None, None], dtype="boolean")
    arr[1] = na
    tm.assert_extension_array_equal(arr, expected)
