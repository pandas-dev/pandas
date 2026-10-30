import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_astype(using_infer_string):
    # with missing values
    arr = pd.array([True, False, None], dtype="boolean")

    with pytest.raises(ValueError, match="cannot convert NA to integer"):
        arr.astype("int64")

    with pytest.raises(ValueError, match="cannot convert float NaN to"):
        arr.astype("bool")

    result = arr.astype("float64")
    expected = np.array([1, 0, np.nan], dtype="float64")
    tm.assert_numpy_array_equal(result, expected)

    result = arr.astype("str")
    if using_infer_string:
        expected = pd.array(
            ["True", "False", None], dtype=pd.StringDtype(na_value=np.nan)
        )
        tm.assert_extension_array_equal(result, expected)
    else:
        expected = np.array(["True", "False", "<NA>"], dtype=f"{tm.ENDIAN}U5")
        tm.assert_numpy_array_equal(result, expected)

    # no missing values
    arr = pd.array([True, False, True], dtype="boolean")
    result = arr.astype("int64")
    expected = np.array([1, 0, 1], dtype="int64")
    tm.assert_numpy_array_equal(result, expected)

    result = arr.astype("bool")
    expected = np.array([True, False, True], dtype="bool")
    tm.assert_numpy_array_equal(result, expected)


def test_astype_to_boolean_array():
    # astype to BooleanArray
    arr = pd.array([True, False, None], dtype="boolean")

    result = arr.astype("boolean")
    tm.assert_extension_array_equal(result, arr)
    result = arr.astype(pd.BooleanDtype())
    tm.assert_extension_array_equal(result, arr)


@pytest.mark.parametrize(
    "values, dtype, expected",
    [
        ([0, 1, 2], "int64", [False, True, True]),
        ([0.0, 0.5, np.nan], "float64", [False, True, None]),
        (["", "False", None], "object", [False, True, None]),
        (["", "False", None], "string", [False, True, None]),
        (["", "False", None], pd.StringDtype(na_value=np.nan), [False, True, None]),
        ([0, 2, None], "category", [False, True, None]),
        (["", "x", None], "category", [False, True, None]),
        ([None, None], "category", [None, None]),
        ([0, 5, None], "timedelta64[ns]", [False, True, None]),
        (["1970-01-01", "2000-01-01", None], "M8[ns]", [False, True, None]),
        (["1970-01-01", "2000-01-01", None], "M8[ns, UTC]", [False, True, None]),
        (["1970-01-01", None], "period[D]", [True, None]),
        ([0, 2, None], "int64[pyarrow]", [False, True, None]),
    ],
)
def test_astype_from_other_dtypes_matches_astype_bool(values, dtype, expected):
    # GH#37614, GH#40566
    if dtype == "int64[pyarrow]":
        pytest.importorskip("pyarrow")
    ser = pd.Series(values, dtype=dtype)
    result = ser.astype("boolean")
    expected = pd.Series(expected, dtype="boolean")
    tm.assert_series_equal(result, expected)

    notna = ser.notna()
    tm.assert_series_equal(result[notna], ser[notna].astype(bool).astype("boolean"))


def test_astype_to_integer_array():
    # astype to IntegerArray
    arr = pd.array([True, False, None], dtype="boolean")

    result = arr.astype("Int64")
    expected = pd.array([1, 0, None], dtype="Int64")
    tm.assert_extension_array_equal(result, expected)
