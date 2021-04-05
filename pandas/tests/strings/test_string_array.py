import numpy as np
import pytest

from pandas._libs import lib

import pandas as pd
from pandas import (
    DataFrame,
    Series,
    _testing as tm,
)


def test_string_array(any_string_method):
    method_name, args, kwargs = any_string_method
    if method_name == "decode":
        pytest.skip("decode requires bytes.")

    data = ["a", "bb", np.nan, "ccc"]
    a = Series(data, dtype=object)
    b = Series(data, dtype="string")

    expected = getattr(a.str, method_name)(*args, **kwargs)
    result = getattr(b.str, method_name)(*args, **kwargs)

    if isinstance(expected, Series):
        if expected.dtype == "object" and lib.is_string_array(
            expected.dropna().values,
        ):
            assert result.dtype == "string"
            result = result.astype(object)

        elif expected.dtype == "object" and lib.is_bool_array(
            expected.values, skipna=True
        ):
            assert result.dtype == "boolean"
            result = result.astype(object)

        elif expected.dtype == "bool":
            assert result.dtype == "boolean"
            result = result.astype("bool")

        elif expected.dtype == "float" and expected.isna().any():
            assert result.dtype == "Int64"
            result = result.astype("float")

    elif isinstance(expected, DataFrame):
        columns = expected.select_dtypes(include="object").columns
        assert all(result[columns].dtypes == "string")
        result[columns] = result[columns].astype(object)
    tm.assert_equal(result, expected)


@pytest.mark.parametrize(
    "method,expected",
    [
        ("count", [2, None]),
        ("find", [0, None]),
        ("index", [0, None]),
        ("rindex", [2, None]),
    ],
)
def test_string_array_numeric_integer_array(method, expected):
    s = Series(["aba", None], dtype="string")
    result = getattr(s.str, method)("a")
    expected = Series(expected, dtype="Int64")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method,expected",
    [
        ("isdigit", [False, None, True]),
        ("isalpha", [True, None, False]),
        ("isalnum", [True, None, True]),
        ("isdigit", [False, None, True]),
    ],
)
def test_string_array_boolean_array(method, expected):
    s = Series(["a", None, "1"], dtype="string")
    result = getattr(s.str, method)()
    expected = Series(expected, dtype="boolean")
    tm.assert_series_equal(result, expected)


def test_string_array_extract():
    # https://github.com/pandas-dev/pandas/issues/30969
    # Only expand=False & multiple groups was failing
    a = Series(["a1", "b2", "cc"], dtype="string")
    b = Series(["a1", "b2", "cc"], dtype="object")
    pat = r"(\w)(\d)"

    result = a.str.extract(pat, expand=False)
    expected = b.str.extract(pat, expand=False)
    assert all(result.dtypes == "string")

    result = result.astype(object)
    tm.assert_equal(result, expected)


def test_str_get_stringarray_multiple_nans():
    s = Series(pd.array(["a", "ab", pd.NA, "abc"]))
    result = s.str.get(2)
    expected = Series(pd.array([pd.NA, pd.NA, pd.NA, "c"]))
    tm.assert_series_equal(result, expected)
