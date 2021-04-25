import operator

import numpy as np
import pytest

from pandas._libs import lib

import pandas as pd
from pandas import (
    DataFrame,
    Series,
    _testing as tm,
)


def test_string_array(nullable_string_dtype, any_string_method):
    method_name, args, kwargs = any_string_method
    if method_name == "decode":
        pytest.skip("decode requires bytes.")

    data = ["a", "bb", np.nan, "ccc"]
    a = Series(data, dtype=object)
    b = Series(data, dtype=nullable_string_dtype)

    expected = getattr(a.str, method_name)(*args, **kwargs)
    result = getattr(b.str, method_name)(*args, **kwargs)

    if isinstance(expected, Series):
        if expected.dtype == "object" and lib.is_string_array(
            expected.dropna().values,
        ):
            assert result.dtype == nullable_string_dtype
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
        assert all(result[columns].dtypes == nullable_string_dtype)
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
def test_string_array_numeric_integer_array(nullable_string_dtype, method, expected):
    s = Series(["aba", None], dtype=nullable_string_dtype)
    result = getattr(s.str, method)("a")
    expected = Series(expected, dtype="Int64")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method,expected",
    [
        ("isdigit", [False, None, True]),
        ("isalpha", [True, None, False]),
        ("isalnum", [True, None, True]),
        ("isnumeric", [False, None, True]),
    ],
)
def test_string_array_boolean_array(nullable_string_dtype, method, expected):
    s = Series(["a", None, "1"], dtype=nullable_string_dtype)
    result = getattr(s.str, method)()
    expected = Series(expected, dtype="boolean")
    tm.assert_series_equal(result, expected)


def test_string_array_extract(nullable_string_dtype):
    # https://github.com/pandas-dev/pandas/issues/30969
    # Only expand=False & multiple groups was failing

    a = Series(["a1", "b2", "cc"], dtype=nullable_string_dtype)
    b = Series(["a1", "b2", "cc"], dtype="object")
    pat = r"(\w)(\d)"

    result = a.str.extract(pat, expand=False)
    expected = b.str.extract(pat, expand=False)
    assert all(result.dtypes == nullable_string_dtype)

    result = result.astype(object)
    tm.assert_equal(result, expected)


def test_str_get_stringarray_multiple_nans(nullable_string_dtype):
    s = Series(pd.array(["a", "ab", pd.NA, "abc"], dtype=nullable_string_dtype))
    result = s.str.get(2)
    expected = Series(pd.array([pd.NA, pd.NA, pd.NA, "c"], dtype=nullable_string_dtype))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "input, method",
    [
        (["a", "b", "c"], operator.methodcaller("capitalize")),
        (["a b", "a bc. de"], operator.methodcaller("capitalize")),
    ],
)
def test_capitalize(input, method, nullable_string_dtype):
    a = Series(input, dtype=nullable_string_dtype)
    b = Series(input, dtype="object")
    result = method(a.str)
    expected = method(b.str)

    assert result.dtype.name == nullable_string_dtype
    tm.assert_series_equal(result.astype(object), expected)
