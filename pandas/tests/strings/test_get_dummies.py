import numpy as np
import pytest

from pandas._config import using_string_dtype

import pandas.util._test_decorators as td

from pandas import (
    NA,
    CategoricalDtype,
    DataFrame,
    Index,
    MultiIndex,
    Series,
    _testing as tm,
)

try:
    import pyarrow as pa
except ImportError:
    pa = None


def test_get_dummies(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|")
    exp_dtype = (
        "boolean"
        if any_string_dtype == "string" and any_string_dtype.na_value is NA
        else "bool"
    )
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=list("abc"), dtype=exp_dtype
    )
    tm.assert_frame_equal(result, expected)

    s = Series(["a;b", "a", 7], dtype=any_string_dtype)
    result = s.str.get_dummies(";")
    expected = DataFrame(
        [[0, 1, 1], [0, 1, 0], [1, 0, 0]], columns=list("7ab"), dtype=exp_dtype
    )
    tm.assert_frame_equal(result, expected)


def test_get_dummies_index():
    # GH9980, GH8028
    idx = Index(["a|b", "a|c", "b|c"])
    result = idx.str.get_dummies("|", dtype=np.int64)

    expected = MultiIndex.from_tuples(
        [(1, 1, 0), (1, 0, 1), (0, 1, 1)], names=("a", "b", "c")
    )
    tm.assert_index_equal(result, expected)


# GH#47872
@pytest.mark.parametrize(
    "dtype",
    [
        np.uint8,
        np.int16,
        np.uint16,
        np.int32,
        np.uint32,
        np.int64,
        np.uint64,
        bool,
        "Int8",
        "Int16",
        "Int32",
        "Int64",
        "boolean",
    ],
)
def test_get_dummies_with_dtype(any_string_dtype, dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=dtype)
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=list("abc"), dtype=dtype
    )
    tm.assert_frame_equal(result, expected)


# GH#47872
@td.skip_if_no("pyarrow")
@pytest.mark.parametrize(
    "dtype",
    [
        "int8[pyarrow]",
        "uint8[pyarrow]",
        "int16[pyarrow]",
        "uint16[pyarrow]",
        "int32[pyarrow]",
        "uint32[pyarrow]",
        "int64[pyarrow]",
        "uint64[pyarrow]",
        "bool[pyarrow]",
    ],
)
def test_get_dummies_with_pyarrow_dtype(any_string_dtype, dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=dtype)
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=list("abc"),
        dtype=dtype,
    )
    tm.assert_frame_equal(result, expected)


# GH#47872
@pytest.mark.xfail(using_string_dtype(), reason="TODO(infer_string)")
def test_get_dummies_with_str_dtype(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=str)
    expected = DataFrame(
        [["T", "T", "F"], ["T", "F", "T"], ["F", "F", "F"]],
        columns=list("abc"),
        dtype=str,
    )
    tm.assert_frame_equal(result, expected)


# GH#47872
@td.skip_if_no("pyarrow")
def test_get_dummies_with_pa_str_dtype(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype="str[pyarrow]")
    expected = DataFrame(
        [
            ["true", "true", "false"],
            ["true", "false", "true"],
            ["false", "false", "false"],
        ],
        columns=list("abc"),
        dtype="str[pyarrow]",
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("dtype_type", ["string", "category"])
def test_get_dummies_ea_dtype(dtype_type, string_dtype_no_object):
    dtype = string_dtype_no_object
    exp_dtype = "boolean" if dtype.na_value is NA else "bool"
    if dtype_type == "category":
        dtype = CategoricalDtype(Index(["a", "b"], dtype))
    s = Series(["a", "b"], dtype=dtype)
    result = s.str.get_dummies()
    expected = DataFrame([[1, 0], [0, 1]], columns=list("ab"), dtype=exp_dtype)
    tm.assert_frame_equal(result, expected)
