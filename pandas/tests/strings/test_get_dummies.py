import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    Series,
    _testing as tm,
)


def test_get_dummies(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|")
    expected = DataFrame([[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=list("abc"))
    tm.assert_frame_equal(result, expected)

    s = Series(["a;b", "a", 7], dtype=any_string_dtype)
    result = s.str.get_dummies(";")
    expected = DataFrame([[0, 1, 1], [0, 1, 0], [1, 0, 0]], columns=list("7ab"))
    tm.assert_frame_equal(result, expected)


def test_get_dummies_index():
    # GH9980, GH8028
    idx = Index(["a|b", "a|c", "b|c"])
    result = idx.str.get_dummies("|")

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
def test_get_dummies_with_str_dtype(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)

    msg = "Only numeric or boolean dtypes are supported for 'dtype'"
    with pytest.raises(ValueError, match=msg):
        s.str.get_dummies("|", dtype=str)

    with pytest.raises(ValueError, match=msg):
        s.str.get_dummies("|", dtype="datetime64[ns]")


def test_get_dummies_empty_and_missing_entries(any_string_dtype):
    # GH#XXXXX an entry that yields no tag must not become a column
    s = Series(["a|b", "", "b", None, "|a|", "a||b"], dtype=any_string_dtype)
    result = s.str.get_dummies("|")
    expected = DataFrame(
        [[1, 1], [0, 0], [0, 1], [0, 0], [1, 0], [1, 1]], columns=list("ab")
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("data", [["", None], [None, None], ["", "|"]])
def test_get_dummies_all_entries_empty(any_string_dtype, data):
    # GH#XXXXX
    s = Series(data, dtype=any_string_dtype)
    result = s.str.get_dummies("|")
    expected = DataFrame(np.empty((2, 0), dtype=np.int64), columns=Index([]))
    tm.assert_frame_equal(result, expected)


@td.skip_if_no("pyarrow")
@pytest.mark.parametrize("pa_type", ["string", "large_string"])
def test_get_dummies_arrow_dtype(pa_type):
    # GH#XXXXX
    import pyarrow as pa

    dtype = pd.ArrowDtype(getattr(pa, pa_type)())
    s = Series(["a|b", "", "b", None, "|a|", "a||b"], dtype=dtype)
    result = s.str.get_dummies("|")
    expected = DataFrame(
        [[1, 1], [0, 0], [0, 1], [0, 0], [1, 0], [1, 1]],
        columns=list("ab"),
        dtype="bool[pyarrow]",
    )
    tm.assert_frame_equal(result, expected)


def test_get_dummies_categorical():
    # GH#XXXXX missing values must not be encoded as a literal "NaN" tag
    s = Series(["a|NaN", "b", None], dtype="category")
    result = s.str.get_dummies("|")
    expected = DataFrame([[1, 1, 0], [0, 0, 1], [0, 0, 0]], columns=["NaN", "a", "b"])
    tm.assert_frame_equal(result, expected)
