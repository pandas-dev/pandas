import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import (
    ArrowDtype,
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


@td.skip_if_no("pyarrow")
def test_get_dummies_with_pyarrow_dtype_int8(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=ArrowDtype(pa.int8()))
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=list("abc"),
        dtype=ArrowDtype(pa.int8()),
    )
    tm.assert_frame_equal(result, expected)


@td.skip_if_no("pyarrow")
def test_get_dummies_with_pyarrow_dtype_uint8(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=ArrowDtype(pa.uint8()))
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=list("abc"),
        dtype=ArrowDtype(pa.uint8()),
    )
    tm.assert_frame_equal(result, expected)


@td.skip_if_no("pyarrow")
def test_get_dummies_with_pyarrow_dtype_int16(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=ArrowDtype(pa.int16()))
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=list("abc"),
        dtype=ArrowDtype(pa.int16()),
    )
    tm.assert_frame_equal(result, expected)


@td.skip_if_no("pyarrow")
def test_get_dummies_with_pyarrow_dtype_uint16(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=ArrowDtype(pa.uint16()))
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=list("abc"),
        dtype=ArrowDtype(pa.uint16()),
    )
    tm.assert_frame_equal(result, expected)


@td.skip_if_no("pyarrow")
def test_get_dummies_with_pyarrow_dtype_int32(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=ArrowDtype(pa.int32()))
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=list("abc"),
        dtype=ArrowDtype(pa.int32()),
    )
    tm.assert_frame_equal(result, expected)


@td.skip_if_no("pyarrow")
def test_get_dummies_with_pyarrow_dtype_uint32(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=ArrowDtype(pa.uint32()))
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=list("abc"),
        dtype=ArrowDtype(pa.uint32()),
    )
    tm.assert_frame_equal(result, expected)


@td.skip_if_no("pyarrow")
def test_get_dummies_with_pyarrow_dtype_int64(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=ArrowDtype(pa.int64()))
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=list("abc"),
        dtype=ArrowDtype(pa.int64()),
    )
    tm.assert_frame_equal(result, expected)


@td.skip_if_no("pyarrow")
def test_get_dummies_with_pyarrow_dtype_uint64(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=ArrowDtype(pa.uint64()))
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=list("abc"),
        dtype=ArrowDtype(pa.uint64()),
    )
    tm.assert_frame_equal(result, expected)


@td.skip_if_no("pyarrow")
def test_get_dummies_with_pyarrow_dtype_bool(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|", dtype=ArrowDtype(pa.bool_()))
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=list("abc"),
        dtype=ArrowDtype(pa.bool_()),
    )
    tm.assert_frame_equal(result, expected)
