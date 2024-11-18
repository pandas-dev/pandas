import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import (
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
@pytest.mark.parametrize("use_string_repr", [True, False])
def test_get_dummies_with_any_string_dtype(
    request, any_string_dtype, any_string_dtype2, use_string_repr, using_infer_string
):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    test_ids = request.node.callspec.id.split("-")
    series_dtype_id = test_ids[0][7:]
    expected_dtype_id = test_ids[1][7:]
    if expected_dtype_id == "object":
        if "pyarrow" in series_dtype_id:
            request.applymarker(
                pytest.mark.xfail(
                    reason=("pyarrow.lib.ArrowTypeError: Expected integer, got bool"),
                    strict=True,
                )
            )
        expected = DataFrame(
            [
                [True, True, False],
                [True, False, True],
                [False, False, False],
            ],
            columns=list("abc"),
            dtype=np.bool_,
        )
    elif expected_dtype_id == "str[pyarrow]" and use_string_repr:
        # data type 'str[pyarrow]' uses pandas.ArrowDtype instead
        expected = DataFrame(
            [
                ["true", "true", "false"],
                ["true", "false", "true"],
                ["false", "false", "false"],
            ],
            columns=list("abc"),
            dtype="str[pyarrow]",
        )
    elif expected_dtype_id == "str[python]" and use_string_repr:
        # data type 'str[python]' not understood"
        expected_dtype_id = str
        if using_infer_string:
            expected = DataFrame(
                [
                    ["True", "True", "False"],
                    ["True", "False", "True"],
                    ["False", "False", "False"],
                ],
                columns=list("abc"),
                dtype=expected_dtype_id,
            )
        else:
            expected = DataFrame(
                [
                    ["T", "T", "F"],
                    ["T", "F", "T"],
                    ["F", "F", "F"],
                ],
                columns=list("abc"),
                dtype=expected_dtype_id,
            )
    else:
        expected = DataFrame(
            [
                ["True", "True", "False"],
                ["True", "False", "True"],
                ["False", "False", "False"],
            ],
            columns=list("abc"),
            dtype=any_string_dtype2,
        )
    if use_string_repr:
        result = s.str.get_dummies("|", dtype=expected_dtype_id)
    else:
        result = s.str.get_dummies("|", dtype=any_string_dtype2)
    tm.assert_frame_equal(result, expected)
