"""
This module tests the functionality of StringArray.
Tests for the str accessors are in pandas/tests/strings/test_string_array.py
"""

import re

import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas.core.dtypes.common import is_dtype_equal

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays.string_ import (
    StringArray,
    StringDtype,
)

skip_if_no_pyarrow = td.skip_if_no("pyarrow", min_version="1.0.0")


def test_repr(string_storage):
    with pd.option_context("string_storage", string_storage):
        df = pd.DataFrame({"A": pd.array(["a", pd.NA, "b"], dtype="string")})
    expected = "      A\n0     a\n1  <NA>\n2     b"
    assert repr(df) == expected

    expected = "0       a\n1    <NA>\n2       b\nName: A, dtype: string"
    assert repr(df.A) == expected

    expected = "<StringArray>\n['a', <NA>, 'b']\nLength: 3, dtype: string"
    assert repr(df.A.array) == expected


def test_none_to_nan(string_storage):
    with pd.option_context("string_storage", string_storage):
        a = StringArray._from_sequence(["a", None, "b"])
    assert a[1] is not None
    assert a[1] is pd.NA


def test_setitem_validates(string_storage):
    with pd.option_context("string_storage", string_storage):
        arr = StringArray._from_sequence(["a", "b"])

    if string_storage == "python":
        msg = "Cannot set non-string value '10' into a StringArray."
    else:
        msg = "Scalar must be NA or str"
    with pytest.raises(ValueError, match=msg):
        arr[0] = 10

    if string_storage == "python":
        msg = "Must provide strings."
    else:
        msg = "Scalar must be NA or str"
    with pytest.raises(ValueError, match=msg):
        arr[:] = np.array([1, 2])


def test_setitem_with_scalar_string(string_storage):
    # is_float_dtype considers some strings, like 'd', to be floats
    # which can cause issues.
    with pd.option_context("string_storage", string_storage):
        arr = pd.array(["a", "c"], dtype="string")
    arr[0] = "d"
    expected = pd.array(["d", "c"], dtype="string")
    tm.assert_extension_array_equal(arr, expected)


def test_astype_roundtrip(string_storage, request):
    if string_storage == "pyarrow":
        reason = "ValueError: Could not convert object to NumPy datetime"
        mark = pytest.mark.xfail(reason=reason, raises=ValueError)
        request.node.add_marker(mark)
    else:
        mark = pytest.mark.xfail(
            reason="GH#36153 casting from StringArray to dt64 fails", raises=ValueError
        )
        request.node.add_marker(mark)

    ser = pd.Series(pd.date_range("2000", periods=12))
    ser[0] = None

    with pd.option_context("string_storage", string_storage):
        casted = ser.astype("string")
    assert is_dtype_equal(casted.dtype, "string")

    result = casted.astype("datetime64[ns]")
    tm.assert_series_equal(result, ser)


def test_add(string_storage, string_storage2, string_storage3, request):
    if string_storage == "pyarrow" or string_storage2 == "pyarrow":
        reason = (
            "unsupported operand type(s) for +: 'ArrowStringArray' and "
            "'ArrowStringArray'"
        )
        mark = pytest.mark.xfail(raises=TypeError, reason=reason)
        request.node.add_marker(mark)

    with pd.option_context("string_storage", string_storage):
        a = pd.Series(["a", "b", "c", None, None], dtype="string")

    with pd.option_context("string_storage", string_storage2):
        b = pd.Series(["x", "y", None, "z", None], dtype="string")

    result = a + b
    with pd.option_context("string_storage", string_storage3):
        expected = pd.Series(["ax", "by", None, None, None], dtype="string")
    tm.assert_series_equal(result, expected)

    result = a.add(b)
    tm.assert_series_equal(result, expected)

    result = a.radd(b)
    with pd.option_context("string_storage", string_storage3):
        expected = pd.Series(["xa", "yb", None, None, None], dtype="string")
    tm.assert_series_equal(result, expected)

    result = a.add(b, fill_value="-")
    with pd.option_context("string_storage", string_storage3):
        expected = pd.Series(["ax", "by", "c-", "-z", None], dtype="string")
    tm.assert_series_equal(result, expected)


def test_add_2d(string_storage, request):
    if string_storage == "pyarrow":
        reason = "Failed: DID NOT RAISE <class 'ValueError'>"
        mark = pytest.mark.xfail(raises=None, reason=reason)
        request.node.add_marker(mark)

    with pd.option_context("string_storage", string_storage):
        a = pd.array(["a", "b", "c"], dtype="string")
    b = np.array([["a", "b", "c"]], dtype=object)
    with pytest.raises(ValueError, match="3 != 1"):
        a + b

    s = pd.Series(a)
    with pytest.raises(ValueError, match="3 != 1"):
        s + b


def test_add_sequence(string_storage, string_storage2, request):
    if string_storage == "pyarrow":
        reason = "unsupported operand type(s) for +: 'ArrowStringArray' and 'list'"
        mark = pytest.mark.xfail(raises=TypeError, reason=reason)
        request.node.add_marker(mark)

    with pd.option_context("string_storage", string_storage):
        a = pd.array(["a", "b", None, None], dtype="string")
    other = ["x", None, "y", None]

    result = a + other
    with pd.option_context("string_storage", string_storage2):
        expected = pd.array(["ax", None, None, None], dtype="string")
    tm.assert_extension_array_equal(result, expected)

    result = other + a
    with pd.option_context("string_storage", string_storage2):
        expected = pd.array(["xa", None, None, None], dtype="string")
    tm.assert_extension_array_equal(result, expected)


def test_mul(string_storage, string_storage2, request):
    if string_storage == "pyarrow":
        reason = "unsupported operand type(s) for *: 'ArrowStringArray' and 'int'"
        mark = pytest.mark.xfail(raises=TypeError, reason=reason)
        request.node.add_marker(mark)

    with pd.option_context("string_storage", string_storage):
        a = pd.array(["a", "b", None], dtype="string")
    result = a * 2

    with pd.option_context("string_storage", string_storage2):
        expected = pd.array(["aa", "bb", None], dtype="string")
    tm.assert_extension_array_equal(result, expected)

    result = 2 * a
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.xfail(reason="GH-28527")
def test_add_strings(string_storage, string_storage2):
    with pd.option_context("string_storage", string_storage):
        arr = pd.array(["a", "b", "c", "d"], dtype="string")
    df = pd.DataFrame([["t", "u", "v", "w"]])
    assert arr.__add__(df) is NotImplemented

    result = arr + df
    with pd.option_context("string_storage", string_storage2):
        expected = pd.DataFrame([["at", "bu", "cv", "dw"]]).astype("string")
    tm.assert_frame_equal(result, expected)

    result = df + arr
    with pd.option_context("string_storage", string_storage2):
        expected = pd.DataFrame([["ta", "ub", "vc", "wd"]]).astype("string")
    tm.assert_frame_equal(result, expected)


@pytest.mark.xfail(reason="GH-28527")
def test_add_frame(string_storage, string_storage2):
    with pd.option_context("string_storage", string_storage):
        arr = pd.array(["a", "b", np.nan, np.nan], dtype="string")
    df = pd.DataFrame([["x", np.nan, "y", np.nan]])

    assert arr.__add__(df) is NotImplemented

    result = arr + df
    with pd.option_context("string_storage", string_storage2):
        expected = pd.DataFrame([["ax", np.nan, np.nan, np.nan]]).astype("string")
    tm.assert_frame_equal(result, expected)

    result = df + arr
    with pd.option_context("string_storage", string_storage2):
        expected = pd.DataFrame([["xa", np.nan, np.nan, np.nan]]).astype("string")
    tm.assert_frame_equal(result, expected)


def test_comparison_methods_scalar(all_compare_operators, string_storage):
    op_name = all_compare_operators
    with pd.option_context("string_storage", string_storage):
        a = pd.array(["a", None, "c"], dtype="string")
    other = "a"
    result = getattr(a, op_name)(other)
    expected = np.array([getattr(item, op_name)(other) for item in a], dtype=object)
    expected = pd.array(expected, dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_comparison_methods_scalar_pd_na(all_compare_operators, string_storage):
    op_name = all_compare_operators
    with pd.option_context("string_storage", string_storage):
        a = pd.array(["a", None, "c"], dtype="string")
    result = getattr(a, op_name)(pd.NA)
    expected = pd.array([None, None, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_comparison_methods_scalar_not_string(
    all_compare_operators, string_storage, request
):
    if all_compare_operators not in ["__eq__", "__ne__"]:
        reason = "comparison op not supported between instances of 'str' and 'int'"
        mark = pytest.mark.xfail(raises=TypeError, reason=reason)
        request.node.add_marker(mark)

    op_name = all_compare_operators
    with pd.option_context("string_storage", string_storage):
        a = pd.array(["a", None, "c"], dtype="string")
    other = 42
    result = getattr(a, op_name)(other)
    expected_data = {"__eq__": [False, None, False], "__ne__": [True, None, True]}[
        op_name
    ]
    expected = pd.array(expected_data, dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_comparison_methods_array(all_compare_operators, string_storage, request):
    if string_storage == "pyarrow":
        mark = pytest.mark.xfail(
            raises=AssertionError, reason="left is not an ExtensionArray"
        )
        request.node.add_marker(mark)

    op_name = all_compare_operators

    with pd.option_context("string_storage", string_storage):
        a = pd.array(["a", None, "c"], dtype="string")
    other = [None, None, "c"]
    result = getattr(a, op_name)(other)
    expected = np.empty_like(a, dtype="object")
    expected[-1] = getattr(other[-1], op_name)(a[-1])
    expected = pd.array(expected, dtype="boolean")
    tm.assert_extension_array_equal(result, expected)

    result = getattr(a, op_name)(pd.NA)
    expected = pd.array([None, None, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_constructor_raises(string_storage):
    if string_storage == "python":
        msg = "StringArray requires a sequence of strings or pandas.NA"
    else:
        msg = "Unsupported type '<class 'numpy.ndarray'>' for ArrowStringArray"

    with pytest.raises(ValueError, match=msg):
        StringArray(np.array(["a", "b"], dtype="S1"), storage=string_storage)

    with pytest.raises(ValueError, match=msg):
        StringArray(np.array([]), storage=string_storage)

    with pytest.raises(ValueError, match=msg):
        StringArray(np.array(["a", np.nan], dtype=object), storage=string_storage)

    with pytest.raises(ValueError, match=msg):
        StringArray(np.array(["a", None], dtype=object), storage=string_storage)

    with pytest.raises(ValueError, match=msg):
        StringArray(np.array(["a", pd.NaT], dtype=object), storage=string_storage)


@pytest.mark.parametrize("copy", [True, False])
def test_from_sequence_no_mutate(copy, string_storage, request):
    if string_storage == "pyarrow" and copy is False:
        mark = pytest.mark.xfail(
            raises=AssertionError, reason="numpy array are different"
        )
        request.node.add_marker(mark)

    nan_arr = np.array(["a", np.nan], dtype=object)
    na_arr = np.array(["a", pd.NA], dtype=object)

    with pd.option_context("string_storage", string_storage):
        result = StringArray._from_sequence(nan_arr, copy=copy)

    if string_storage == "pyarrow":
        import pyarrow as pa

        expected = StringArray(
            pa.array(na_arr, type=pa.string(), from_pandas=True), storage="pyarrow"
        )
    else:
        expected = StringArray(na_arr, storage="python")

    tm.assert_extension_array_equal(result, expected)

    expected = nan_arr if copy else na_arr
    tm.assert_numpy_array_equal(nan_arr, expected)


def test_astype_int(string_storage):
    with pd.option_context("string_storage", string_storage):
        arr = pd.array(["1", "2", "3"], dtype="string")
    result = arr.astype("int64")
    expected = np.array([1, 2, 3], dtype="int64")
    tm.assert_numpy_array_equal(result, expected)

    with pd.option_context("string_storage", string_storage):
        arr = pd.array(["1", pd.NA, "3"], dtype="string")
    msg = re.escape("int() argument must be a string, a bytes-like object or a number")
    with pytest.raises(TypeError, match=msg):
        arr.astype("int64")


def test_astype_nullable_int(string_storage):
    with pd.option_context("string_storage", string_storage):
        arr = pd.array(["1", pd.NA, "3"], dtype="string")

    result = arr.astype("Int64")
    expected = pd.array([1, pd.NA, 3], dtype="Int64")
    tm.assert_extension_array_equal(result, expected)


def test_astype_float(string_storage, any_float_allowed_nullable_dtype):
    # Don't compare arrays (37974)
    with pd.option_context("string_storage", string_storage):
        ser = pd.Series(["1.1", pd.NA, "3.3"], dtype="string")
    result = ser.astype(any_float_allowed_nullable_dtype)
    expected = pd.Series([1.1, np.nan, 3.3], dtype=any_float_allowed_nullable_dtype)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.xfail(reason="Not implemented StringArray.sum")
def test_reduce(skipna, string_storage):
    with pd.option_context("string_storage", string_storage):
        arr = pd.Series(["a", "b", "c"], dtype="string")
    result = arr.sum(skipna=skipna)
    assert result == "abc"


@pytest.mark.parametrize("method", ["min", "max"])
@pytest.mark.parametrize("skipna", [True, False])
def test_min_max(method, skipna, string_storage, request):
    if string_storage == "pyarrow":
        mark = pytest.mark.xfail(raises=NotImplementedError, reason="not implemented")
        request.node.add_marker(mark)

    with pd.option_context("string_storage", string_storage):
        arr = pd.Series(["a", "b", "c", None], dtype="string")
    result = getattr(arr, method)(skipna=skipna)
    if skipna:
        expected = "a" if method == "min" else "c"
        assert result == expected
    else:
        assert result is pd.NA


@pytest.mark.parametrize("method", ["min", "max"])
@pytest.mark.parametrize("box", [pd.Series, pd.array])
def test_min_max_numpy(method, box, string_storage, request):
    if string_storage == "pyarrow":
        mark = pytest.mark.xfail(raises=NotImplementedError, reason="not implemented")
        request.node.add_marker(mark)

    with pd.option_context("string_storage", string_storage):
        arr = box(["a", "b", "c", None], dtype="string")
    result = getattr(np, method)(arr)
    expected = "a" if method == "min" else "c"
    assert result == expected


@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.xfail(reason="Not implemented StringArray.sum")
def test_reduce_missing(skipna, string_storage):
    with pd.option_context("string_storage", string_storage):
        arr = pd.Series([None, "a", None, "b", "c", None], dtype="string")
    result = arr.sum(skipna=skipna)
    if skipna:
        assert result == "abc"
    else:
        assert pd.isna(result)


def test_fillna_args(string_storage, string_storage2, request):
    # GH 37987

    if string_storage == "pyarrow":
        reason = (
            "Regex pattern \"Cannot set non-string value '1' into "
            "a StringArray.\" does not match 'Scalar must be NA or str'"
        )
        mark = pytest.mark.xfail(raises=AssertionError, reason=reason)
        request.node.add_marker(mark)

    with pd.option_context("string_storage", string_storage):
        arr = pd.array(["a", pd.NA], dtype="string")

    res = arr.fillna(value="b")
    with pd.option_context("string_storage", string_storage2):
        expected = pd.array(["a", "b"], dtype="string")
    tm.assert_extension_array_equal(res, expected)

    res = arr.fillna(value=np.str_("b"))
    with pd.option_context("string_storage", string_storage2):
        expected = pd.array(["a", "b"], dtype="string")
    tm.assert_extension_array_equal(res, expected)

    msg = "Cannot set non-string value '1' into a StringArray."
    with pytest.raises(ValueError, match=msg):
        arr.fillna(value=1)


@td.skip_if_no("pyarrow")
def test_arrow_array(string_storage):
    # protocol added in 0.15.0
    import pyarrow as pa

    with pd.option_context("string_storage", string_storage):
        data = pd.array(["a", "b", "c"], dtype="string")
    arr = pa.array(data)
    expected = pa.array(list(data), type=pa.string(), from_pandas=True)
    if string_storage == "pyarrow":
        expected = pa.chunked_array(expected)

    assert arr.equals(expected)


@td.skip_if_no("pyarrow")
def test_arrow_roundtrip(string_storage):
    # roundtrip possible from arrow 1.0.0
    import pyarrow as pa

    with pd.option_context("string_storage", string_storage):
        data = pd.array(["a", "b", None], dtype="string")
    df = pd.DataFrame({"a": data})
    table = pa.table(df)
    assert table.field("a").type == "string"
    result = table.to_pandas()
    assert isinstance(result["a"].dtype, StringDtype)
    tm.assert_frame_equal(result, df)
    # ensure the missing value is represented by NA and not np.nan or None
    assert result.loc[2, "a"] is pd.NA


@td.skip_if_no("pyarrow")
def test_arrow_load_from_zero_chunks(string_storage):
    # GH-41040
    import pyarrow as pa

    with pd.option_context("string_storage", string_storage):
        data = pd.array([], dtype="string")
    df = pd.DataFrame({"a": data})
    table = pa.table(df)
    assert table.field("a").type == "string"
    # Instantiate the same table with no chunks at all
    table = pa.table([pa.chunked_array([], type=pa.string())], schema=table.schema)
    result = table.to_pandas()
    assert isinstance(result["a"].dtype, StringDtype)
    tm.assert_frame_equal(result, df)


def test_value_counts_na(string_storage):
    with pd.option_context("string_storage", string_storage):
        arr = pd.array(["a", "b", "a", pd.NA], dtype="string")
    result = arr.value_counts(dropna=False)
    expected = pd.Series([2, 1, 1], index=["a", "b", pd.NA], dtype="Int64")
    tm.assert_series_equal(result, expected)

    result = arr.value_counts(dropna=True)
    expected = pd.Series([2, 1], index=["a", "b"], dtype="Int64")
    tm.assert_series_equal(result, expected)


def test_value_counts_with_normalize(string_storage):
    with pd.option_context("string_storage", string_storage):
        s = pd.Series(["a", "b", "a", pd.NA], dtype="string")
    result = s.value_counts(normalize=True)
    expected = pd.Series([2, 1], index=["a", "b"], dtype="Float64") / 3
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "values, expected",
    [
        (["a", "b", "c"], np.array([False, False, False])),
        (["a", "b", None], np.array([False, False, True])),
    ],
)
def test_use_inf_as_na(values, expected, string_storage):
    # https://github.com/pandas-dev/pandas/issues/33655
    with pd.option_context("string_storage", string_storage):
        values = pd.array(values, dtype="string")
    with pd.option_context("mode.use_inf_as_na", True):
        result = values.isna()
        tm.assert_numpy_array_equal(result, expected)

        result = pd.Series(values).isna()
        expected = pd.Series(expected)
        tm.assert_series_equal(result, expected)

        result = pd.DataFrame(values).isna()
        expected = pd.DataFrame(expected)
        tm.assert_frame_equal(result, expected)


def test_memory_usage(string_storage, request):
    # GH 33963

    if string_storage == "pyarrow":
        pytest.skip("not applicable")

    with pd.option_context("string_storage", string_storage):
        series = pd.Series(["a", "b", "c"], dtype="string")

    assert 0 < series.nbytes <= series.memory_usage() < series.memory_usage(deep=True)


@pytest.mark.parametrize("float_dtype", [np.float16, np.float32, np.float64])
def test_astype_from_float_dtype(float_dtype, string_storage, string_storage2):
    # https://github.com/pandas-dev/pandas/issues/36451
    s = pd.Series([0.1], dtype=float_dtype)
    with pd.option_context("string_storage", string_storage):
        result = s.astype("string")
    with pd.option_context("string_storage", string_storage2):
        expected = pd.Series(["0.1"], dtype="string")
    tm.assert_series_equal(result, expected)


def test_to_numpy_returns_pdna_default(string_storage):
    with pd.option_context("string_storage", string_storage):
        arr = pd.array(["a", pd.NA, "b"], dtype="string")
    result = np.array(arr)
    expected = np.array(["a", pd.NA, "b"], dtype=object)
    tm.assert_numpy_array_equal(result, expected)


def test_to_numpy_na_value(string_storage, nulls_fixture):
    na_value = nulls_fixture
    with pd.option_context("string_storage", string_storage):
        arr = pd.array(["a", pd.NA, "b"], dtype="string")
    result = arr.to_numpy(na_value=na_value)
    expected = np.array(["a", na_value, "b"], dtype=object)
    tm.assert_numpy_array_equal(result, expected)


def test_isin(string_storage, request):
    with pd.option_context("string_storage", string_storage):
        s = pd.Series(["a", "b", None], dtype="string")

    result = s.isin(["a", "c"])
    expected = pd.Series([True, False, False])
    tm.assert_series_equal(result, expected)

    result = s.isin(["a", pd.NA])
    expected = pd.Series([True, False, True])
    tm.assert_series_equal(result, expected)

    result = s.isin([])
    expected = pd.Series([False, False, False])
    tm.assert_series_equal(result, expected)

    result = s.isin(["a", pd.Timestamp.now()])
    expected = pd.Series([True, False, False])
    tm.assert_series_equal(result, expected)
