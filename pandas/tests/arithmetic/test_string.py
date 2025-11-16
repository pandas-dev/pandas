import operator
from pathlib import Path

import numpy as np
import pytest

from pandas.compat import HAS_PYARROW
from pandas.errors import Pandas4Warning
import pandas.util._test_decorators as td

import pandas as pd
from pandas import (
    NA,
    ArrowDtype,
    Series,
    StringDtype,
)
import pandas._testing as tm
from pandas.core.construction import extract_array


def string_dtype_highest_priority(dtype1, dtype2):
    if HAS_PYARROW:
        DTYPE_HIERARCHY = [
            StringDtype("python", na_value=np.nan),
            StringDtype("pyarrow", na_value=np.nan),
            StringDtype("python", na_value=NA),
            StringDtype("pyarrow", na_value=NA),
        ]
    else:
        DTYPE_HIERARCHY = [
            StringDtype("python", na_value=np.nan),
            StringDtype("python", na_value=NA),
        ]

    h1 = DTYPE_HIERARCHY.index(dtype1)
    h2 = DTYPE_HIERARCHY.index(dtype2)
    return DTYPE_HIERARCHY[max(h1, h2)]


def test_eq_all_na():
    pytest.importorskip("pyarrow")
    a = pd.array([NA, NA], dtype=StringDtype("pyarrow"))
    result = a == a
    expected = pd.array([NA, NA], dtype="boolean[pyarrow]")
    tm.assert_extension_array_equal(result, expected)


def test_reversed_logical_ops(any_string_dtype):
    # GH#60234
    dtype = any_string_dtype
    warn = None if dtype == object else Pandas4Warning
    left = Series([True, False, False, True])
    right = Series(["", "", "b", "c"], dtype=dtype)

    msg = "operations between boolean dtype and"
    with tm.assert_produces_warning(warn, match=msg):
        result = left | right
    expected = left | right.astype(bool)
    tm.assert_series_equal(result, expected)

    with tm.assert_produces_warning(warn, match=msg):
        result = left & right
    expected = left & right.astype(bool)
    tm.assert_series_equal(result, expected)

    with tm.assert_produces_warning(warn, match=msg):
        result = left ^ right
    expected = left ^ right.astype(bool)
    tm.assert_series_equal(result, expected)


def test_pathlib_path_division(any_string_dtype, request):
    # GH#61940
    if any_string_dtype == object:
        mark = pytest.mark.xfail(
            reason="with NA present we go through _masked_arith_op which "
            "raises TypeError bc Path is not recognized by lib.is_scalar."
        )
        request.applymarker(mark)

    item = Path("/Users/Irv/")
    ser = Series(["A", "B", NA], dtype=any_string_dtype)

    result = item / ser
    expected = Series([item / "A", item / "B", ser.dtype.na_value], dtype=object)
    tm.assert_series_equal(result, expected)

    result = ser / item
    expected = Series(["A" / item, "B" / item, ser.dtype.na_value], dtype=object)
    tm.assert_series_equal(result, expected)


def test_mixed_object_comparison(any_string_dtype):
    # GH#60228
    dtype = any_string_dtype
    ser = Series(["a", "b"], dtype=dtype)

    mixed = Series([1, "b"], dtype=object)

    result = ser == mixed
    expected = Series([False, True], dtype=bool)
    if dtype == object:
        pass
    elif dtype.storage == "python" and dtype.na_value is NA:
        expected = expected.astype("boolean")
    elif dtype.storage == "pyarrow" and dtype.na_value is NA:
        expected = expected.astype("bool[pyarrow]")

    tm.assert_series_equal(result, expected)


def test_pyarrow_numpy_string_invalid():
    # GH#56008
    pa = pytest.importorskip("pyarrow")
    ser = Series([False, True])
    ser2 = Series(["a", "b"], dtype=StringDtype(na_value=np.nan))
    result = ser == ser2
    expected_eq = Series(False, index=ser.index)
    tm.assert_series_equal(result, expected_eq)

    result = ser != ser2
    expected_ne = Series(True, index=ser.index)
    tm.assert_series_equal(result, expected_ne)

    with pytest.raises(TypeError, match="Invalid comparison"):
        ser > ser2

    # GH#59505
    ser3 = ser2.astype("string[pyarrow]")
    result3_eq = ser3 == ser
    tm.assert_series_equal(result3_eq, expected_eq.astype("bool[pyarrow]"))
    result3_ne = ser3 != ser
    tm.assert_series_equal(result3_ne, expected_ne.astype("bool[pyarrow]"))

    with pytest.raises(TypeError, match="Invalid comparison"):
        ser > ser3

    ser4 = ser2.astype(ArrowDtype(pa.string()))
    result4_eq = ser4 == ser
    tm.assert_series_equal(result4_eq, expected_eq.astype("bool[pyarrow]"))
    result4_ne = ser4 != ser
    tm.assert_series_equal(result4_ne, expected_ne.astype("bool[pyarrow]"))

    with pytest.raises(TypeError, match="Invalid comparison"):
        ser > ser4


def test_mul_bool_invalid(any_string_dtype):
    # GH#62595
    dtype = any_string_dtype
    ser = Series(["a", "b", "c"], dtype=dtype)

    if dtype == object:
        pytest.skip("This is not expect to raise")
    elif dtype.storage == "python":
        msg = "Cannot multiply StringArray by bools. Explicitly cast to integers"
    else:
        msg = "Can only string multiply by an integer"

    with pytest.raises(TypeError, match=msg):
        False * ser
    with pytest.raises(TypeError, match=msg):
        ser * True
    with pytest.raises(TypeError, match=msg):
        ser * np.array([True, False, True], dtype=bool)
    with pytest.raises(TypeError, match=msg):
        np.array([True, False, True], dtype=bool) * ser


def test_add(any_string_dtype, request):
    dtype = any_string_dtype
    if dtype == object:
        mark = pytest.mark.xfail(
            reason="Need to update expected for numpy object dtype"
        )
        request.applymarker(mark)

    a = Series(["a", "b", "c", None, None], dtype=dtype)
    b = Series(["x", "y", None, "z", None], dtype=dtype)

    result = a + b
    expected = Series(["ax", "by", None, None, None], dtype=dtype)
    tm.assert_series_equal(result, expected)

    result = a.add(b)
    tm.assert_series_equal(result, expected)

    result = a.radd(b)
    expected = Series(["xa", "yb", None, None, None], dtype=dtype)
    tm.assert_series_equal(result, expected)

    result = a.add(b, fill_value="-")
    expected = Series(["ax", "by", "c-", "-z", None], dtype=dtype)
    tm.assert_series_equal(result, expected)


def test_add_2d(any_string_dtype, request):
    dtype = any_string_dtype

    if dtype == object or dtype.storage == "pyarrow":
        reason = "Failed: DID NOT RAISE <class 'ValueError'>"
        mark = pytest.mark.xfail(raises=None, reason=reason)
        request.applymarker(mark)

    a = pd.array(["a", "b", "c"], dtype=dtype)
    b = np.array([["a", "b", "c"]], dtype=object)
    with pytest.raises(ValueError, match="3 != 1"):
        a + b

    s = Series(a)
    with pytest.raises(ValueError, match="3 != 1"):
        s + b


def test_add_sequence(any_string_dtype, request, using_infer_string):
    dtype = any_string_dtype
    if (
        dtype != object
        and dtype.storage == "python"
        and dtype.na_value is np.nan
        and HAS_PYARROW
        and using_infer_string
    ):
        mark = pytest.mark.xfail(
            reason="As of GH#62522, the list gets wrapped with sanitize_array, "
            "which casts to a higher-priority StringArray, so we get "
            "NotImplemented."
        )
        request.applymarker(mark)
    if dtype == np.dtype(object) and using_infer_string:
        mark = pytest.mark.xfail(reason="Cannot broadcast list")
        request.applymarker(mark)

    a = pd.array(["a", "b", None, None], dtype=dtype)
    other = ["x", None, "y", None]

    result = a + other
    expected = pd.array(["ax", None, None, None], dtype=dtype)
    tm.assert_extension_array_equal(result, expected)

    result = other + a
    expected = pd.array(["xa", None, None, None], dtype=dtype)
    tm.assert_extension_array_equal(result, expected)


def test_mul(any_string_dtype):
    dtype = any_string_dtype
    a = pd.array(["a", "b", None], dtype=dtype)
    result = a * 2
    expected = pd.array(["aa", "bb", None], dtype=dtype)
    tm.assert_extension_array_equal(result, expected)

    result = 2 * a
    tm.assert_extension_array_equal(result, expected)


def test_add_strings(any_string_dtype, request):
    dtype = any_string_dtype
    if dtype != np.dtype(object):
        mark = pytest.mark.xfail(reason="GH-28527")
        request.applymarker(mark)
    arr = pd.array(["a", "b", "c", "d"], dtype=dtype)
    df = pd.DataFrame([["t", "y", "v", "w"]], dtype=object)
    assert arr.__add__(df) is NotImplemented

    result = arr + df
    expected = pd.DataFrame([["at", "by", "cv", "dw"]]).astype(dtype)
    tm.assert_frame_equal(result, expected)

    result = df + arr
    expected = pd.DataFrame([["ta", "yb", "vc", "wd"]]).astype(dtype)
    tm.assert_frame_equal(result, expected)


@pytest.mark.xfail(reason="GH-28527")
def test_add_frame(dtype):
    arr = pd.array(["a", "b", np.nan, np.nan], dtype=dtype)
    df = pd.DataFrame([["x", np.nan, "y", np.nan]])

    assert arr.__add__(df) is NotImplemented

    result = arr + df
    expected = pd.DataFrame([["ax", np.nan, np.nan, np.nan]]).astype(dtype)
    tm.assert_frame_equal(result, expected)

    result = df + arr
    expected = pd.DataFrame([["xa", np.nan, np.nan, np.nan]]).astype(dtype)
    tm.assert_frame_equal(result, expected)


def test_comparison_methods_scalar(comparison_op, any_string_dtype):
    dtype = any_string_dtype
    op_name = f"__{comparison_op.__name__}__"
    a = pd.array(["a", None, "c"], dtype=dtype)
    other = "a"
    result = getattr(a, op_name)(other)
    if dtype == object or dtype.na_value is np.nan:
        expected = np.array([getattr(item, op_name)(other) for item in a])
        if comparison_op == operator.ne:
            expected[1] = True
        else:
            expected[1] = False
        result = extract_array(result, extract_numpy=True)
        tm.assert_numpy_array_equal(result, expected.astype(np.bool_))
    else:
        expected_dtype = "boolean[pyarrow]" if dtype.storage == "pyarrow" else "boolean"
        expected = np.array([getattr(item, op_name)(other) for item in a], dtype=object)
        expected = pd.array(expected, dtype=expected_dtype)
        tm.assert_extension_array_equal(result, expected)


def test_comparison_methods_scalar_pd_na(comparison_op, any_string_dtype):
    dtype = any_string_dtype
    op_name = f"__{comparison_op.__name__}__"
    a = pd.array(["a", None, "c"], dtype=dtype)
    result = getattr(a, op_name)(NA)

    if dtype == np.dtype(object) or dtype.na_value is np.nan:
        if operator.ne == comparison_op:
            expected = np.array([True, True, True])
        else:
            expected = np.array([False, False, False])
        result = extract_array(result, extract_numpy=True)
        tm.assert_numpy_array_equal(result, expected)
    else:
        expected_dtype = "boolean[pyarrow]" if dtype.storage == "pyarrow" else "boolean"
        expected = pd.array([None, None, None], dtype=expected_dtype)
        tm.assert_extension_array_equal(result, expected)
        tm.assert_extension_array_equal(result, expected)


def test_comparison_methods_scalar_not_string(comparison_op, any_string_dtype):
    op_name = f"__{comparison_op.__name__}__"
    dtype = any_string_dtype

    a = pd.array(["a", None, "c"], dtype=dtype)
    other = 42

    if op_name not in ["__eq__", "__ne__"]:
        with pytest.raises(TypeError, match="Invalid comparison|not supported between"):
            getattr(a, op_name)(other)

        return

    result = getattr(a, op_name)(other)
    result = extract_array(result, extract_numpy=True)

    if dtype == np.dtype(object) or dtype.na_value is np.nan:
        expected_data = {
            "__eq__": [False, False, False],
            "__ne__": [True, True, True],
        }[op_name]
        expected = np.array(expected_data)
        tm.assert_numpy_array_equal(result, expected)
    else:
        expected_data = {"__eq__": [False, None, False], "__ne__": [True, None, True]}[
            op_name
        ]
        expected_dtype = "boolean[pyarrow]" if dtype.storage == "pyarrow" else "boolean"
        expected = pd.array(expected_data, dtype=expected_dtype)
        tm.assert_extension_array_equal(result, expected)


def test_comparison_methods_array(comparison_op, any_string_dtype, any_string_dtype2):
    op_name = f"__{comparison_op.__name__}__"
    dtype = any_string_dtype
    dtype2 = any_string_dtype2

    a = pd.array(["a", None, "c"], dtype=dtype)
    other = pd.array([None, None, "c"], dtype=dtype2)
    result = comparison_op(a, other)
    result = extract_array(result, extract_numpy=True)

    # ensure operation is commutative
    result2 = comparison_op(other, a)
    result2 = extract_array(result2, extract_numpy=True)
    tm.assert_equal(result, result2)

    if (dtype == object or dtype.na_value is np.nan) and (
        dtype2 == object or dtype2.na_value is np.nan
    ):
        if operator.ne == comparison_op:
            expected = np.array([True, True, False])
        else:
            expected = np.array([False, False, False])
            expected[-1] = getattr(other[-1], op_name)(a[-1])
        result = extract_array(result, extract_numpy=True)
        tm.assert_numpy_array_equal(result, expected)

    else:
        if dtype == object:
            max_dtype = dtype2
        elif dtype2 == object:
            max_dtype = dtype
        else:
            max_dtype = string_dtype_highest_priority(dtype, dtype2)
        if max_dtype.storage == "python":
            expected_dtype = "boolean"
        else:
            expected_dtype = "bool[pyarrow]"

        expected = np.full(len(a), fill_value=None, dtype="object")
        expected[-1] = getattr(other[-1], op_name)(a[-1])
        expected = pd.array(expected, dtype=expected_dtype)
        tm.assert_equal(result, expected)


@td.skip_if_no("pyarrow")
def test_comparison_methods_array_arrow_extension(comparison_op, any_string_dtype):
    # Test pd.ArrowDtype(pa.string()) against other string arrays
    import pyarrow as pa

    dtype2 = any_string_dtype

    op_name = f"__{comparison_op.__name__}__"
    dtype = ArrowDtype(pa.string())
    a = pd.array(["a", None, "c"], dtype=dtype)
    other = pd.array([None, None, "c"], dtype=dtype2)
    result = comparison_op(a, other)

    # ensure operation is commutative
    result2 = comparison_op(other, a)
    tm.assert_equal(result, result2)

    expected = pd.array([None, None, True], dtype="bool[pyarrow]")
    expected[-1] = getattr(other[-1], op_name)(a[-1])
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize("box", [pd.array, pd.Index, Series])
def test_comparison_methods_list(comparison_op, any_string_dtype, box, request):
    dtype = any_string_dtype

    if box is pd.array and dtype != object and dtype.na_value is np.nan:
        mark = pytest.mark.xfail(
            reason="After wrapping list, op returns NotImplemented, see GH#62522"
        )
        request.applymarker(mark)

    op_name = f"__{comparison_op.__name__}__"

    a = box(pd.array(["a", None, "c"], dtype=dtype))
    item = "c"
    other = [None, None, "c"]
    result = comparison_op(a, other)

    # ensure operation is commutative
    result2 = comparison_op(other, a)
    tm.assert_equal(result, result2)

    if dtype == np.dtype(object) or dtype.na_value is np.nan:
        if operator.ne == comparison_op:
            expected = np.array([True, True, False])
        else:
            expected = np.array([False, False, False])
            expected[-1] = getattr(item, op_name)(item)
        if box is not pd.Index:
            # if GH#62766 is addressed this check can be removed
            expected = box(expected, dtype=expected.dtype)
        tm.assert_equal(result, expected)

    else:
        expected_dtype = "boolean[pyarrow]" if dtype.storage == "pyarrow" else "boolean"
        expected = np.full(len(a), fill_value=None, dtype="object")
        expected[-1] = getattr(item, op_name)(item)
        expected = pd.array(expected, dtype=expected_dtype)
        expected = extract_array(expected, extract_numpy=True)
        if box is not pd.Index:
            # if GH#62766 is addressed this check can be removed
            expected = tm.box_expected(expected, box)
        tm.assert_equal(result, expected)
