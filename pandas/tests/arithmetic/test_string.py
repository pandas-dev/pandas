from pathlib import Path

import numpy as np
import pytest

from pandas.errors import Pandas4Warning

from pandas import (
    NA,
    ArrowDtype,
    Series,
    StringDtype,
)
import pandas._testing as tm


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
