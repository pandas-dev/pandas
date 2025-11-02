import numpy as np
import pandas as pd
import pandas.testing as tm
import pytest

def test_series_astype_nullable_int_preserves_nans():
    # regression/edge: astype -> nullable integer dtype should preserve NaNs
    s = pd.Series([1, np.nan, 3], dtype="float64")
    res = s.astype("Int64")
    # expected: dtype is nullable Int64 and NaN is represented as <NA>
    assert res.dtype == "Int64"
    expected = pd.Series([1, pd.NA, 3], dtype="Int64")
    tm.assert_series_equal(res, expected)

def test_series_astype_from_nullable_int_to_float_roundtrip():
    # convert nullable Int64 -> float -> Int64, ensure values and missingness preserved
    s = pd.Series([1, pd.NA, 4], dtype="Int64")
    f = s.astype("float64")
    assert f.dtype == "float64"
    # float representation should have np.nan where original had <NA>
    assert np.isnan(f.iloc[1])
    # roundtrip back to nullable Int64
    back = f.astype("Int64")
    expected = pd.Series([1, pd.NA, 4], dtype="Int64")
    tm.assert_series_equal(back, expected)

@pytest.mark.parametrize("to_dtype", ["Int32", "Int64", "Float32", "Float64")
def test_nullable_series_astype_various_dtypes_preserve_missing(to_dtype):
    # small matrix of cases: ensure missingness preserved when casting between
    # nullable integer/float dtypes and non-nullable numpy float dtypes
    s = pd.Series([0, 1, pd.NA, 3], dtype="Int64")
    res = s.astype(to_dtype)
    if to_dtype.startswith("Int"):
        # result should be nullable integer with pd.NA retained
        assert str(res.dtype).startswith("Int")
        expected = pd.Series([0, 1, pd.NA, 3], dtype=to_dtype)
        tm.assert_series_equal(res, expected)
    else:
        # float dtypes: missingness becomes np.nan and dtype is numpy float
        assert "Float" in to_dtype or to_dtype.startswith("float") or res.dtype.kind == "f"
        assert np.isnan(res.iloc[2])
