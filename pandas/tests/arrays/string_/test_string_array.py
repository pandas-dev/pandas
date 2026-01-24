import pandas as pd


def test_string_array_none_dtype_consistency():
    arr = pd.array([1, None], dtype=str)
    ser_arr = pd.Series([1, None], dtype=str).array

    assert arr.dtype == ser_arr.dtype

