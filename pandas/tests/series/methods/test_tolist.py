import pytest

from pandas import (
    Interval,
    Period,
    Series,
    Timedelta,
    Timestamp,
)


@pytest.mark.parametrize(
    "values, dtype, expected_dtype",
    (
        ([1], "int64", int),
        ([1], "Int64", int),
        ([1], "int64[pyarrow]", int),
        ([1.0], "float64", float),
        ([1.0], "Float64", float),
        ([1.0], "float64[pyarrow]", float),
        (["abc"], "object", str),
        (["abc"], "string", str),
        (["abc"], "string[pyarrow]", str),
        ([Interval(1, 3)], "interval", Interval),
        ([Period("2000-01-01", "D")], "period[D]", Period),
        ([Timedelta(days=1)], "timedelta64[ns]", Timedelta),
        ([Timestamp("2000-01-01")], "datetime64[ns]", Timestamp),
    ),
)
def test_tolist_scalar_dtype(values, dtype, expected_dtype):
    # GH49890
    ser = Series(values, dtype=dtype)
    result_dtype = type(ser.tolist()[0])
    assert result_dtype == expected_dtype
