import pytest

import pandas.util._test_decorators as td

import pandas as pd


@pytest.mark.parametrize(
    "values, dtype, expected_dtype",
    (
        ([1], "int64", int),
        ([1], "Int64", int),
        ([1.0], "float64", float),
        ([1.0], "Float64", float),
        (["abc"], "object", str),
        (["abc"], "string", str),
        ([pd.Interval(1, 3)], "interval", pd.Interval),
        ([pd.Period("2000-01-01", "D")], "period[D]", pd.Period),
        ([pd.Timedelta(days=1)], "timedelta64[ns]", pd.Timedelta),
        ([pd.Timestamp("2000-01-01")], "datetime64[ns]", pd.Timestamp),
        pytest.param([1], "int64[pyarrow]", int, marks=td.skip_if_no("pyarrow")),
        pytest.param([1.0], "float64[pyarrow]", float, marks=td.skip_if_no("pyarrow")),
        pytest.param(["abc"], "string[pyarrow]", str, marks=td.skip_if_no("pyarrow")),
    ),
)
def test_tolist_scalar_dtype(values, dtype, expected_dtype):
    # GH49890
    ser = pd.Series(values, dtype=dtype)
    result_dtype = type(ser.tolist()[0])
    assert result_dtype == expected_dtype
