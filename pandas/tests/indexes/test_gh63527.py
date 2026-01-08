import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_pyarrow_index_datetime_properties_delegation():
    """
    Test that PyArrow-backed Index delegates datetime properties
    (like .year, .month) correctly.
    GH#63527
    """
    pytest.importorskip("pyarrow")
    dates = ["2023-01-01", "2023-06-15", "2023-12-31"]
    dti = pd.to_datetime(dates)

    df = pd.DataFrame({"a": [1, 2, 3]}, index=dti)
    df_pa = df.convert_dtypes(dtype_backend="pyarrow")

    if not isinstance(df_pa.index.dtype, pd.ArrowDtype):
        df_pa.index = df_pa.index.astype("timestamp[ns][pyarrow]")

    idx_pa = df_pa.index

    properties = ["year", "month", "day", "is_month_start", "quarter"]

    for prop in properties:
        expected = getattr(dti, prop)
        result = getattr(idx_pa, prop)

        res_array = np.asarray(result)
        exp_array = np.asarray(expected)

        tm.assert_numpy_array_equal(res_array, exp_array, check_dtype=False)

        if hasattr(result, "dtype"):
            assert isinstance(result.dtype, pd.ArrowDtype)


def test_pyarrow_index_groupby_functionality():
    """
    Ensure the properties work inside a groupby operation.
    """
    pytest.importorskip("pyarrow")
    dates = pd.to_datetime(["2021-01-01", "2021-01-02", "2021-02-01"])
    df = pd.DataFrame({"val": [10, 20, 30]}, index=dates)

    df_pa = df.convert_dtypes(dtype_backend="pyarrow")
    if not isinstance(df_pa.index.dtype, pd.ArrowDtype):
        df_pa.index = df_pa.index.astype("timestamp[ns][pyarrow]")

    result = df_pa.groupby(df_pa.index.month)["val"].sum()

    expected_index = pd.Index([1, 2])
    expected_data = [30, 30]

    tm.assert_numpy_array_equal(
        np.asarray(result.index),
        np.asarray(expected_index),
        check_dtype=False,
    )

    assert result.tolist() == expected_data
