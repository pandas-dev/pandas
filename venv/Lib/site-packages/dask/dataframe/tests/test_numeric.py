from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from dask.array import Array, from_array
from dask.dataframe import Series, _dask_expr_enabled, from_pandas, to_numeric
from dask.dataframe.utils import pyarrow_strings_enabled
from dask.delayed import Delayed

DASK_EXPR_ENABLED = _dask_expr_enabled()


@pytest.mark.parametrize("arg", ["5", 5, "5 "])
def test_to_numeric_on_scalars(arg):
    output = to_numeric(arg)
    assert isinstance(output, Delayed)
    assert output.compute() == 5


def test_to_numeric_on_dask_array():
    arg = from_array(["1.0", "2", "-3", "5.1"])
    expected = np.array([1.0, 2.0, -3.0, 5.1])
    output = to_numeric(arg)
    expected_dtype = "int64"
    assert output.dtype == expected_dtype
    assert isinstance(output, Array)
    assert list(output.compute()) == list(expected)


def test_to_numeric_on_dask_array_with_meta():
    arg = from_array(["1.0", "2", "-3", "5.1"])
    expected = np.array([1.0, 2.0, -3.0, 5.1])
    output = to_numeric(arg, meta=np.array((), dtype="float64"))
    expected_dtype = "float64"
    assert output.dtype == expected_dtype
    assert isinstance(output, Array)
    assert list(output.compute()) == list(expected)


def test_to_numeric_on_dask_dataframe_series():
    s = pd.Series(["1.0", "2", -3, -5.1])
    arg = from_pandas(s, npartitions=2)
    expected = pd.to_numeric(s)
    output = to_numeric(arg)
    expected_dtype = "int64"
    if pyarrow_strings_enabled():
        # `to_numeric` output depends on input dtype
        expected_dtype = "Int64"
    assert output.dtype == expected_dtype
    assert isinstance(output, Series)
    assert list(output.compute()) == list(expected)


def test_to_numeric_on_dask_dataframe_series_with_meta():
    s = pd.Series(["1.0", "2", -3, -5.1])
    arg = from_pandas(s, npartitions=2)
    expected = pd.to_numeric(s)
    output = to_numeric(arg, meta=pd.Series([], dtype="float64"))
    assert output.dtype == "float64"
    assert isinstance(output, Series)
    assert list(output.compute()) == list(expected)


def test_to_numeric_on_dask_dataframe_dataframe_raises_error():
    s = pd.Series(["1.0", "2", -3, -5.1])
    df = pd.DataFrame({"a": s, "b": s})
    arg = from_pandas(df, npartitions=2)
    with pytest.raises(TypeError, match="arg must be a list, tuple, dask."):
        to_numeric(arg)


def test_to_numeric_raises():
    with pytest.raises(ValueError, match="invalid error value"):
        to_numeric("10", errors="invalid")
    with pytest.raises(KeyError, match="``meta`` is not allowed"):
        to_numeric("10", meta=pd.Series([], dtype="float64"))
