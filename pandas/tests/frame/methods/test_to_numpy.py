import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame
import pandas._testing as tm


@pytest.mark.parametrize(
    "data",
    [
        {"a": [1, 2, 3], "b": [1, 2, None]},
        {"a": np.array([1, 2, 3]), "b": np.array([1, 2, np.nan])},
        {"a": pd.array([1, 2, 3]), "b": pd.array([1, 2, None])},
    ],
)
@pytest.mark.parametrize("dtype, na_value", [(float, np.nan), (object, None)])
def test_to_numpy_dataframe_na_value(data, dtype, na_value):
    # https://github.com/pandas-dev/pandas/issues/33820
    df = DataFrame(data)
    result = df.to_numpy(dtype=dtype, na_value=na_value)
    expected = np.array([[1, 1], [2, 2], [3, na_value]], dtype=dtype)
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize(
    "data, expected",
    [
        (
            {"a": pd.array([1, 2, None])},
            np.array([[1.0], [2.0], [np.nan]], dtype=float),
        ),
        (
            {"a": [1, 2, 3], "b": [1, 2, 3]},
            np.array([[1, 1], [2, 2], [3, 3]], dtype=float),
        ),
    ],
)
def test_to_numpy_dataframe_single_block(data, expected):
    # https://github.com/pandas-dev/pandas/issues/33820
    df = DataFrame(data)
    result = df.to_numpy(dtype=float, na_value=np.nan)
    tm.assert_numpy_array_equal(result, expected)


def test_to_numpy_dataframe_single_block_no_mutate():
    # https://github.com/pandas-dev/pandas/issues/33820
    result = DataFrame(np.array([1.0, 2.0, np.nan]))
    expected = DataFrame(np.array([1.0, 2.0, np.nan]))
    result.to_numpy(na_value=0.0)
    tm.assert_frame_equal(result, expected)
