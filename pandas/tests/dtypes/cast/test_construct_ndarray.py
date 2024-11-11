import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.construction import sanitize_array


@pytest.mark.parametrize(
    "values, dtype, expected",
    [
        ([1, 2, 3], None, np.array([1, 2, 3], dtype=np.int64)),
        (np.array([1, 2, 3]), None, np.array([1, 2, 3])),
        (["1", "2", None], None, np.array(["1", "2", None])),
        (["1", "2", None], np.dtype("str"), np.array(["1", "2", None])),
        ([1, 2, None], np.dtype("str"), np.array(["1", "2", None])),
        (
            np.ma.masked_array([1, 2, 3], mask=[False, True, False]),
            "int64",
            np.array([1, 2, 3], dtype=np.int64),
        ),
        (
            np.ma.masked_array([1, 2, 3], mask=[False, True, False]),
            "float64",
            np.array([1, 2, 3], dtype=np.float64),
        ),
        (
            np.ma.masked_array([1, 2, 3], mask=[False, True, False]),
            "UInt64",
            np.array([1, 2, 3], dtype=np.uint64),
        ),
        (
            np.ma.masked_array([1.0, 2.0, 3.0], mask=[False, True, False]),
            "float64",
            np.array([1.0, 2.0, 3.0], dtype=np.float64),
        ),
        (
            np.ma.masked_array([1.0, 2.0, 3.0], mask=[False, True, False]),
            "Int64",
            np.array([1, 2, 3], dtype=np.int64),
        ),
    ],
)
def test_construct_1d_ndarray_preserving_na(
    values, dtype, expected, using_infer_string
):
    result = sanitize_array(values, index=None, dtype=dtype)
    if using_infer_string and expected.dtype == object and dtype is None:
        tm.assert_extension_array_equal(result, pd.array(expected, dtype="str"))
    else:
        tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize("dtype", ["m8[ns]", "M8[ns]"])
def test_construct_1d_ndarray_preserving_na_datetimelike(dtype):
    arr = np.arange(5, dtype=np.int64).view(dtype)
    expected = np.array(list(arr), dtype=object)
    assert all(isinstance(x, type(arr[0])) for x in expected)

    result = sanitize_array(arr, index=None, dtype=np.dtype(object))
    tm.assert_numpy_array_equal(result, expected)


def test_sanitize_masked_array_with_masked_ea(values, dtype, expected):
    result = sanitize_array(values, dtype=dtype)
    tm.assert_numpy_array_equal(result, expected)
