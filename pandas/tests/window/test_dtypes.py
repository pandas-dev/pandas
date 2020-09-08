import numpy as np
import pytest

from pandas import DataFrame, Series
import pandas._testing as tm
from pandas.core.base import DataError

# gh-12373 : rolling functions error on float32 data
# make sure rolling functions works for different dtypes
#
# further note that we are only checking rolling for fully dtype
# compliance (though both expanding and ewm inherit)


def get_dtype(dtype, coerce_int=None):
    if coerce_int is False and "int" in dtype:
        return None
    if dtype != "category":
        return np.dtype(dtype)
    return dtype


@pytest.mark.parametrize(
    "method, data, expected_data, coerce_int",
    [
        ("count", np.arange(5), [1, 2, 2, 2, 2], True),
        ("count", np.arange(10, 0, -2), [1, 2, 2, 2, 2], True),
        ("count", [0, 1, 2, np.nan, 4], [1, 2, 2, 1, 1], False),
        ("max", np.arange(5), [np.nan, 1, 2, 3, 4], True),
        ("max", np.arange(10, 0, -2), [np.nan, 10, 8, 6, 4], True),
        ("max", [0, 1, 2, np.nan, 4], [np.nan, 1, 2, np.nan, np.nan], False),
        ("min", np.arange(5), [np.nan, 0, 1, 2, 3], True),
        ("min", np.arange(10, 0, -2), [np.nan, 8, 6, 4, 2], True),
        ("min", [0, 1, 2, np.nan, 4], [np.nan, 0, 1, np.nan, np.nan], False),
        ("sum", np.arange(5), [np.nan, 1, 3, 5, 7], True),
        ("sum", np.arange(10, 0, -2), [np.nan, 18, 14, 10, 6], True),
        ("sum", [0, 1, 2, np.nan, 4], [np.nan, 1, 3, np.nan, np.nan], False),
        ("mean", np.arange(5), [np.nan, 0.5, 1.5, 2.5, 3.5], True),
        ("mean", np.arange(10, 0, -2), [np.nan, 9, 7, 5, 3], True),
        ("mean", [0, 1, 2, np.nan, 4], [np.nan, 0.5, 1.5, np.nan, np.nan], False),
        ("std", np.arange(5), [np.nan] + [np.sqrt(0.5)] * 4, True),
        ("std", np.arange(10, 0, -2), [np.nan] + [np.sqrt(2)] * 4, True),
        (
            "std",
            [0, 1, 2, np.nan, 4],
            [np.nan] + [np.sqrt(0.5)] * 2 + [np.nan] * 2,
            False,
        ),
        ("var", np.arange(5), [np.nan, 0.5, 0.5, 0.5, 0.5], True),
        ("var", np.arange(10, 0, -2), [np.nan, 2, 2, 2, 2], True),
        ("var", [0, 1, 2, np.nan, 4], [np.nan, 0.5, 0.5, np.nan, np.nan], False),
        ("median", np.arange(5), [np.nan, 0.5, 1.5, 2.5, 3.5], True),
        ("median", np.arange(10, 0, -2), [np.nan, 9, 7, 5, 3], True),
        ("median", [0, 1, 2, np.nan, 4], [np.nan, 0.5, 1.5, np.nan, np.nan], False),
    ],
)
def test_series_dtypes(method, data, expected_data, coerce_int, dtypes):
    s = Series(data, dtype=get_dtype(dtypes, coerce_int=coerce_int))
    if dtypes in ("m8[ns]", "M8[ns]") and method != "count":
        msg = "No numeric types to aggregate"
        with pytest.raises(DataError, match=msg):
            getattr(s.rolling(2), method)()
    else:
        result = getattr(s.rolling(2), method)()
        expected = Series(expected_data, dtype="float64")
        tm.assert_almost_equal(result, expected)


@pytest.mark.parametrize(
    "method, expected_data",
    [
        ("count", {0: Series([1, 2, 2, 2, 2]), 1: Series([1, 2, 2, 2, 2])}),
        ("max", {0: Series([np.nan, 2, 4, 6, 8]), 1: Series([np.nan, 3, 5, 7, 9])}),
        ("min", {0: Series([np.nan, 0, 2, 4, 6]), 1: Series([np.nan, 1, 3, 5, 7])}),
        (
            "sum",
            {0: Series([np.nan, 2, 6, 10, 14]), 1: Series([np.nan, 4, 8, 12, 16])},
        ),
        ("mean", {0: Series([np.nan, 1, 3, 5, 7]), 1: Series([np.nan, 2, 4, 6, 8])}),
        (
            "std",
            {
                0: Series([np.nan] + [np.sqrt(2)] * 4),
                1: Series([np.nan] + [np.sqrt(2)] * 4),
            },
        ),
        ("var", {0: Series([np.nan, 2, 2, 2, 2]), 1: Series([np.nan, 2, 2, 2, 2])}),
        ("median", {0: Series([np.nan, 1, 3, 5, 7]), 1: Series([np.nan, 2, 4, 6, 8])}),
    ],
)
def test_dataframe_dtypes(method, expected_data, dtypes):
    if dtypes == "category":
        pytest.skip("Category dataframe testing not implemented.")
    df = DataFrame(np.arange(10).reshape((5, 2)), dtype=get_dtype(dtypes))
    if dtypes in ("m8[ns]", "M8[ns]") and method != "count":
        msg = "No numeric types to aggregate"
        with pytest.raises(DataError, match=msg):
            getattr(df.rolling(2), method)()
    else:
        result = getattr(df.rolling(2), method)()
        expected = DataFrame(expected_data, dtype="float64")
        tm.assert_frame_equal(result, expected)
