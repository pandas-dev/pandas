import pytest

import pandas.core.dtypes.concat as _concat

import pandas as pd
from pandas import Series
import pandas._testing as tm


def test_concat_mismatched_categoricals_with_empty():
    # concat_compat behavior on series._values should match pd.concat on series
    ser1 = Series(["a", "b", "c"], dtype="category")
    ser2 = Series([], dtype="category")

    result = _concat.concat_compat([ser1._values, ser2._values])
    expected = pd.concat([ser1, ser2])._values
    tm.assert_numpy_array_equal(result, expected)


def test_concat_single_dataframe_tz_aware():
    # https://github.com/pandas-dev/pandas/issues/25257
    df = pd.DataFrame(
        {"timestamp": [pd.Timestamp("2020-04-08 09:00:00.709949+0000", tz="UTC")]}
    )
    expected = df.copy()
    result = pd.concat([df])
    tm.assert_frame_equal(result, expected)


def test_concat_periodarray_2d():
    pi = pd.period_range("2016-01-01", periods=36, freq="D")
    arr = pi._data.reshape(6, 6)

    result = _concat.concat_compat([arr[:2], arr[2:]], axis=0)
    tm.assert_period_array_equal(result, arr)

    result = _concat.concat_compat([arr[:, :2], arr[:, 2:]], axis=1)
    tm.assert_period_array_equal(result, arr)

    msg = (
        "all the input array dimensions.* for the concatenation axis must match exactly"
    )
    with pytest.raises(ValueError, match=msg):
        _concat.concat_compat([arr[:, :2], arr[:, 2:]], axis=0)

    with pytest.raises(ValueError, match=msg):
        _concat.concat_compat([arr[:2], arr[2:]], axis=1)


def test_concat_series_between_empty_and_tzaware_series(using_infer_string):
    tzaware_time = pd.Timestamp("2020-01-01T00:00:00+00:00")
    ser1 = Series(index=[tzaware_time], data=0, dtype=float)
    ser2 = Series(dtype=float)

    result = pd.concat([ser1, ser2], axis=1)
    expected = pd.DataFrame(
        data=[
            (0.0, None),
        ],
        index=[tzaware_time]
        if using_infer_string
        else pd.Index([tzaware_time], dtype=object),
        columns=[0, 1],
        dtype=float,
    )
    tm.assert_frame_equal(result, expected)
