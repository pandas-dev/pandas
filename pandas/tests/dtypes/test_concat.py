import pytest

import pandas.core.dtypes.concat as _concat

import pandas as pd
import pandas._testing as tm


def test_concat_mismatched_categoricals_with_empty():
    # concat_compat behavior on series._values should match pd.concat on series
    ser1 = pd.Series(["a", "b", "c"], dtype="category")
    ser2 = pd.Series([], dtype="category")

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
    ser1 = pd.Series(index=[tzaware_time], data=0, dtype=float)
    ser2 = pd.Series(dtype=float)

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


@pytest.mark.parametrize(
    "value, dtype",
    [(2**60 + 1, "uint64[pyarrow]"), (1.5, "double[pyarrow]"), (1.5, "float64")],
)
def test_concat_arrow_int_with_float_common_type_rounds(value, dtype):
    # GH#58819 ints beyond 2**53 round like numpy instead of raising
    pytest.importorskip("pyarrow")
    ser = pd.Series([2**60 + 1, None], dtype="int64[pyarrow]")
    other = pd.Series([value], dtype=dtype)
    result = pd.concat([ser, other], ignore_index=True)
    expected = pd.Series([float(2**60), None, float(value)], dtype="double[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_concat_numpy_int_with_arrow_float_rounds():
    # GH#58819
    pytest.importorskip("pyarrow")
    ser = pd.Series([2**60 + 1])
    other = pd.Series([1.5], dtype="double[pyarrow]")
    result = pd.concat([ser, other], ignore_index=True)
    expected = pd.Series([float(2**60), 1.5], dtype="double[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_index_append_arrow_uint64_with_int64_rounds():
    # GH#58819
    pytest.importorskip("pyarrow")
    idx = pd.Index([2**60 + 1], dtype="uint64[pyarrow]")
    other = pd.Index([2**60 + 1], dtype="int64[pyarrow]")
    result = idx.append(other)
    expected = pd.Index([float(2**60)] * 2, dtype="double[pyarrow]")
    tm.assert_index_equal(result, expected)
