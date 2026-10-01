from __future__ import annotations

from datetime import (
    date,
    datetime,
    time,
    timedelta,
)
from decimal import Decimal
import locale

import numpy as np
import pytest

from pandas._libs.tslibs import timezones
from pandas.compat import (
    is_platform_windows,
)
from pandas.compat.pyarrow import pa_version_under22p0
from pandas.errors import (
    OutOfBoundsDatetime,
    OutOfBoundsTimedelta,
    Pandas4Warning,
)

from pandas.core.dtypes.dtypes import (
    ArrowDtype,
)

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")

from pandas.core.arrays.arrow.array import ArrowExtensionArray


def _require_timezone_database(request):
    if is_platform_windows() and pa_version_under22p0:
        mark = pytest.mark.xfail(
            raises=pa.ArrowInvalid,
            reason=(
                "TODO: Set ARROW_TIMEZONE_DATABASE environment variable "
                "on CI to path to the tzdata for pyarrow."
            ),
        )
        request.applymarker(mark)


@pytest.mark.parametrize("pa_type", [pa.date32(), pa.date64()], ids=str)
def test_date_std_sem(pa_type):
    # GH#69752 the result was read as seconds without converting it from days
    # (date32) or milliseconds (date64)
    dates = [date(2020, 1, 1), date(2020, 1, 3), date(2020, 1, 2), date(2020, 1, 6)]
    ser = pd.Series(dates, dtype=ArrowDtype(pa_type))
    unit = "s" if pa.types.is_date32(pa_type) else "ms"
    expected = pd.Series(pd.to_datetime(dates).as_unit(unit))

    result = ser.std()
    assert result == expected.std()
    assert result.unit == unit
    assert ser.std(ddof=0) == expected.std(ddof=0)
    expected_sem = pd.Timedelta(days=1, hours=1, minutes=55, seconds=22)
    if unit == "ms":
        expected_sem += pd.Timedelta(milliseconds=666)
    assert ser.sem() == expected_sem

    frame_result = pd.DataFrame({"a": ser}).std()
    assert frame_result.dtype == ArrowDtype(pa.duration(unit))
    assert frame_result["a"] == result

    with_na = pd.Series([*dates, None], dtype=ArrowDtype(pa_type))
    assert with_na.std() == result
    assert with_na.std(skipna=False) is pd.NA


def test_date64_std_keeps_milliseconds():
    # GH#69752
    ser = pd.Series(ArrowExtensionArray(pa.array([0, 1, 2, 3], pa.date64())))
    expected = pd.Series(np.array([0, 1, 2, 3], dtype="M8[ms]")).std()
    assert ser.std() == expected == pd.Timedelta(milliseconds=1)


@pytest.mark.parametrize("unit", ["ns", "us", "ms", "s"])
def test_duration_from_strings_with_nat(unit):
    # GH51175
    strings = ["1000", "NaT"]
    pa_type = pa.duration(unit)
    dtype = ArrowDtype(pa_type)
    result = ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)
    expected = ArrowExtensionArray(pa.array([1000, None], type=pa_type))
    tm.assert_extension_array_equal(result, expected)


def test_unsupported_dt(data):
    pa_dtype = data.dtype.pyarrow_dtype
    if not pa.types.is_temporal(pa_dtype):
        with pytest.raises(
            AttributeError, match="Can only use .dt accessor with datetimelike values"
        ):
            pd.Series(data).dt


@pytest.mark.parametrize(
    "prop, expected",
    [
        ["year", 2023],
        ["day", 2],
        ["day_of_week", 0],
        ["weekday", 0],
        ["day_of_year", 2],
        ["hour", 3],
        ["minute", 4],
        ["is_leap_year", False],
        ["microsecond", 2000],
        ["month", 1],
        ["nanosecond", 6],
        ["quarter", 1],
        ["second", 7],
        ["date", date(2023, 1, 2)],
        ["time", time(3, 4, 7, 2000)],
    ],
)
def test_dt_properties(prop, expected):
    ser = pd.Series(
        [
            pd.Timestamp(
                year=2023,
                month=1,
                day=2,
                hour=3,
                minute=4,
                second=7,
                microsecond=2000,
                nanosecond=6,
            ),
            None,
        ],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    if prop == "weekday":
        # GH#12816
        warn = Pandas4Warning
    else:
        warn = None
    with tm.assert_produces_warning(warn, match="weekday"):
        result = getattr(ser.dt, prop)
    exp_type = None
    if isinstance(expected, date):
        exp_type = pa.date32()
    elif isinstance(expected, time):
        exp_type = pa.time64("ns")
    expected = pd.Series(ArrowExtensionArray(pa.array([expected, None], type=exp_type)))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("microsecond", [2000, 5, 0])
def test_dt_microsecond(microsecond):
    # GH 59183
    ser = pd.Series(
        [
            pd.Timestamp(
                year=2024,
                month=7,
                day=7,
                second=5,
                microsecond=microsecond,
                nanosecond=6,
            ),
            None,
        ],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    result = ser.dt.microsecond
    expected = pd.Series([microsecond, None], dtype="int64[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_dt_is_month_start_end():
    ser = pd.Series(
        [
            datetime(year=2023, month=12, day=2, hour=3),
            datetime(year=2023, month=1, day=1, hour=3),
            datetime(year=2023, month=3, day=31, hour=3),
            None,
        ],
        dtype=ArrowDtype(pa.timestamp("us")),
    )
    result = ser.dt.is_month_start
    expected = pd.Series([False, True, False, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)

    result = ser.dt.is_month_end
    expected = pd.Series([False, False, True, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


def test_dt_is_year_start_end():
    ser = pd.Series(
        [
            datetime(year=2023, month=12, day=31, hour=3),
            datetime(year=2023, month=1, day=1, hour=3),
            datetime(year=2023, month=3, day=31, hour=3),
            None,
        ],
        dtype=ArrowDtype(pa.timestamp("us")),
    )
    result = ser.dt.is_year_start
    expected = pd.Series([False, True, False, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)

    result = ser.dt.is_year_end
    expected = pd.Series([True, False, False, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


def test_dt_is_quarter_start_end():
    ser = pd.Series(
        [
            datetime(year=2023, month=11, day=30, hour=3),
            datetime(year=2023, month=1, day=1, hour=3),
            datetime(year=2023, month=3, day=31, hour=3),
            None,
        ],
        dtype=ArrowDtype(pa.timestamp("us")),
    )
    result = ser.dt.is_quarter_start
    expected = pd.Series([False, True, False, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)

    result = ser.dt.is_quarter_end
    expected = pd.Series([False, False, True, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("method", ["days_in_month"])
def test_dt_days_in_month(method):
    ser = pd.Series(
        [
            datetime(year=2023, month=3, day=30, hour=3),
            datetime(year=2023, month=4, day=1, hour=3),
            datetime(year=2023, month=2, day=3, hour=3),
            None,
        ],
        dtype=ArrowDtype(pa.timestamp("us")),
    )
    result = getattr(ser.dt, method)
    expected = pd.Series([31, 30, 28, None], dtype=ArrowDtype(pa.int64()))
    tm.assert_series_equal(result, expected)


def test_dt_normalize():
    ser = pd.Series(
        [
            datetime(year=2023, month=3, day=30),
            datetime(year=2023, month=4, day=1, hour=3),
            datetime(year=2023, month=2, day=3, hour=23, minute=59, second=59),
            None,
        ],
        dtype=ArrowDtype(pa.timestamp("us")),
    )
    result = ser.dt.normalize()
    expected = pd.Series(
        [
            datetime(year=2023, month=3, day=30),
            datetime(year=2023, month=4, day=1),
            datetime(year=2023, month=2, day=3),
            None,
        ],
        dtype=ArrowDtype(pa.timestamp("us")),
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("unit", ["us", "ns"])
def test_dt_time_preserve_unit(unit):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp(unit)),
    )
    assert ser.dt.unit == unit

    result = ser.dt.time
    expected = pd.Series(
        ArrowExtensionArray(pa.array([time(3, 0), None], type=pa.time64(unit)))
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("tz", [None, "UTC", "US/Pacific"])
def test_dt_tz(tz):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns", tz=tz)),
    )
    result = ser.dt.tz
    assert result == timezones.maybe_get_tz(tz)


def test_dt_isocalendar():
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    result = ser.dt.isocalendar()
    expected = pd.DataFrame(
        [[2023, 1, 1], [0, 0, 0]],
        columns=["year", "week", "day"],
        dtype="int64[pyarrow]",
    )
    tm.assert_frame_equal(result, expected)


def test_dt_isocalendar_preserves_index():
    # GH#65894
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
        index=["a", "b"],
    )
    result = ser.dt.isocalendar()
    expected = pd.DataFrame(
        [[2023, 1, 1], [0, 0, 0]],
        columns=["year", "week", "day"],
        index=["a", "b"],
        dtype="int64[pyarrow]",
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "method, exp", [["day_name", "Sunday"], ["month_name", "January"]]
)
def test_dt_day_month_name(method, exp, request):
    # GH 52388
    _require_timezone_database(request)

    ser = pd.Series([datetime(2023, 1, 1), None], dtype=ArrowDtype(pa.timestamp("ms")))
    result = getattr(ser.dt, method)()
    expected = pd.Series([exp, None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


def test_dt_strftime(request):
    _require_timezone_database(request)

    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    result = ser.dt.strftime("%Y-%m-%dT%H:%M:%S")
    expected = pd.Series(
        ["2023-01-02T03:00:00.000000000", None], dtype=ArrowDtype(pa.string())
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("method", ["ceil", "floor", "round"])
def test_dt_roundlike_tz_options_not_supported(method):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    with pytest.raises(NotImplementedError, match="ambiguous is not supported."):
        getattr(ser.dt, method)("1h", ambiguous="NaT")

    with pytest.raises(NotImplementedError, match="nonexistent is not supported."):
        getattr(ser.dt, method)("1h", nonexistent="NaT")


@pytest.mark.parametrize("method", ["ceil", "floor", "round"])
def test_dt_roundlike_unsupported_freq(method):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    with pytest.raises(ValueError, match="freq='1B' is not supported"):
        getattr(ser.dt, method)("1B")

    with pytest.raises(ValueError, match="Must specify a valid frequency: None"):
        getattr(ser.dt, method)(None)


@pytest.mark.parametrize("freq", ["D", "h", "min", "s", "ms", "us", "ns"])
@pytest.mark.parametrize("method", ["ceil", "floor", "round"])
def test_dt_ceil_year_floor(freq, method):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=1), None],
    )
    pa_dtype = ArrowDtype(pa.timestamp("ns"))
    expected = getattr(ser.dt, method)(f"1{freq}").astype(pa_dtype)
    result = getattr(ser.astype(pa_dtype).dt, method)(f"1{freq}")
    tm.assert_series_equal(result, expected)


def test_dt_to_pydatetime():
    # GH 51859
    data = [datetime(2022, 1, 1), datetime(2023, 1, 1)]
    ser = pd.Series(data, dtype=ArrowDtype(pa.timestamp("ns")))
    result = ser.dt.to_pydatetime()
    expected = pd.Series(data, dtype=object)
    tm.assert_series_equal(result, expected)
    assert all(type(expected.iloc[i]) is datetime for i in range(len(expected)))

    expected = ser.astype("datetime64[ns]").dt.to_pydatetime()
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("date_type", [32, 64])
def test_dt_to_pydatetime_date_error(date_type):
    # GH 52812
    ser = pd.Series(
        [date(2022, 12, 31)],
        dtype=ArrowDtype(getattr(pa, f"date{date_type}")()),
    )
    with pytest.raises(ValueError, match="to_pydatetime cannot be called with"):
        ser.dt.to_pydatetime()


def test_dt_tz_localize_unsupported_tz_options():
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    with pytest.raises(NotImplementedError, match="ambiguous='NaT' is not supported"):
        ser.dt.tz_localize("UTC", ambiguous="NaT")

    with pytest.raises(NotImplementedError, match="nonexistent='NaT' is not supported"):
        ser.dt.tz_localize("UTC", nonexistent="NaT")


def test_dt_tz_localize_none(request):
    _require_timezone_database(request)

    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns", tz="US/Pacific")),
    )
    result = ser.dt.tz_localize(None)
    expected = pd.Series(
        [ser[0].tz_localize(None), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("unit", ["us", "ns"])
def test_dt_tz_localize(unit, request):
    _require_timezone_database(request)

    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp(unit)),
    )
    result = ser.dt.tz_localize("US/Pacific")
    exp_data = pa.array(
        [datetime(year=2023, month=1, day=2, hour=3), None], type=pa.timestamp(unit)
    )
    exp_data = pa.compute.assume_timezone(exp_data, "US/Pacific")
    expected = pd.Series(ArrowExtensionArray(exp_data))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "nonexistent, exp_date",
    [
        ["shift_forward", datetime(year=2023, month=3, day=12, hour=3)],
        ["shift_backward", pd.Timestamp("2023-03-12 01:59:59.999999999")],
    ],
)
def test_dt_tz_localize_nonexistent(nonexistent, exp_date, request):
    _require_timezone_database(request)

    ser = pd.Series(
        [datetime(year=2023, month=3, day=12, hour=2, minute=30), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    result = ser.dt.tz_localize("US/Pacific", nonexistent=nonexistent)
    exp_data = pa.array([exp_date, None], type=pa.timestamp("ns"))
    exp_data = pa.compute.assume_timezone(exp_data, "US/Pacific")
    expected = pd.Series(ArrowExtensionArray(exp_data))
    tm.assert_series_equal(result, expected)


def test_dt_tz_convert_not_tz_raises():
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    with pytest.raises(TypeError, match="Cannot convert tz-naive timestamps"):
        ser.dt.tz_convert("UTC")


def test_dt_tz_convert_none():
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns", "US/Pacific")),
    )
    result = ser.dt.tz_convert(None)
    expected = pd.Series(
        [ser[0].tz_convert(None), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("unit", ["us", "ns"])
def test_dt_tz_convert(unit):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp(unit, "US/Pacific")),
    )
    result = ser.dt.tz_convert("US/Eastern")
    expected = pd.Series(
        [ser[0].tz_convert("US/Eastern"), None],
        dtype=ArrowDtype(pa.timestamp(unit, "US/Eastern")),
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("dtype", ["timestamp[ms][pyarrow]", "duration[ms][pyarrow]"])
def test_as_unit(dtype):
    # GH 52284
    ser = pd.Series([1000, None], dtype=dtype)
    result = ser.dt.as_unit("ns")
    expected = ser.astype(dtype.replace("ms", "ns"))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "from_unit,to_unit",
    [
        ("ns", "us"),
        ("ns", "ms"),
        ("ns", "s"),
        ("us", "ms"),
        ("us", "s"),
        ("ms", "s"),
        ("s", "ms"),
        ("s", "us"),
        ("s", "ns"),
        ("ms", "us"),
        ("ms", "ns"),
        ("us", "ns"),
    ],
)
def test_as_unit_duration_truncation(from_unit, to_unit):
    # Test that as_unit truncates correctly (matches NumPy behavior)
    # Value with sub-unit precision to test truncation
    ser_numpy = pd.Series(
        pd.to_timedelta([93784567890123, None], unit="ns").as_unit(from_unit)
    )
    ser_arrow = ser_numpy.astype(f"duration[{from_unit}][pyarrow]")

    result = ser_arrow.dt.as_unit(to_unit)
    expected = ser_numpy.dt.as_unit(to_unit).astype(f"duration[{to_unit}][pyarrow]")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "from_unit,to_unit",
    [
        ("ns", "us"),
        ("ns", "ms"),
        ("ns", "s"),
        ("s", "ns"),
        ("ms", "ns"),
        ("us", "ns"),
    ],
)
def test_as_unit_timestamp(from_unit, to_unit):
    # Test timestamp as_unit matches NumPy behavior
    # Create Arrow series directly to preserve nulls correctly
    ser_arrow = pd.Series(
        [pd.Timestamp("2024-01-15 12:30:45.123456789"), None],
        dtype=f"timestamp[{from_unit}][pyarrow]",
    )
    ser_numpy = ser_arrow.astype(f"datetime64[{from_unit}]")

    result = ser_arrow.dt.as_unit(to_unit)
    expected_numpy = ser_numpy.dt.as_unit(to_unit)
    # Compare values (excluding null handling differences)
    tm.assert_almost_equal(
        result.dropna().to_numpy(dtype=f"datetime64[{to_unit}]"),
        expected_numpy.dropna().to_numpy(),
    )
    # Verify nulls are preserved
    assert result.isna().sum() == ser_arrow.isna().sum()


@pytest.mark.parametrize("to_unit", ["s", "ms", "us", "ns"])
def test_as_unit_timestamp_with_timezone(to_unit):
    # Test that timezone is preserved
    ser_numpy = pd.Series(
        pd.to_datetime(["2024-01-15 12:30:45.123456789"])
        .tz_localize("US/Eastern")
        .as_unit("ns")
    )
    ser_arrow = ser_numpy.astype("timestamp[ns, US/Eastern][pyarrow]")

    result = ser_arrow.dt.as_unit(to_unit)
    expected = ser_numpy.dt.as_unit(to_unit).astype(
        f"timestamp[{to_unit}, US/Eastern][pyarrow]"
    )
    tm.assert_series_equal(result, expected)
    assert str(result.dtype) == f"timestamp[{to_unit}, tz=US/Eastern][pyarrow]"


def test_as_unit_date_raises():
    # as_unit should raise for date types
    ser = pd.Series([1, 2], dtype=ArrowDtype(pa.date32()))
    with pytest.raises(NotImplementedError, match="as_unit not implemented"):
        ser.dt.as_unit("ns")


@pytest.mark.parametrize(
    "from_unit, to_unit",
    [("ns", "s"), ("ns", "ms"), ("ns", "us"), ("us", "s"), ("ms", "s")],
)
def test_as_unit_duration_negative_floors(from_unit, to_unit):
    # GH#63573 downcasting a negative duration must floor toward -inf like
    # numpy, not truncate toward zero
    values = [93784567890123, -93784567890123, None]
    ser_arrow = pd.Series(pd.to_timedelta(values, unit="ns").as_unit(from_unit)).astype(
        f"duration[{from_unit}][pyarrow]"
    )

    result = ser_arrow.dt.as_unit(to_unit)
    expected = pd.Series(
        pd.to_timedelta(values, unit="ns").as_unit(from_unit).as_unit(to_unit)
    ).astype(f"duration[{to_unit}][pyarrow]")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("to_unit", ["s", "ms", "us"])
def test_as_unit_timestamp_pre_epoch_floors(to_unit):
    # GH#63573 downcasting a pre-epoch timestamp must floor toward -inf like
    # numpy, not truncate toward zero (which would cross the epoch)
    ser_arrow = pd.Series(
        [pd.Timestamp("1969-12-31 23:59:59.123456789"), None],
        dtype="timestamp[ns][pyarrow]",
    )
    result = ser_arrow.dt.as_unit(to_unit)
    expected = pd.Series(ser_arrow.astype("datetime64[ns]").dt.as_unit(to_unit)).astype(
        f"timestamp[{to_unit}][pyarrow]"
    )
    tm.assert_series_equal(result, expected)


def test_as_unit_timestamp_overflow_raises():
    # GH#63573 upcasting out-of-bounds values must raise instead of silently
    # wrapping via int64 overflow
    ser = pd.Series(
        [pd.Timestamp("2600-01-01").as_unit("s"), None],
        dtype="timestamp[s][pyarrow]",
    )
    with pytest.raises(OutOfBoundsDatetime, match="overflow"):
        ser.dt.as_unit("ns")


def test_as_unit_duration_overflow_raises():
    # GH#63573 upcasting out-of-bounds durations must raise instead of silently
    # wrapping via int64 overflow
    ser = pd.Series([19880899200, None], dtype="duration[s][pyarrow]")
    with pytest.raises(OutOfBoundsTimedelta, match="overflow"):
        ser.dt.as_unit("ns")


@pytest.mark.parametrize(
    "prop, expected",
    [
        ["days", 1],
        ["seconds", 2],
        ["microseconds", 3],
        ["nanoseconds", 4],
    ],
)
def test_dt_timedelta_properties(prop, expected):
    # GH 52284
    ser = pd.Series(
        [
            pd.Timedelta(
                days=1,
                seconds=2,
                microseconds=3,
                nanoseconds=4,
            ),
            None,
        ],
        dtype=ArrowDtype(pa.duration("ns")),
    )
    result = getattr(ser.dt, prop)
    expected = pd.Series(
        ArrowExtensionArray(pa.array([expected, None], type=pa.int32()))
    )
    tm.assert_series_equal(result, expected)


def test_dt_timedelta_total_seconds():
    # GH 52284
    ser = pd.Series(
        [
            pd.Timedelta(
                days=1,
                seconds=2,
                microseconds=3,
                nanoseconds=4,
            ),
            None,
        ],
        dtype=ArrowDtype(pa.duration("ns")),
    )
    result = ser.dt.total_seconds()
    expected = pd.Series(
        ArrowExtensionArray(pa.array([86402.000003, None], type=pa.float64()))
    )
    tm.assert_series_equal(result, expected)


def test_dt_timedelta_components():
    # Test that .dt.components returns correct values for Arrow duration arrays
    ser = pd.Series(
        [
            pd.Timedelta(days=5, hours=3, minutes=2, seconds=1, milliseconds=4),
            pd.Timedelta(days=0, hours=23, minutes=59, seconds=59),
            None,
        ],
        dtype=ArrowDtype(pa.duration("ns")),
    )
    result = ser.dt.components
    expected = pd.DataFrame(
        {
            "days": pd.array([5, 0, None], dtype="Int32[pyarrow]"),
            "hours": pd.array([3, 23, None], dtype="Int32[pyarrow]"),
            "minutes": pd.array([2, 59, None], dtype="Int32[pyarrow]"),
            "seconds": pd.array([1, 59, None], dtype="Int32[pyarrow]"),
            "milliseconds": pd.array([4, 0, None], dtype="Int32[pyarrow]"),
            "microseconds": pd.array([0, 0, None], dtype="Int32[pyarrow]"),
            "nanoseconds": pd.array([0, 0, None], dtype="Int32[pyarrow]"),
        }
    )
    tm.assert_frame_equal(result, expected)


def test_dt_timedelta_components_preserves_index():
    # Test that .dt.components preserves the parent Series index
    idx = pd.Index(["a", "b", "c"])
    ser = pd.Series(
        [
            pd.Timedelta(days=1, hours=2),
            pd.Timedelta(days=3, hours=4),
            None,
        ],
        dtype=ArrowDtype(pa.duration("ns")),
        index=idx,
    )
    result = ser.dt.components
    tm.assert_index_equal(result.index, idx)


def test_dt_timedelta_components_negative_and_nulls():
    # Negative durations, unit boundaries, and interspersed nulls should all
    # match the NumPy-backed result. pandas normalizes negatives so that the
    # sub-day components stay non-negative (e.g. -22h57m57s -> -1 day + 1h2m3s).
    values = [
        pd.Timedelta(days=1, hours=2, minutes=3),  # ordinary positive
        pd.Timedelta(days=0),  # zero
        pd.NaT,  # null interspersed
        pd.Timedelta(hours=-22, minutes=-57, seconds=-57),  # -1 day + 1h2m3s
        pd.Timedelta(days=-2, hours=5),  # multi-day negative
        pd.Timedelta(days=-1),  # exactly -1 day, sub-day fields all 0
        pd.Timedelta(days=-100, hours=12),  # large negative
        pd.Timedelta(nanoseconds=-1),  # -1 day + 23:59:59.999999999
        pd.Timedelta(microseconds=-1),  # -1 of each unit below exercises borrow
        pd.Timedelta(milliseconds=-1),
        pd.Timedelta(seconds=-1),
        pd.Timedelta(minutes=-1),
        pd.Timedelta(hours=-1),
        pd.Timedelta(days=-1, nanoseconds=-1),  # -1 day - 1ns
        pd.NaT,
    ]
    ser_arrow = pd.Series(values, dtype=ArrowDtype(pa.duration("ns")))
    ser_numpy = pd.Series(values)

    # The direct .dt.<component> accessors should match the NumPy-backed result
    for attr in ["days", "seconds", "microseconds", "nanoseconds"]:
        expected = pd.Series(
            getattr(ser_numpy.dt, attr).values, dtype="Int32[pyarrow]"
        ).where(ser_arrow.notna(), pd.NA)
        tm.assert_series_equal(getattr(ser_arrow.dt, attr), expected)

    # .dt.components
    result = ser_arrow.dt.components
    expected = ser_numpy.dt.components.astype("Int32[pyarrow]")
    tm.assert_frame_equal(result, expected)


def test_dt_timedelta_accessors_match_python_timedelta():
    # GH 63470, GH#63283: .dt.seconds/.dt.microseconds previously returned the
    # .dt.components field values (0-59 / 0-999) instead of Python timedelta
    # semantics (total sub-day seconds, total sub-second microseconds)
    td = timedelta(
        days=1, hours=2, minutes=3, seconds=45, milliseconds=678, microseconds=123
    )
    ser = pd.Series([pd.Timedelta(td), None], dtype=ArrowDtype(pa.duration("ns")))

    # .dt.days should match Python timedelta.days
    assert ser.dt.days.iloc[0] == td.days

    # .dt.seconds should be total seconds in sub-day portion (0-86399)
    # = hours*3600 + minutes*60 + seconds = 2*3600 + 3*60 + 45 = 7425
    assert ser.dt.seconds.iloc[0] == td.seconds

    # .dt.microseconds should be total microseconds in sub-second portion (0-999999)
    # = milliseconds*1000 + microseconds = 678*1000 + 123 = 678123
    assert ser.dt.microseconds.iloc[0] == td.microseconds

    # Null handling
    assert pd.isna(ser.dt.days.iloc[1])
    assert pd.isna(ser.dt.seconds.iloc[1])
    assert pd.isna(ser.dt.microseconds.iloc[1])


@pytest.mark.parametrize("unit", ["s", "ms", "us", "ns"])
def test_dt_timedelta_components_different_units(unit):
    # Test that components work correctly across all duration units,
    # including coarser units where finer components should be zero.
    td = pd.Timedelta(days=1, hours=2, minutes=3, seconds=4, milliseconds=5)
    ser_arrow = pd.Series([td, None], dtype=ArrowDtype(pa.duration(unit)))

    result = ser_arrow.dt.components

    # Days, hours, minutes, seconds are representable at every unit
    assert result["days"].iloc[0] == 1
    assert result["hours"].iloc[0] == 2
    assert result["minutes"].iloc[0] == 3
    assert result["seconds"].iloc[0] == 4
    # Milliseconds are 0 for the coarser "s" unit, otherwise the stored value.
    # Micro/nanoseconds are 0 here regardless of unit (td has none).
    assert result["milliseconds"].iloc[0] == (5 if unit != "s" else 0)
    assert result["microseconds"].iloc[0] == 0
    assert result["nanoseconds"].iloc[0] == 0

    # Null handling across all units
    for col in result.columns:
        assert pd.isna(result[col].iloc[1])

    # GH 63470, GH#63283: the direct .dt.<component> accessors must match the
    # NumPy-backed result at every unit. In particular .dt.microseconds on a
    # coarser unit (e.g. "ms") must scale up the sub-second portion rather
    # than returning 0.
    ser_numpy = pd.Series([td.as_unit(unit), None])
    for attr in ["days", "seconds", "microseconds", "nanoseconds"]:
        expected = pd.Series(
            getattr(ser_numpy.dt, attr).values, dtype="Int32[pyarrow]"
        ).where(ser_arrow.notna(), pd.NA)
        tm.assert_series_equal(getattr(ser_arrow.dt, attr), expected)


def test_dt_timedelta_components_all_null():
    # Test all-null array hits the min_scalar.is_valid fast path
    ser = pd.Series(
        [None, None, None],
        dtype=ArrowDtype(pa.duration("ns")),
    )
    result = ser.dt.components
    for col in result.columns:
        assert result[col].isna().all()


@pytest.mark.parametrize("attr", ["days", "seconds", "microseconds", "nanoseconds"])
def test_dt_duration_components_reject_timestamp(attr):
    # GH 63470: the duration component accessors must not silently treat a
    # timestamp's underlying int64 as a duration; they should raise like the
    # NumPy-backed datetime accessors do.
    ser = pd.Series(
        pd.to_datetime(["2020-01-01 01:02:03"]), dtype="timestamp[ns][pyarrow]"
    )
    with pytest.raises(NotImplementedError, match=f"dt.{attr} is not supported"):
        getattr(ser.dt, attr)


def test_dt_to_pytimedelta():
    # GH 52284
    data = [timedelta(1, 2, 3), timedelta(1, 2, 4)]
    ser = pd.Series(data, dtype=ArrowDtype(pa.duration("ns")))

    msg = "The behavior of ArrowTemporalProperties.to_pytimedelta is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = ser.dt.to_pytimedelta()
    expected = np.array(data, dtype=object)
    tm.assert_numpy_array_equal(result, expected)
    assert all(type(res) is timedelta for res in result)

    msg = "The behavior of TimedeltaProperties.to_pytimedelta is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        expected = ser.astype("timedelta64[ns]").dt.to_pytimedelta()
    tm.assert_numpy_array_equal(result, expected)


def test_dt_components():
    # GH 52284
    ser = pd.Series(
        [
            pd.Timedelta(
                days=1,
                seconds=2,
                microseconds=3,
                nanoseconds=4,
            ),
            None,
        ],
        dtype=ArrowDtype(pa.duration("ns")),
    )
    result = ser.dt.components
    expected = pd.DataFrame(
        [[1, 0, 0, 2, 0, 3, 4], [pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA]],
        columns=[
            "days",
            "hours",
            "minutes",
            "seconds",
            "milliseconds",
            "microseconds",
            "nanoseconds",
        ],
        dtype="int32[pyarrow]",
    )
    tm.assert_frame_equal(result, expected)


def test_dt_components_large_values():
    ser = pd.Series(
        [
            pd.Timedelta("365 days 23:59:59.999000"),
            None,
        ],
        dtype=ArrowDtype(pa.duration("ns")),
    )
    result = ser.dt.components
    expected = pd.DataFrame(
        [
            [365, 23, 59, 59, 999, 0, 0],
            [pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA],
        ],
        columns=[
            "days",
            "hours",
            "minutes",
            "seconds",
            "milliseconds",
            "microseconds",
            "nanoseconds",
        ],
        dtype="int32[pyarrow]",
    )
    tm.assert_frame_equal(result, expected)


def test_dt_day_remainder_cache_invalidation():
    # Test that _dt_day_remainder cache is invalidated after __setitem__
    ser = pd.Series(
        [pd.Timedelta("1 days 1:02:03"), pd.Timedelta("2 days 4:05:06")],
        dtype=ArrowDtype(pa.duration("ns")),
    )
    arr = ser.array

    # Access to cache the day remainder and verify initial values
    result_before = pd.Series(arr._dt_seconds, dtype="int32[pyarrow]")
    expected_before = pd.Series(
        [1 * 3600 + 2 * 60 + 3, 4 * 3600 + 5 * 60 + 6], dtype="int32[pyarrow]"
    )
    tm.assert_series_equal(result_before, expected_before)
    assert "_dt_day_remainder" in arr._cache

    # Modify the array
    arr[0] = pd.Timedelta("3 days 7:08:09")

    # Cache should be invalidated
    assert "_dt_day_remainder" not in arr._cache

    # Accessing again should give correct (recomputed) values, not stale cached values
    result_after = pd.Series(arr._dt_seconds, dtype="int32[pyarrow]")
    expected_after = pd.Series(
        [7 * 3600 + 8 * 60 + 9, 4 * 3600 + 5 * 60 + 6], dtype="int32[pyarrow]"
    )
    tm.assert_series_equal(result_after, expected_after)

    # Verify the first value actually changed (not returning stale cache)
    assert result_after.iloc[0] != result_before.iloc[0]


def test_dt_day_remainder_cache_not_pickled(tmp_path):
    # The _dt_day_remainder cache can be recomputed, so it should not be
    # serialized with the array
    ser = pd.Series(
        [pd.Timedelta("1 days 1:02:03"), pd.Timedelta("2 days 4:05:06")],
        dtype=ArrowDtype(pa.duration("ns")),
    )
    arr = ser.array
    arr._dt_day_remainder  # populate the cache
    assert "_dt_day_remainder" in arr._cache

    result = tm.round_trip_pickle(arr, tmp_path / "cache.pkl")
    assert result._cache == {}
    tm.assert_extension_array_equal(result, arr)


def test_dt_day_remainder_cache_invalidated_on_sort():
    # GH 63470: sort() mutates _pa_array in place, so the _dt_day_remainder
    # cache must be cleared or the duration components would be stale
    arr = pd.array(
        [pd.Timedelta("2 days 4:05:06"), pd.Timedelta("1 days 1:02:03")],
        dtype=ArrowDtype(pa.duration("ns")),
    )
    arr._dt_day_remainder  # populate the cache
    assert "_dt_day_remainder" in arr._cache

    arr.sort()
    assert arr._cache == {}

    # Components reflect the sorted order, not the pre-sort cache
    assert arr._dt_seconds.tolist() == [3723, 14706]


def test_date32_repr():
    # GH48238
    arrow_dt = pa.array([date.fromisoformat("2020-01-01")], type=pa.date32())
    ser = pd.Series(arrow_dt, dtype=ArrowDtype(arrow_dt.type))
    assert repr(ser) == "0    2020-01-01\ndtype: date32[day][pyarrow]"


def test_duration_overflow_from_ndarray_containing_nat():
    # GH52843
    data_ts = pd.to_datetime([1, None])
    data_td = pd.to_timedelta([1, None])
    ser_ts = pd.Series(data_ts, dtype=ArrowDtype(pa.timestamp("ns")))
    ser_td = pd.Series(data_td, dtype=ArrowDtype(pa.duration("ns")))
    result = ser_ts + ser_td
    expected = pd.Series([2, None], dtype=ArrowDtype(pa.timestamp("ns")))
    tm.assert_series_equal(result, expected)


def test_timestamp_dtype_disallows_decimal():
    # GH#61773 constructing with pyarrow timestamp dtype should disallow
    #  Decimal NaN, just like pd.to_datetime
    vals = [pd.Timestamp("2016-01-02 03:04:05"), Decimal("NaN")]

    msg = "<class 'decimal.Decimal'> is not convertible to datetime"
    with pytest.raises(TypeError, match=msg):
        # Check that the non-pyarrow version raises as expected
        pd.to_datetime(vals)

    with pytest.raises(TypeError, match=msg):
        pd.array(vals, dtype=ArrowDtype(pa.timestamp("us")))


def test_timestamp_dtype_matches_to_datetime():
    # GH#61775
    dtype1 = "datetime64[ns, US/Eastern]"
    dtype2 = "timestamp[ns, US/Eastern][pyarrow]"

    ts = pd.Timestamp("2025-07-03 18:10")

    result = pd.Series([ts], dtype=dtype2)
    expected = pd.Series([ts], dtype=dtype1).convert_dtypes(dtype_backend="pyarrow")

    tm.assert_series_equal(result, expected)


def test_timestamp_vs_dt64_comparison():
    # GH#60937
    left = pd.Series(["2016-01-01"], dtype="timestamp[ns][pyarrow]")
    right = left.astype("datetime64[ns]")

    result = left == right
    expected = pd.Series([True], dtype="bool[pyarrow]")
    tm.assert_series_equal(result, expected)

    result = right == left
    tm.assert_series_equal(result, expected)


# TODO: reuse assert_invalid_comparison?
def test_date_vs_timestamp_scalar_comparison():
    # GH#62157 match non-pyarrow behavior
    ser = pd.Series(["2016-01-01"], dtype="date32[pyarrow]")
    ser2 = ser.astype("timestamp[ns][pyarrow]")

    ts = ser2[0]
    dt = ser[0]

    # date dtype don't match a Timestamp object
    assert not (ser == ts).any()
    assert not (ts == ser).any()

    # timestamp dtype doesn't match date object
    assert not (ser2 == dt).any()
    assert not (dt == ser2).any()


# TODO: reuse assert_invalid_comparison?
def test_date_vs_timestamp_array_comparison():
    # GH#62157 match non-pyarrow behavior
    # GH#60937
    ser = pd.Series(["2016-01-01"], dtype="date32[pyarrow]")
    ser2 = ser.astype("timestamp[ns][pyarrow]")
    ser3 = ser.astype("datetime64[ns]")

    assert not (ser == ser2).any()
    assert not (ser2 == ser).any()
    assert (ser != ser2).all()
    assert (ser2 != ser).all()

    assert not (ser == ser3).any()
    assert not (ser3 == ser).any()
    assert (ser != ser3).all()
    assert (ser3 != ser).all()


@pytest.mark.parametrize(
    "offset",
    [
        pd.offsets.Hour(),
        pd.offsets.Minute(),
        pd.offsets.Second(),
        pd.offsets.Milli(),
        pd.offsets.Micro(),
        pd.offsets.Nano(),
    ],
)
@pytest.mark.parametrize("dtype", ["date32[pyarrow]", "date64[pyarrow]"])
def test_date32_pyarrow_intraday_offset_raises(offset, dtype):
    ser = pd.Series([date(2022, 12, 30)], dtype=dtype)
    with pytest.raises(TypeError, match="intra-day"):
        ser + offset
    with pytest.raises(TypeError, match="intra-day"):
        ser - offset
    with pytest.raises(TypeError, match="intra-day"):
        offset + ser


@pytest.mark.parametrize(
    "offset",
    [
        pd.offsets.MonthEnd(),
        pd.offsets.MonthBegin(),
        pd.offsets.Day(5),
        pd.DateOffset(years=1),
    ],
)
@pytest.mark.parametrize("dtype", ["date32[pyarrow]", "date64[pyarrow]"])
def test_date32_pyarrow_dateoffset_add(offset, dtype):
    ser = pd.Series([date(2022, 12, 30)], dtype=dtype)

    result = ser + offset
    expected = offset + date(2022, 12, 30)
    if isinstance(expected, pd.Timestamp):
        expected = expected.date()
    assert result[0] == expected

    result = ser - offset
    expected = date(2022, 12, 30) - offset
    if isinstance(expected, pd.Timestamp):
        expected = expected.date()
    assert result[0] == expected

    result = offset + ser
    expected = offset + date(2022, 12, 30)
    if isinstance(expected, pd.Timestamp):
        expected = expected.date()
    assert result[0] == expected


@pytest.mark.parametrize("dtype", ["date32[pyarrow]", "date64[pyarrow]"])
def test_date32_pyarrow_dateoffset_with_nulls(dtype):
    ser = pd.Series([date(2022, 12, 30), None], dtype=dtype)
    result = ser + pd.offsets.MonthEnd()
    assert result[0] == date(2022, 12, 31)
    assert pd.isna(result[1])  # handles NA, NaT, None uniformly


@pytest.mark.parametrize("method", ["sum", "min", "max", "mean", "median"])
def test_duration_reduction_consistency(unit, method):
    # GH#63170
    dtype = f"duration[{unit}][pyarrow]"
    ser = pd.Series([timedelta(seconds=1), timedelta(seconds=2)], dtype=dtype)
    result = getattr(ser, method)()
    assert isinstance(result, pd.Timedelta), (
        f"{method} for {unit} returned {type(result)}"
    )
    assert result.unit == unit


@pytest.mark.parametrize("method", ["min", "max", "median"])
def test_timestamp_reduction_consistency(unit, method):
    # GH#63170
    dtype = f"timestamp[{unit}][pyarrow]"
    ser = pd.Series([datetime(2024, 1, 1), datetime(2024, 1, 3)], dtype=dtype)
    result = getattr(ser, method)()
    assert isinstance(result, pd.Timestamp), (
        f"{method} for {unit} returned {type(result)}"
    )
    assert result.unit == unit


@pytest.mark.parametrize(
    "target_tz, value_tz", [(None, "US/Eastern"), ("US/Eastern", None)]
)
@pytest.mark.parametrize("method", ["where", "putmask"])
def test_index_tz_mismatch_casts_to_object(target_tz, value_tz, method):
    # GH#69029 Index catches the TypeError and widens, as DatetimeIndex does
    idx = pd.Index(pd.array([1, 2], dtype=ArrowDtype(pa.timestamp("ns", tz=target_tz))))
    value = pd.Timestamp("2016-01-01", tz=value_tz)
    result = getattr(idx, method)(
        [False, True] if method == "where" else [True, False], value
    )
    assert result.dtype == object
    assert result[0] == value
    # the entry the mask did not select keeps its timezone; without the
    #  _get_common_dtype fix the whole index is relabelled tz-naive
    assert result[1] == idx[1]


@pytest.mark.parametrize(
    "pa_type", tm.DATETIME_PYARROW_DTYPES + tm.TIMEDELTA_PYARROW_DTYPES
)
def test_quantile_temporal(pa_type):
    # GH52678
    data = [1, 2, 3]
    ser = pd.Series(data, dtype=ArrowDtype(pa_type))
    result = ser.quantile(0.1)
    expected = ser[0]
    assert result == expected


@pytest.mark.parametrize(
    "pa_type", tm.DATETIME_PYARROW_DTYPES + tm.TIMEDELTA_PYARROW_DTYPES
)
def test_from_sequence_temporal(pa_type):
    # GH 53171
    val = 3
    unit = pa_type.unit
    if pa.types.is_duration(pa_type):
        seq = [pd.Timedelta(val, unit=unit).as_unit(unit)]
    else:
        seq = [pd.Timestamp(val, unit=unit, tz=pa_type.tz).as_unit(unit)]

    result = ArrowExtensionArray._from_sequence(seq, dtype=pa_type)
    expected = ArrowExtensionArray(pa.array([val], type=pa_type))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "pa_type", tm.DATETIME_PYARROW_DTYPES + tm.TIMEDELTA_PYARROW_DTYPES
)
def test_setitem_temporal(pa_type):
    # GH 53171
    unit = pa_type.unit
    if pa.types.is_duration(pa_type):
        val = pd.Timedelta(1, unit=unit).as_unit(unit)
    else:
        val = pd.Timestamp(1, unit=unit, tz=pa_type.tz).as_unit(unit)

    arr = ArrowExtensionArray(pa.array([1, 2, 3], type=pa_type))

    result = arr.copy()
    result[:] = val
    expected = ArrowExtensionArray(pa.array([1, 1, 1], type=pa_type))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "pa_type", tm.DATETIME_PYARROW_DTYPES + tm.TIMEDELTA_PYARROW_DTYPES
)
def test_arithmetic_temporal(pa_type, request):
    # GH 53171
    arr = ArrowExtensionArray(pa.array([1, 2, 3], type=pa_type))
    unit = pa_type.unit
    result = arr - pd.Timedelta(1, unit=unit).as_unit(unit)
    expected = ArrowExtensionArray(pa.array([0, 1, 2], type=pa_type))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "pa_type", tm.DATETIME_PYARROW_DTYPES + tm.TIMEDELTA_PYARROW_DTYPES
)
def test_comparison_temporal(pa_type):
    # GH 53171
    unit = pa_type.unit
    if pa.types.is_duration(pa_type):
        val = pd.Timedelta(1, unit=unit).as_unit(unit)
    else:
        val = pd.Timestamp(1, unit=unit, tz=pa_type.tz).as_unit(unit)

    arr = ArrowExtensionArray(pa.array([1, 2, 3], type=pa_type))

    result = arr > val
    expected = ArrowExtensionArray(pa.array([False, True, True], type=pa.bool_()))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "pa_type", tm.DATETIME_PYARROW_DTYPES + tm.TIMEDELTA_PYARROW_DTYPES
)
def test_getitem_temporal(pa_type):
    # GH 53326
    arr = ArrowExtensionArray(pa.array([1, 2, 3], type=pa_type))
    result = arr[1]
    if pa.types.is_duration(pa_type):
        expected = pd.Timedelta(2, unit=pa_type.unit).as_unit(pa_type.unit)
        assert isinstance(result, pd.Timedelta)
    else:
        expected = pd.Timestamp(2, unit=pa_type.unit, tz=pa_type.tz).as_unit(
            pa_type.unit
        )
        assert isinstance(result, pd.Timestamp)
    assert result.unit == expected.unit
    assert result == expected


@pytest.mark.parametrize(
    "pa_type", tm.DATETIME_PYARROW_DTYPES + tm.TIMEDELTA_PYARROW_DTYPES
)
def test_iter_temporal(pa_type):
    # GH 53326
    arr = ArrowExtensionArray(pa.array([1, None], type=pa_type))
    result = list(arr)
    if pa.types.is_duration(pa_type):
        expected = [
            pd.Timedelta(1, unit=pa_type.unit).as_unit(pa_type.unit),
            pd.NA,
        ]
        assert isinstance(result[0], pd.Timedelta)
    else:
        expected = [
            pd.Timestamp(1, unit=pa_type.unit, tz=pa_type.tz).as_unit(pa_type.unit),
            pd.NA,
        ]
        assert isinstance(result[0], pd.Timestamp)
    assert result[0].unit == expected[0].unit
    assert result == expected


@pytest.mark.parametrize(
    "pa_type", tm.DATETIME_PYARROW_DTYPES + tm.TIMEDELTA_PYARROW_DTYPES, ids=repr
)
@pytest.mark.parametrize("dtype", [None, object])
def test_to_numpy_temporal(pa_type, dtype):
    # GH 53326
    # GH 55997: Return datetime64/timedelta64 types with NaT if possible
    arr = ArrowExtensionArray(pa.array([1, None], type=pa_type))
    result = arr.to_numpy(dtype=dtype)
    if pa.types.is_duration(pa_type):
        value = pd.Timedelta(1, unit=pa_type.unit).as_unit(pa_type.unit)
    else:
        value = pd.Timestamp(1, unit=pa_type.unit, tz=pa_type.tz).as_unit(pa_type.unit)

    if dtype == object or (pa.types.is_timestamp(pa_type) and pa_type.tz is not None):
        if dtype == object:
            na = pd.NA
        else:
            na = pd.NaT
        expected = np.array([value, na], dtype=object)
        assert result[0].unit == value.unit
    else:
        na = pa_type.to_pandas_dtype().type("nat", pa_type.unit)
        value = value.to_numpy()
        expected = np.array([value, na])
        assert np.datetime_data(result[0])[0] == pa_type.unit
    tm.assert_numpy_array_equal(result, expected)


def test_string_to_datetime_parsing_cast():
    # GH 56266
    string_dates = ["2020-01-01 04:30:00", "2020-01-02 00:00:00", "2020-01-03 00:00:00"]
    result = pd.Series(string_dates, dtype="timestamp[s][pyarrow]")

    pd_res = pd.to_datetime(string_dates).as_unit("s")
    expected = pd.Series(ArrowExtensionArray(pa.array(pd_res, from_pandas=True)))
    tm.assert_series_equal(result, expected)


def test_string_to_time_parsing_cast():
    # GH 56463
    string_times = ["11:41:43.076160"]
    result = pd.Series(string_times, dtype="time64[us][pyarrow]")
    expected = pd.Series(
        ArrowExtensionArray(pa.array([time(11, 41, 43, 76160)], from_pandas=True))
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("dtype", ["time32[s][pyarrow]", "time64[us][pyarrow]"])
def test_string_to_time_parsing_cast_meridiem(dtype):
    with tm.set_locale("C", locale.LC_TIME):
        # GH#18793 the space before AM/PM used to make these coerce to null
        result = pd.Series(["3:25:00 PM"], dtype=dtype)
        expected = pd.Series(["15:25:00"], dtype=dtype)
        tm.assert_series_equal(result, expected)
