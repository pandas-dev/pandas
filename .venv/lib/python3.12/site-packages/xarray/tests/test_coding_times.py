from __future__ import annotations

import warnings
from datetime import datetime, timedelta
from itertools import product, starmap
from typing import Literal

import numpy as np
import pandas as pd
import pytest
from pandas.errors import OutOfBoundsDatetime, OutOfBoundsTimedelta

from xarray import (
    DataArray,
    Dataset,
    Variable,
    conventions,
    date_range,
    decode_cf,
)
from xarray.coders import CFDatetimeCoder, CFTimedeltaCoder
from xarray.coding.times import (
    _encode_datetime_with_cftime,
    _netcdf_to_numpy_timeunit,
    _numpy_to_netcdf_timeunit,
    _should_cftime_be_used,
    cftime_to_nptime,
    decode_cf_datetime,
    decode_cf_timedelta,
    encode_cf_datetime,
    encode_cf_timedelta,
    format_cftime_datetime,
    infer_datetime_units,
    infer_timedelta_units,
)
from xarray.coding.variables import SerializationWarning
from xarray.conventions import _update_bounds_attributes, cf_encoder
from xarray.core.common import contains_cftime_datetimes
from xarray.core.types import PDDatetimeUnitOptions
from xarray.core.utils import is_duck_dask_array
from xarray.testing import assert_equal, assert_identical
from xarray.tests import (
    _ALL_CALENDARS,
    _NON_STANDARD_CALENDARS,
    _STANDARD_CALENDAR_NAMES,
    _STANDARD_CALENDARS,
    DuckArrayWrapper,
    FirstElementAccessibleArray,
    _all_cftime_date_types,
    arm_xfail,
    assert_array_equal,
    assert_duckarray_allclose,
    assert_duckarray_equal,
    assert_no_warnings,
    has_cftime,
    requires_cftime,
    requires_dask,
)

_CF_DATETIME_NUM_DATES_UNITS = [
    (np.arange(10), "days since 2000-01-01", "s"),
    (np.arange(10).astype("float64"), "days since 2000-01-01", "s"),
    (np.arange(10).astype("float32"), "days since 2000-01-01", "s"),
    (np.arange(10).reshape(2, 5), "days since 2000-01-01", "s"),
    (12300 + np.arange(5), "hours since 1680-01-01 00:00:00", "s"),
    # here we add a couple minor formatting errors to test
    # the robustness of the parsing algorithm.
    (12300 + np.arange(5), "hour since 1680-01-01  00:00:00", "s"),
    (12300 + np.arange(5), "Hour  since 1680-01-01 00:00:00", "s"),
    (12300 + np.arange(5), " Hour  since  1680-01-01 00:00:00 ", "s"),
    (10, "days since 2000-01-01", "s"),
    ([10], "daYs  since 2000-01-01", "s"),
    ([[10]], "days since 2000-01-01", "s"),
    ([10, 10], "days since 2000-01-01", "s"),
    (np.array(10), "days since 2000-01-01", "s"),
    (0, "days since 1000-01-01", "s"),
    ([0], "days since 1000-01-01", "s"),
    ([[0]], "days since 1000-01-01", "s"),
    (np.arange(2), "days since 1000-01-01", "s"),
    (np.arange(0, 100000, 20000), "days since 1900-01-01", "s"),
    (np.arange(0, 100000, 20000), "days since 1-01-01", "s"),
    (17093352.0, "hours since 1-1-1 00:00:0.0", "s"),
    ([0.5, 1.5], "hours since 1900-01-01T00:00:00", "s"),
    (0, "milliseconds since 2000-01-01T00:00:00", "s"),
    (0, "microseconds since 2000-01-01T00:00:00", "s"),
    (np.int32(788961600), "seconds since 1981-01-01", "s"),  # GH2002
    (12300 + np.arange(5), "hour since 1680-01-01 00:00:00.500000", "us"),
    (164375, "days since 1850-01-01 00:00:00", "s"),
    (164374.5, "days since 1850-01-01 00:00:00", "s"),
    ([164374.5, 168360.5], "days since 1850-01-01 00:00:00", "s"),
]
_CF_DATETIME_TESTS = [
    num_dates_units + (calendar,)
    for num_dates_units, calendar in product(
        _CF_DATETIME_NUM_DATES_UNITS, _STANDARD_CALENDAR_NAMES
    )
]


@requires_cftime
@pytest.mark.filterwarnings("ignore:Ambiguous reference date string")
@pytest.mark.filterwarnings("ignore:Times can't be serialized faithfully")
@pytest.mark.parametrize(
    ["num_dates", "units", "minimum_resolution", "calendar"], _CF_DATETIME_TESTS
)
def test_cf_datetime(
    num_dates,
    units: str,
    minimum_resolution: PDDatetimeUnitOptions,
    calendar: str,
    time_unit: PDDatetimeUnitOptions,
) -> None:
    import cftime

    expected = cftime.num2date(
        num_dates, units, calendar, only_use_cftime_datetimes=True
    )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Unable to decode time axis")
        actual = decode_cf_datetime(num_dates, units, calendar, time_unit=time_unit)

    if actual.dtype.kind != "O":
        if np.timedelta64(1, time_unit) > np.timedelta64(1, minimum_resolution):
            expected_unit = minimum_resolution
        else:
            expected_unit = time_unit
        expected = cftime_to_nptime(expected, time_unit=expected_unit)

    assert_array_equal(actual, expected)
    encoded1, _, _ = encode_cf_datetime(actual, units, calendar)

    assert_array_equal(num_dates, encoded1)

    if hasattr(num_dates, "ndim") and num_dates.ndim == 1 and "1000" not in units:
        # verify that wrapping with a pandas.Index works
        # note that it *does not* currently work to put
        # non-datetime64 compatible dates into a pandas.Index
        encoded2, _, _ = encode_cf_datetime(pd.Index(actual), units, calendar)
        assert_array_equal(num_dates, encoded2)


@requires_cftime
def test_decode_cf_datetime_overflow(time_unit: PDDatetimeUnitOptions) -> None:
    # checks for
    # https://github.com/pydata/pandas/issues/14068
    # https://github.com/pydata/xarray/issues/975
    from cftime import DatetimeGregorian

    datetime = DatetimeGregorian
    units = "days since 2000-01-01 00:00:00"

    # date after 2262 and before 1678
    days = (-117710, 95795)
    expected = (datetime(1677, 9, 20), datetime(2262, 4, 12))
    for i, day in enumerate(days):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Unable to decode time axis")
            result = decode_cf_datetime(
                day, units, calendar="standard", time_unit=time_unit
            )
        assert result == expected[i]
        # additional check to see if type/dtypes are correct
        if time_unit == "ns":
            assert isinstance(result.item(), datetime)
        else:
            assert result.dtype == np.dtype(f"=M8[{time_unit}]")


def test_decode_cf_datetime_non_standard_units() -> None:
    expected = pd.date_range(periods=100, start="1970-01-01", freq="h")
    # netCDFs from madis.noaa.gov use this format for their time units
    # they cannot be parsed by cftime, but pd.Timestamp works
    units = "hours since 1-1-1970"
    actual = decode_cf_datetime(np.arange(100), units)
    assert_array_equal(actual, expected)


@requires_cftime
def test_decode_cf_datetime_non_iso_strings() -> None:
    # datetime strings that are _almost_ ISO compliant but not quite,
    # but which cftime.num2date can still parse correctly
    expected = pd.date_range(periods=100, start="2000-01-01", freq="h")
    cases = [
        (np.arange(100), "hours since 2000-01-01 0"),
        (np.arange(100), "hours since 2000-1-1 0"),
        (np.arange(100), "hours since 2000-01-01 0:00"),
    ]
    for num_dates, units in cases:
        actual = decode_cf_datetime(num_dates, units)
        assert_array_equal(actual, expected)


@requires_cftime
@pytest.mark.parametrize("calendar", _STANDARD_CALENDARS)
def test_decode_standard_calendar_inside_timestamp_range(
    calendar, time_unit: PDDatetimeUnitOptions
) -> None:
    import cftime

    units = "hours since 0001-01-01"
    times = pd.date_range(
        "2001-04-01-00", end="2001-04-30-23", unit=time_unit, freq="h"
    )
    # to_pydatetime() will return microsecond
    time = cftime.date2num(times.to_pydatetime(), units, calendar=calendar)
    expected = times.values
    # for cftime we get "us" resolution
    # ns resolution is handled by cftime due to the reference date
    # being out of bounds, but the times themselves are
    # representable with nanosecond resolution.
    actual = decode_cf_datetime(time, units, calendar=calendar, time_unit=time_unit)
    assert actual.dtype == np.dtype(f"=M8[{time_unit}]")
    assert_array_equal(actual, expected)


@requires_cftime
@pytest.mark.parametrize("calendar", _NON_STANDARD_CALENDARS)
def test_decode_non_standard_calendar_inside_timestamp_range(calendar) -> None:
    import cftime

    units = "days since 0001-01-01"
    times = pd.date_range("2001-04-01-00", end="2001-04-30-23", freq="h")
    non_standard_time = cftime.date2num(times.to_pydatetime(), units, calendar=calendar)

    expected = cftime.num2date(
        non_standard_time, units, calendar=calendar, only_use_cftime_datetimes=True
    )
    expected_dtype = np.dtype("O")

    actual = decode_cf_datetime(non_standard_time, units, calendar=calendar)
    assert actual.dtype == expected_dtype
    assert_array_equal(actual, expected)


@requires_cftime
@pytest.mark.parametrize("calendar", _ALL_CALENDARS)
def test_decode_dates_outside_timestamp_range(
    calendar, time_unit: PDDatetimeUnitOptions
) -> None:
    import cftime

    units = "days since 0001-01-01"
    times = [datetime(1, 4, 1, h) for h in range(1, 5)]
    time = cftime.date2num(times, units, calendar=calendar)

    expected = cftime.num2date(
        time, units, calendar=calendar, only_use_cftime_datetimes=True
    )
    if calendar == "proleptic_gregorian" and time_unit != "ns":
        expected = cftime_to_nptime(expected, time_unit=time_unit)
    expected_date_type = type(expected[0])

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Unable to decode time axis")
        actual = decode_cf_datetime(time, units, calendar=calendar, time_unit=time_unit)
    assert all(isinstance(value, expected_date_type) for value in actual)
    assert_array_equal(actual, expected)


@requires_cftime
@pytest.mark.parametrize("calendar", _STANDARD_CALENDARS)
@pytest.mark.parametrize("num_time", [735368, [735368], [[735368]]])
def test_decode_standard_calendar_single_element_inside_timestamp_range(
    calendar,
    time_unit: PDDatetimeUnitOptions,
    num_time,
) -> None:
    units = "days since 0001-01-01"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Unable to decode time axis")
        actual = decode_cf_datetime(
            num_time, units, calendar=calendar, time_unit=time_unit
        )

    assert actual.dtype == np.dtype(f"=M8[{time_unit}]")


@requires_cftime
@pytest.mark.parametrize("calendar", _NON_STANDARD_CALENDARS)
def test_decode_non_standard_calendar_single_element_inside_timestamp_range(
    calendar,
) -> None:
    units = "days since 0001-01-01"
    for num_time in [735368, [735368], [[735368]]]:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Unable to decode time axis")
            actual = decode_cf_datetime(num_time, units, calendar=calendar)
        assert actual.dtype == np.dtype("O")


@requires_cftime
@pytest.mark.parametrize("calendar", _NON_STANDARD_CALENDARS)
def test_decode_single_element_outside_timestamp_range(calendar) -> None:
    import cftime

    units = "days since 0001-01-01"
    for days in [1, 1470376]:
        for num_time in [days, [days], [[days]]]:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", "Unable to decode time axis")
                actual = decode_cf_datetime(num_time, units, calendar=calendar)

            expected = cftime.num2date(
                days, units, calendar, only_use_cftime_datetimes=True
            )
            assert isinstance(actual.item(), type(expected))


@requires_cftime
@pytest.mark.parametrize("calendar", _STANDARD_CALENDARS)
def test_decode_standard_calendar_multidim_time_inside_timestamp_range(
    calendar,
    time_unit: PDDatetimeUnitOptions,
) -> None:
    import cftime

    units = "days since 0001-01-01"
    times1 = pd.date_range("2001-04-01", end="2001-04-05", freq="D")
    times2 = pd.date_range("2001-05-01", end="2001-05-05", freq="D")
    time1 = cftime.date2num(times1.to_pydatetime(), units, calendar=calendar)
    time2 = cftime.date2num(times2.to_pydatetime(), units, calendar=calendar)
    mdim_time = np.empty((len(time1), 2))
    mdim_time[:, 0] = time1
    mdim_time[:, 1] = time2

    expected1 = times1.values
    expected2 = times2.values

    actual = decode_cf_datetime(
        mdim_time, units, calendar=calendar, time_unit=time_unit
    )
    assert actual.dtype == np.dtype(f"=M8[{time_unit}]")
    assert_array_equal(actual[:, 0], expected1)
    assert_array_equal(actual[:, 1], expected2)


@requires_cftime
@pytest.mark.parametrize("calendar", _NON_STANDARD_CALENDARS)
def test_decode_nonstandard_calendar_multidim_time_inside_timestamp_range(
    calendar,
) -> None:
    import cftime

    units = "days since 0001-01-01"
    times1 = pd.date_range("2001-04-01", end="2001-04-05", freq="D")
    times2 = pd.date_range("2001-05-01", end="2001-05-05", freq="D")
    time1 = cftime.date2num(times1.to_pydatetime(), units, calendar=calendar)
    time2 = cftime.date2num(times2.to_pydatetime(), units, calendar=calendar)
    mdim_time = np.empty((len(time1), 2))
    mdim_time[:, 0] = time1
    mdim_time[:, 1] = time2

    if cftime.__name__ == "cftime":
        expected1 = cftime.num2date(
            time1, units, calendar, only_use_cftime_datetimes=True
        )
        expected2 = cftime.num2date(
            time2, units, calendar, only_use_cftime_datetimes=True
        )
    else:
        expected1 = cftime.num2date(time1, units, calendar)
        expected2 = cftime.num2date(time2, units, calendar)

    expected_dtype = np.dtype("O")

    actual = decode_cf_datetime(mdim_time, units, calendar=calendar)

    assert actual.dtype == expected_dtype
    assert_array_equal(actual[:, 0], expected1)
    assert_array_equal(actual[:, 1], expected2)


@requires_cftime
@pytest.mark.parametrize("calendar", _ALL_CALENDARS)
def test_decode_multidim_time_outside_timestamp_range(
    calendar, time_unit: PDDatetimeUnitOptions
) -> None:
    import cftime

    units = "days since 0001-01-01"
    times1 = [datetime(1, 4, day) for day in range(1, 6)]
    times2 = [datetime(1, 5, day) for day in range(1, 6)]
    time1 = cftime.date2num(times1, units, calendar=calendar)
    time2 = cftime.date2num(times2, units, calendar=calendar)
    mdim_time = np.empty((len(time1), 2))
    mdim_time[:, 0] = time1
    mdim_time[:, 1] = time2

    expected1 = cftime.num2date(time1, units, calendar, only_use_cftime_datetimes=True)
    expected2 = cftime.num2date(time2, units, calendar, only_use_cftime_datetimes=True)

    if calendar == "proleptic_gregorian" and time_unit != "ns":
        expected1 = cftime_to_nptime(expected1, time_unit=time_unit)
        expected2 = cftime_to_nptime(expected2, time_unit=time_unit)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Unable to decode time axis")
        actual = decode_cf_datetime(
            mdim_time, units, calendar=calendar, time_unit=time_unit
        )

    dtype: np.dtype
    dtype = np.dtype("O")
    if calendar == "proleptic_gregorian" and time_unit != "ns":
        dtype = np.dtype(f"=M8[{time_unit}]")

    assert actual.dtype == dtype
    assert_array_equal(actual[:, 0], expected1)
    assert_array_equal(actual[:, 1], expected2)


@requires_cftime
@pytest.mark.parametrize(
    ("calendar", "num_time"),
    [("360_day", 720058.0), ("all_leap", 732059.0), ("366_day", 732059.0)],
)
def test_decode_non_standard_calendar_single_element(calendar, num_time) -> None:
    import cftime

    units = "days since 0001-01-01"

    actual = decode_cf_datetime(num_time, units, calendar=calendar)

    expected = np.asarray(
        cftime.num2date(num_time, units, calendar, only_use_cftime_datetimes=True)
    )
    assert actual.dtype == np.dtype("O")
    assert expected == actual


@requires_cftime
def test_decode_360_day_calendar() -> None:
    import cftime

    calendar = "360_day"
    # ensure leap year doesn't matter
    for year in [2010, 2011, 2012, 2013, 2014]:
        units = f"days since {year}-01-01"
        num_times = np.arange(100)

        expected = cftime.num2date(
            num_times, units, calendar, only_use_cftime_datetimes=True
        )

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            actual = decode_cf_datetime(num_times, units, calendar=calendar)
            assert len(w) == 0

        assert actual.dtype == np.dtype("O")
        assert_array_equal(actual, expected)


@requires_cftime
def test_decode_abbreviation() -> None:
    """Test making sure we properly fall back to cftime on abbreviated units."""
    import cftime

    val = np.array([1586628000000.0])
    units = "msecs since 1970-01-01T00:00:00Z"
    actual = decode_cf_datetime(val, units)
    expected = cftime_to_nptime(cftime.num2date(val, units))
    assert_array_equal(actual, expected)


@arm_xfail
@requires_cftime
@pytest.mark.parametrize(
    ["num_dates", "units", "expected_list"],
    [
        ([np.nan], "days since 2000-01-01", ["NaT"]),
        ([np.nan, 0], "days since 2000-01-01", ["NaT", "2000-01-01T00:00:00Z"]),
        (
            [np.nan, 0, 1],
            "days since 2000-01-01",
            ["NaT", "2000-01-01T00:00:00Z", "2000-01-02T00:00:00Z"],
        ),
    ],
)
def test_cf_datetime_nan(num_dates, units, expected_list) -> None:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "All-NaN")
        actual = decode_cf_datetime(num_dates, units)
    # use pandas because numpy will deprecate timezone-aware conversions
    expected = pd.to_datetime(expected_list).to_numpy(dtype="datetime64[ns]")
    assert_array_equal(expected, actual)


@requires_cftime
def test_decoded_cf_datetime_array_2d(time_unit: PDDatetimeUnitOptions) -> None:
    # regression test for GH1229
    variable = Variable(
        ("x", "y"), np.array([[0, 1], [2, 3]]), {"units": "days since 2000-01-01"}
    )
    result = CFDatetimeCoder(time_unit=time_unit).decode(variable)
    assert result.dtype == f"datetime64[{time_unit}]"
    expected = pd.date_range("2000-01-01", periods=4).values.reshape(2, 2)
    assert_array_equal(np.asarray(result), expected)


@pytest.mark.parametrize("decode_times", [True, False])
@pytest.mark.parametrize("mask_and_scale", [True, False])
def test_decode_datetime_mask_and_scale(
    decode_times: bool, mask_and_scale: bool
) -> None:
    attrs = {
        "units": "nanoseconds since 1970-01-01",
        "calendar": "proleptic_gregorian",
        "_FillValue": np.int16(-1),
        "add_offset": 100000.0,
    }
    encoded = Variable(["time"], np.array([0, -1, 1], "int16"), attrs=attrs)
    decoded = conventions.decode_cf_variable(
        "foo", encoded, mask_and_scale=mask_and_scale, decode_times=decode_times
    )
    result = conventions.encode_cf_variable(decoded, name="foo")
    assert_identical(encoded, result)
    assert encoded.dtype == result.dtype


FREQUENCIES_TO_ENCODING_UNITS = {
    "ns": "nanoseconds",
    "us": "microseconds",
    "ms": "milliseconds",
    "s": "seconds",
    "min": "minutes",
    "h": "hours",
    "D": "days",
}


@pytest.mark.parametrize(("freq", "units"), FREQUENCIES_TO_ENCODING_UNITS.items())
def test_infer_datetime_units(freq, units) -> None:
    dates = pd.date_range("2000", periods=2, freq=freq)
    expected = f"{units} since 2000-01-01 00:00:00"
    assert expected == infer_datetime_units(dates)


@pytest.mark.parametrize(
    ["dates", "expected"],
    [
        (
            pd.to_datetime(["1900-01-01", "1900-01-02", "NaT"], unit="ns"),
            "days since 1900-01-01 00:00:00",
        ),
        (
            pd.to_datetime(["NaT", "1900-01-01"], unit="ns"),
            "days since 1900-01-01 00:00:00",
        ),
        (pd.to_datetime(["NaT"], unit="ns"), "days since 1970-01-01 00:00:00"),
    ],
)
def test_infer_datetime_units_with_NaT(dates, expected) -> None:
    assert expected == infer_datetime_units(dates)


_CFTIME_DATETIME_UNITS_TESTS = [
    ([(1900, 1, 1), (1900, 1, 1)], "days since 1900-01-01 00:00:00.000000"),
    (
        [(1900, 1, 1), (1900, 1, 2), (1900, 1, 2, 0, 0, 1)],
        "seconds since 1900-01-01 00:00:00.000000",
    ),
    (
        [(1900, 1, 1), (1900, 1, 8), (1900, 1, 16)],
        "days since 1900-01-01 00:00:00.000000",
    ),
]


@requires_cftime
@pytest.mark.parametrize(
    "calendar", _NON_STANDARD_CALENDARS + ["gregorian", "proleptic_gregorian"]
)
@pytest.mark.parametrize(("date_args", "expected"), _CFTIME_DATETIME_UNITS_TESTS)
def test_infer_cftime_datetime_units(calendar, date_args, expected) -> None:
    date_type = _all_cftime_date_types()[calendar]
    dates = list(starmap(date_type, date_args))
    assert expected == infer_datetime_units(dates)


@pytest.mark.filterwarnings("ignore:Timedeltas can't be serialized faithfully")
@pytest.mark.parametrize(
    ["timedeltas", "units", "numbers"],
    [
        ("1D", "days", np.int64(1)),
        (["1D", "2D", "3D"], "days", np.array([1, 2, 3], "int64")),
        ("1h", "hours", np.int64(1)),
        ("1ms", "milliseconds", np.int64(1)),
        ("1us", "microseconds", np.int64(1)),
        ("1ns", "nanoseconds", np.int64(1)),
        (["NaT", "0s", "1s"], None, [np.iinfo(np.int64).min, 0, 1]),
        (["30m", "60m"], "hours", [0.5, 1.0]),
        ("NaT", "days", np.iinfo(np.int64).min),
        (["NaT", "NaT"], "days", [np.iinfo(np.int64).min, np.iinfo(np.int64).min]),
    ],
)
def test_cf_timedelta(timedeltas, units, numbers) -> None:
    if timedeltas == "NaT":
        timedeltas = np.timedelta64("NaT", "ns")
    else:
        timedeltas = pd.to_timedelta(timedeltas).as_unit("ns").to_numpy()
    numbers = np.array(numbers)

    expected = numbers
    actual, _ = encode_cf_timedelta(timedeltas, units)
    assert_array_equal(expected, actual)
    assert expected.dtype == actual.dtype

    if units is not None:
        expected = timedeltas
        actual = decode_cf_timedelta(numbers, units)
        assert_array_equal(expected, actual)
        assert expected.dtype == actual.dtype

    expected = np.timedelta64("NaT", "ns")
    actual = decode_cf_timedelta(np.array(np.nan), "days")
    assert_array_equal(expected, actual)
    assert expected.dtype == actual.dtype


def test_cf_timedelta_2d() -> None:
    units = "days"
    numbers = np.atleast_2d([1, 2, 3])

    timedeltas = pd.to_timedelta(["1D", "2D", "3D"]).as_unit("ns")
    timedeltas_2d = np.atleast_2d(timedeltas.to_numpy())
    expected = timedeltas_2d

    actual = decode_cf_timedelta(numbers, units)
    assert_array_equal(expected, actual)
    assert expected.dtype == actual.dtype


@pytest.mark.parametrize("encoding_unit", FREQUENCIES_TO_ENCODING_UNITS.values())
def test_decode_cf_timedelta_time_unit(
    time_unit: PDDatetimeUnitOptions, encoding_unit
) -> None:
    encoded = 1
    encoding_unit_as_numpy = _netcdf_to_numpy_timeunit(encoding_unit)
    if np.timedelta64(1, time_unit) > np.timedelta64(1, encoding_unit_as_numpy):
        expected = np.timedelta64(encoded, encoding_unit_as_numpy)
    else:
        expected = np.timedelta64(encoded, encoding_unit_as_numpy).astype(
            f"timedelta64[{time_unit}]"
        )
    result = decode_cf_timedelta(encoded, encoding_unit, time_unit)
    assert result == expected
    assert result.dtype == expected.dtype


def test_decode_cf_timedelta_time_unit_out_of_bounds(
    time_unit: PDDatetimeUnitOptions,
) -> None:
    # Define a scale factor that will guarantee overflow with the given
    # time_unit.
    scale_factor = np.timedelta64(1, time_unit) // np.timedelta64(1, "ns")
    encoded = scale_factor * 300 * 365
    with pytest.raises(OutOfBoundsTimedelta):
        decode_cf_timedelta(encoded, "days", time_unit)


def test_cf_timedelta_roundtrip_large_value(time_unit: PDDatetimeUnitOptions) -> None:
    value = np.timedelta64(np.iinfo(np.int64).max, time_unit)
    encoded, units = encode_cf_timedelta(value)
    decoded = decode_cf_timedelta(encoded, units, time_unit=time_unit)
    assert value == decoded
    assert value.dtype == decoded.dtype


@pytest.mark.parametrize(
    ["deltas", "expected"],
    [
        (pd.to_timedelta(["1 day", "2 days"]), "days"),
        (pd.to_timedelta(["1h", "1 day 1 hour"]), "hours"),
        (pd.to_timedelta(["1m", "2m", np.nan]), "minutes"),
        (pd.to_timedelta(["1m3s", "1m4s"]), "seconds"),
    ],
)
def test_infer_timedelta_units(deltas, expected) -> None:
    assert expected == infer_timedelta_units(deltas)


@requires_cftime
@pytest.mark.parametrize(
    ["date_args", "expected"],
    [
        ((1, 2, 3, 4, 5, 6), "0001-02-03 04:05:06.000000"),
        ((10, 2, 3, 4, 5, 6), "0010-02-03 04:05:06.000000"),
        ((100, 2, 3, 4, 5, 6), "0100-02-03 04:05:06.000000"),
        ((1000, 2, 3, 4, 5, 6), "1000-02-03 04:05:06.000000"),
    ],
)
def test_format_cftime_datetime(date_args, expected) -> None:
    date_types = _all_cftime_date_types()
    for date_type in date_types.values():
        result = format_cftime_datetime(date_type(*date_args))
        assert result == expected


@pytest.mark.parametrize("calendar", _ALL_CALENDARS)
def test_decode_cf(calendar, time_unit: PDDatetimeUnitOptions) -> None:
    days = [1.0, 2.0, 3.0]
    # TODO: GH5690 — do we want to allow this type for `coords`?
    da = DataArray(days, coords=[days], dims=["time"], name="test")
    ds = da.to_dataset()

    for v in ["test", "time"]:
        ds[v].attrs["units"] = "days since 2001-01-01"
        ds[v].attrs["calendar"] = calendar

    if not has_cftime and calendar not in _STANDARD_CALENDAR_NAMES:
        with pytest.raises(ValueError):
            ds = decode_cf(ds)
    else:
        ds = decode_cf(ds, decode_times=CFDatetimeCoder(time_unit=time_unit))

        if calendar not in _STANDARD_CALENDAR_NAMES:
            assert ds.test.dtype == np.dtype("O")
        else:
            assert ds.test.dtype == np.dtype(f"=M8[{time_unit}]")


def test_decode_cf_time_bounds(time_unit: PDDatetimeUnitOptions) -> None:
    da = DataArray(
        np.arange(6, dtype="int64").reshape((3, 2)),
        coords={"time": [1, 2, 3]},
        dims=("time", "nbnd"),
        name="time_bnds",
    )

    attrs = {
        "units": "days since 2001-01",
        "calendar": "standard",
        "bounds": "time_bnds",
    }

    ds = da.to_dataset()
    ds["time"].attrs.update(attrs)
    _update_bounds_attributes(ds.variables)
    assert ds.variables["time_bnds"].attrs == {
        "units": "days since 2001-01",
        "calendar": "standard",
    }
    dsc = decode_cf(ds, decode_times=CFDatetimeCoder(time_unit=time_unit))
    assert dsc.time_bnds.dtype == np.dtype(f"=M8[{time_unit}]")
    dsc = decode_cf(ds, decode_times=False)
    assert dsc.time_bnds.dtype == np.dtype("int64")

    # Do not overwrite existing attrs
    ds = da.to_dataset()
    ds["time"].attrs.update(attrs)
    bnd_attr = {"units": "hours since 2001-01", "calendar": "noleap"}
    ds["time_bnds"].attrs.update(bnd_attr)
    _update_bounds_attributes(ds.variables)
    assert ds.variables["time_bnds"].attrs == bnd_attr

    # If bounds variable not available do not complain
    ds = da.to_dataset()
    ds["time"].attrs.update(attrs)
    ds["time"].attrs["bounds"] = "fake_var"
    _update_bounds_attributes(ds.variables)


@requires_cftime
def test_encode_time_bounds() -> None:
    time = pd.date_range("2000-01-16", periods=1)
    time_bounds = pd.date_range("2000-01-01", periods=2, freq="MS")
    ds = Dataset(dict(time=time, time_bounds=time_bounds))
    ds.time.attrs = {"bounds": "time_bounds"}
    ds.time.encoding = {"calendar": "noleap", "units": "days since 2000-01-01"}

    expected = {}
    # expected['time'] = Variable(data=np.array([15]), dims=['time'])
    expected["time_bounds"] = Variable(data=np.array([0, 31]), dims=["time_bounds"])

    encoded, _ = cf_encoder(ds.variables, ds.attrs)
    assert_equal(encoded["time_bounds"], expected["time_bounds"])
    assert "calendar" not in encoded["time_bounds"].attrs
    assert "units" not in encoded["time_bounds"].attrs

    # if time_bounds attrs are same as time attrs, it doesn't matter
    ds.time_bounds.encoding = {"calendar": "noleap", "units": "days since 2000-01-01"}
    encoded, _ = cf_encoder(dict(ds.variables.items()), ds.attrs)
    assert_equal(encoded["time_bounds"], expected["time_bounds"])
    assert "calendar" not in encoded["time_bounds"].attrs
    assert "units" not in encoded["time_bounds"].attrs

    # for CF-noncompliant case of time_bounds attrs being different from
    # time attrs; preserve them for faithful roundtrip
    ds.time_bounds.encoding = {"calendar": "noleap", "units": "days since 1849-01-01"}
    encoded, _ = cf_encoder(dict(ds.variables.items()), ds.attrs)
    with pytest.raises(AssertionError):
        assert_equal(encoded["time_bounds"], expected["time_bounds"])
    assert "calendar" not in encoded["time_bounds"].attrs
    assert encoded["time_bounds"].attrs["units"] == ds.time_bounds.encoding["units"]

    ds.time.encoding = {}
    with pytest.warns(UserWarning):
        cf_encoder(ds.variables, ds.attrs)


@pytest.fixture(params=_ALL_CALENDARS)
def calendar(request):
    return request.param


@pytest.fixture
def times(calendar):
    import cftime

    return cftime.num2date(
        np.arange(4),
        units="hours since 2000-01-01",
        calendar=calendar,
        only_use_cftime_datetimes=True,
    )


@pytest.fixture
def data(times):
    data = np.random.rand(2, 2, 4)
    lons = np.linspace(0, 11, 2)
    lats = np.linspace(0, 20, 2)
    return DataArray(
        data, coords=[lons, lats, times], dims=["lon", "lat", "time"], name="data"
    )


@pytest.fixture
def times_3d(times):
    lons = np.linspace(0, 11, 2)
    lats = np.linspace(0, 20, 2)
    times_arr = np.random.choice(times, size=(2, 2, 4))
    return DataArray(
        times_arr, coords=[lons, lats, times], dims=["lon", "lat", "time"], name="data"
    )


@requires_cftime
def test_contains_cftime_datetimes_1d(data) -> None:
    assert contains_cftime_datetimes(data.time.variable)


@requires_cftime
@requires_dask
def test_contains_cftime_datetimes_dask_1d(data) -> None:
    assert contains_cftime_datetimes(data.time.variable.chunk())


@requires_cftime
def test_contains_cftime_datetimes_3d(times_3d) -> None:
    assert contains_cftime_datetimes(times_3d.variable)


@requires_cftime
@requires_dask
def test_contains_cftime_datetimes_dask_3d(times_3d) -> None:
    assert contains_cftime_datetimes(times_3d.variable.chunk())


@pytest.mark.parametrize("non_cftime_data", [DataArray([]), DataArray([1, 2])])
def test_contains_cftime_datetimes_non_cftimes(non_cftime_data) -> None:
    assert not contains_cftime_datetimes(non_cftime_data.variable)


@requires_dask
@pytest.mark.parametrize("non_cftime_data", [DataArray([]), DataArray([1, 2])])
def test_contains_cftime_datetimes_non_cftimes_dask(non_cftime_data) -> None:
    assert not contains_cftime_datetimes(non_cftime_data.variable.chunk())


@requires_cftime
@pytest.mark.parametrize("shape", [(24,), (8, 3), (2, 4, 3)])
def test_encode_cf_datetime_overflow(shape) -> None:
    # Test for fix to GH 2272
    dates = pd.date_range("2100", periods=24).values.reshape(shape)
    units = "days since 1800-01-01"
    calendar = "standard"

    num, _, _ = encode_cf_datetime(dates, units, calendar)
    roundtrip = decode_cf_datetime(num, units, calendar)
    np.testing.assert_array_equal(dates, roundtrip)


def test_encode_expected_failures() -> None:
    dates = pd.date_range("2000", periods=3)
    with pytest.raises(ValueError, match="invalid time units"):
        encode_cf_datetime(dates, units="days after 2000-01-01")
    with pytest.raises(ValueError, match="invalid reference date"):
        encode_cf_datetime(dates, units="days since NO_YEAR")


def test_encode_cf_datetime_pandas_min() -> None:
    # GH 2623
    dates = pd.date_range("2000", periods=3)
    num, units, calendar = encode_cf_datetime(dates)
    expected_num = np.array([0.0, 1.0, 2.0])
    expected_units = "days since 2000-01-01 00:00:00"
    expected_calendar = "proleptic_gregorian"
    np.testing.assert_array_equal(num, expected_num)
    assert units == expected_units
    assert calendar == expected_calendar


@requires_cftime
def test_encode_cf_datetime_invalid_pandas_valid_cftime() -> None:
    num, units, calendar = encode_cf_datetime(
        pd.date_range("2000", periods=3),
        # Pandas fails to parse this unit, but cftime is quite happy with it
        "days since 1970-01-01 00:00:00 00",
        "standard",
    )

    expected_num = [10957, 10958, 10959]
    expected_units = "days since 1970-01-01 00:00:00 00"
    expected_calendar = "standard"
    assert_array_equal(num, expected_num)
    assert units == expected_units
    assert calendar == expected_calendar


@requires_cftime
def test_time_units_with_timezone_roundtrip(calendar) -> None:
    # Regression test for GH 2649
    expected_units = "days since 2000-01-01T00:00:00-05:00"
    expected_num_dates = np.array([1, 2, 3])
    dates = decode_cf_datetime(expected_num_dates, expected_units, calendar)

    # Check that dates were decoded to UTC; here the hours should all
    # equal 5.
    result_hours = DataArray(dates).dt.hour
    expected_hours = DataArray([5, 5, 5])
    assert_equal(result_hours, expected_hours)

    # Check that the encoded values are accurately roundtripped.
    result_num_dates, result_units, result_calendar = encode_cf_datetime(
        dates, expected_units, calendar
    )

    if calendar in _STANDARD_CALENDARS:
        assert_duckarray_equal(result_num_dates, expected_num_dates)
    else:
        # cftime datetime arithmetic is not quite exact.
        assert_duckarray_allclose(result_num_dates, expected_num_dates)

    assert result_units == expected_units
    assert result_calendar == calendar


@pytest.mark.parametrize("calendar", _STANDARD_CALENDARS)
def test_use_cftime_default_standard_calendar_in_range(calendar) -> None:
    numerical_dates = [0, 1]
    units = "days since 2000-01-01"
    expected = pd.date_range("2000", periods=2)

    with assert_no_warnings():
        result = decode_cf_datetime(numerical_dates, units, calendar)
        np.testing.assert_array_equal(result, expected)


@requires_cftime
@pytest.mark.parametrize("calendar", ["standard", "gregorian"])
@pytest.mark.parametrize("units_year", [1500, 1580])
def test_use_cftime_default_standard_calendar_out_of_range(
    calendar, units_year
) -> None:
    from cftime import num2date

    numerical_dates = [0, 1]
    units = f"days since {units_year}-01-01"
    expected = num2date(
        numerical_dates, units, calendar, only_use_cftime_datetimes=True
    )

    with pytest.warns(SerializationWarning):
        result = decode_cf_datetime(numerical_dates, units, calendar)
        np.testing.assert_array_equal(result, expected)


@requires_cftime
@pytest.mark.parametrize("calendar", _NON_STANDARD_CALENDARS)
@pytest.mark.parametrize("units_year", [1500, 2000, 2500])
def test_use_cftime_default_non_standard_calendar(
    calendar, units_year, time_unit: PDDatetimeUnitOptions
) -> None:
    from cftime import num2date

    numerical_dates = [0, 1]
    units = f"days since {units_year}-01-01"
    expected = num2date(
        numerical_dates, units, calendar, only_use_cftime_datetimes=True
    )

    if time_unit == "ns" and units_year == 2500:
        with pytest.warns(SerializationWarning, match="Unable to decode time axis"):
            result = decode_cf_datetime(
                numerical_dates, units, calendar, time_unit=time_unit
            )
    else:
        with assert_no_warnings():
            result = decode_cf_datetime(
                numerical_dates, units, calendar, time_unit=time_unit
            )

    np.testing.assert_array_equal(result, expected)


@requires_cftime
@pytest.mark.parametrize("calendar", _ALL_CALENDARS)
@pytest.mark.parametrize("units_year", [1500, 2000, 2500])
def test_use_cftime_true(calendar, units_year) -> None:
    from cftime import num2date

    numerical_dates = [0, 1]
    units = f"days since {units_year}-01-01"
    expected = num2date(
        numerical_dates, units, calendar, only_use_cftime_datetimes=True
    )

    with assert_no_warnings():
        result = decode_cf_datetime(numerical_dates, units, calendar, use_cftime=True)
        np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize("calendar", _STANDARD_CALENDARS)
def test_use_cftime_false_standard_calendar_in_range(calendar) -> None:
    numerical_dates = [0, 1]
    units = "days since 2000-01-01"
    expected = pd.date_range("2000", periods=2)

    with assert_no_warnings():
        result = decode_cf_datetime(numerical_dates, units, calendar, use_cftime=False)
        np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize("calendar", ["standard", "gregorian"])
@pytest.mark.parametrize("units_year", [1500, 1582])
def test_use_cftime_false_standard_calendar_out_of_range(calendar, units_year) -> None:
    numerical_dates = [0, 1]
    units = f"days since {units_year}-01-01"
    with pytest.raises(OutOfBoundsDatetime):
        decode_cf_datetime(numerical_dates, units, calendar, use_cftime=False)


@pytest.mark.parametrize("calendar", _NON_STANDARD_CALENDARS)
@pytest.mark.parametrize("units_year", [1500, 2000, 2500])
def test_use_cftime_false_non_standard_calendar(calendar, units_year) -> None:
    numerical_dates = [0, 1]
    units = f"days since {units_year}-01-01"
    with pytest.raises(OutOfBoundsDatetime):
        decode_cf_datetime(numerical_dates, units, calendar, use_cftime=False)


@requires_cftime
@pytest.mark.parametrize("calendar", _ALL_CALENDARS)
def test_decode_ambiguous_time_warns(calendar) -> None:
    # GH 4422, 4506
    from cftime import num2date

    # we don't decode non-standard calendards with
    # pandas so expect no warning will be emitted
    is_standard_calendar = calendar in _STANDARD_CALENDAR_NAMES

    dates = [1, 2, 3]
    units = "days since 1-1-1"
    expected = num2date(dates, units, calendar=calendar, only_use_cftime_datetimes=True)

    if is_standard_calendar:
        with pytest.warns(SerializationWarning) as record:
            result = decode_cf_datetime(dates, units, calendar=calendar)
        relevant_warnings = [
            r
            for r in record.list
            if str(r.message).startswith("Ambiguous reference date string: 1-1-1")
        ]
        assert len(relevant_warnings) == 1
    else:
        with assert_no_warnings():
            result = decode_cf_datetime(dates, units, calendar=calendar)

    np.testing.assert_array_equal(result, expected)


@pytest.mark.filterwarnings("ignore:Times can't be serialized faithfully")
@pytest.mark.parametrize("encoding_units", FREQUENCIES_TO_ENCODING_UNITS.values())
@pytest.mark.parametrize("freq", FREQUENCIES_TO_ENCODING_UNITS.keys())
@pytest.mark.parametrize("use_cftime", [True, False])
def test_encode_cf_datetime_defaults_to_correct_dtype(
    encoding_units, freq, use_cftime
) -> None:
    if not has_cftime and use_cftime:
        pytest.skip("Test requires cftime")
    if (freq == "ns" or encoding_units == "nanoseconds") and use_cftime:
        pytest.skip("Nanosecond frequency is not valid for cftime dates.")
    times = date_range("2000", periods=3, freq=freq, use_cftime=use_cftime)
    units = f"{encoding_units} since 2000-01-01"
    encoded, _units, _ = encode_cf_datetime(times, units)

    numpy_timeunit = _netcdf_to_numpy_timeunit(encoding_units)
    encoding_units_as_timedelta = np.timedelta64(1, numpy_timeunit)
    if pd.to_timedelta(1, freq) >= encoding_units_as_timedelta:
        assert encoded.dtype == np.int64
    else:
        assert encoded.dtype == np.float64


@pytest.mark.parametrize("freq", FREQUENCIES_TO_ENCODING_UNITS.keys())
def test_encode_decode_roundtrip_datetime64(
    freq, time_unit: PDDatetimeUnitOptions
) -> None:
    # See GH 4045. Prior to GH 4684 this test would fail for frequencies of
    # "s", "ms", "us", and "ns".
    initial_time = pd.date_range("1678-01-01", periods=1)
    times = initial_time.append(pd.date_range("1968", periods=2, freq=freq))
    variable = Variable(["time"], times)
    encoded = conventions.encode_cf_variable(variable)
    decoded = conventions.decode_cf_variable(
        "time", encoded, decode_times=CFDatetimeCoder(time_unit=time_unit)
    )
    assert_equal(variable, decoded)


@requires_cftime
@pytest.mark.parametrize("freq", ["us", "ms", "s", "min", "h", "D"])
def test_encode_decode_roundtrip_cftime(freq) -> None:
    initial_time = date_range("0001", periods=1, use_cftime=True)
    times = initial_time.append(
        date_range("0001", periods=2, freq=freq, use_cftime=True)
        + timedelta(days=291000 * 365)
    )
    variable = Variable(["time"], times)
    encoded = conventions.encode_cf_variable(variable)
    decoder = CFDatetimeCoder(use_cftime=True)
    decoded = conventions.decode_cf_variable("time", encoded, decode_times=decoder)
    assert_equal(variable, decoded)


@requires_cftime
def test__encode_datetime_with_cftime() -> None:
    # See GH 4870. cftime versions > 1.4.0 required us to adapt the
    # way _encode_datetime_with_cftime was written.
    import cftime

    calendar = "gregorian"
    times = cftime.num2date([0, 1], "hours since 2000-01-01", calendar)

    encoding_units = "days since 2000-01-01"
    # Since netCDF files do not support storing float128 values, we ensure that
    # float64 values are used by setting longdouble=False in num2date.  This try
    # except logic can be removed when xarray's minimum version of cftime is at
    # least 1.6.2.
    try:
        expected = cftime.date2num(times, encoding_units, calendar, longdouble=False)
    except TypeError:
        expected = cftime.date2num(times, encoding_units, calendar)
    result = _encode_datetime_with_cftime(times, encoding_units, calendar)
    np.testing.assert_equal(result, expected)


@requires_cftime
def test_round_trip_standard_calendar_cftime_datetimes_pre_reform() -> None:
    from cftime import DatetimeGregorian

    dates = np.array([DatetimeGregorian(1, 1, 1), DatetimeGregorian(2000, 1, 1)])
    encoded = encode_cf_datetime(dates, "seconds since 2000-01-01", "standard")
    with pytest.warns(SerializationWarning, match="Unable to decode time axis"):
        decoded = decode_cf_datetime(*encoded)
    np.testing.assert_equal(decoded, dates)


@pytest.mark.parametrize("calendar", ["standard", "gregorian"])
def test_encode_cf_datetime_gregorian_proleptic_gregorian_mismatch_error(
    calendar: str,
    time_unit: PDDatetimeUnitOptions,
) -> None:
    if time_unit == "ns":
        pytest.skip("datetime64[ns] values can only be defined post reform")
    dates = np.array(["0001-01-01", "2001-01-01"], dtype=f"datetime64[{time_unit}]")
    with pytest.raises(ValueError, match="proleptic_gregorian"):
        encode_cf_datetime(dates, "seconds since 2000-01-01", calendar)


@pytest.mark.parametrize("calendar", ["gregorian", "Gregorian", "GREGORIAN"])
def test_decode_encode_roundtrip_with_non_lowercase_letters(
    calendar, time_unit: PDDatetimeUnitOptions
) -> None:
    # See GH 5093.
    times = [0, 1]
    units = "days since 2000-01-01"
    attrs = {"calendar": calendar, "units": units}
    variable = Variable(["time"], times, attrs)
    decoded = conventions.decode_cf_variable(
        "time", variable, decode_times=CFDatetimeCoder(time_unit=time_unit)
    )
    encoded = conventions.encode_cf_variable(decoded)

    # Previously this would erroneously be an array of cftime.datetime
    # objects.  We check here that it is decoded properly to np.datetime64.
    assert np.issubdtype(decoded.dtype, np.datetime64)

    # Use assert_identical to ensure that the calendar attribute maintained its
    # original form throughout the roundtripping process, uppercase letters and
    # all.
    assert_identical(variable, encoded)


@requires_cftime
def test_should_cftime_be_used_source_outside_range():
    src = date_range(
        "1000-01-01", periods=100, freq="MS", calendar="noleap", use_cftime=True
    )
    with pytest.raises(
        ValueError, match=r"Source time range is not valid for numpy datetimes."
    ):
        _should_cftime_be_used(src, "standard", False)


@requires_cftime
def test_should_cftime_be_used_target_not_npable():
    src = date_range(
        "2000-01-01", periods=100, freq="MS", calendar="noleap", use_cftime=True
    )
    with pytest.raises(
        ValueError, match=r"Calendar 'noleap' is only valid with cftime."
    ):
        _should_cftime_be_used(src, "noleap", False)


@pytest.mark.parametrize(
    "dtype",
    [np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64],
)
def test_decode_cf_datetime_varied_integer_dtypes(dtype):
    units = "seconds since 2018-08-22T03:23:03Z"
    num_dates = dtype(50)
    # Set use_cftime=False to ensure we cannot mask a failure by falling back
    # to cftime.
    result = decode_cf_datetime(num_dates, units, use_cftime=False)
    expected = np.asarray(np.datetime64("2018-08-22T03:23:53", "ns"))
    np.testing.assert_equal(result, expected)


@requires_cftime
def test_decode_cf_datetime_uint64_with_cftime():
    units = "days since 1700-01-01"
    num_dates = np.uint64(182621)
    result = decode_cf_datetime(num_dates, units)
    expected = np.asarray(np.datetime64("2200-01-01", "ns"))
    np.testing.assert_equal(result, expected)


def test_decode_cf_datetime_uint64_with_pandas_overflow_error():
    units = "nanoseconds since 1970-01-01"
    calendar = "standard"
    num_dates = np.uint64(1_000_000 * 86_400 * 360 * 500_000)
    with pytest.raises(OutOfBoundsTimedelta):
        decode_cf_datetime(num_dates, units, calendar, use_cftime=False)


@requires_cftime
def test_decode_cf_datetime_uint64_with_cftime_overflow_error():
    units = "microseconds since 1700-01-01"
    calendar = "360_day"
    num_dates = np.uint64(1_000_000 * 86_400 * 360 * 500_000)
    with pytest.raises(OverflowError):
        decode_cf_datetime(num_dates, units, calendar)


@pytest.mark.parametrize("use_cftime", [True, False])
def test_decode_0size_datetime(use_cftime):
    # GH1329
    if use_cftime and not has_cftime:
        pytest.skip()

    dtype = object if use_cftime else "=M8[ns]"
    expected = np.array([], dtype=dtype)
    actual = decode_cf_datetime(
        np.zeros(shape=0, dtype=np.int64),
        units="days since 1970-01-01 00:00:00",
        calendar="proleptic_gregorian",
        use_cftime=use_cftime,
    )
    np.testing.assert_equal(expected, actual)


def test_decode_float_datetime():
    num_dates = np.array([1867128, 1867134, 1867140], dtype="float32")
    units = "hours since 1800-01-01"
    calendar = "standard"

    expected = np.array(
        ["2013-01-01T00:00:00", "2013-01-01T06:00:00", "2013-01-01T12:00:00"],
        dtype="datetime64[ns]",
    )

    actual = decode_cf_datetime(
        num_dates, units=units, calendar=calendar, use_cftime=False
    )
    np.testing.assert_equal(actual, expected)


@pytest.mark.parametrize("time_unit", ["ms", "us", "ns"])
def test_decode_float_datetime_with_decimals(
    time_unit: PDDatetimeUnitOptions,
) -> None:
    # test resolution enhancement for floats
    values = np.array([0, 0.125, 0.25, 0.375, 0.75, 1.0], dtype="float32")
    expected = np.array(
        [
            "2000-01-01T00:00:00.000",
            "2000-01-01T00:00:00.125",
            "2000-01-01T00:00:00.250",
            "2000-01-01T00:00:00.375",
            "2000-01-01T00:00:00.750",
            "2000-01-01T00:00:01.000",
        ],
        dtype=f"=M8[{time_unit}]",
    )

    units = "seconds since 2000-01-01"
    calendar = "standard"
    actual = decode_cf_datetime(values, units, calendar, time_unit=time_unit)
    assert actual.dtype == expected.dtype
    np.testing.assert_equal(actual, expected)


@pytest.mark.parametrize(
    "time_unit, num", [("s", 0.123), ("ms", 0.1234), ("us", 0.1234567)]
)
def test_coding_float_datetime_warning(
    time_unit: PDDatetimeUnitOptions, num: float
) -> None:
    units = "seconds since 2000-01-01"
    calendar = "standard"
    values = np.array([num], dtype="float32")
    with pytest.warns(
        SerializationWarning,
        match=f"Can't decode floating point datetimes to {time_unit!r}",
    ):
        decode_cf_datetime(values, units, calendar, time_unit=time_unit)


@requires_cftime
def test_scalar_unit() -> None:
    # test that a scalar units (often NaN when using to_netcdf) does not raise an error
    variable = Variable(("x", "y"), np.array([[0, 1], [2, 3]]), {"units": np.nan})
    result = CFDatetimeCoder().decode(variable)
    assert np.isnan(result.attrs["units"])


@requires_cftime
def test_contains_cftime_lazy() -> None:
    import cftime

    from xarray.core.common import _contains_cftime_datetimes

    times = np.array(
        [cftime.DatetimeGregorian(1, 1, 2, 0), cftime.DatetimeGregorian(1, 1, 2, 0)],
        dtype=object,
    )
    array = FirstElementAccessibleArray(times)
    assert _contains_cftime_datetimes(array)


@pytest.mark.parametrize(
    "timestr, format, dtype, fill_value, use_encoding",
    [
        ("1677-09-21T00:12:43.145224193", "ns", np.int64, 20, True),
        ("1970-09-21T00:12:44.145224808", "ns", np.float64, 1e30, True),
        (
            "1677-09-21T00:12:43.145225216",
            "ns",
            np.float64,
            -9.223372036854776e18,
            True,
        ),
        ("1677-09-21T00:12:43.145224193", "ns", np.int64, None, False),
        ("1677-09-21T00:12:43.145225", "us", np.int64, None, False),
        ("1970-01-01T00:00:01.000001", "us", np.int64, None, False),
        ("1677-09-21T00:21:52.901038080", "ns", np.float32, 20.0, True),
    ],
)
def test_roundtrip_datetime64_nanosecond_precision(
    timestr: str,
    format: Literal["ns", "us"],
    dtype: np.typing.DTypeLike | None,
    fill_value: int | float | None,
    use_encoding: bool,
    time_unit: PDDatetimeUnitOptions,
) -> None:
    # test for GH7817
    time = np.datetime64(timestr, format)
    times = [np.datetime64("1970-01-01T00:00:00", format), np.datetime64("NaT"), time]

    if use_encoding:
        encoding = dict(dtype=dtype, _FillValue=fill_value)
    else:
        encoding = {}

    var = Variable(["time"], times, encoding=encoding)
    assert var.dtype == np.dtype(f"=M8[{format}]")

    encoded_var = conventions.encode_cf_variable(var)
    assert (
        encoded_var.attrs["units"]
        == f"{_numpy_to_netcdf_timeunit(format)} since 1970-01-01 00:00:00"
    )
    assert encoded_var.attrs["calendar"] == "proleptic_gregorian"
    assert encoded_var.data.dtype == dtype
    decoded_var = conventions.decode_cf_variable(
        "foo", encoded_var, decode_times=CFDatetimeCoder(time_unit=time_unit)
    )

    result_unit = (
        format
        if np.timedelta64(1, format) <= np.timedelta64(1, time_unit)
        else time_unit
    )
    assert decoded_var.dtype == np.dtype(f"=M8[{result_unit}]")
    assert (
        decoded_var.encoding["units"]
        == f"{_numpy_to_netcdf_timeunit(format)} since 1970-01-01 00:00:00"
    )
    assert decoded_var.encoding["dtype"] == dtype
    assert decoded_var.encoding["calendar"] == "proleptic_gregorian"
    assert_identical(var, decoded_var)


def test_roundtrip_datetime64_nanosecond_precision_warning(
    time_unit: PDDatetimeUnitOptions,
) -> None:
    # test warning if times can't be serialized faithfully
    times = [
        np.datetime64("1970-01-01T00:01:00", time_unit),
        np.datetime64("NaT", time_unit),
        np.datetime64("1970-01-02T00:01:00", time_unit),
    ]
    units = "days since 1970-01-10T01:01:00"
    needed_units = "hours"
    new_units = f"{needed_units} since 1970-01-10T01:01:00"

    encoding = dict(dtype=None, _FillValue=20, units=units)
    var = Variable(["time"], times, encoding=encoding)
    with pytest.warns(UserWarning, match=f"Resolution of {needed_units!r} needed."):
        encoded_var = conventions.encode_cf_variable(var)
    assert encoded_var.dtype == np.float64
    assert encoded_var.attrs["units"] == units
    assert encoded_var.attrs["_FillValue"] == 20.0

    decoded_var = conventions.decode_cf_variable("foo", encoded_var)
    assert_identical(var, decoded_var)

    encoding = dict(dtype="int64", _FillValue=20, units=units)
    var = Variable(["time"], times, encoding=encoding)
    with pytest.warns(
        UserWarning, match=f"Serializing with units {new_units!r} instead."
    ):
        encoded_var = conventions.encode_cf_variable(var)
    assert encoded_var.dtype == np.int64
    assert encoded_var.attrs["units"] == new_units
    assert encoded_var.attrs["_FillValue"] == 20

    decoded_var = conventions.decode_cf_variable(
        "foo", encoded_var, decode_times=CFDatetimeCoder(time_unit=time_unit)
    )
    assert_identical(var, decoded_var)

    encoding = dict(dtype="float64", _FillValue=20, units=units)
    var = Variable(["time"], times, encoding=encoding)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        encoded_var = conventions.encode_cf_variable(var)
    assert encoded_var.dtype == np.float64
    assert encoded_var.attrs["units"] == units
    assert encoded_var.attrs["_FillValue"] == 20.0

    decoded_var = conventions.decode_cf_variable(
        "foo", encoded_var, decode_times=CFDatetimeCoder(time_unit=time_unit)
    )
    assert_identical(var, decoded_var)

    encoding = dict(dtype="int64", _FillValue=20, units=new_units)
    var = Variable(["time"], times, encoding=encoding)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        encoded_var = conventions.encode_cf_variable(var)
    assert encoded_var.dtype == np.int64
    assert encoded_var.attrs["units"] == new_units
    assert encoded_var.attrs["_FillValue"] == 20

    decoded_var = conventions.decode_cf_variable(
        "foo", encoded_var, decode_times=CFDatetimeCoder(time_unit=time_unit)
    )
    assert_identical(var, decoded_var)


@pytest.mark.parametrize(
    "dtype, fill_value",
    [(np.int64, 20), (np.int64, np.iinfo(np.int64).min), (np.float64, 1e30)],
)
def test_roundtrip_timedelta64_nanosecond_precision(
    dtype: np.typing.DTypeLike | None,
    fill_value: int | float,
    time_unit: PDDatetimeUnitOptions,
) -> None:
    # test for GH7942
    one_day = np.timedelta64(1, "ns")
    nat = np.timedelta64("nat", "ns")
    timedelta_values = (np.arange(5) * one_day).astype("timedelta64[ns]")
    timedelta_values[2] = nat
    timedelta_values[4] = nat

    encoding = dict(dtype=dtype, _FillValue=fill_value, units="nanoseconds")
    var = Variable(["time"], timedelta_values, encoding=encoding)

    encoded_var = conventions.encode_cf_variable(var)
    decoded_var = conventions.decode_cf_variable(
        "foo",
        encoded_var,
        decode_times=CFDatetimeCoder(time_unit=time_unit),
        decode_timedelta=CFTimedeltaCoder(time_unit=time_unit),
    )

    assert_identical(var, decoded_var)


def test_roundtrip_timedelta64_nanosecond_precision_warning() -> None:
    # test warning if timedeltas can't be serialized faithfully
    one_day = np.timedelta64(1, "D")
    nat = np.timedelta64("nat", "ns")
    timedelta_values = (np.arange(5) * one_day).astype("timedelta64[ns]")
    timedelta_values[2] = nat
    timedelta_values[4] = np.timedelta64(12, "h").astype("timedelta64[ns]")

    units = "days"
    needed_units = "hours"
    wmsg = (
        f"Timedeltas can't be serialized faithfully with requested units {units!r}. "
        f"Serializing with units {needed_units!r} instead."
    )
    encoding = dict(dtype=np.int64, _FillValue=20, units=units)
    var = Variable(["time"], timedelta_values, encoding=encoding)
    with pytest.warns(UserWarning, match=wmsg):
        encoded_var = conventions.encode_cf_variable(var)
    assert encoded_var.dtype == np.int64
    assert encoded_var.attrs["units"] == needed_units
    assert encoded_var.attrs["_FillValue"] == 20
    decoded_var = conventions.decode_cf_variable(
        "foo", encoded_var, decode_timedelta=CFTimedeltaCoder(time_unit="ns")
    )
    assert_identical(var, decoded_var)
    assert decoded_var.encoding["dtype"] == np.int64


_TEST_ROUNDTRIP_FLOAT_TIMES_TESTS = {
    "GH-8271": (
        20.0,
        np.array(
            ["1970-01-01 00:00:00", "1970-01-01 06:00:00", "NaT"],
            dtype="datetime64[ns]",
        ),
        "days since 1960-01-01",
        np.array([3653, 3653.25, 20.0]),
    ),
    "GH-9488-datetime64[ns]": (
        1.0e20,
        np.array(["2010-01-01 12:00:00", "NaT"], dtype="datetime64[ns]"),
        "seconds since 2010-01-01",
        np.array([43200, 1.0e20]),
    ),
    "GH-9488-timedelta64[ns]": (
        1.0e20,
        np.array([1_000_000_000, "NaT"], dtype="timedelta64[ns]"),
        "seconds",
        np.array([1.0, 1.0e20]),
    ),
}


@pytest.mark.parametrize(
    ("fill_value", "times", "units", "encoded_values"),
    _TEST_ROUNDTRIP_FLOAT_TIMES_TESTS.values(),
    ids=_TEST_ROUNDTRIP_FLOAT_TIMES_TESTS.keys(),
)
def test_roundtrip_float_times(fill_value, times, units, encoded_values) -> None:
    # Regression test for GitHub issues #8271 and #9488
    var = Variable(
        ["time"],
        times,
        encoding=dict(dtype=np.float64, _FillValue=fill_value, units=units),
    )

    encoded_var = conventions.encode_cf_variable(var)
    np.testing.assert_array_equal(encoded_var, encoded_values)
    assert encoded_var.attrs["units"] == units
    assert encoded_var.attrs["_FillValue"] == fill_value

    decoded_var = conventions.decode_cf_variable(
        "foo", encoded_var, decode_timedelta=CFTimedeltaCoder(time_unit="ns")
    )
    assert_identical(var, decoded_var)
    assert decoded_var.encoding["units"] == units
    assert decoded_var.encoding["_FillValue"] == fill_value


_ENCODE_DATETIME64_VIA_DASK_TESTS = {
    "pandas-encoding-with-prescribed-units-and-dtype": (
        "D",
        "days since 1700-01-01",
        np.dtype("int32"),
    ),
    "mixed-cftime-pandas-encoding-with-prescribed-units-and-dtype": pytest.param(
        "250YS", "days since 1700-01-01", np.dtype("int32"), marks=requires_cftime
    ),
    "pandas-encoding-with-default-units-and-dtype": ("250YS", None, None),
}


@requires_dask
@pytest.mark.parametrize(
    ("freq", "units", "dtype"),
    _ENCODE_DATETIME64_VIA_DASK_TESTS.values(),
    ids=_ENCODE_DATETIME64_VIA_DASK_TESTS.keys(),
)
def test_encode_cf_datetime_datetime64_via_dask(
    freq, units, dtype, time_unit: PDDatetimeUnitOptions
) -> None:
    import dask.array

    times_pd = pd.date_range(start="1700", freq=freq, periods=3, unit=time_unit)
    times = dask.array.from_array(times_pd, chunks=1)
    encoded_times, encoding_units, encoding_calendar = encode_cf_datetime(
        times, units, None, dtype
    )

    assert is_duck_dask_array(encoded_times)
    assert encoded_times.chunks == times.chunks

    if units is not None and dtype is not None:
        assert encoding_units == units
        assert encoded_times.dtype == dtype
    else:
        expected_netcdf_time_unit = _numpy_to_netcdf_timeunit(time_unit)
        assert encoding_units == f"{expected_netcdf_time_unit} since 1970-01-01"
        assert encoded_times.dtype == np.dtype("int64")

    assert encoding_calendar == "proleptic_gregorian"

    decoded_times = decode_cf_datetime(
        encoded_times, encoding_units, encoding_calendar, time_unit=time_unit
    )
    np.testing.assert_equal(decoded_times, times)
    assert decoded_times.dtype == times.dtype


@requires_dask
@pytest.mark.parametrize(
    ("range_function", "start", "units", "dtype"),
    [
        (pd.date_range, "2000", None, np.dtype("int32")),
        (pd.date_range, "2000", "days since 2000-01-01", None),
        (pd.timedelta_range, "0D", None, np.dtype("int32")),
        (pd.timedelta_range, "0D", "days", None),
    ],
)
def test_encode_via_dask_cannot_infer_error(
    range_function, start, units, dtype
) -> None:
    values = range_function(start=start, freq="D", periods=3)
    encoding = dict(units=units, dtype=dtype)
    variable = Variable(["time"], values, encoding=encoding).chunk({"time": 1})
    with pytest.raises(ValueError, match="When encoding chunked arrays"):
        conventions.encode_cf_variable(variable)


@requires_cftime
@requires_dask
@pytest.mark.parametrize(
    ("units", "dtype"), [("days since 1700-01-01", np.dtype("int32")), (None, None)]
)
def test_encode_cf_datetime_cftime_datetime_via_dask(units, dtype) -> None:
    import dask.array

    calendar = "standard"
    times_idx = date_range(
        start="1700", freq="D", periods=3, calendar=calendar, use_cftime=True
    )
    times = dask.array.from_array(times_idx, chunks=1)
    encoded_times, encoding_units, encoding_calendar = encode_cf_datetime(
        times, units, None, dtype
    )

    assert is_duck_dask_array(encoded_times)
    assert encoded_times.chunks == times.chunks

    if units is not None and dtype is not None:
        assert encoding_units == units
        assert encoded_times.dtype == dtype
    else:
        assert encoding_units == "microseconds since 1970-01-01"
        assert encoded_times.dtype == np.int64

    assert encoding_calendar == calendar

    decoded_times = decode_cf_datetime(
        encoded_times, encoding_units, encoding_calendar, use_cftime=True
    )
    np.testing.assert_equal(decoded_times, times)


@pytest.mark.parametrize(
    "use_cftime", [False, pytest.param(True, marks=requires_cftime)]
)
@pytest.mark.parametrize("use_dask", [False, pytest.param(True, marks=requires_dask)])
def test_encode_cf_datetime_units_change(use_cftime, use_dask) -> None:
    times = date_range(start="2000", freq="12h", periods=3, use_cftime=use_cftime)
    encoding = dict(units="days since 2000-01-01", dtype=np.dtype("int64"))
    variable = Variable(["time"], times, encoding=encoding)

    if use_dask:
        variable = variable.chunk({"time": 1})
        with pytest.raises(ValueError, match="Times can't be serialized"):
            conventions.encode_cf_variable(variable).compute()
    else:
        with pytest.warns(UserWarning, match="Times can't be serialized"):
            encoded = conventions.encode_cf_variable(variable)
        if use_cftime:
            expected_units = "hours since 2000-01-01 00:00:00.000000"
        else:
            expected_units = "hours since 2000-01-01"
        assert encoded.attrs["units"] == expected_units
        decoded = conventions.decode_cf_variable(
            "name", encoded, decode_times=CFDatetimeCoder(use_cftime=use_cftime)
        )
        assert_equal(variable, decoded)


@pytest.mark.parametrize("use_dask", [False, pytest.param(True, marks=requires_dask)])
def test_encode_cf_datetime_precision_loss_regression_test(use_dask) -> None:
    # Regression test for
    # https://github.com/pydata/xarray/issues/9134#issuecomment-2191446463
    times = date_range("2000", periods=5, freq="ns")
    encoding = dict(units="seconds since 1970-01-01", dtype=np.dtype("int64"))
    variable = Variable(["time"], times, encoding=encoding)

    if use_dask:
        variable = variable.chunk({"time": 1})
        with pytest.raises(ValueError, match="Times can't be serialized"):
            conventions.encode_cf_variable(variable).compute()
    else:
        with pytest.warns(UserWarning, match="Times can't be serialized"):
            encoded = conventions.encode_cf_variable(variable)
        decoded = conventions.decode_cf_variable("name", encoded)
        assert_equal(variable, decoded)


@requires_dask
@pytest.mark.parametrize(
    ("units", "dtype"), [("days", np.dtype("int32")), (None, None)]
)
def test_encode_cf_timedelta_via_dask(
    units: str | None, dtype: np.dtype | None, time_unit: PDDatetimeUnitOptions
) -> None:
    import dask.array

    times_pd = pd.timedelta_range(start="0D", freq="D", periods=3, unit=time_unit)  # type: ignore[call-arg,unused-ignore]
    times = dask.array.from_array(times_pd, chunks=1)
    encoded_times, encoding_units = encode_cf_timedelta(times, units, dtype)

    assert is_duck_dask_array(encoded_times)
    assert encoded_times.chunks == times.chunks

    if units is not None and dtype is not None:
        assert encoding_units == units
        assert encoded_times.dtype == dtype
    else:
        assert encoding_units == _numpy_to_netcdf_timeunit(time_unit)
        assert encoded_times.dtype == np.dtype("int64")

    decoded_times = decode_cf_timedelta(
        encoded_times, encoding_units, time_unit=time_unit
    )
    np.testing.assert_equal(decoded_times, times)
    assert decoded_times.dtype == times.dtype


@pytest.mark.parametrize("use_dask", [False, pytest.param(True, marks=requires_dask)])
def test_encode_cf_timedelta_units_change(use_dask) -> None:
    timedeltas = pd.timedelta_range(start="0h", freq="12h", periods=3)
    encoding = dict(units="days", dtype=np.dtype("int64"))
    variable = Variable(["time"], timedeltas, encoding=encoding)

    if use_dask:
        variable = variable.chunk({"time": 1})
        with pytest.raises(ValueError, match="Timedeltas can't be serialized"):
            conventions.encode_cf_variable(variable).compute()
    else:
        # In this case we automatically modify the encoding units to continue
        # encoding with integer values.
        with pytest.warns(UserWarning, match="Timedeltas can't be serialized"):
            encoded = conventions.encode_cf_variable(variable)
        assert encoded.attrs["units"] == "hours"
        decoded = conventions.decode_cf_variable(
            "name", encoded, decode_timedelta=CFTimedeltaCoder(time_unit="ns")
        )
        assert_equal(variable, decoded)


@pytest.mark.parametrize("use_dask", [False, pytest.param(True, marks=requires_dask)])
def test_encode_cf_timedelta_small_dtype_missing_value(use_dask) -> None:
    # Regression test for GitHub issue #9134
    timedeltas = np.array([1, 2, "NaT", 4], dtype="timedelta64[D]").astype(
        "timedelta64[ns]"
    )
    encoding = dict(units="days", dtype=np.dtype("int16"), _FillValue=np.int16(-1))
    variable = Variable(["time"], timedeltas, encoding=encoding)

    if use_dask:
        variable = variable.chunk({"time": 1})

    encoded = conventions.encode_cf_variable(variable)
    decoded = conventions.decode_cf_variable("name", encoded, decode_timedelta=True)
    assert_equal(variable, decoded)


_DECODE_TIMEDELTA_VIA_UNITS_TESTS = {
    "default": (True, None, np.dtype("timedelta64[ns]"), True),
    "decode_timedelta=True": (True, True, np.dtype("timedelta64[ns]"), False),
    "decode_timedelta=False": (True, False, np.dtype("int64"), False),
    "inherit-time_unit-from-decode_times": (
        CFDatetimeCoder(time_unit="s"),
        None,
        np.dtype("timedelta64[s]"),
        True,
    ),
    "set-time_unit-via-CFTimedeltaCoder-decode_times=True": (
        True,
        CFTimedeltaCoder(time_unit="s"),
        np.dtype("timedelta64[s]"),
        False,
    ),
    "set-time_unit-via-CFTimedeltaCoder-decode_times=False": (
        False,
        CFTimedeltaCoder(time_unit="s"),
        np.dtype("timedelta64[s]"),
        False,
    ),
    "override-time_unit-from-decode_times": (
        CFDatetimeCoder(time_unit="ns"),
        CFTimedeltaCoder(time_unit="s"),
        np.dtype("timedelta64[s]"),
        False,
    ),
}


@pytest.mark.parametrize(
    ("decode_times", "decode_timedelta", "expected_dtype", "warns"),
    list(_DECODE_TIMEDELTA_VIA_UNITS_TESTS.values()),
    ids=list(_DECODE_TIMEDELTA_VIA_UNITS_TESTS.keys()),
)
def test_decode_timedelta_via_units(
    decode_times, decode_timedelta, expected_dtype, warns
) -> None:
    timedeltas = pd.timedelta_range(0, freq="D", periods=3)
    attrs = {"units": "days"}
    var = Variable(["time"], timedeltas, encoding=attrs)
    encoded = Variable(["time"], np.array([0, 1, 2]), attrs=attrs)
    if warns:
        with pytest.warns(
            FutureWarning,
            match="xarray will not decode the variable 'foo' into a timedelta64 dtype",
        ):
            decoded = conventions.decode_cf_variable(
                "foo",
                encoded,
                decode_times=decode_times,
                decode_timedelta=decode_timedelta,
            )
    else:
        decoded = conventions.decode_cf_variable(
            "foo", encoded, decode_times=decode_times, decode_timedelta=decode_timedelta
        )
    if decode_timedelta is False:
        assert_equal(encoded, decoded)
    else:
        assert_equal(var, decoded)
    assert decoded.dtype == expected_dtype


_DECODE_TIMEDELTA_VIA_DTYPE_TESTS = {
    "default": (True, None, "ns", np.dtype("timedelta64[ns]")),
    "decode_timedelta=False": (True, False, "ns", np.dtype("int64")),
    "decode_timedelta=True": (True, True, "ns", np.dtype("timedelta64[ns]")),
    "use-original-units": (True, True, "s", np.dtype("timedelta64[s]")),
    "inherit-time_unit-from-decode_times": (
        CFDatetimeCoder(time_unit="s"),
        None,
        "ns",
        np.dtype("timedelta64[s]"),
    ),
    "set-time_unit-via-CFTimedeltaCoder-decode_times=True": (
        True,
        CFTimedeltaCoder(time_unit="s"),
        "ns",
        np.dtype("timedelta64[s]"),
    ),
    "set-time_unit-via-CFTimedeltaCoder-decode_times=False": (
        False,
        CFTimedeltaCoder(time_unit="s"),
        "ns",
        np.dtype("timedelta64[s]"),
    ),
    "override-time_unit-from-decode_times": (
        CFDatetimeCoder(time_unit="ns"),
        CFTimedeltaCoder(time_unit="s"),
        "ns",
        np.dtype("timedelta64[s]"),
    ),
    "decode-different-units": (
        True,
        CFTimedeltaCoder(time_unit="us"),
        "s",
        np.dtype("timedelta64[us]"),
    ),
}


@pytest.mark.parametrize(
    ("decode_times", "decode_timedelta", "original_unit", "expected_dtype"),
    list(_DECODE_TIMEDELTA_VIA_DTYPE_TESTS.values()),
    ids=list(_DECODE_TIMEDELTA_VIA_DTYPE_TESTS.keys()),
)
def test_decode_timedelta_via_dtype(
    decode_times, decode_timedelta, original_unit, expected_dtype
) -> None:
    timedeltas = pd.timedelta_range(0, freq="D", periods=3, unit=original_unit)  # type: ignore[call-arg,unused-ignore]
    encoding = {"units": "days"}
    var = Variable(["time"], timedeltas, encoding=encoding)
    encoded = conventions.encode_cf_variable(var)
    assert encoded.attrs["dtype"] == f"timedelta64[{original_unit}]"
    assert encoded.attrs["units"] == encoding["units"]
    decoded = conventions.decode_cf_variable(
        "foo", encoded, decode_times=decode_times, decode_timedelta=decode_timedelta
    )
    if decode_timedelta is False:
        assert_equal(encoded, decoded)
    else:
        assert_equal(var, decoded)
    assert decoded.dtype == expected_dtype


@pytest.mark.parametrize("dtype", [np.uint64, np.int64, np.float64])
def test_decode_timedelta_dtypes(dtype) -> None:
    encoded = Variable(["time"], np.arange(10), {"units": "seconds"})
    coder = CFTimedeltaCoder(time_unit="s")
    decoded = coder.decode(encoded)
    assert decoded.dtype.kind == "m"
    assert_equal(coder.encode(decoded), encoded)


def test_lazy_decode_timedelta_unexpected_dtype() -> None:
    attrs = {"units": "seconds"}
    encoded = Variable(["time"], [0, 0.5, 1], attrs=attrs)
    decoded = conventions.decode_cf_variable(
        "foo", encoded, decode_timedelta=CFTimedeltaCoder(time_unit="s")
    )

    expected_dtype_upon_lazy_decoding = np.dtype("timedelta64[s]")
    assert decoded.dtype == expected_dtype_upon_lazy_decoding

    expected_dtype_upon_loading = np.dtype("timedelta64[ms]")
    with pytest.warns(SerializationWarning, match="Can't decode floating"):
        assert decoded.load().dtype == expected_dtype_upon_loading


def test_lazy_decode_timedelta_error() -> None:
    attrs = {"units": "seconds"}
    encoded = Variable(["time"], [0, np.iinfo(np.int64).max, 1], attrs=attrs)
    decoded = conventions.decode_cf_variable(
        "foo", encoded, decode_timedelta=CFTimedeltaCoder(time_unit="ms")
    )
    with pytest.raises(OutOfBoundsTimedelta, match="overflow"):
        decoded.load()


@pytest.mark.parametrize(
    "calendar",
    [
        "standard",
        pytest.param(
            "360_day", marks=pytest.mark.skipif(not has_cftime, reason="no cftime")
        ),
    ],
)
def test_duck_array_decode_times(calendar) -> None:
    from xarray.core.indexing import LazilyIndexedArray

    days = LazilyIndexedArray(DuckArrayWrapper(np.array([1.0, 2.0, 3.0])))
    var = Variable(
        ["time"], days, {"units": "days since 2001-01-01", "calendar": calendar}
    )
    decoded = conventions.decode_cf_variable(
        "foo", var, decode_times=CFDatetimeCoder(use_cftime=None)
    )
    if calendar not in _STANDARD_CALENDAR_NAMES:
        assert decoded.dtype == np.dtype("O")
    else:
        assert decoded.dtype == np.dtype("=M8[ns]")


@pytest.mark.parametrize("decode_timedelta", [True, False])
@pytest.mark.parametrize("mask_and_scale", [True, False])
def test_decode_timedelta_mask_and_scale(
    decode_timedelta: bool, mask_and_scale: bool
) -> None:
    attrs = {
        "dtype": "timedelta64[ns]",
        "units": "nanoseconds",
        "_FillValue": np.int16(-1),
        "add_offset": 100000.0,
    }
    encoded = Variable(["time"], np.array([0, -1, 1], "int16"), attrs=attrs)
    decoded = conventions.decode_cf_variable(
        "foo", encoded, mask_and_scale=mask_and_scale, decode_timedelta=decode_timedelta
    )
    result = conventions.encode_cf_variable(decoded, name="foo")
    assert_identical(encoded, result)
    assert encoded.dtype == result.dtype


def test_decode_floating_point_timedelta_no_serialization_warning() -> None:
    attrs = {"units": "seconds"}
    encoded = Variable(["time"], [0, 0.1, 0.2], attrs=attrs)
    decoded = conventions.decode_cf_variable("foo", encoded, decode_timedelta=True)
    with assert_no_warnings():
        decoded.load()


def test_timedelta64_coding_via_dtype(time_unit: PDDatetimeUnitOptions) -> None:
    timedeltas = np.array([0, 1, "NaT"], dtype=f"timedelta64[{time_unit}]")
    variable = Variable(["time"], timedeltas)
    expected_units = _numpy_to_netcdf_timeunit(time_unit)

    encoded = conventions.encode_cf_variable(variable)
    assert encoded.attrs["dtype"] == f"timedelta64[{time_unit}]"
    assert encoded.attrs["units"] == expected_units

    decoded = conventions.decode_cf_variable("timedeltas", encoded)
    assert decoded.encoding["dtype"] == np.dtype("int64")
    assert decoded.encoding["units"] == expected_units

    assert_identical(decoded, variable)
    assert decoded.dtype == variable.dtype

    reencoded = conventions.encode_cf_variable(decoded)
    assert_identical(reencoded, encoded)
    assert reencoded.dtype == encoded.dtype


def test_timedelta_coding_via_dtype_non_pandas_coarse_resolution_warning() -> None:
    attrs = {"dtype": "timedelta64[D]", "units": "days"}
    encoded = Variable(["time"], [0, 1, 2], attrs=attrs)
    with pytest.warns(UserWarning, match="xarray only supports"):
        decoded = conventions.decode_cf_variable("timedeltas", encoded)
    expected_array = np.array([0, 1, 2], dtype="timedelta64[D]")
    expected_array = expected_array.astype("timedelta64[s]")
    expected = Variable(["time"], expected_array)
    assert_identical(decoded, expected)
    assert decoded.dtype == np.dtype("timedelta64[s]")


@pytest.mark.xfail(reason="xarray does not recognize picoseconds as time-like")
def test_timedelta_coding_via_dtype_non_pandas_fine_resolution_warning() -> None:
    attrs = {"dtype": "timedelta64[ps]", "units": "picoseconds"}
    encoded = Variable(["time"], [0, 1000, 2000], attrs=attrs)
    with pytest.warns(UserWarning, match="xarray only supports"):
        decoded = conventions.decode_cf_variable("timedeltas", encoded)
    expected_array = np.array([0, 1000, 2000], dtype="timedelta64[ps]")
    expected_array = expected_array.astype("timedelta64[ns]")
    expected = Variable(["time"], expected_array)
    assert_identical(decoded, expected)
    assert decoded.dtype == np.dtype("timedelta64[ns]")


def test_timedelta_decode_via_dtype_invalid_encoding() -> None:
    attrs = {"dtype": "timedelta64[s]", "units": "seconds"}
    encoding = {"units": "foo"}
    encoded = Variable(["time"], [0, 1, 2], attrs=attrs, encoding=encoding)
    with pytest.raises(ValueError, match=r"Key .* already exists"):
        conventions.decode_cf_variable("timedeltas", encoded)


@pytest.mark.parametrize("attribute", ["dtype", "units"])
def test_timedelta_encode_via_dtype_invalid_attribute(attribute) -> None:
    timedeltas = pd.timedelta_range(0, freq="D", periods=3)
    attrs = {attribute: "foo"}
    variable = Variable(["time"], timedeltas, attrs=attrs)
    with pytest.raises(ValueError, match=r"Key .* already exists"):
        conventions.encode_cf_variable(variable)


@pytest.mark.parametrize(
    ("decode_via_units", "decode_via_dtype", "attrs", "expect_timedelta64"),
    [
        (True, True, {"units": "seconds"}, True),
        (True, False, {"units": "seconds"}, True),
        (False, True, {"units": "seconds"}, False),
        (False, False, {"units": "seconds"}, False),
        (True, True, {"dtype": "timedelta64[s]", "units": "seconds"}, True),
        (True, False, {"dtype": "timedelta64[s]", "units": "seconds"}, True),
        (False, True, {"dtype": "timedelta64[s]", "units": "seconds"}, True),
        (False, False, {"dtype": "timedelta64[s]", "units": "seconds"}, False),
    ],
    ids=lambda x: f"{x!r}",
)
def test_timedelta_decoding_options(
    decode_via_units, decode_via_dtype, attrs, expect_timedelta64
) -> None:
    array = np.array([0, 1, 2], dtype=np.dtype("int64"))
    encoded = Variable(["time"], array, attrs=attrs)

    # Confirm we decode to the expected dtype.
    decode_timedelta = CFTimedeltaCoder(
        time_unit="s",
        decode_via_units=decode_via_units,
        decode_via_dtype=decode_via_dtype,
    )
    decoded = conventions.decode_cf_variable(
        "foo", encoded, decode_timedelta=decode_timedelta
    )
    if expect_timedelta64:
        assert decoded.dtype == np.dtype("timedelta64[s]")
    else:
        assert decoded.dtype == np.dtype("int64")

    # Confirm we exactly roundtrip.
    reencoded = conventions.encode_cf_variable(decoded)

    expected = encoded.copy()
    if "dtype" not in attrs and decode_via_units:
        expected.attrs["dtype"] = "timedelta64[s]"
    assert_identical(reencoded, expected)


def test_timedelta_encoding_explicit_non_timedelta64_dtype() -> None:
    encoding = {"dtype": np.dtype("int32")}
    timedeltas = pd.timedelta_range(0, freq="D", periods=3)
    variable = Variable(["time"], timedeltas, encoding=encoding)

    encoded = conventions.encode_cf_variable(variable)
    assert encoded.attrs["units"] == "days"
    assert encoded.attrs["dtype"] == "timedelta64[ns]"
    assert encoded.dtype == np.dtype("int32")

    decoded = conventions.decode_cf_variable("foo", encoded)
    assert_identical(decoded, variable)

    reencoded = conventions.encode_cf_variable(decoded)
    assert_identical(reencoded, encoded)
    assert encoded.attrs["units"] == "days"
    assert encoded.attrs["dtype"] == "timedelta64[ns]"
    assert encoded.dtype == np.dtype("int32")


@pytest.mark.parametrize("mask_attribute", ["_FillValue", "missing_value"])
def test_timedelta64_coding_via_dtype_with_mask(
    time_unit: PDDatetimeUnitOptions, mask_attribute: str
) -> None:
    timedeltas = np.array([0, 1, "NaT"], dtype=f"timedelta64[{time_unit}]")
    mask = 10
    variable = Variable(["time"], timedeltas, encoding={mask_attribute: mask})
    expected_dtype = f"timedelta64[{time_unit}]"
    expected_units = _numpy_to_netcdf_timeunit(time_unit)

    encoded = conventions.encode_cf_variable(variable)
    assert encoded.attrs["dtype"] == expected_dtype
    assert encoded.attrs["units"] == expected_units
    assert encoded.attrs[mask_attribute] == mask
    assert encoded[-1] == mask

    decoded = conventions.decode_cf_variable("timedeltas", encoded)
    assert decoded.encoding["dtype"] == np.dtype("int64")
    assert decoded.encoding["units"] == expected_units
    assert decoded.encoding[mask_attribute] == mask
    assert np.isnat(decoded[-1])

    assert_identical(decoded, variable)
    assert decoded.dtype == variable.dtype

    reencoded = conventions.encode_cf_variable(decoded)
    assert_identical(reencoded, encoded)
    assert reencoded.dtype == encoded.dtype


def test_roundtrip_0size_timedelta(time_unit: PDDatetimeUnitOptions) -> None:
    # regression test for GitHub issue #10310
    encoding = {"units": "days", "dtype": np.dtype("int64")}
    data = np.array([], dtype=f"=m8[{time_unit}]")
    decoded = Variable(["time"], data, encoding=encoding)
    encoded = conventions.encode_cf_variable(decoded, name="foo")
    assert encoded.dtype == encoding["dtype"]
    assert encoded.attrs["units"] == encoding["units"]
    decoded = conventions.decode_cf_variable("foo", encoded, decode_timedelta=True)
    assert decoded.dtype == np.dtype(f"=m8[{time_unit}]")
    with assert_no_warnings():
        decoded.load()
    assert decoded.dtype == np.dtype("=m8[s]")
    assert decoded.encoding == encoding


def test_roundtrip_empty_datetime64_array(time_unit: PDDatetimeUnitOptions) -> None:
    # Regression test for GitHub issue #10722.
    encoding = {
        "units": "days since 1990-1-1",
        "dtype": np.dtype("float64"),
        "calendar": "standard",
    }
    times = date_range("2000", periods=0, unit=time_unit)
    variable = Variable(["time"], times, encoding=encoding)

    encoded = conventions.encode_cf_variable(variable, name="foo")
    assert encoded.dtype == np.dtype("float64")

    decode_times = CFDatetimeCoder(time_unit=time_unit)
    roundtripped = conventions.decode_cf_variable(
        "foo", encoded, decode_times=decode_times
    )
    assert_identical(variable, roundtripped)
    assert roundtripped.dtype == variable.dtype
