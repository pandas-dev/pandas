from __future__ import annotations

import pickle
from datetime import timedelta
from textwrap import dedent

import numpy as np
import pandas as pd
import pytest

import xarray as xr
from xarray.coding.cftimeindex import (
    CFTimeIndex,
    _parse_array_of_cftime_strings,
    _parsed_string_to_bounds,
    assert_all_valid_date_type,
)
from xarray.coding.times import (
    _parse_iso8601,
    parse_iso8601_like,
)
from xarray.core.types import PDDatetimeUnitOptions
from xarray.tests import (
    _ALL_CALENDARS,
    _NON_STANDARD_CALENDAR_NAMES,
    _all_cftime_date_types,
    assert_array_equal,
    assert_identical,
    has_cftime,
    has_pandas_3,
    requires_cftime,
)

# cftime 1.5.2 renames "gregorian" to "standard"
standard_or_gregorian = ""
if has_cftime:
    standard_or_gregorian = "standard"


def date_dict(
    year=None,
    month=None,
    day=None,
    hour=None,
    minute=None,
    second=None,
    microsecond=None,
):
    return dict(
        year=year,
        month=month,
        day=day,
        hour=hour,
        minute=minute,
        second=second,
        microsecond=microsecond,
    )


ISO8601_LIKE_STRING_TESTS = {
    "year": ("1999", date_dict(year="1999")),
    "month": ("199901", date_dict(year="1999", month="01")),
    "month-dash": ("1999-01", date_dict(year="1999", month="01")),
    "day": ("19990101", date_dict(year="1999", month="01", day="01")),
    "day-dash": ("1999-01-01", date_dict(year="1999", month="01", day="01")),
    "hour": ("19990101T12", date_dict(year="1999", month="01", day="01", hour="12")),
    "hour-dash": (
        "1999-01-01T12",
        date_dict(year="1999", month="01", day="01", hour="12"),
    ),
    "hour-space-separator": (
        "1999-01-01 12",
        date_dict(year="1999", month="01", day="01", hour="12"),
    ),
    "minute": (
        "19990101T1234",
        date_dict(year="1999", month="01", day="01", hour="12", minute="34"),
    ),
    "minute-dash": (
        "1999-01-01T12:34",
        date_dict(year="1999", month="01", day="01", hour="12", minute="34"),
    ),
    "minute-space-separator": (
        "1999-01-01 12:34",
        date_dict(year="1999", month="01", day="01", hour="12", minute="34"),
    ),
    "second": (
        "19990101T123456",
        date_dict(
            year="1999", month="01", day="01", hour="12", minute="34", second="56"
        ),
    ),
    "second-dash": (
        "1999-01-01T12:34:56",
        date_dict(
            year="1999", month="01", day="01", hour="12", minute="34", second="56"
        ),
    ),
    "second-space-separator": (
        "1999-01-01 12:34:56",
        date_dict(
            year="1999", month="01", day="01", hour="12", minute="34", second="56"
        ),
    ),
    "microsecond-1": (
        "19990101T123456.123456",
        date_dict(
            year="1999",
            month="01",
            day="01",
            hour="12",
            minute="34",
            second="56",
            microsecond="123456",
        ),
    ),
    "microsecond-2": (
        "19990101T123456.1",
        date_dict(
            year="1999",
            month="01",
            day="01",
            hour="12",
            minute="34",
            second="56",
            microsecond="1",
        ),
    ),
}


@pytest.mark.parametrize(
    ("string", "expected"),
    list(ISO8601_LIKE_STRING_TESTS.values()),
    ids=list(ISO8601_LIKE_STRING_TESTS.keys()),
)
@pytest.mark.parametrize(
    "five_digit_year", [False, True], ids=["four-digit-year", "five-digit-year"]
)
@pytest.mark.parametrize("sign", ["", "+", "-"], ids=["None", "plus", "minus"])
def test_parse_iso8601_like(
    five_digit_year: bool, sign: str, string: str, expected: dict
) -> None:
    pre = "1" if five_digit_year else ""
    datestring = sign + pre + string
    result = parse_iso8601_like(datestring)
    expected = expected.copy()
    expected.update(year=sign + pre + expected["year"])
    assert result == expected

    # check malformed single digit addendum
    # this check is only performed when we have at least "hour" given
    # like "1999010101", where a single added digit should raise
    # for "1999" (year), "199901" (month) and "19990101" (day)
    # and a single added digit the string would just be interpreted
    # as having a 5-digit year.
    if result["microsecond"] is None and result["hour"] is not None:
        with pytest.raises(ValueError):
            parse_iso8601_like(datestring + "3")

    # check malformed floating point addendum
    if result["second"] is None or result["microsecond"] is not None:
        with pytest.raises(ValueError):
            parse_iso8601_like(datestring + ".3")


_CFTIME_CALENDARS = [
    "365_day",
    "360_day",
    "julian",
    "all_leap",
    "366_day",
    "gregorian",
    "proleptic_gregorian",
]


@pytest.fixture(params=_CFTIME_CALENDARS)
def date_type(request):
    return _all_cftime_date_types()[request.param]


@pytest.fixture
def index(date_type):
    dates = [
        date_type(1, 1, 1),
        date_type(1, 2, 1),
        date_type(2, 1, 1),
        date_type(2, 2, 1),
    ]
    return CFTimeIndex(dates)


@pytest.fixture
def monotonic_decreasing_index(date_type):
    dates = [
        date_type(2, 2, 1),
        date_type(2, 1, 1),
        date_type(1, 2, 1),
        date_type(1, 1, 1),
    ]
    return CFTimeIndex(dates)


@pytest.fixture
def length_one_index(date_type):
    dates = [date_type(1, 1, 1)]
    return CFTimeIndex(dates)


@pytest.fixture
def da(index):
    return xr.DataArray([1, 2, 3, 4], coords=[index], dims=["time"])


@pytest.fixture
def series(index):
    return pd.Series([1, 2, 3, 4], index=index)


@pytest.fixture
def df(index):
    return pd.DataFrame([1, 2, 3, 4], index=index)


@pytest.fixture
def feb_days(date_type):
    import cftime

    if date_type is cftime.DatetimeAllLeap:
        return 29
    elif date_type is cftime.Datetime360Day:
        return 30
    else:
        return 28


@pytest.fixture
def dec_days(date_type):
    import cftime

    if date_type is cftime.Datetime360Day:
        return 30
    else:
        return 31


@pytest.fixture
def index_with_name(date_type):
    dates = [
        date_type(1, 1, 1),
        date_type(1, 2, 1),
        date_type(2, 1, 1),
        date_type(2, 2, 1),
    ]
    return CFTimeIndex(dates, name="foo")


@requires_cftime
@pytest.mark.parametrize(("name", "expected_name"), [("bar", "bar"), (None, "foo")])
def test_constructor_with_name(index_with_name, name, expected_name):
    result = CFTimeIndex(index_with_name, name=name).name
    assert result == expected_name


@requires_cftime
def test_assert_all_valid_date_type(date_type, index):
    import cftime

    if date_type is cftime.DatetimeNoLeap:
        mixed_date_types = np.array(
            [date_type(1, 1, 1), cftime.DatetimeAllLeap(1, 2, 1)]
        )
    else:
        mixed_date_types = np.array(
            [date_type(1, 1, 1), cftime.DatetimeNoLeap(1, 2, 1)]
        )
    with pytest.raises(TypeError):
        assert_all_valid_date_type(mixed_date_types)

    with pytest.raises(TypeError):
        assert_all_valid_date_type(np.array([1, date_type(1, 1, 1)]))

    assert_all_valid_date_type(np.array([date_type(1, 1, 1), date_type(1, 2, 1)]))


@requires_cftime
@pytest.mark.parametrize(
    ("field", "expected"),
    [
        ("year", [1, 1, 2, 2]),
        ("month", [1, 2, 1, 2]),
        ("day", [1, 1, 1, 1]),
        ("hour", [0, 0, 0, 0]),
        ("minute", [0, 0, 0, 0]),
        ("second", [0, 0, 0, 0]),
        ("microsecond", [0, 0, 0, 0]),
    ],
)
def test_cftimeindex_field_accessors(index, field, expected):
    result = getattr(index, field)
    expected = np.array(expected, dtype=np.int64)
    assert_array_equal(result, expected)
    assert result.dtype == expected.dtype


@requires_cftime
@pytest.mark.parametrize(
    ("field"),
    [
        "year",
        "month",
        "day",
        "hour",
        "minute",
        "second",
        "microsecond",
        "dayofyear",
        "dayofweek",
        "days_in_month",
    ],
)
def test_empty_cftimeindex_field_accessors(field):
    index = CFTimeIndex([])
    result = getattr(index, field)
    expected = np.array([], dtype=np.int64)
    assert_array_equal(result, expected)
    assert result.dtype == expected.dtype


@requires_cftime
def test_cftimeindex_dayofyear_accessor(index):
    result = index.dayofyear
    expected = np.array([date.dayofyr for date in index], dtype=np.int64)
    assert_array_equal(result, expected)
    assert result.dtype == expected.dtype


@requires_cftime
def test_cftimeindex_dayofweek_accessor(index):
    result = index.dayofweek
    expected = np.array([date.dayofwk for date in index], dtype=np.int64)
    assert_array_equal(result, expected)
    assert result.dtype == expected.dtype


@requires_cftime
def test_cftimeindex_days_in_month_accessor(index):
    result = index.days_in_month
    expected = np.array([date.daysinmonth for date in index], dtype=np.int64)
    assert_array_equal(result, expected)
    assert result.dtype == expected.dtype


@requires_cftime
@pytest.mark.parametrize(
    ("string", "date_args", "reso"),
    [
        ("1999", (1999, 1, 1), "year"),
        ("199902", (1999, 2, 1), "month"),
        ("19990202", (1999, 2, 2), "day"),
        ("19990202T01", (1999, 2, 2, 1), "hour"),
        ("19990202T0101", (1999, 2, 2, 1, 1), "minute"),
        ("19990202T010156", (1999, 2, 2, 1, 1, 56), "second"),
        ("19990202T010156.123456", (1999, 2, 2, 1, 1, 56, 123456), "microsecond"),
    ],
)
def test_parse_iso8601_with_reso(date_type, string, date_args, reso):
    expected_date = date_type(*date_args)
    expected_reso = reso
    result_date, result_reso = _parse_iso8601(date_type, string)
    assert result_date == expected_date
    assert result_reso == expected_reso


@requires_cftime
def test_parse_string_to_bounds_year(date_type, dec_days):
    parsed = date_type(2, 2, 10, 6, 2, 8, 1)
    expected_start = date_type(2, 1, 1)
    expected_end = date_type(2, 12, dec_days, 23, 59, 59, 999999)
    result_start, result_end = _parsed_string_to_bounds(date_type, "year", parsed)
    assert result_start == expected_start
    assert result_end == expected_end


@requires_cftime
def test_parse_string_to_bounds_month_feb(date_type, feb_days):
    parsed = date_type(2, 2, 10, 6, 2, 8, 1)
    expected_start = date_type(2, 2, 1)
    expected_end = date_type(2, 2, feb_days, 23, 59, 59, 999999)
    result_start, result_end = _parsed_string_to_bounds(date_type, "month", parsed)
    assert result_start == expected_start
    assert result_end == expected_end


@requires_cftime
def test_parse_string_to_bounds_month_dec(date_type, dec_days):
    parsed = date_type(2, 12, 1)
    expected_start = date_type(2, 12, 1)
    expected_end = date_type(2, 12, dec_days, 23, 59, 59, 999999)
    result_start, result_end = _parsed_string_to_bounds(date_type, "month", parsed)
    assert result_start == expected_start
    assert result_end == expected_end


@requires_cftime
@pytest.mark.parametrize(
    ("reso", "ex_start_args", "ex_end_args"),
    [
        ("day", (2, 2, 10), (2, 2, 10, 23, 59, 59, 999999)),
        ("hour", (2, 2, 10, 6), (2, 2, 10, 6, 59, 59, 999999)),
        ("minute", (2, 2, 10, 6, 2), (2, 2, 10, 6, 2, 59, 999999)),
        ("second", (2, 2, 10, 6, 2, 8), (2, 2, 10, 6, 2, 8, 999999)),
    ],
)
def test_parsed_string_to_bounds_sub_monthly(
    date_type, reso, ex_start_args, ex_end_args
):
    parsed = date_type(2, 2, 10, 6, 2, 8, 123456)
    expected_start = date_type(*ex_start_args)
    expected_end = date_type(*ex_end_args)

    result_start, result_end = _parsed_string_to_bounds(date_type, reso, parsed)
    assert result_start == expected_start
    assert result_end == expected_end


@requires_cftime
def test_parsed_string_to_bounds_raises(date_type):
    with pytest.raises(KeyError):
        _parsed_string_to_bounds(date_type, "a", date_type(1, 1, 1))


@requires_cftime
def test_get_loc(date_type, index):
    result = index.get_loc("0001")
    assert result == slice(0, 2)

    result = index.get_loc(date_type(1, 2, 1))
    assert result == 1

    result = index.get_loc("0001-02-01")
    assert result == slice(1, 2)

    with pytest.raises(KeyError, match=r"1234"):
        index.get_loc("1234")


@requires_cftime
def test_get_slice_bound(date_type, index):
    result = index.get_slice_bound("0001", "left")
    expected = 0
    assert result == expected

    result = index.get_slice_bound("0001", "right")
    expected = 2
    assert result == expected

    result = index.get_slice_bound(date_type(1, 3, 1), "left")
    expected = 2
    assert result == expected

    result = index.get_slice_bound(date_type(1, 3, 1), "right")
    expected = 2
    assert result == expected


@requires_cftime
def test_get_slice_bound_decreasing_index(date_type, monotonic_decreasing_index):
    result = monotonic_decreasing_index.get_slice_bound("0001", "left")
    expected = 2
    assert result == expected

    result = monotonic_decreasing_index.get_slice_bound("0001", "right")
    expected = 4
    assert result == expected

    result = monotonic_decreasing_index.get_slice_bound(date_type(1, 3, 1), "left")
    expected = 2
    assert result == expected

    result = monotonic_decreasing_index.get_slice_bound(date_type(1, 3, 1), "right")
    expected = 2
    assert result == expected


@requires_cftime
def test_get_slice_bound_length_one_index(date_type, length_one_index):
    result = length_one_index.get_slice_bound("0001", "left")
    expected = 0
    assert result == expected

    result = length_one_index.get_slice_bound("0001", "right")
    expected = 1
    assert result == expected

    result = length_one_index.get_slice_bound(date_type(1, 3, 1), "left")
    expected = 1
    assert result == expected

    result = length_one_index.get_slice_bound(date_type(1, 3, 1), "right")
    expected = 1
    assert result == expected


@requires_cftime
def test_string_slice_length_one_index(length_one_index):
    da = xr.DataArray([1], coords=[length_one_index], dims=["time"])
    result = da.sel(time=slice("0001", "0001"))
    assert_identical(result, da)


@requires_cftime
def test_date_type_property(date_type, index):
    assert index.date_type is date_type


@requires_cftime
def test_contains(date_type, index):
    assert "0001-01-01" in index
    assert "0001" in index
    assert "0003" not in index
    assert date_type(1, 1, 1) in index
    assert date_type(3, 1, 1) not in index


@requires_cftime
def test_groupby(da):
    result = da.groupby("time.month").sum("time")
    expected = xr.DataArray([4, 6], coords=[[1, 2]], dims=["month"])
    assert_identical(result, expected)


SEL_STRING_OR_LIST_TESTS = {
    "string": "0001",
    "string-slice": slice("0001-01-01", "0001-12-30"),
    "bool-list": [True, True, False, False],
}


@requires_cftime
@pytest.mark.parametrize(
    "sel_arg",
    list(SEL_STRING_OR_LIST_TESTS.values()),
    ids=list(SEL_STRING_OR_LIST_TESTS.keys()),
)
def test_sel_string_or_list(da, index, sel_arg):
    expected = xr.DataArray([1, 2], coords=[index[:2]], dims=["time"])
    result = da.sel(time=sel_arg)
    assert_identical(result, expected)


@requires_cftime
def test_sel_date_slice_or_list(da, index, date_type):
    expected = xr.DataArray([1, 2], coords=[index[:2]], dims=["time"])
    result = da.sel(time=slice(date_type(1, 1, 1), date_type(1, 12, 30)))
    assert_identical(result, expected)

    result = da.sel(time=[date_type(1, 1, 1), date_type(1, 2, 1)])
    assert_identical(result, expected)


@requires_cftime
def test_sel_date_scalar(da, date_type, index):
    expected = xr.DataArray(1).assign_coords(time=index[0])
    result = da.sel(time=date_type(1, 1, 1))
    assert_identical(result, expected)


@requires_cftime
def test_sel_date_distant_date(da, date_type, index):
    expected = xr.DataArray(4).assign_coords(time=index[3])
    result = da.sel(time=date_type(2000, 1, 1), method="nearest")
    assert_identical(result, expected)


@requires_cftime
@pytest.mark.parametrize(
    "sel_kwargs",
    [
        {"method": "nearest"},
        {"method": "nearest", "tolerance": timedelta(days=70)},
        {"method": "nearest", "tolerance": timedelta(days=1800000)},
    ],
)
def test_sel_date_scalar_nearest(da, date_type, index, sel_kwargs):
    expected = xr.DataArray(2).assign_coords(time=index[1])
    result = da.sel(time=date_type(1, 4, 1), **sel_kwargs)
    assert_identical(result, expected)

    expected = xr.DataArray(3).assign_coords(time=index[2])
    result = da.sel(time=date_type(1, 11, 1), **sel_kwargs)
    assert_identical(result, expected)


@requires_cftime
@pytest.mark.parametrize(
    "sel_kwargs",
    [{"method": "pad"}, {"method": "pad", "tolerance": timedelta(days=365)}],
)
def test_sel_date_scalar_pad(da, date_type, index, sel_kwargs):
    expected = xr.DataArray(2).assign_coords(time=index[1])
    result = da.sel(time=date_type(1, 4, 1), **sel_kwargs)
    assert_identical(result, expected)

    expected = xr.DataArray(2).assign_coords(time=index[1])
    result = da.sel(time=date_type(1, 11, 1), **sel_kwargs)
    assert_identical(result, expected)


@requires_cftime
@pytest.mark.parametrize(
    "sel_kwargs",
    [{"method": "backfill"}, {"method": "backfill", "tolerance": timedelta(days=365)}],
)
def test_sel_date_scalar_backfill(da, date_type, index, sel_kwargs):
    expected = xr.DataArray(3).assign_coords(time=index[2])
    result = da.sel(time=date_type(1, 4, 1), **sel_kwargs)
    assert_identical(result, expected)

    expected = xr.DataArray(3).assign_coords(time=index[2])
    result = da.sel(time=date_type(1, 11, 1), **sel_kwargs)
    assert_identical(result, expected)


@requires_cftime
@pytest.mark.parametrize(
    "sel_kwargs",
    [
        {"method": "pad", "tolerance": timedelta(days=20)},
        {"method": "backfill", "tolerance": timedelta(days=20)},
        {"method": "nearest", "tolerance": timedelta(days=20)},
    ],
)
def test_sel_date_scalar_tolerance_raises(da, date_type, sel_kwargs):
    with pytest.raises(KeyError):
        da.sel(time=date_type(1, 5, 1), **sel_kwargs)


@requires_cftime
@pytest.mark.parametrize(
    "sel_kwargs",
    [{"method": "nearest"}, {"method": "nearest", "tolerance": timedelta(days=70)}],
)
def test_sel_date_list_nearest(da, date_type, index, sel_kwargs):
    expected = xr.DataArray([2, 2], coords=[[index[1], index[1]]], dims=["time"])
    result = da.sel(time=[date_type(1, 3, 1), date_type(1, 4, 1)], **sel_kwargs)
    assert_identical(result, expected)

    expected = xr.DataArray([2, 3], coords=[[index[1], index[2]]], dims=["time"])
    result = da.sel(time=[date_type(1, 3, 1), date_type(1, 12, 1)], **sel_kwargs)
    assert_identical(result, expected)

    expected = xr.DataArray([3, 3], coords=[[index[2], index[2]]], dims=["time"])
    result = da.sel(time=[date_type(1, 11, 1), date_type(1, 12, 1)], **sel_kwargs)
    assert_identical(result, expected)


@requires_cftime
@pytest.mark.parametrize(
    "sel_kwargs",
    [{"method": "pad"}, {"method": "pad", "tolerance": timedelta(days=365)}],
)
def test_sel_date_list_pad(da, date_type, index, sel_kwargs):
    expected = xr.DataArray([2, 2], coords=[[index[1], index[1]]], dims=["time"])
    result = da.sel(time=[date_type(1, 3, 1), date_type(1, 4, 1)], **sel_kwargs)
    assert_identical(result, expected)


@requires_cftime
@pytest.mark.parametrize(
    "sel_kwargs",
    [{"method": "backfill"}, {"method": "backfill", "tolerance": timedelta(days=365)}],
)
def test_sel_date_list_backfill(da, date_type, index, sel_kwargs):
    expected = xr.DataArray([3, 3], coords=[[index[2], index[2]]], dims=["time"])
    result = da.sel(time=[date_type(1, 3, 1), date_type(1, 4, 1)], **sel_kwargs)
    assert_identical(result, expected)


@requires_cftime
@pytest.mark.parametrize(
    "sel_kwargs",
    [
        {"method": "pad", "tolerance": timedelta(days=20)},
        {"method": "backfill", "tolerance": timedelta(days=20)},
        {"method": "nearest", "tolerance": timedelta(days=20)},
    ],
)
def test_sel_date_list_tolerance_raises(da, date_type, sel_kwargs):
    with pytest.raises(KeyError):
        da.sel(time=[date_type(1, 2, 1), date_type(1, 5, 1)], **sel_kwargs)


@requires_cftime
def test_isel(da, index):
    expected = xr.DataArray(1).assign_coords(time=index[0])
    result = da.isel(time=0)
    assert_identical(result, expected)

    expected = xr.DataArray([1, 2], coords=[index[:2]], dims=["time"])
    result = da.isel(time=[0, 1])
    assert_identical(result, expected)


@pytest.fixture
def scalar_args(date_type):
    return [date_type(1, 1, 1)]


@pytest.fixture
def range_args(date_type):
    return [
        "0001",
        slice("0001-01-01", "0001-12-30"),
        slice(None, "0001-12-30"),
        slice(date_type(1, 1, 1), date_type(1, 12, 30)),
        slice(None, date_type(1, 12, 30)),
    ]


@requires_cftime
def test_indexing_in_series_getitem(series, index, scalar_args, range_args):
    for arg in scalar_args:
        assert series[arg] == 1

    expected = pd.Series([1, 2], index=index[:2])
    for arg in range_args:
        assert series[arg].equals(expected)


@requires_cftime
def test_indexing_in_series_loc(series, index, scalar_args, range_args):
    for arg in scalar_args:
        assert series.loc[arg] == 1

    expected = pd.Series([1, 2], index=index[:2])
    for arg in range_args:
        assert series.loc[arg].equals(expected)


@requires_cftime
def test_indexing_in_series_iloc(series, index):
    expected1 = 1
    assert series.iloc[0] == expected1

    expected2 = pd.Series([1, 2], index=index[:2])
    assert series.iloc[:2].equals(expected2)


@requires_cftime
def test_series_dropna(index):
    series = pd.Series([0.0, 1.0, np.nan, np.nan], index=index)
    expected = series.iloc[:2]
    result = series.dropna()
    assert result.equals(expected)


@requires_cftime
def test_indexing_in_dataframe_loc(df, index, scalar_args, range_args):
    expected_s = pd.Series([1], name=index[0])
    for arg in scalar_args:
        result_s = df.loc[arg]
        assert result_s.equals(expected_s)

    expected_df = pd.DataFrame([1, 2], index=index[:2])
    for arg in range_args:
        result_df = df.loc[arg]
        assert result_df.equals(expected_df)


@requires_cftime
def test_indexing_in_dataframe_iloc(df, index):
    expected_s = pd.Series([1], name=index[0])
    result_s = df.iloc[0]
    assert result_s.equals(expected_s)
    assert result_s.equals(expected_s)

    expected_df = pd.DataFrame([1, 2], index=index[:2])
    result_df = df.iloc[:2]
    assert result_df.equals(expected_df)


@requires_cftime
def test_concat_cftimeindex(date_type):
    da1 = xr.DataArray(
        [1.0, 2.0], coords=[[date_type(1, 1, 1), date_type(1, 2, 1)]], dims=["time"]
    )
    da2 = xr.DataArray(
        [3.0, 4.0], coords=[[date_type(1, 3, 1), date_type(1, 4, 1)]], dims=["time"]
    )
    da = xr.concat([da1, da2], dim="time")

    assert isinstance(da.xindexes["time"].to_pandas_index(), CFTimeIndex)


@requires_cftime
def test_empty_cftimeindex():
    index = CFTimeIndex([])
    assert index.date_type is None


@requires_cftime
def test_cftimeindex_add(index):
    date_type = index.date_type
    expected_dates = [
        date_type(1, 1, 2),
        date_type(1, 2, 2),
        date_type(2, 1, 2),
        date_type(2, 2, 2),
    ]
    expected = CFTimeIndex(expected_dates)
    result = index + timedelta(days=1)
    assert result.equals(expected)
    assert isinstance(result, CFTimeIndex)


@requires_cftime
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_cftimeindex_add_timedeltaindex(calendar) -> None:
    a = xr.date_range("2000", periods=5, calendar=calendar, use_cftime=True)
    deltas = pd.TimedeltaIndex([timedelta(days=2) for _ in range(5)])
    result = a + deltas
    expected = a.shift(2, "D")
    assert result.equals(expected)
    assert isinstance(result, CFTimeIndex)


@requires_cftime
@pytest.mark.parametrize("n", [2.0, 1.5])
@pytest.mark.parametrize(
    "freq,units",
    [
        ("h", "h"),
        ("min", "min"),
        ("s", "s"),
        ("ms", "ms"),
    ],
)
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_cftimeindex_shift_float(n, freq, units, calendar) -> None:
    a = xr.date_range("2000", periods=3, calendar=calendar, freq="D", use_cftime=True)
    result = a + pd.Timedelta(n, units)
    expected = a.shift(n, freq)
    assert result.equals(expected)
    assert isinstance(result, CFTimeIndex)


@requires_cftime
def test_cftimeindex_shift_float_us() -> None:
    a = xr.date_range("2000", periods=3, freq="D", use_cftime=True)
    with pytest.raises(
        ValueError, match="Could not convert to integer offset at any resolution"
    ):
        a.shift(2.5, "us")


@requires_cftime
@pytest.mark.parametrize("freq", ["YS", "YE", "QS", "QE", "MS", "ME", "D"])
def test_cftimeindex_shift_float_fails_for_non_tick_freqs(freq) -> None:
    a = xr.date_range("2000", periods=3, freq="D", use_cftime=True)
    with pytest.raises(TypeError, match="unsupported operand type"):
        a.shift(2.5, freq)


@requires_cftime
def test_cftimeindex_radd(index):
    date_type = index.date_type
    expected_dates = [
        date_type(1, 1, 2),
        date_type(1, 2, 2),
        date_type(2, 1, 2),
        date_type(2, 2, 2),
    ]
    expected = CFTimeIndex(expected_dates)
    result = timedelta(days=1) + index
    assert result.equals(expected)
    assert isinstance(result, CFTimeIndex)


@requires_cftime
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_timedeltaindex_add_cftimeindex(calendar) -> None:
    a = xr.date_range("2000", periods=5, calendar=calendar, use_cftime=True)
    deltas = pd.TimedeltaIndex([timedelta(days=2) for _ in range(5)])
    result = deltas + a
    expected = a.shift(2, "D")
    assert result.equals(expected)
    assert isinstance(result, CFTimeIndex)


@requires_cftime
def test_cftimeindex_sub_timedelta(index):
    date_type = index.date_type
    expected_dates = [
        date_type(1, 1, 2),
        date_type(1, 2, 2),
        date_type(2, 1, 2),
        date_type(2, 2, 2),
    ]
    expected = CFTimeIndex(expected_dates)
    result = index + timedelta(days=2)
    result = result - timedelta(days=1)
    assert result.equals(expected)
    assert isinstance(result, CFTimeIndex)


@requires_cftime
@pytest.mark.parametrize(
    "other",
    [np.array(4 * [timedelta(days=1)]), np.array(timedelta(days=1))],
    ids=["1d-array", "scalar-array"],
)
def test_cftimeindex_sub_timedelta_array(index, other):
    date_type = index.date_type
    expected_dates = [
        date_type(1, 1, 2),
        date_type(1, 2, 2),
        date_type(2, 1, 2),
        date_type(2, 2, 2),
    ]
    expected = CFTimeIndex(expected_dates)
    result = index + timedelta(days=2)
    result = result - other
    assert result.equals(expected)
    assert isinstance(result, CFTimeIndex)


@requires_cftime
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_cftimeindex_sub_cftimeindex(calendar) -> None:
    a = xr.date_range("2000", periods=5, calendar=calendar, use_cftime=True)
    b = a.shift(2, "D")
    result = b - a
    expected = pd.TimedeltaIndex([timedelta(days=2) for _ in range(5)])
    assert result.equals(expected)
    assert isinstance(result, pd.TimedeltaIndex)


@requires_cftime
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_cftimeindex_sub_cftime_datetime(calendar):
    a = xr.date_range("2000", periods=5, calendar=calendar, use_cftime=True)
    result = a - a[0]
    expected = pd.TimedeltaIndex([timedelta(days=i) for i in range(5)])
    assert result.equals(expected)
    assert isinstance(result, pd.TimedeltaIndex)


@requires_cftime
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_cftime_datetime_sub_cftimeindex(calendar):
    a = xr.date_range("2000", periods=5, calendar=calendar, use_cftime=True)
    result = a[0] - a
    expected = pd.TimedeltaIndex([timedelta(days=-i) for i in range(5)])
    assert result.equals(expected)
    assert isinstance(result, pd.TimedeltaIndex)


@requires_cftime
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_distant_cftime_datetime_sub_cftimeindex(calendar):
    a = xr.date_range("2000", periods=5, calendar=calendar, use_cftime=True)
    if not has_pandas_3:
        with pytest.raises(ValueError, match="difference exceeds"):
            a.date_type(1, 1, 1) - a
    else:
        result = a.date_type(1, 1, 1) - a
        assert isinstance(result, pd.TimedeltaIndex)
        assert result.unit == "us"

        # Check that we can recover original index from subtracting timedeltas
        roundtrip = CFTimeIndex(a.date_type(1, 1, 1) - result.to_pytimedelta())
        assert roundtrip.equals(a)


@requires_cftime
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_cftimeindex_sub_timedeltaindex(calendar) -> None:
    a = xr.date_range("2000", periods=5, calendar=calendar, use_cftime=True)
    deltas = pd.TimedeltaIndex([timedelta(days=2) for _ in range(5)])
    result = a - deltas
    expected = a.shift(-2, "D")
    assert result.equals(expected)
    assert isinstance(result, CFTimeIndex)


@requires_cftime
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_cftimeindex_sub_index_of_cftime_datetimes(calendar):
    a = xr.date_range("2000", periods=5, calendar=calendar, use_cftime=True)
    b = pd.Index(a.values)
    expected = a - a
    result = a - b
    assert result.equals(expected)
    assert isinstance(result, pd.TimedeltaIndex)


@requires_cftime
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_cftimeindex_sub_not_implemented(calendar):
    a = xr.date_range("2000", periods=5, calendar=calendar, use_cftime=True)
    with pytest.raises(TypeError, match="unsupported operand"):
        a - 1


@requires_cftime
def test_cftimeindex_rsub(index):
    with pytest.raises(TypeError):
        timedelta(days=1) - index


@requires_cftime
@pytest.mark.parametrize("freq", ["D", timedelta(days=1)])
def test_cftimeindex_shift(index, freq) -> None:
    date_type = index.date_type
    expected_dates = [
        date_type(1, 1, 3),
        date_type(1, 2, 3),
        date_type(2, 1, 3),
        date_type(2, 2, 3),
    ]
    expected = CFTimeIndex(expected_dates)
    result = index.shift(2, freq)
    assert result.equals(expected)
    assert isinstance(result, CFTimeIndex)


@requires_cftime
def test_cftimeindex_shift_invalid_periods() -> None:
    index = xr.date_range("2000", periods=3, use_cftime=True)
    with pytest.raises(TypeError):
        index.shift("a", "D")


@requires_cftime
def test_cftimeindex_shift_invalid_freq() -> None:
    index = xr.date_range("2000", periods=3, use_cftime=True)
    with pytest.raises(TypeError):
        index.shift(1, 1)


@requires_cftime
@pytest.mark.parametrize(
    ("calendar", "expected"),
    [
        ("noleap", "noleap"),
        ("365_day", "noleap"),
        ("360_day", "360_day"),
        ("julian", "julian"),
        ("gregorian", standard_or_gregorian),
        ("standard", standard_or_gregorian),
        ("proleptic_gregorian", "proleptic_gregorian"),
    ],
)
def test_cftimeindex_calendar_property(calendar, expected):
    index = xr.date_range(start="2000", periods=3, calendar=calendar, use_cftime=True)
    assert index.calendar == expected


@requires_cftime
def test_empty_cftimeindex_calendar_property():
    index = CFTimeIndex([])
    assert index.calendar is None


@requires_cftime
@pytest.mark.parametrize(
    "calendar",
    [
        "noleap",
        "365_day",
        "360_day",
        "julian",
        "gregorian",
        "standard",
        "proleptic_gregorian",
    ],
)
def test_cftimeindex_freq_property_none_size_lt_3(calendar):
    for periods in range(3):
        index = xr.date_range(
            start="2000", periods=periods, calendar=calendar, use_cftime=True
        )
        assert index.freq is None


@requires_cftime
@pytest.mark.parametrize(
    ("calendar", "expected"),
    [
        ("noleap", "noleap"),
        ("365_day", "noleap"),
        ("360_day", "360_day"),
        ("julian", "julian"),
        ("gregorian", standard_or_gregorian),
        ("standard", standard_or_gregorian),
        ("proleptic_gregorian", "proleptic_gregorian"),
    ],
)
def test_cftimeindex_calendar_repr(calendar, expected):
    """Test that cftimeindex has calendar property in repr."""
    index = xr.date_range(start="2000", periods=3, calendar=calendar, use_cftime=True)
    repr_str = index.__repr__()
    assert f" calendar='{expected}'" in repr_str
    assert "2000-01-01 00:00:00, 2000-01-02 00:00:00" in repr_str


@requires_cftime
@pytest.mark.parametrize("periods", [2, 40])
def test_cftimeindex_periods_repr(periods):
    """Test that cftimeindex has periods property in repr."""
    index = xr.date_range(start="2000", periods=periods, use_cftime=True)
    repr_str = index.__repr__()
    assert f" length={periods}" in repr_str


@requires_cftime
@pytest.mark.parametrize("calendar", ["noleap", "360_day", "standard"])
@pytest.mark.parametrize("freq", ["D", "h"])
def test_cftimeindex_freq_in_repr(freq, calendar):
    """Test that cftimeindex has frequency property in repr."""
    index = xr.date_range(
        start="2000", periods=3, freq=freq, calendar=calendar, use_cftime=True
    )
    repr_str = index.__repr__()
    assert f", freq='{freq}'" in repr_str


@requires_cftime
@pytest.mark.parametrize(
    "periods,expected",
    [
        (
            2,
            f"""\
CFTimeIndex([2000-01-01 00:00:00, 2000-01-02 00:00:00],
            dtype='object', length=2, calendar='{standard_or_gregorian}', freq=None)""",
        ),
        (
            4,
            f"""\
CFTimeIndex([2000-01-01 00:00:00, 2000-01-02 00:00:00, 2000-01-03 00:00:00,
             2000-01-04 00:00:00],
            dtype='object', length=4, calendar='{standard_or_gregorian}', freq='D')""",
        ),
        (
            101,
            f"""\
CFTimeIndex([2000-01-01 00:00:00, 2000-01-02 00:00:00, 2000-01-03 00:00:00,
             2000-01-04 00:00:00, 2000-01-05 00:00:00, 2000-01-06 00:00:00,
             2000-01-07 00:00:00, 2000-01-08 00:00:00, 2000-01-09 00:00:00,
             2000-01-10 00:00:00,
             ...
             2000-04-01 00:00:00, 2000-04-02 00:00:00, 2000-04-03 00:00:00,
             2000-04-04 00:00:00, 2000-04-05 00:00:00, 2000-04-06 00:00:00,
             2000-04-07 00:00:00, 2000-04-08 00:00:00, 2000-04-09 00:00:00,
             2000-04-10 00:00:00],
            dtype='object', length=101, calendar='{standard_or_gregorian}', freq='D')""",
        ),
    ],
)
def test_cftimeindex_repr_formatting(periods, expected):
    """Test that cftimeindex.__repr__ is formatted similar to pd.Index.__repr__."""
    index = xr.date_range(start="2000", periods=periods, freq="D", use_cftime=True)
    expected = dedent(expected)
    assert expected == repr(index)


@requires_cftime
@pytest.mark.parametrize("display_width", [40, 80, 100])
@pytest.mark.parametrize("periods", [2, 3, 4, 100, 101])
def test_cftimeindex_repr_formatting_width(periods, display_width):
    """Test that cftimeindex is sensitive to OPTIONS['display_width']."""
    index = xr.date_range(start="2000", periods=periods, use_cftime=True)
    len_intro_str = len("CFTimeIndex(")
    with xr.set_options(display_width=display_width):
        repr_str = index.__repr__()
        splitted = repr_str.split("\n")
        for i, s in enumerate(splitted):
            # check that lines not longer than OPTIONS['display_width']
            assert len(s) <= display_width, f"{len(s)} {s} {display_width}"
            if i > 0:
                # check for initial spaces
                assert s[:len_intro_str] == " " * len_intro_str


@requires_cftime
@pytest.mark.parametrize("periods", [22, 50, 100])
def test_cftimeindex_repr_101_shorter(periods):
    index_101 = xr.date_range(start="2000", periods=101, use_cftime=True)
    index_periods = xr.date_range(start="2000", periods=periods, use_cftime=True)
    index_101_repr_str = index_101.__repr__()
    index_periods_repr_str = index_periods.__repr__()
    assert len(index_101_repr_str) < len(index_periods_repr_str)


@requires_cftime
def test_parse_array_of_cftime_strings():
    from cftime import DatetimeNoLeap

    strings = np.array([["2000-01-01", "2000-01-02"], ["2000-01-03", "2000-01-04"]])
    expected = np.array(
        [
            [DatetimeNoLeap(2000, 1, 1), DatetimeNoLeap(2000, 1, 2)],
            [DatetimeNoLeap(2000, 1, 3), DatetimeNoLeap(2000, 1, 4)],
        ]
    )

    result = _parse_array_of_cftime_strings(strings, DatetimeNoLeap)
    np.testing.assert_array_equal(result, expected)

    # Test scalar array case
    strings = np.array("2000-01-01")
    expected = np.array(DatetimeNoLeap(2000, 1, 1))
    result = _parse_array_of_cftime_strings(strings, DatetimeNoLeap)
    np.testing.assert_array_equal(result, expected)


@requires_cftime
@pytest.mark.parametrize("calendar", _ALL_CALENDARS)
def test_strftime_of_cftime_array(calendar):
    date_format = "%Y%m%d%H%M"
    cf_values = xr.date_range("2000", periods=5, calendar=calendar, use_cftime=True)
    dt_values = pd.date_range("2000", periods=5)
    expected = pd.Index(dt_values.strftime(date_format))
    result = cf_values.strftime(date_format)
    assert result.equals(expected)


@requires_cftime
@pytest.mark.parametrize("calendar", _ALL_CALENDARS)
@pytest.mark.parametrize("unsafe", [False, True])
def test_to_datetimeindex(calendar, unsafe) -> None:
    index = xr.date_range("2000", periods=5, calendar=calendar, use_cftime=True)
    expected = pd.date_range("2000", periods=5, unit="ns")

    if calendar in _NON_STANDARD_CALENDAR_NAMES and not unsafe:
        with pytest.warns(RuntimeWarning, match="non-standard"):
            result = index.to_datetimeindex(time_unit="ns")
    else:
        result = index.to_datetimeindex(unsafe=unsafe, time_unit="ns")

    assert result.equals(expected)
    np.testing.assert_array_equal(result, expected)
    assert isinstance(result, pd.DatetimeIndex)


@requires_cftime
def test_to_datetimeindex_future_warning() -> None:
    index = xr.date_range("2000", periods=5, use_cftime=True)
    expected = pd.date_range("2000", periods=5, unit="ns")
    with pytest.warns(FutureWarning, match="In a future version"):
        result = index.to_datetimeindex()
    assert result.equals(expected)
    assert result.dtype == expected.dtype


@requires_cftime
@pytest.mark.parametrize("calendar", _ALL_CALENDARS)
def test_to_datetimeindex_out_of_range(calendar) -> None:
    index = xr.date_range("0001", periods=5, calendar=calendar, use_cftime=True)
    with pytest.raises(ValueError, match="0001"):
        index.to_datetimeindex(time_unit="ns")


@requires_cftime
@pytest.mark.parametrize("unsafe", [False, True])
def test_to_datetimeindex_gregorian_pre_reform(unsafe) -> None:
    index = xr.date_range("1582", periods=5, calendar="gregorian", use_cftime=True)
    if unsafe:
        result = index.to_datetimeindex(time_unit="us", unsafe=unsafe)
    else:
        with pytest.warns(RuntimeWarning, match="reform"):
            result = index.to_datetimeindex(time_unit="us", unsafe=unsafe)
    expected = pd.date_range("1582", periods=5, unit="us")
    assert result.equals(expected)
    assert result.dtype == expected.dtype


@requires_cftime
@pytest.mark.parametrize("calendar", ["all_leap", "360_day"])
def test_to_datetimeindex_feb_29(calendar) -> None:
    index = xr.date_range("2001-02-28", periods=2, calendar=calendar, use_cftime=True)
    with pytest.raises(ValueError, match="29"):
        index.to_datetimeindex(time_unit="ns")


@pytest.mark.xfail(reason="fails on pandas main branch")
@requires_cftime
def test_multiindex():
    index = xr.date_range(
        "2001-01-01", periods=100, calendar="360_day", use_cftime=True
    )
    mindex = pd.MultiIndex.from_arrays([index])
    assert mindex.get_loc("2001-01") == slice(0, 30)


@requires_cftime
@pytest.mark.parametrize("freq", ["3663s", "33min", "2h"])
@pytest.mark.parametrize("method", ["floor", "ceil", "round"])
def test_rounding_methods_against_datetimeindex(freq, method) -> None:
    # for now unit="us" seems good enough
    expected = pd.date_range("2000-01-02T01:03:51", periods=10, freq="1777s", unit="ns")
    expected = getattr(expected, method)(freq)
    result = xr.date_range(
        "2000-01-02T01:03:51", periods=10, freq="1777s", use_cftime=True
    )
    result = getattr(result, method)(freq).to_datetimeindex(time_unit="ns")
    assert result.equals(expected)


@requires_cftime
@pytest.mark.parametrize("method", ["floor", "ceil", "round"])
def test_rounding_methods_empty_cftimindex(method):
    index = CFTimeIndex([])
    result = getattr(index, method)("2s")

    expected = CFTimeIndex([])

    assert result.equals(expected)
    assert result is not index


@requires_cftime
@pytest.mark.parametrize("method", ["floor", "ceil", "round"])
def test_rounding_methods_invalid_freq(method):
    index = xr.date_range(
        "2000-01-02T01:03:51", periods=10, freq="1777s", use_cftime=True
    )
    with pytest.raises(ValueError, match="fixed"):
        getattr(index, method)("MS")


@pytest.fixture
def rounding_index(date_type):
    return xr.CFTimeIndex(
        [
            date_type(1, 1, 1, 1, 59, 59, 999512),
            date_type(1, 1, 1, 3, 0, 1, 500001),
            date_type(1, 1, 1, 7, 0, 6, 499999),
        ]
    )


@requires_cftime
def test_ceil(rounding_index, date_type):
    result = rounding_index.ceil("s")
    expected = xr.CFTimeIndex(
        [
            date_type(1, 1, 1, 2, 0, 0, 0),
            date_type(1, 1, 1, 3, 0, 2, 0),
            date_type(1, 1, 1, 7, 0, 7, 0),
        ]
    )
    assert result.equals(expected)


@requires_cftime
def test_floor(rounding_index, date_type):
    result = rounding_index.floor("s")
    expected = xr.CFTimeIndex(
        [
            date_type(1, 1, 1, 1, 59, 59, 0),
            date_type(1, 1, 1, 3, 0, 1, 0),
            date_type(1, 1, 1, 7, 0, 6, 0),
        ]
    )
    assert result.equals(expected)


@requires_cftime
def test_round(rounding_index, date_type):
    result = rounding_index.round("s")
    expected = xr.CFTimeIndex(
        [
            date_type(1, 1, 1, 2, 0, 0, 0),
            date_type(1, 1, 1, 3, 0, 2, 0),
            date_type(1, 1, 1, 7, 0, 6, 0),
        ]
    )
    assert result.equals(expected)


@requires_cftime
def test_asi8(date_type):
    index = xr.CFTimeIndex([date_type(1970, 1, 1), date_type(1970, 1, 2)])
    result = index.asi8
    expected = 1000000 * 86400 * np.array([0, 1])
    np.testing.assert_array_equal(result, expected)


@requires_cftime
def test_asi8_distant_date():
    """Test that asi8 conversion is truly exact."""
    import cftime

    date_type = cftime.DatetimeProlepticGregorian
    index = xr.CFTimeIndex([date_type(10731, 4, 22, 3, 25, 45, 123456)])
    result = index.asi8
    expected = np.array([1000000 * 86400 * 400 * 8000 + 12345 * 1000000 + 123456])
    np.testing.assert_array_equal(result, expected)


@requires_cftime
def test_asi8_empty_cftimeindex():
    index = xr.CFTimeIndex([])
    result = index.asi8
    expected = np.array([], dtype=np.int64)
    np.testing.assert_array_equal(result, expected)


@requires_cftime
def test_infer_freq_valid_types(time_unit: PDDatetimeUnitOptions) -> None:
    cf_index = xr.date_range("2000-01-01", periods=3, freq="D", use_cftime=True)
    assert xr.infer_freq(cf_index) == "D"
    assert xr.infer_freq(xr.DataArray(cf_index)) == "D"

    pd_index = pd.date_range("2000-01-01", periods=3, freq="D").as_unit(time_unit)
    assert xr.infer_freq(pd_index) == "D"
    assert xr.infer_freq(xr.DataArray(pd_index)) == "D"

    pd_td_index = pd.timedelta_range(start="1D", periods=3, freq="D").as_unit(time_unit)
    assert xr.infer_freq(pd_td_index) == "D"
    assert xr.infer_freq(xr.DataArray(pd_td_index)) == "D"


@requires_cftime
def test_infer_freq_invalid_inputs():
    # Non-datetime DataArray
    with pytest.raises(ValueError, match="must contain datetime-like objects"):
        xr.infer_freq(xr.DataArray([0, 1, 2]))

    index = xr.date_range("1990-02-03", periods=4, freq="MS", use_cftime=True)
    # 2D DataArray
    with pytest.raises(ValueError, match="must be 1D"):
        xr.infer_freq(xr.DataArray([index, index]))

    # CFTimeIndex too short
    with pytest.raises(ValueError, match="Need at least 3 dates to infer frequency"):
        xr.infer_freq(index[:2])

    # Non-monotonic input
    assert xr.infer_freq(index[np.array([0, 2, 1, 3])]) is None

    # Non-unique input
    assert xr.infer_freq(index[np.array([0, 1, 1, 2])]) is None

    # No unique frequency (here 1st step is MS, second is 2MS)
    assert xr.infer_freq(index[np.array([0, 1, 3])]) is None

    # Same, but for QS
    index = xr.date_range("1990-02-03", periods=4, freq="QS", use_cftime=True)
    assert xr.infer_freq(index[np.array([0, 1, 3])]) is None


@requires_cftime
@pytest.mark.parametrize(
    "freq",
    [
        "300YS-JAN",
        "YE-DEC",
        "YS-JUL",
        "2YS-FEB",
        "QE-NOV",
        "3QS-DEC",
        "MS",
        "4ME",
        "7D",
        "D",
        "30h",
        "5min",
        "40s",
    ],
)
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_infer_freq(freq, calendar):
    index = xr.date_range(
        "2000-01-01", periods=3, freq=freq, calendar=calendar, use_cftime=True
    )
    out = xr.infer_freq(index)
    assert out == freq


@requires_cftime
@pytest.mark.parametrize("calendar", _CFTIME_CALENDARS)
def test_pickle_cftimeindex(calendar):
    idx = xr.date_range(
        "2000-01-01", periods=3, freq="D", calendar=calendar, use_cftime=True
    )
    idx_pkl = pickle.loads(pickle.dumps(idx))
    assert (idx == idx_pkl).all()
