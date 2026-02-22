from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import xarray as xr
from xarray.tests import (
    _CFTIME_CALENDARS,
    _all_cftime_date_types,
    assert_allclose,
    assert_array_equal,
    assert_chunks_equal,
    assert_equal,
    assert_identical,
    raise_if_dask_computes,
    requires_cftime,
    requires_dask,
)


class TestDatetimeAccessor:
    @pytest.fixture(autouse=True)
    def setup(self):
        nt = 100
        data = np.random.rand(10, 10, nt)
        lons = np.linspace(0, 11, 10)
        lats = np.linspace(0, 20, 10)
        self.times = pd.date_range(start="2000/01/01", freq="h", periods=nt)

        self.data = xr.DataArray(
            data,
            coords=[lons, lats, self.times],
            dims=["lon", "lat", "time"],
            name="data",
        )

        self.times_arr = np.random.choice(self.times, size=(10, 10, nt))
        self.times_data = xr.DataArray(
            self.times_arr,
            coords=[lons, lats, self.times],
            dims=["lon", "lat", "time"],
            name="data",
        )

    @pytest.mark.parametrize(
        "field",
        [
            "year",
            "month",
            "day",
            "hour",
            "minute",
            "second",
            "microsecond",
            "nanosecond",
            "week",
            "weekofyear",
            "dayofweek",
            "weekday",
            "dayofyear",
            "quarter",
            "date",
            "time",
            "daysinmonth",
            "days_in_month",
            "is_month_start",
            "is_month_end",
            "is_quarter_start",
            "is_quarter_end",
            "is_year_start",
            "is_year_end",
            "is_leap_year",
        ],
    )
    def test_field_access(self, field) -> None:
        if field in ["week", "weekofyear"]:
            data = self.times.isocalendar()["week"]
        else:
            data = getattr(self.times, field)

        if data.dtype.kind != "b" and field not in ("date", "time"):
            # pandas 2.0 returns int32 for integer fields now
            data = data.astype("int64")

        translations = {
            "weekday": "dayofweek",
            "daysinmonth": "days_in_month",
            "weekofyear": "week",
        }
        name = translations.get(field, field)

        expected = xr.DataArray(data, name=name, coords=[self.times], dims=["time"])

        if field in ["week", "weekofyear"]:
            with pytest.warns(
                FutureWarning, match="dt.weekofyear and dt.week have been deprecated"
            ):
                actual = getattr(self.data.time.dt, field)
        else:
            actual = getattr(self.data.time.dt, field)
        assert not isinstance(actual.variable, xr.IndexVariable)

        assert expected.dtype == actual.dtype
        assert_identical(expected, actual)

    def test_total_seconds(self) -> None:
        # Subtract a value in the middle of the range to ensure that some values
        # are negative
        delta = self.data.time - np.datetime64("2000-01-03")
        actual = delta.dt.total_seconds()
        expected = xr.DataArray(
            np.arange(-48, 52, dtype=np.float64) * 3600,
            name="total_seconds",
            coords=[self.data.time],
        )
        # This works with assert_identical when pandas is >=1.5.0.
        assert_allclose(expected, actual)

    @pytest.mark.parametrize(
        "field, pandas_field",
        [
            ("year", "year"),
            ("week", "week"),
            ("weekday", "day"),
        ],
    )
    def test_isocalendar(self, field, pandas_field) -> None:
        # pandas isocalendar has dtypy UInt32Dtype, convert to Int64
        expected = pd.Index(getattr(self.times.isocalendar(), pandas_field).astype(int))
        expected = xr.DataArray(
            expected, name=field, coords=[self.times], dims=["time"]
        )

        actual = self.data.time.dt.isocalendar()[field]
        assert_equal(expected, actual)

    def test_calendar(self) -> None:
        cal = self.data.time.dt.calendar
        assert cal == "proleptic_gregorian"

    def test_strftime(self) -> None:
        assert (
            "2000-01-01 01:00:00" == self.data.time.dt.strftime("%Y-%m-%d %H:%M:%S")[1]
        )

    @requires_cftime
    @pytest.mark.parametrize(
        "calendar,expected",
        [("standard", 366), ("noleap", 365), ("360_day", 360), ("all_leap", 366)],
    )
    def test_days_in_year(self, calendar, expected) -> None:
        assert (
            self.data.convert_calendar(calendar, align_on="year").time.dt.days_in_year
            == expected
        ).all()

    def test_not_datetime_type(self) -> None:
        nontime_data = self.data.copy()
        int_data = np.arange(len(self.data.time)).astype("int8")
        nontime_data = nontime_data.assign_coords(time=int_data)
        with pytest.raises(AttributeError, match=r"dt"):
            _ = nontime_data.time.dt

    @pytest.mark.filterwarnings("ignore:dt.weekofyear and dt.week have been deprecated")
    @requires_dask
    @pytest.mark.parametrize(
        "field",
        [
            "year",
            "month",
            "day",
            "hour",
            "minute",
            "second",
            "microsecond",
            "nanosecond",
            "week",
            "weekofyear",
            "dayofweek",
            "weekday",
            "dayofyear",
            "quarter",
            "date",
            "time",
            "is_month_start",
            "is_month_end",
            "is_quarter_start",
            "is_quarter_end",
            "is_year_start",
            "is_year_end",
            "is_leap_year",
            "days_in_year",
        ],
    )
    def test_dask_field_access(self, field) -> None:
        import dask.array as da

        expected = getattr(self.times_data.dt, field)

        dask_times_arr = da.from_array(self.times_arr, chunks=(5, 5, 50))
        dask_times_2d = xr.DataArray(
            dask_times_arr, coords=self.data.coords, dims=self.data.dims, name="data"
        )

        with raise_if_dask_computes():
            actual = getattr(dask_times_2d.dt, field)

        assert isinstance(actual.data, da.Array)
        assert_chunks_equal(actual, dask_times_2d)
        assert_equal(actual.compute(), expected.compute())

    @requires_dask
    @pytest.mark.parametrize(
        "field",
        [
            "year",
            "week",
            "weekday",
        ],
    )
    def test_isocalendar_dask(self, field) -> None:
        import dask.array as da

        expected = getattr(self.times_data.dt.isocalendar(), field)

        dask_times_arr = da.from_array(self.times_arr, chunks=(5, 5, 50))
        dask_times_2d = xr.DataArray(
            dask_times_arr, coords=self.data.coords, dims=self.data.dims, name="data"
        )

        with raise_if_dask_computes():
            actual = dask_times_2d.dt.isocalendar()[field]

        assert isinstance(actual.data, da.Array)
        assert_chunks_equal(actual, dask_times_2d)
        assert_equal(actual.compute(), expected.compute())

    @requires_dask
    @pytest.mark.parametrize(
        "method, parameters",
        [
            ("floor", "D"),
            ("ceil", "D"),
            ("round", "D"),
            ("strftime", "%Y-%m-%d %H:%M:%S"),
        ],
    )
    def test_dask_accessor_method(self, method, parameters) -> None:
        import dask.array as da

        expected = getattr(self.times_data.dt, method)(parameters)
        dask_times_arr = da.from_array(self.times_arr, chunks=(5, 5, 50))
        dask_times_2d = xr.DataArray(
            dask_times_arr, coords=self.data.coords, dims=self.data.dims, name="data"
        )

        with raise_if_dask_computes():
            actual = getattr(dask_times_2d.dt, method)(parameters)

        assert isinstance(actual.data, da.Array)
        assert_chunks_equal(actual, dask_times_2d)
        assert_equal(actual.compute(), expected.compute())

    def test_seasons(self) -> None:
        dates = xr.date_range(
            start="2000/01/01", freq="ME", periods=12, use_cftime=False
        )
        dates = dates.append(pd.Index([np.datetime64("NaT")]))
        dates = xr.DataArray(dates)
        seasons = xr.DataArray(
            [
                "DJF",
                "DJF",
                "MAM",
                "MAM",
                "MAM",
                "JJA",
                "JJA",
                "JJA",
                "SON",
                "SON",
                "SON",
                "DJF",
                "nan",
            ]
        )

        assert_array_equal(seasons.values, dates.dt.season.values)

    @pytest.mark.parametrize(
        "method, parameters", [("floor", "D"), ("ceil", "D"), ("round", "D")]
    )
    def test_accessor_method(self, method, parameters) -> None:
        dates = pd.date_range("2014-01-01", "2014-05-01", freq="h")
        xdates = xr.DataArray(dates, dims=["time"])
        expected = getattr(dates, method)(parameters)
        actual = getattr(xdates.dt, method)(parameters)
        assert_array_equal(expected, actual)


class TestTimedeltaAccessor:
    @pytest.fixture(autouse=True)
    def setup(self):
        nt = 100
        data = np.random.rand(10, 10, nt)
        lons = np.linspace(0, 11, 10)
        lats = np.linspace(0, 20, 10)
        self.times = pd.timedelta_range(start="1 day", freq="6h", periods=nt)

        self.data = xr.DataArray(
            data,
            coords=[lons, lats, self.times],
            dims=["lon", "lat", "time"],
            name="data",
        )

        self.times_arr = np.random.choice(self.times, size=(10, 10, nt))
        self.times_data = xr.DataArray(
            self.times_arr,
            coords=[lons, lats, self.times],
            dims=["lon", "lat", "time"],
            name="data",
        )

    def test_not_datetime_type(self) -> None:
        nontime_data = self.data.copy()
        int_data = np.arange(len(self.data.time)).astype("int8")
        nontime_data = nontime_data.assign_coords(time=int_data)
        with pytest.raises(AttributeError, match=r"dt"):
            _ = nontime_data.time.dt

    @pytest.mark.parametrize(
        "field", ["days", "seconds", "microseconds", "nanoseconds"]
    )
    def test_field_access(self, field) -> None:
        expected = xr.DataArray(
            getattr(self.times, field), name=field, coords=[self.times], dims=["time"]
        )
        actual = getattr(self.data.time.dt, field)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "method, parameters", [("floor", "D"), ("ceil", "D"), ("round", "D")]
    )
    def test_accessor_methods(self, method, parameters) -> None:
        dates = pd.timedelta_range(start="1 day", end="30 days", freq="6h")
        xdates = xr.DataArray(dates, dims=["time"])
        expected = getattr(dates, method)(parameters)
        actual = getattr(xdates.dt, method)(parameters)
        assert_array_equal(expected, actual)

    @requires_dask
    @pytest.mark.parametrize(
        "field", ["days", "seconds", "microseconds", "nanoseconds"]
    )
    def test_dask_field_access(self, field) -> None:
        import dask.array as da

        expected = getattr(self.times_data.dt, field)

        dask_times_arr = da.from_array(self.times_arr, chunks=(5, 5, 50))
        dask_times_2d = xr.DataArray(
            dask_times_arr, coords=self.data.coords, dims=self.data.dims, name="data"
        )

        with raise_if_dask_computes():
            actual = getattr(dask_times_2d.dt, field)

        assert isinstance(actual.data, da.Array)
        assert_chunks_equal(actual, dask_times_2d)
        assert_equal(actual, expected)

    @requires_dask
    @pytest.mark.parametrize(
        "method, parameters", [("floor", "D"), ("ceil", "D"), ("round", "D")]
    )
    def test_dask_accessor_method(self, method, parameters) -> None:
        import dask.array as da

        expected = getattr(self.times_data.dt, method)(parameters)
        dask_times_arr = da.from_array(self.times_arr, chunks=(5, 5, 50))
        dask_times_2d = xr.DataArray(
            dask_times_arr, coords=self.data.coords, dims=self.data.dims, name="data"
        )

        with raise_if_dask_computes():
            actual = getattr(dask_times_2d.dt, method)(parameters)

        assert isinstance(actual.data, da.Array)
        assert_chunks_equal(actual, dask_times_2d)
        assert_equal(actual.compute(), expected.compute())


_NT = 100


@pytest.fixture(params=_CFTIME_CALENDARS)
def calendar(request):
    return request.param


@pytest.fixture
def cftime_date_type(calendar):
    if calendar == "standard":
        calendar = "proleptic_gregorian"
    return _all_cftime_date_types()[calendar]


@pytest.fixture
def times(calendar):
    import cftime

    return cftime.num2date(
        np.arange(_NT),
        units="hours since 2000-01-01",
        calendar=calendar,
        only_use_cftime_datetimes=True,
    )


@pytest.fixture
def data(times):
    data = np.random.rand(10, 10, _NT)
    lons = np.linspace(0, 11, 10)
    lats = np.linspace(0, 20, 10)
    return xr.DataArray(
        data, coords=[lons, lats, times], dims=["lon", "lat", "time"], name="data"
    )


@pytest.fixture
def times_3d(times):
    lons = np.linspace(0, 11, 10)
    lats = np.linspace(0, 20, 10)
    times_arr = np.random.choice(times, size=(10, 10, _NT))
    return xr.DataArray(
        times_arr, coords=[lons, lats, times], dims=["lon", "lat", "time"], name="data"
    )


@requires_cftime
@pytest.mark.parametrize(
    "field", ["year", "month", "day", "hour", "dayofyear", "dayofweek"]
)
def test_field_access(data, field) -> None:
    result = getattr(data.time.dt, field)
    expected = xr.DataArray(
        getattr(xr.coding.cftimeindex.CFTimeIndex(data.time.values), field),
        name=field,
        coords=data.time.coords,
        dims=data.time.dims,
    )

    assert_equal(result, expected)


@requires_cftime
def test_calendar_cftime(data) -> None:
    expected = data.time.values[0].calendar
    assert data.time.dt.calendar == expected


def test_calendar_datetime64_2d() -> None:
    data = xr.DataArray(np.zeros((4, 5), dtype="datetime64[ns]"), dims=("x", "y"))
    assert data.dt.calendar == "proleptic_gregorian"


@requires_dask
def test_calendar_datetime64_3d_dask() -> None:
    import dask.array as da

    data = xr.DataArray(
        da.zeros((4, 5, 6), dtype="datetime64[ns]"), dims=("x", "y", "z")
    )
    with raise_if_dask_computes():
        assert data.dt.calendar == "proleptic_gregorian"


@requires_dask
@requires_cftime
def test_calendar_dask_cftime() -> None:
    from cftime import num2date

    # 3D lazy dask
    data = xr.DataArray(
        num2date(
            np.random.randint(1, 1000000, size=(4, 5, 6)),
            "hours since 1970-01-01T00:00",
            calendar="noleap",
        ),
        dims=("x", "y", "z"),
    ).chunk()
    with raise_if_dask_computes(max_computes=2):
        assert data.dt.calendar == "noleap"


@requires_cftime
def test_isocalendar_cftime(data) -> None:
    with pytest.raises(
        AttributeError, match=r"'CFTimeIndex' object has no attribute 'isocalendar'"
    ):
        data.time.dt.isocalendar()


@requires_cftime
def test_date_cftime(data) -> None:
    with pytest.raises(
        AttributeError,
        match=r"'CFTimeIndex' object has no attribute `date`. Consider using the floor method instead, for instance: `.time.dt.floor\('D'\)`.",
    ):
        data.time.dt.date()


@requires_cftime
@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_cftime_strftime_access(data) -> None:
    """compare cftime formatting against datetime formatting"""
    date_format = "%Y%m%d%H"
    result = data.time.dt.strftime(date_format)
    datetime_array = xr.DataArray(
        xr.coding.cftimeindex.CFTimeIndex(data.time.values).to_datetimeindex(
            time_unit="ns"
        ),
        name="stftime",
        coords=data.time.coords,
        dims=data.time.dims,
    )
    expected = datetime_array.dt.strftime(date_format)
    assert_equal(result, expected)


@requires_cftime
@requires_dask
@pytest.mark.parametrize(
    "field", ["year", "month", "day", "hour", "dayofyear", "dayofweek"]
)
def test_dask_field_access_1d(data, field) -> None:
    import dask.array as da

    expected = xr.DataArray(
        getattr(xr.coding.cftimeindex.CFTimeIndex(data.time.values), field),
        name=field,
        dims=["time"],
    )
    times = xr.DataArray(data.time.values, dims=["time"]).chunk({"time": 50})
    result = getattr(times.dt, field)
    assert isinstance(result.data, da.Array)
    assert result.chunks == times.chunks
    assert_equal(result.compute(), expected)


@requires_cftime
@requires_dask
@pytest.mark.parametrize(
    "field", ["year", "month", "day", "hour", "dayofyear", "dayofweek"]
)
def test_dask_field_access(times_3d, data, field) -> None:
    import dask.array as da

    expected = xr.DataArray(
        getattr(
            xr.coding.cftimeindex.CFTimeIndex(times_3d.values.ravel()), field
        ).reshape(times_3d.shape),
        name=field,
        coords=times_3d.coords,
        dims=times_3d.dims,
    )
    times_3d = times_3d.chunk({"lon": 5, "lat": 5, "time": 50})
    result = getattr(times_3d.dt, field)
    assert isinstance(result.data, da.Array)
    assert result.chunks == times_3d.chunks
    assert_equal(result.compute(), expected)


@requires_cftime
def test_seasons(cftime_date_type) -> None:
    dates = xr.DataArray(
        np.array([cftime_date_type(2000, month, 15) for month in range(1, 13)])
    )
    seasons = xr.DataArray(
        [
            "DJF",
            "DJF",
            "MAM",
            "MAM",
            "MAM",
            "JJA",
            "JJA",
            "JJA",
            "SON",
            "SON",
            "SON",
            "DJF",
        ]
    )

    assert_array_equal(seasons.values, dates.dt.season.values)


@pytest.fixture
def cftime_rounding_dataarray(cftime_date_type):
    return xr.DataArray(
        [
            [cftime_date_type(1, 1, 1, 1), cftime_date_type(1, 1, 1, 15)],
            [cftime_date_type(1, 1, 1, 23), cftime_date_type(1, 1, 2, 1)],
        ]
    )


@requires_cftime
@requires_dask
@pytest.mark.parametrize("use_dask", [False, True])
def test_cftime_floor_accessor(
    cftime_rounding_dataarray, cftime_date_type, use_dask
) -> None:
    import dask.array as da

    freq = "D"
    expected = xr.DataArray(
        [
            [cftime_date_type(1, 1, 1, 0), cftime_date_type(1, 1, 1, 0)],
            [cftime_date_type(1, 1, 1, 0), cftime_date_type(1, 1, 2, 0)],
        ],
        name="floor",
    )

    if use_dask:
        chunks = {"dim_0": 1}
        # Currently a compute is done to inspect a single value of the array
        # if it is of object dtype to check if it is a cftime.datetime (if not
        # we raise an error when using the dt accessor).
        with raise_if_dask_computes(max_computes=1):
            result = cftime_rounding_dataarray.chunk(chunks).dt.floor(freq)
        expected = expected.chunk(chunks)
        assert isinstance(result.data, da.Array)
        assert result.chunks == expected.chunks
    else:
        result = cftime_rounding_dataarray.dt.floor(freq)

    assert_identical(result, expected)


@requires_cftime
@requires_dask
@pytest.mark.parametrize("use_dask", [False, True])
def test_cftime_ceil_accessor(
    cftime_rounding_dataarray, cftime_date_type, use_dask
) -> None:
    import dask.array as da

    freq = "D"
    expected = xr.DataArray(
        [
            [cftime_date_type(1, 1, 2, 0), cftime_date_type(1, 1, 2, 0)],
            [cftime_date_type(1, 1, 2, 0), cftime_date_type(1, 1, 3, 0)],
        ],
        name="ceil",
    )

    if use_dask:
        chunks = {"dim_0": 1}
        # Currently a compute is done to inspect a single value of the array
        # if it is of object dtype to check if it is a cftime.datetime (if not
        # we raise an error when using the dt accessor).
        with raise_if_dask_computes(max_computes=1):
            result = cftime_rounding_dataarray.chunk(chunks).dt.ceil(freq)
        expected = expected.chunk(chunks)
        assert isinstance(result.data, da.Array)
        assert result.chunks == expected.chunks
    else:
        result = cftime_rounding_dataarray.dt.ceil(freq)

    assert_identical(result, expected)


@requires_cftime
@requires_dask
@pytest.mark.parametrize("use_dask", [False, True])
def test_cftime_round_accessor(
    cftime_rounding_dataarray, cftime_date_type, use_dask
) -> None:
    import dask.array as da

    freq = "D"
    expected = xr.DataArray(
        [
            [cftime_date_type(1, 1, 1, 0), cftime_date_type(1, 1, 2, 0)],
            [cftime_date_type(1, 1, 2, 0), cftime_date_type(1, 1, 2, 0)],
        ],
        name="round",
    )

    if use_dask:
        chunks = {"dim_0": 1}
        # Currently a compute is done to inspect a single value of the array
        # if it is of object dtype to check if it is a cftime.datetime (if not
        # we raise an error when using the dt accessor).
        with raise_if_dask_computes(max_computes=1):
            result = cftime_rounding_dataarray.chunk(chunks).dt.round(freq)
        expected = expected.chunk(chunks)
        assert isinstance(result.data, da.Array)
        assert result.chunks == expected.chunks
    else:
        result = cftime_rounding_dataarray.dt.round(freq)

    assert_identical(result, expected)


@pytest.mark.parametrize(
    "use_cftime",
    [False, pytest.param(True, marks=requires_cftime)],
    ids=lambda x: f"use_cftime={x}",
)
@pytest.mark.parametrize(
    "use_dask",
    [False, pytest.param(True, marks=requires_dask)],
    ids=lambda x: f"use_dask={x}",
)
def test_decimal_year(use_cftime, use_dask) -> None:
    year = 2000
    periods = 10
    freq = "h"

    shape = (2, 5)
    dims = ["x", "y"]
    hours_in_year = 24 * 366

    times = xr.date_range(f"{year}", periods=periods, freq=freq, use_cftime=use_cftime)

    da = xr.DataArray(times.values.reshape(shape), dims=dims)

    if use_dask:
        da = da.chunk({"y": 2})
        # Computing the decimal year for a cftime datetime array requires a
        # number of small computes (6):
        # - 4x one compute per .dt accessor call (requires inspecting one
        #   object-dtype array element to see if it is time-like)
        # - 2x one compute per calendar inference (requires inspecting one
        #   array element to read off the calendar)
        max_computes = 6 * use_cftime
        with raise_if_dask_computes(max_computes=max_computes):
            result = da.dt.decimal_year
    else:
        result = da.dt.decimal_year

    expected = xr.DataArray(
        year + np.arange(periods).reshape(shape) / hours_in_year, dims=dims
    )
    xr.testing.assert_equal(result, expected)
