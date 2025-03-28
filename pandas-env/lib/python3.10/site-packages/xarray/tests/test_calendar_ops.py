from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from xarray import CFTimeIndex, DataArray, infer_freq
from xarray.coding.calendar_ops import convert_calendar, interp_calendar
from xarray.coding.cftime_offsets import date_range
from xarray.testing import assert_identical
from xarray.tests import requires_cftime

cftime = pytest.importorskip("cftime")


@pytest.mark.parametrize(
    "source, target, use_cftime, freq",
    [
        ("standard", "noleap", None, "D"),
        ("noleap", "proleptic_gregorian", True, "D"),
        ("noleap", "all_leap", None, "D"),
        ("all_leap", "proleptic_gregorian", False, "4h"),
    ],
)
def test_convert_calendar(source, target, use_cftime, freq):
    src = DataArray(
        date_range("2004-01-01", "2004-12-31", freq=freq, calendar=source),
        dims=("time",),
        name="time",
    )
    da_src = DataArray(
        np.linspace(0, 1, src.size), dims=("time",), coords={"time": src}
    )

    conv = convert_calendar(da_src, target, use_cftime=use_cftime)

    assert conv.time.dt.calendar == target

    if source != "noleap":
        expected_times = date_range(
            "2004-01-01",
            "2004-12-31",
            freq=freq,
            use_cftime=use_cftime,
            calendar=target,
        )
    else:
        expected_times_pre_leap = date_range(
            "2004-01-01",
            "2004-02-28",
            freq=freq,
            use_cftime=use_cftime,
            calendar=target,
        )
        expected_times_post_leap = date_range(
            "2004-03-01",
            "2004-12-31",
            freq=freq,
            use_cftime=use_cftime,
            calendar=target,
        )
        expected_times = expected_times_pre_leap.append(expected_times_post_leap)
    np.testing.assert_array_equal(conv.time, expected_times)


@pytest.mark.parametrize(
    "source,target,freq",
    [
        ("standard", "360_day", "D"),
        ("360_day", "proleptic_gregorian", "D"),
        ("proleptic_gregorian", "360_day", "4h"),
    ],
)
@pytest.mark.parametrize("align_on", ["date", "year"])
def test_convert_calendar_360_days(source, target, freq, align_on):
    src = DataArray(
        date_range("2004-01-01", "2004-12-30", freq=freq, calendar=source),
        dims=("time",),
        name="time",
    )
    da_src = DataArray(
        np.linspace(0, 1, src.size), dims=("time",), coords={"time": src}
    )

    conv = convert_calendar(da_src, target, align_on=align_on)

    assert conv.time.dt.calendar == target

    if align_on == "date":
        np.testing.assert_array_equal(
            conv.time.resample(time="ME").last().dt.day,
            [30, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30],
        )
    elif target == "360_day":
        np.testing.assert_array_equal(
            conv.time.resample(time="ME").last().dt.day,
            [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 29],
        )
    else:
        np.testing.assert_array_equal(
            conv.time.resample(time="ME").last().dt.day,
            [30, 29, 30, 30, 31, 30, 30, 31, 30, 31, 29, 31],
        )
    if source == "360_day" and align_on == "year":
        assert conv.size == 360 if freq == "D" else 360 * 4
    else:
        assert conv.size == 359 if freq == "D" else 359 * 4


def test_convert_calendar_360_days_random():
    da_std = DataArray(
        np.linspace(0, 1, 366),
        dims=("time",),
        coords={
            "time": date_range(
                "2004-01-01",
                "2004-12-31",
                freq="D",
                calendar="standard",
                use_cftime=False,
            )
        },
    )
    da_360 = DataArray(
        np.linspace(0, 1, 360),
        dims=("time",),
        coords={
            "time": date_range("2004-01-01", "2004-12-30", freq="D", calendar="360_day")
        },
    )

    conv = convert_calendar(da_std, "360_day", align_on="random")
    conv2 = convert_calendar(da_std, "360_day", align_on="random")
    assert (conv != conv2).any()

    conv = convert_calendar(da_360, "standard", use_cftime=False, align_on="random")
    assert np.datetime64("2004-02-29") not in conv.time
    conv2 = convert_calendar(da_360, "standard", use_cftime=False, align_on="random")
    assert (conv2 != conv).any()

    # Ensure that added days are evenly distributed in the 5 fifths of each year
    conv = convert_calendar(da_360, "noleap", align_on="random", missing=np.nan)
    conv = conv.where(conv.isnull(), drop=True)
    nandoys = conv.time.dt.dayofyear[:366]
    assert all(nandoys < np.array([74, 147, 220, 293, 366]))
    assert all(nandoys > np.array([0, 73, 146, 219, 292]))


@requires_cftime
@pytest.mark.parametrize(
    "source,target,freq",
    [
        ("standard", "noleap", "D"),
        ("noleap", "proleptic_gregorian", "4h"),
        ("noleap", "all_leap", "ME"),
        ("360_day", "noleap", "D"),
        ("noleap", "360_day", "D"),
    ],
)
def test_convert_calendar_missing(source, target, freq):
    src = DataArray(
        date_range(
            "2004-01-01",
            "2004-12-31" if source != "360_day" else "2004-12-30",
            freq=freq,
            calendar=source,
        ),
        dims=("time",),
        name="time",
    )
    da_src = DataArray(
        np.linspace(0, 1, src.size), dims=("time",), coords={"time": src}
    )
    out = convert_calendar(da_src, target, missing=np.nan, align_on="date")

    expected_freq = freq
    assert infer_freq(out.time) == expected_freq

    expected = date_range(
        "2004-01-01",
        "2004-12-31" if target != "360_day" else "2004-12-30",
        freq=freq,
        calendar=target,
    )
    np.testing.assert_array_equal(out.time, expected)

    if freq != "ME":
        out_without_missing = convert_calendar(da_src, target, align_on="date")
        expected_nan = out.isel(time=~out.time.isin(out_without_missing.time))
        assert expected_nan.isnull().all()

        expected_not_nan = out.sel(time=out_without_missing.time)
        assert_identical(expected_not_nan, out_without_missing)


@requires_cftime
def test_convert_calendar_errors():
    src_nl = DataArray(
        date_range("0000-01-01", "0000-12-31", freq="D", calendar="noleap"),
        dims=("time",),
        name="time",
    )
    # no align_on for conversion to 360_day
    with pytest.raises(ValueError, match="Argument `align_on` must be specified"):
        convert_calendar(src_nl, "360_day")

    # Standard doesn't support year 0
    with pytest.raises(
        ValueError, match="Source time coordinate contains dates with year 0"
    ):
        convert_calendar(src_nl, "standard")

    # no align_on for conversion from 360 day
    src_360 = convert_calendar(src_nl, "360_day", align_on="year")
    with pytest.raises(ValueError, match="Argument `align_on` must be specified"):
        convert_calendar(src_360, "noleap")

    # Datetime objects
    da = DataArray([0, 1, 2], dims=("x",), name="x")
    with pytest.raises(ValueError, match="Coordinate x must contain datetime objects."):
        convert_calendar(da, "standard", dim="x")


def test_convert_calendar_same_calendar():
    src = DataArray(
        date_range("2000-01-01", periods=12, freq="6h", use_cftime=False),
        dims=("time",),
        name="time",
    )
    out = convert_calendar(src, "proleptic_gregorian")
    assert src is out


@pytest.mark.parametrize(
    "source,target",
    [
        ("standard", "noleap"),
        ("noleap", "proleptic_gregorian"),
        ("standard", "360_day"),
        ("360_day", "proleptic_gregorian"),
        ("noleap", "all_leap"),
        ("360_day", "noleap"),
    ],
)
def test_interp_calendar(source, target):
    src = DataArray(
        date_range("2004-01-01", "2004-07-30", freq="D", calendar=source),
        dims=("time",),
        name="time",
    )
    tgt = DataArray(
        date_range("2004-01-01", "2004-07-30", freq="D", calendar=target),
        dims=("time",),
        name="time",
    )
    da_src = DataArray(
        np.linspace(0, 1, src.size), dims=("time",), coords={"time": src}
    )
    conv = interp_calendar(da_src, tgt)

    assert_identical(tgt.time, conv.time)

    np.testing.assert_almost_equal(conv.max(), 1, 2)
    assert conv.min() == 0


@requires_cftime
def test_interp_calendar_errors():
    src_nl = DataArray(
        [1] * 100,
        dims=("time",),
        coords={
            "time": date_range("0000-01-01", periods=100, freq="MS", calendar="noleap")
        },
    )
    tgt_360 = date_range("0001-01-01", "0001-12-30", freq="MS", calendar="standard")

    with pytest.raises(
        ValueError, match="Source time coordinate contains dates with year 0"
    ):
        interp_calendar(src_nl, tgt_360)

    da1 = DataArray([0, 1, 2], dims=("x",), name="x")
    da2 = da1 + 1

    with pytest.raises(
        ValueError, match="Both 'source.x' and 'target' must contain datetime objects."
    ):
        interp_calendar(da1, da2, dim="x")


@requires_cftime
@pytest.mark.parametrize(
    ("source_calendar", "target_calendar", "expected_index"),
    [("standard", "noleap", CFTimeIndex), ("all_leap", "standard", pd.DatetimeIndex)],
)
def test_convert_calendar_produces_time_index(
    source_calendar, target_calendar, expected_index
):
    # https://github.com/pydata/xarray/issues/9138
    time = date_range("2000-01-01", "2002-01-01", freq="D", calendar=source_calendar)
    temperature = np.ones(len(time))
    da = DataArray(
        data=temperature,
        dims=["time"],
        coords=dict(
            time=time,
        ),
    )
    converted = da.convert_calendar(target_calendar)
    assert isinstance(converted.indexes["time"], expected_index)
