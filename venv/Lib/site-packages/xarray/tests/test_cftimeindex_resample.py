from __future__ import annotations

import datetime
from typing import TypedDict

import numpy as np
import pandas as pd
import pytest
from packaging.version import Version

import xarray as xr
from xarray.coding.cftime_offsets import _new_to_legacy_freq
from xarray.coding.cftimeindex import CFTimeIndex
from xarray.core.resample_cftime import CFTimeGrouper

cftime = pytest.importorskip("cftime")


# Create a list of pairs of similar-length initial and resample frequencies
# that cover:
# - Resampling from shorter to longer frequencies
# - Resampling from longer to shorter frequencies
# - Resampling from one initial frequency to another.
# These are used to test the cftime version of resample against pandas
# with a standard calendar.
FREQS = [
    ("8003D", "4001D"),
    ("8003D", "16006D"),
    ("8003D", "21YS"),
    ("6h", "3h"),
    ("6h", "12h"),
    ("6h", "400min"),
    ("3D", "D"),
    ("3D", "6D"),
    ("11D", "MS"),
    ("3MS", "MS"),
    ("3MS", "6MS"),
    ("3MS", "85D"),
    ("7ME", "3ME"),
    ("7ME", "14ME"),
    ("7ME", "2QS-APR"),
    ("43QS-AUG", "21QS-AUG"),
    ("43QS-AUG", "86QS-AUG"),
    ("43QS-AUG", "11YE-JUN"),
    ("11QE-JUN", "5QE-JUN"),
    ("11QE-JUN", "22QE-JUN"),
    ("11QE-JUN", "51MS"),
    ("3YS-MAR", "YS-MAR"),
    ("3YS-MAR", "6YS-MAR"),
    ("3YS-MAR", "14QE-FEB"),
    ("7YE-MAY", "3YE-MAY"),
    ("7YE-MAY", "14YE-MAY"),
    ("7YE-MAY", "85ME"),
]


def compare_against_pandas(
    da_datetimeindex,
    da_cftimeindex,
    freq,
    closed=None,
    label=None,
    offset=None,
    origin=None,
) -> None:
    if isinstance(origin, tuple):
        origin_pandas = pd.Timestamp(datetime.datetime(*origin))
        origin_cftime = cftime.DatetimeGregorian(*origin)
    else:
        origin_pandas = origin
        origin_cftime = origin

    try:
        result_datetimeindex = da_datetimeindex.resample(
            time=freq,
            closed=closed,
            label=label,
            offset=offset,
            origin=origin_pandas,
        ).mean()
    except ValueError:
        with pytest.raises(ValueError):
            da_cftimeindex.resample(
                time=freq,
                closed=closed,
                label=label,
                origin=origin_cftime,
                offset=offset,
            ).mean()
    else:
        result_cftimeindex = da_cftimeindex.resample(
            time=freq,
            closed=closed,
            label=label,
            origin=origin_cftime,
            offset=offset,
        ).mean()
    # TODO (benbovy - flexible indexes): update when CFTimeIndex is a xarray Index subclass
    result_cftimeindex["time"] = (
        result_cftimeindex.xindexes["time"].to_pandas_index().to_datetimeindex()
    )
    xr.testing.assert_identical(result_cftimeindex, result_datetimeindex)


def da(index) -> xr.DataArray:
    return xr.DataArray(
        np.arange(100.0, 100.0 + index.size), coords=[index], dims=["time"]
    )


@pytest.mark.parametrize("freqs", FREQS, ids=lambda x: "{}->{}".format(*x))
@pytest.mark.parametrize("closed", [None, "left", "right"])
@pytest.mark.parametrize("label", [None, "left", "right"])
@pytest.mark.parametrize("offset", [None, "5s"], ids=lambda x: f"{x}")
def test_resample(freqs, closed, label, offset) -> None:
    initial_freq, resample_freq = freqs
    if (
        resample_freq == "4001D"
        and closed == "right"
        and Version(pd.__version__) < Version("2.2")
    ):
        pytest.skip(
            "Pandas fixed a bug in this test case in version 2.2, which we "
            "ported to xarray, so this test no longer produces the same "
            "result as pandas for earlier pandas versions."
        )
    start = "2000-01-01T12:07:01"
    origin = "start"

    datetime_index = pd.date_range(
        start=start, periods=5, freq=_new_to_legacy_freq(initial_freq)
    )
    cftime_index = xr.cftime_range(start=start, periods=5, freq=initial_freq)
    da_datetimeindex = da(datetime_index)
    da_cftimeindex = da(cftime_index)

    compare_against_pandas(
        da_datetimeindex,
        da_cftimeindex,
        resample_freq,
        closed=closed,
        label=label,
        offset=offset,
        origin=origin,
    )


@pytest.mark.parametrize(
    ("freq", "expected"),
    [
        ("s", "left"),
        ("min", "left"),
        ("h", "left"),
        ("D", "left"),
        ("ME", "right"),
        ("MS", "left"),
        ("QE", "right"),
        ("QS", "left"),
        ("YE", "right"),
        ("YS", "left"),
    ],
)
def test_closed_label_defaults(freq, expected) -> None:
    assert CFTimeGrouper(freq=freq).closed == expected
    assert CFTimeGrouper(freq=freq).label == expected


@pytest.mark.filterwarnings("ignore:Converting a CFTimeIndex")
@pytest.mark.parametrize(
    "calendar", ["gregorian", "noleap", "all_leap", "360_day", "julian"]
)
def test_calendars(calendar: str) -> None:
    # Limited testing for non-standard calendars
    freq, closed, label = "8001min", None, None
    xr_index = xr.cftime_range(
        start="2004-01-01T12:07:01", periods=7, freq="3D", calendar=calendar
    )
    pd_index = pd.date_range(start="2004-01-01T12:07:01", periods=7, freq="3D")
    da_cftime = da(xr_index).resample(time=freq, closed=closed, label=label).mean()
    da_datetime = da(pd_index).resample(time=freq, closed=closed, label=label).mean()
    # TODO (benbovy - flexible indexes): update when CFTimeIndex is a xarray Index subclass
    new_pd_index = da_cftime.xindexes["time"].to_pandas_index()
    assert isinstance(new_pd_index, CFTimeIndex)  # shouldn't that be a pd.Index?
    da_cftime["time"] = new_pd_index.to_datetimeindex()
    xr.testing.assert_identical(da_cftime, da_datetime)


class DateRangeKwargs(TypedDict):
    start: str
    periods: int
    freq: str


@pytest.mark.parametrize("closed", ["left", "right"])
@pytest.mark.parametrize(
    "origin",
    ["start_day", "start", "end", "end_day", "epoch", (1970, 1, 1, 3, 2)],
    ids=lambda x: f"{x}",
)
def test_origin(closed, origin) -> None:
    initial_freq, resample_freq = ("3h", "9h")
    start = "1969-12-31T12:07:01"
    index_kwargs: DateRangeKwargs = dict(start=start, periods=12, freq=initial_freq)
    datetime_index = pd.date_range(**index_kwargs)
    cftime_index = xr.cftime_range(**index_kwargs)
    da_datetimeindex = da(datetime_index)
    da_cftimeindex = da(cftime_index)

    compare_against_pandas(
        da_datetimeindex,
        da_cftimeindex,
        resample_freq,
        closed=closed,
        origin=origin,
    )


@pytest.mark.parametrize("offset", ["foo", "5MS", 10])
def test_invalid_offset_error(offset: str | int) -> None:
    cftime_index = xr.cftime_range("2000", periods=5)
    da_cftime = da(cftime_index)
    with pytest.raises(ValueError, match="offset must be"):
        da_cftime.resample(time="2D", offset=offset)  # type: ignore[arg-type]


def test_timedelta_offset() -> None:
    timedelta = datetime.timedelta(seconds=5)
    string = "5s"

    cftime_index = xr.cftime_range("2000", periods=5)
    da_cftime = da(cftime_index)

    timedelta_result = da_cftime.resample(time="2D", offset=timedelta).mean()
    string_result = da_cftime.resample(time="2D", offset=string).mean()
    xr.testing.assert_identical(timedelta_result, string_result)
