from datetime import datetime

import numpy as np
import pytest

from pandas.errors import Pandas4Warning

from pandas.core.dtypes.common import is_extension_array_dtype

import pandas as pd
from pandas import (
    DataFrame,
    DatetimeIndex,
    Index,
    MultiIndex,
    NaT,
    PeriodIndex,
    Series,
    TimedeltaIndex,
)
import pandas._testing as tm
from pandas.core.groupby.groupby import DataError
from pandas.core.groupby.grouper import Grouper
from pandas.core.indexes.datetimes import date_range
from pandas.core.indexes.period import period_range
from pandas.core.indexes.timedeltas import timedelta_range
from pandas.core.resample import _asfreq_compat


@pytest.fixture(
    params=[
        "linear",
        "time",
        "index",
        "values",
        "nearest",
        "zero",
        "slinear",
        "quadratic",
        "cubic",
        "barycentric",
        "krogh",
        "from_derivatives",
        "piecewise_polynomial",
        "pchip",
        "akima",
    ],
)
def all_1d_no_arg_interpolation_methods(request):
    return request.param


@pytest.mark.parametrize("freq", ["2D", "1h"])
@pytest.mark.parametrize(
    "index",
    [
        timedelta_range("1 day", "10 day", freq="D"),
        date_range(datetime(2005, 1, 1), datetime(2005, 1, 10), freq="D"),
    ],
)
def test_asfreq(frame_or_series, index, freq):
    obj = frame_or_series(range(len(index)), index=index)
    idx_range = date_range if isinstance(index, DatetimeIndex) else timedelta_range

    result = obj.resample(freq).asfreq()
    new_index = idx_range(obj.index[0], obj.index[-1], freq=freq)
    expected = obj.reindex(new_index)
    tm.assert_almost_equal(result, expected)


@pytest.mark.parametrize(
    "index",
    [
        timedelta_range("1 day", "10 day", freq="D"),
        date_range(datetime(2005, 1, 1), datetime(2005, 1, 10), freq="D"),
    ],
)
def test_asfreq_fill_value(index):
    # test for fill value during resampling, issue 3715

    ser = Series(range(len(index)), index=index, name="a")
    idx_range = date_range if isinstance(index, DatetimeIndex) else timedelta_range

    result = ser.resample("1h").asfreq()
    new_index = idx_range(ser.index[0], ser.index[-1], freq="1h")
    expected = ser.reindex(new_index)
    tm.assert_series_equal(result, expected)

    # Explicit cast to float to avoid implicit cast when setting None
    frame = ser.astype("float").to_frame("value")
    frame.iloc[1] = None
    result = frame.resample("1h").asfreq(fill_value=4.0)
    new_index = idx_range(frame.index[0], frame.index[-1], freq="1h")
    expected = frame.reindex(new_index, fill_value=4.0)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "index",
    [
        timedelta_range("1 day", "3 day", freq="D"),
        date_range(datetime(2005, 1, 1), datetime(2005, 1, 3), freq="D"),
        period_range(datetime(2005, 1, 1), datetime(2005, 1, 3), freq="D"),
    ],
)
def test_resample_interpolate(index):
    # GH#12925
    df = DataFrame(range(len(index)), index=index)
    result = df.resample("1min").asfreq().interpolate()
    expected = df.resample("1min").interpolate()
    tm.assert_frame_equal(result, expected)


def test_resample_interpolate_inplace_deprecated():
    # GH#58690
    dti = date_range(datetime(2005, 1, 1), datetime(2005, 1, 10), freq="D")

    df = DataFrame(range(len(dti)), index=dti)
    rs = df.resample("1min")
    msg = "The 'inplace' keyword in DatetimeIndexResampler.interpolate"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        rs.interpolate(inplace=False)

    msg2 = "Cannot interpolate inplace on a resampled object"
    with pytest.raises(ValueError, match=msg2):
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            rs.interpolate(inplace=True)


def test_resample_interpolate_regular_sampling_off_grid(
    all_1d_no_arg_interpolation_methods,
):
    pytest.importorskip("scipy")
    # GH#21351
    index = date_range("2000-01-01 00:01:00", periods=5, freq="2h")
    ser = Series(np.arange(5.0), index)

    method = all_1d_no_arg_interpolation_methods
    result = ser.resample("1h").interpolate(method)

    if method == "linear":
        values = np.repeat(np.arange(0.0, 4.0), 2) + np.tile([1 / 3, 2 / 3], 4)
    elif method == "nearest":
        values = np.repeat(np.arange(0.0, 5.0), 2)[1:-1]
    elif method == "zero":
        values = np.repeat(np.arange(0.0, 4.0), 2)
    else:
        values = 0.491667 + np.arange(0.0, 4.0, 0.5)
    values = np.insert(values, 0, np.nan)
    index = date_range("2000-01-01 00:00:00", periods=9, freq="1h")
    expected = Series(values, index=index)
    tm.assert_series_equal(result, expected)


def test_resample_interpolate_irregular_sampling(all_1d_no_arg_interpolation_methods):
    pytest.importorskip("scipy")
    # GH#21351
    ser = Series(
        np.linspace(0.0, 1.0, 5),
        index=DatetimeIndex(
            [
                "2000-01-01 00:00:03",
                "2000-01-01 00:00:22",
                "2000-01-01 00:00:24",
                "2000-01-01 00:00:31",
                "2000-01-01 00:00:39",
            ]
        ),
    )

    # Resample to 5 second sampling and interpolate with the given method
    ser_resampled = ser.resample("5s").interpolate(all_1d_no_arg_interpolation_methods)

    # Check that none of the resampled values are NaN, except the first one
    # which lies 3 seconds before the first actual data point
    assert np.isnan(ser_resampled.iloc[0])
    assert not ser_resampled.iloc[1:].isna().any()


def test_raises_on_non_datetimelike_index():
    # this is a non datetimelike index
    xp = DataFrame()
    msg = (
        "Only valid with DatetimeIndex, TimedeltaIndex or PeriodIndex, "
        "but got an instance of 'RangeIndex'"
    )
    with pytest.raises(TypeError, match=msg):
        xp.resample("YE")


@pytest.mark.parametrize(
    "index",
    [
        PeriodIndex([], freq="D", name="a"),
        DatetimeIndex([], name="a"),
        TimedeltaIndex([], name="a"),
    ],
)
@pytest.mark.parametrize("freq", ["ME", "D", "h"])
def test_resample_empty_series(freq, index, resample_method):
    # GH12771 & GH12868

    ser = Series(index=index, dtype=float)
    if freq == "ME" and isinstance(ser.index, TimedeltaIndex):
        msg = (
            "Resampling on a TimedeltaIndex requires fixed-duration `freq`, "
            "e.g. '24h' or '3D', not <MonthEnd>"
        )
        with pytest.raises(ValueError, match=msg):
            ser.resample(freq)
        return
    elif freq == "ME" and isinstance(ser.index, PeriodIndex):
        # index is PeriodIndex, so convert to corresponding Period freq
        freq = "M"
    rs = ser.resample(freq)
    result = getattr(rs, resample_method)()

    if resample_method == "ohlc":
        expected = DataFrame(
            [], index=ser.index[:0], columns=["open", "high", "low", "close"]
        )
        expected.index = _asfreq_compat(ser.index, freq)
        tm.assert_frame_equal(result, expected, check_dtype=False)
    else:
        expected = ser.copy()
        expected.index = _asfreq_compat(ser.index, freq)
        tm.assert_series_equal(result, expected, check_dtype=False)

    tm.assert_index_equal(result.index, expected.index)
    assert result.index.freq == expected.index.freq


@pytest.mark.parametrize("min_count", [0, 1])
def test_resample_empty_sum_string(string_dtype_no_object, min_count):
    # https://github.com/pandas-dev/pandas/issues/60229
    dtype = string_dtype_no_object
    ser = Series(
        pd.NA,
        index=DatetimeIndex(
            [
                "2000-01-01 00:00:00",
                "2000-01-01 00:00:10",
                "2000-01-01 00:00:20",
                "2000-01-01 00:00:30",
            ]
        ),
        dtype=dtype,
    )
    rs = ser.resample("20s")
    result = rs.sum(min_count=min_count)

    value = "" if min_count == 0 else pd.NA
    index = date_range(start="2000-01-01", freq="20s", periods=2, unit="us")
    expected = Series(value, index=index, dtype=dtype)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "freq",
    [
        pytest.param("ME", marks=pytest.mark.xfail(reason="Don't know why this fails")),
        "D",
        "h",
    ],
)
def test_resample_nat_index_series(freq, resample_method):
    # GH39227

    ser = Series(range(5), index=PeriodIndex([NaT] * 5, freq=freq))

    rs = ser.resample(freq)
    result = getattr(rs, resample_method)()

    if resample_method == "ohlc":
        expected = DataFrame(
            [], index=ser.index[:0], columns=["open", "high", "low", "close"]
        )
        tm.assert_frame_equal(result, expected, check_dtype=False)
    else:
        expected = ser[:0].copy()
        tm.assert_series_equal(result, expected, check_dtype=False)
    tm.assert_index_equal(result.index, expected.index)
    assert result.index.freq == expected.index.freq


@pytest.mark.parametrize(
    "index",
    [
        PeriodIndex([], freq="D", name="a"),
        DatetimeIndex([], name="a"),
        TimedeltaIndex([], name="a"),
    ],
)
@pytest.mark.parametrize("freq", ["ME", "D", "h"])
@pytest.mark.parametrize("resample_method", ["count", "size"])
def test_resample_count_empty_series(freq, index, resample_method):
    # GH28427
    ser = Series(index=index)
    if freq == "ME" and isinstance(ser.index, TimedeltaIndex):
        msg = (
            "Resampling on a TimedeltaIndex requires fixed-duration `freq`, "
            "e.g. '24h' or '3D', not <MonthEnd>"
        )
        with pytest.raises(ValueError, match=msg):
            ser.resample(freq)
        return
    elif freq == "ME" and isinstance(ser.index, PeriodIndex):
        # index is PeriodIndex, so convert to corresponding Period freq
        freq = "M"
    rs = ser.resample(freq)

    result = getattr(rs, resample_method)()

    index = _asfreq_compat(ser.index, freq)

    expected = Series([], dtype="int64", index=index, name=ser.name)

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "index", [DatetimeIndex([]), TimedeltaIndex([]), PeriodIndex([], freq="D")]
)
@pytest.mark.parametrize("freq", ["ME", "D", "h"])
def test_resample_empty_dataframe(index, freq, resample_method):
    # GH13212
    df = DataFrame(index=index)
    # count retains dimensions too
    if freq == "ME" and isinstance(df.index, TimedeltaIndex):
        msg = (
            "Resampling on a TimedeltaIndex requires fixed-duration `freq`, "
            "e.g. '24h' or '3D', not <MonthEnd>"
        )
        with pytest.raises(ValueError, match=msg):
            df.resample(freq, group_keys=False)
        return
    elif freq == "ME" and isinstance(df.index, PeriodIndex):
        # index is PeriodIndex, so convert to corresponding Period freq
        freq = "M"
    rs = df.resample(freq, group_keys=False)
    result = getattr(rs, resample_method)()
    if resample_method == "ohlc":
        # TODO: no tests with len(df.columns) > 0
        mi = MultiIndex.from_product([df.columns, ["open", "high", "low", "close"]])
        expected = DataFrame([], index=df.index[:0], columns=mi, dtype=np.float64)
        expected.index = _asfreq_compat(df.index, freq)

    elif resample_method != "size":
        expected = df.copy()
    else:
        # GH14962
        expected = Series([], dtype=np.int64)

    expected.index = _asfreq_compat(df.index, freq)

    tm.assert_index_equal(result.index, expected.index)
    assert result.index.freq == expected.index.freq
    tm.assert_almost_equal(result, expected)

    # test size for GH13212 (currently stays as df)


@pytest.mark.parametrize(
    "index", [DatetimeIndex([]), TimedeltaIndex([]), PeriodIndex([], freq="D")]
)
@pytest.mark.parametrize("freq", ["ME", "D", "h"])
def test_resample_count_empty_dataframe(freq, index):
    # GH28427
    empty_frame_dti = DataFrame(index=index, columns=Index(["a"], dtype=object))

    if freq == "ME" and isinstance(empty_frame_dti.index, TimedeltaIndex):
        msg = (
            "Resampling on a TimedeltaIndex requires fixed-duration `freq`, "
            "e.g. '24h' or '3D', not <MonthEnd>"
        )
        with pytest.raises(ValueError, match=msg):
            empty_frame_dti.resample(freq)
        return
    elif freq == "ME" and isinstance(empty_frame_dti.index, PeriodIndex):
        # index is PeriodIndex, so convert to corresponding Period freq
        freq = "M"
    result = empty_frame_dti.resample(freq).count()

    index = _asfreq_compat(empty_frame_dti.index, freq)

    expected = DataFrame(dtype="int64", index=index, columns=Index(["a"], dtype=object))

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "index", [DatetimeIndex([]), TimedeltaIndex([]), PeriodIndex([], freq="D")]
)
@pytest.mark.parametrize("freq", ["ME", "D", "h"])
def test_resample_size_empty_dataframe(freq, index):
    # GH28427

    empty_frame_dti = DataFrame(index=index, columns=Index(["a"], dtype=object))

    if freq == "ME" and isinstance(empty_frame_dti.index, TimedeltaIndex):
        msg = (
            "Resampling on a TimedeltaIndex requires fixed-duration `freq`, "
            "e.g. '24h' or '3D', not <MonthEnd>"
        )
        with pytest.raises(ValueError, match=msg):
            empty_frame_dti.resample(freq)
        return
    elif freq == "ME" and isinstance(empty_frame_dti.index, PeriodIndex):
        # index is PeriodIndex, so convert to corresponding Period freq
        freq = "M"
    result = empty_frame_dti.resample(freq).size()

    index = _asfreq_compat(empty_frame_dti.index, freq)

    expected = Series([], dtype="int64", index=index)

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("index", [DatetimeIndex([]), TimedeltaIndex([])])
@pytest.mark.parametrize("freq", ["D", "h"])
@pytest.mark.parametrize(
    "method", ["ffill", "bfill", "nearest", "asfreq", "interpolate", "mean"]
)
def test_resample_apply_empty_dataframe(index, freq, method):
    # GH#55572
    empty_frame_dti = DataFrame(index=index)

    rs = empty_frame_dti.resample(freq)
    result = rs.apply(getattr(rs, method))

    expected_index = _asfreq_compat(empty_frame_dti.index, freq)
    expected = DataFrame([], index=expected_index)

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "index",
    [
        PeriodIndex([], freq="M", name="a"),
        DatetimeIndex([], name="a"),
        TimedeltaIndex([], name="a"),
    ],
)
@pytest.mark.parametrize("dtype", [float, int, object, "datetime64[ns]"])
def test_resample_empty_dtypes(index, dtype, resample_method):
    # Empty series were sometimes causing a segfault (for the functions
    # with Cython bounds-checking disabled) or an IndexError.  We just run
    # them to ensure they no longer do.  (GH #10228)
    empty_series_dti = Series([], index, dtype)
    rs = empty_series_dti.resample("D", group_keys=False)
    try:
        getattr(rs, resample_method)()
    except DataError:
        # Ignore these since some combinations are invalid
        # (ex: doing mean with dtype of np.object_)
        pass


@pytest.mark.parametrize(
    "index",
    [
        PeriodIndex([], freq="D", name="a"),
        DatetimeIndex([], name="a"),
        TimedeltaIndex([], name="a"),
    ],
)
@pytest.mark.parametrize("freq", ["ME", "D", "h"])
def test_apply_to_empty_series(index, freq):
    # GH 14313
    ser = Series(index=index)

    if freq == "ME" and isinstance(ser.index, TimedeltaIndex):
        msg = (
            "Resampling on a TimedeltaIndex requires fixed-duration `freq`, "
            "e.g. '24h' or '3D', not <MonthEnd>"
        )
        with pytest.raises(ValueError, match=msg):
            ser.resample(freq)
        return
    elif freq == "ME" and isinstance(ser.index, PeriodIndex):
        # index is PeriodIndex, so convert to corresponding Period freq
        freq = "M"
    result = ser.resample(freq, group_keys=False).apply(lambda x: 1)
    expected = ser.resample(freq).apply("sum")

    tm.assert_series_equal(result, expected, check_dtype=False)


@pytest.mark.parametrize(
    "index",
    [
        timedelta_range("1 day", "10 day", freq="D"),
        date_range(datetime(2005, 1, 1), datetime(2005, 1, 10), freq="D"),
        period_range(datetime(2005, 1, 1), datetime(2005, 1, 10), freq="D"),
    ],
)
def test_resampler_is_iterable(index):
    # GH 15314
    series = Series(range(len(index)), index=index)
    freq = "h"
    tg = Grouper(freq=freq, convention="start")
    grouped = series.groupby(tg)
    resampled = series.resample(freq)
    for (rk, rv), (gk, gv) in zip(resampled, grouped):
        assert rk == gk
        tm.assert_series_equal(rv, gv)


@pytest.mark.parametrize(
    "index",
    [
        timedelta_range("1 day", "10 day", freq="D"),
        date_range(datetime(2005, 1, 1), datetime(2005, 1, 10), freq="D"),
        period_range(datetime(2005, 1, 1), datetime(2005, 1, 10), freq="D"),
    ],
)
def test_resample_quantile(index):
    # GH 15023
    ser = Series(range(len(index)), index=index)
    q = 0.75
    freq = "h"

    result = ser.resample(freq).quantile(q)
    expected = ser.resample(freq).agg(lambda x: x.quantile(q)).rename(ser.name)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("how", ["first", "last"])
def test_first_last_skipna(any_real_nullable_dtype, skipna, how):
    # GH#57019
    if is_extension_array_dtype(any_real_nullable_dtype):
        na_value = Series(dtype=any_real_nullable_dtype).dtype.na_value
    else:
        na_value = np.nan
    df = DataFrame(
        {
            "a": [2, 1, 1, 2],
            "b": [na_value, 3.0, na_value, 4.0],
            "c": [na_value, 3.0, na_value, 4.0],
        },
        index=date_range("2020-01-01", periods=4, freq="D", unit="ns"),
        dtype=any_real_nullable_dtype,
    )
    rs = df.resample("ME")
    method = getattr(rs, how)
    result = method(skipna=skipna)

    ts = pd.to_datetime("2020-01-31").as_unit("ns")
    gb = df.groupby(df.shape[0] * [ts])
    expected = getattr(gb, how)(skipna=skipna)
    expected.index.freq = "ME"
    tm.assert_frame_equal(result, expected)
