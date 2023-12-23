from datetime import datetime

import numpy as np
import pytest

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


@pytest.mark.parametrize("freq", ["2D", "1h"])
@pytest.mark.parametrize(
    "index",
    [
        timedelta_range("1 day", "10 day", freq="D"),
        date_range(datetime(2005, 1, 1), datetime(2005, 1, 10), freq="D"),
    ],
)
@pytest.mark.parametrize("klass", [DataFrame, Series])
def test_asfreq(klass, index, freq):
    obj = klass(range(len(index)), index=index)
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
        timedelta_range("1 day", "10 day", freq="D"),
        date_range(datetime(2005, 1, 1), datetime(2005, 1, 10), freq="D"),
        period_range(datetime(2005, 1, 1), datetime(2005, 1, 10), freq="D"),
    ],
)
def test_resample_interpolate(index):
    # GH#12925
    df = DataFrame(range(len(index)), index=index)
    warn = None
    if isinstance(df.index, PeriodIndex):
        warn = FutureWarning
    msg = "Resampling with a PeriodIndex is deprecated"
    with tm.assert_produces_warning(warn, match=msg):
        result = df.resample("1min").asfreq().interpolate()
        expected = df.resample("1min").interpolate()
    tm.assert_frame_equal(result, expected)


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

    warn = None
    if isinstance(ser.index, PeriodIndex):
        warn = FutureWarning
    msg = "Resampling with a PeriodIndex is deprecated"
    with tm.assert_produces_warning(warn, match=msg):
        rs = ser.resample(freq)
    result = getattr(rs, resample_method)()

    if resample_method == "ohlc":
        expected = DataFrame(
            [], index=ser.index[:0].copy(), columns=["open", "high", "low", "close"]
        )
        expected.index = _asfreq_compat(ser.index, freq)
        tm.assert_frame_equal(result, expected, check_dtype=False)
    else:
        expected = ser.copy()
        expected.index = _asfreq_compat(ser.index, freq)
        tm.assert_series_equal(result, expected, check_dtype=False)

    tm.assert_index_equal(result.index, expected.index)
    assert result.index.freq == expected.index.freq


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

    msg = "Resampling with a PeriodIndex is deprecated"
    with tm.assert_produces_warning(FutureWarning, match=msg):
        rs = ser.resample(freq)
    result = getattr(rs, resample_method)()

    if resample_method == "ohlc":
        expected = DataFrame(
            [], index=ser.index[:0].copy(), columns=["open", "high", "low", "close"]
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

    warn = None
    if isinstance(ser.index, PeriodIndex):
        warn = FutureWarning
    msg = "Resampling with a PeriodIndex is deprecated"
    with tm.assert_produces_warning(warn, match=msg):
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

    warn = None
    if isinstance(df.index, PeriodIndex):
        warn = FutureWarning
    msg = "Resampling with a PeriodIndex is deprecated"
    with tm.assert_produces_warning(warn, match=msg):
        rs = df.resample(freq, group_keys=False)
    result = getattr(rs, resample_method)()
    if resample_method == "ohlc":
        # TODO: no tests with len(df.columns) > 0
        mi = MultiIndex.from_product([df.columns, ["open", "high", "low", "close"]])
        expected = DataFrame(
            [], index=df.index[:0].copy(), columns=mi, dtype=np.float64
        )
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

    warn = None
    if isinstance(empty_frame_dti.index, PeriodIndex):
        warn = FutureWarning
    msg = "Resampling with a PeriodIndex is deprecated"
    with tm.assert_produces_warning(warn, match=msg):
        rs = empty_frame_dti.resample(freq)
    result = rs.count()

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

    msg = "Resampling with a PeriodIndex"
    warn = None
    if isinstance(empty_frame_dti.index, PeriodIndex):
        warn = FutureWarning
    with tm.assert_produces_warning(warn, match=msg):
        rs = empty_frame_dti.resample(freq)
    result = rs.size()

    index = _asfreq_compat(empty_frame_dti.index, freq)

    expected = Series([], dtype="int64", index=index)

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "index",
    [
        PeriodIndex([], freq="M", name="a"),
        DatetimeIndex([], name="a"),
        TimedeltaIndex([], name="a"),
    ],
)
@pytest.mark.parametrize("dtype", [float, int, object, "datetime64[ns]"])
@pytest.mark.filterwarnings(r"ignore:PeriodDtype\[B\] is deprecated:FutureWarning")
def test_resample_empty_dtypes(index, dtype, resample_method):
    # Empty series were sometimes causing a segfault (for the functions
    # with Cython bounds-checking disabled) or an IndexError.  We just run
    # them to ensure they no longer do.  (GH #10228)
    warn = None
    if isinstance(index, PeriodIndex):
        # GH#53511
        index = PeriodIndex([], freq="B", name=index.name)
        warn = FutureWarning
    msg = "Resampling with a PeriodIndex is deprecated"

    empty_series_dti = Series([], index, dtype)
    with tm.assert_produces_warning(warn, match=msg):
        rs = empty_series_dti.resample("d", group_keys=False)
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

    msg = "Resampling with a PeriodIndex"
    warn = None
    if isinstance(ser.index, PeriodIndex):
        warn = FutureWarning

    with tm.assert_produces_warning(warn, match=msg):
        rs = ser.resample(freq, group_keys=False)

    result = rs.apply(lambda x: 1)
    with tm.assert_produces_warning(warn, match=msg):
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
    msg = "Resampling with a PeriodIndex"
    warn = None
    if isinstance(series.index, PeriodIndex):
        warn = FutureWarning

    with tm.assert_produces_warning(warn, match=msg):
        grouped = series.groupby(tg)

    with tm.assert_produces_warning(warn, match=msg):
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

    msg = "Resampling with a PeriodIndex"
    warn = None
    if isinstance(ser.index, PeriodIndex):
        warn = FutureWarning
    with tm.assert_produces_warning(warn, match=msg):
        result = ser.resample(freq).quantile(q)
        expected = ser.resample(freq).agg(lambda x: x.quantile(q)).rename(ser.name)
    tm.assert_series_equal(result, expected)
