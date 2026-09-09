from datetime import timedelta

import numpy as np
import pytest

from pandas.core.dtypes.common import is_integer

import pandas as pd
import pandas._testing as tm

from pandas.tseries.offsets import Day


@pytest.fixture(params=[None, "foo"])
def name(request):
    return request.param


class TestIntervalRange:
    @pytest.mark.parametrize("freq, periods", [(1, 100), (2.5, 40), (5, 20), (25, 4)])
    def test_constructor_numeric(self, closed, name, freq, periods):
        start, end = 0, 100
        breaks = np.arange(101, step=freq)
        if breaks.dtype.kind == "f":
            breaks = breaks.astype(np.float64)
        else:
            breaks = breaks.astype(np.int64)
        expected = pd.IntervalIndex.from_breaks(breaks, name=name, closed=closed)

        # defined from start/end/freq
        result = pd.interval_range(
            start=start, end=end, freq=freq, name=name, closed=closed
        )
        tm.assert_index_equal(result, expected)

        # defined from start/periods/freq
        result = pd.interval_range(
            start=start, periods=periods, freq=freq, name=name, closed=closed
        )
        tm.assert_index_equal(result, expected)

        # defined from end/periods/freq
        result = pd.interval_range(
            end=end, periods=periods, freq=freq, name=name, closed=closed
        )
        tm.assert_index_equal(result, expected)

        # GH 20976: linspace behavior defined from start/end/periods
        result = pd.interval_range(
            start=start, end=end, periods=periods, name=name, closed=closed
        )
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("tz", [None, "US/Eastern"])
    @pytest.mark.parametrize(
        "freq, periods", [("D", 364), ("2D", 182), ("22D18h", 16), ("ME", 11)]
    )
    def test_constructor_timestamp(self, closed, name, freq, periods, tz):
        start, end = pd.Timestamp("20180101", tz=tz), pd.Timestamp("20181231", tz=tz)
        breaks = pd.date_range(start=start, end=end, freq=freq)
        expected = pd.IntervalIndex.from_breaks(breaks, name=name, closed=closed)

        # defined from start/end/freq
        result = pd.interval_range(
            start=start, end=end, freq=freq, name=name, closed=closed
        )
        tm.assert_index_equal(result, expected)

        # defined from start/periods/freq
        result = pd.interval_range(
            start=start, periods=periods, freq=freq, name=name, closed=closed
        )
        tm.assert_index_equal(result, expected)

        # defined from end/periods/freq
        result = pd.interval_range(
            end=end, periods=periods, freq=freq, name=name, closed=closed
        )
        tm.assert_index_equal(result, expected)

        # GH 20976: linspace behavior defined from start/end/periods
        if not breaks.freq.n == 1 and tz is None:
            result = pd.interval_range(
                start=start, end=end, periods=periods, name=name, closed=closed
            )
            tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "freq, periods", [("D", 100), ("2D12h", 40), ("5D", 20), ("25D", 4)]
    )
    def test_constructor_timedelta(self, closed, name, freq, periods):
        start, end = pd.Timedelta("0 days"), pd.Timedelta("100 days")
        breaks = pd.timedelta_range(start=start, end=end, freq=freq)
        expected = pd.IntervalIndex.from_breaks(breaks, name=name, closed=closed)

        # defined from start/end/freq
        result = pd.interval_range(
            start=start, end=end, freq=freq, name=name, closed=closed
        )
        tm.assert_index_equal(result, expected)

        # defined from start/periods/freq
        result = pd.interval_range(
            start=start, periods=periods, freq=freq, name=name, closed=closed
        )
        tm.assert_index_equal(result, expected)

        # defined from end/periods/freq
        result = pd.interval_range(
            end=end, periods=periods, freq=freq, name=name, closed=closed
        )
        tm.assert_index_equal(result, expected)

        # GH 20976: linspace behavior defined from start/end/periods
        result = pd.interval_range(
            start=start, end=end, periods=periods, name=name, closed=closed
        )
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "start, end, freq, expected_endpoint",
        [
            (0, 10, 3, 9),
            (0, 10, 1.5, 9),
            (0.5, 10, 3, 9.5),
            (pd.Timedelta("0D"), pd.Timedelta("10D"), "2D4h", pd.Timedelta("8D16h")),
            (
                pd.Timestamp("2018-01-01"),
                pd.Timestamp("2018-02-09"),
                "MS",
                pd.Timestamp("2018-02-01"),
            ),
            (
                pd.Timestamp("2018-01-01", tz="US/Eastern"),
                pd.Timestamp("2018-01-20", tz="US/Eastern"),
                "5D12h",
                pd.Timestamp("2018-01-17 12:00:00", tz="US/Eastern"),
            ),
        ],
    )
    def test_early_truncation(self, start, end, freq, expected_endpoint):
        # index truncates early if freq causes end to be skipped
        result = pd.interval_range(start=start, end=end, freq=freq)
        result_endpoint = result.right[-1]
        assert result_endpoint == expected_endpoint

    @pytest.mark.parametrize(
        "start, end, freq",
        [(0.5, None, None), (None, 4.5, None), (0.5, None, 1.5), (None, 6.5, 1.5)],
    )
    def test_no_invalid_float_truncation(self, start, end, freq):
        # GH 21161
        if freq is None:
            breaks = [0.5, 1.5, 2.5, 3.5, 4.5]
        else:
            breaks = [0.5, 2.0, 3.5, 5.0, 6.5]
        expected = pd.IntervalIndex.from_breaks(breaks)

        result = pd.interval_range(start=start, end=end, periods=4, freq=freq)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "start, mid, end",
        [
            (
                pd.Timestamp("2018-03-10", tz="US/Eastern"),
                pd.Timestamp("2018-03-10 23:30:00", tz="US/Eastern"),
                pd.Timestamp("2018-03-12", tz="US/Eastern"),
            ),
            (
                pd.Timestamp("2018-11-03", tz="US/Eastern"),
                pd.Timestamp("2018-11-04 00:30:00", tz="US/Eastern"),
                pd.Timestamp("2018-11-05", tz="US/Eastern"),
            ),
        ],
    )
    def test_linspace_dst_transition(self, start, mid, end):
        # GH 20976: linspace behavior defined from start/end/periods
        # accounts for the hour gained/lost during DST transition
        start = start.as_unit("ns")
        mid = mid.as_unit("ns")
        end = end.as_unit("ns")
        result = pd.interval_range(start=start, end=end, periods=2)
        expected = pd.IntervalIndex.from_breaks([start, mid, end])
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("freq", [2, 2.0])
    @pytest.mark.parametrize("end", [10, 10.0])
    @pytest.mark.parametrize("start", [0, 0.0])
    def test_float_subtype(self, start, end, freq):
        # Has float subtype if any of start/end/freq are float, even if all
        # resulting endpoints can safely be upcast to integers

        # defined from start/end/freq
        index = pd.interval_range(start=start, end=end, freq=freq)
        result = index.dtype.subtype
        expected = "int64" if is_integer(start + end + freq) else "float64"
        assert result == expected

        # defined from start/periods/freq
        index = pd.interval_range(start=start, periods=5, freq=freq)
        result = index.dtype.subtype
        expected = "int64" if is_integer(start + freq) else "float64"
        assert result == expected

        # defined from end/periods/freq
        index = pd.interval_range(end=end, periods=5, freq=freq)
        result = index.dtype.subtype
        expected = "int64" if is_integer(end + freq) else "float64"
        assert result == expected

        # GH 20976: linspace behavior defined from start/end/periods
        index = pd.interval_range(start=start, end=end, periods=5)
        result = index.dtype.subtype
        expected = "int64" if is_integer(start + end) else "float64"
        assert result == expected

    @pytest.mark.parametrize(
        "start, end, expected",
        [
            (np.int8(1), np.int8(10), np.dtype("int8")),
            (np.int8(1), np.float16(10), np.dtype("float64")),
            (np.float32(1), np.float32(10), np.dtype("float32")),
            (1, 10, np.dtype("int64")),
            (1, 10.0, np.dtype("float64")),
        ],
    )
    def test_interval_dtype(self, start, end, expected):
        result = pd.interval_range(start=start, end=end).dtype.subtype
        assert result == expected

    def test_interval_range_fractional_period(self):
        # float value for periods
        msg = "periods must be an integer, got 10.5"
        ts = pd.Timestamp("2024-03-25")
        with pytest.raises(TypeError, match=msg):
            pd.interval_range(ts, periods=10.5)

    def test_constructor_coverage(self):
        # equivalent timestamp-like start/end
        start, end = pd.Timestamp("2017-01-01"), pd.Timestamp("2017-01-15")
        expected = pd.interval_range(start=start, end=end)

        result = pd.interval_range(start=start.to_pydatetime(), end=end.to_pydatetime())
        tm.assert_index_equal(result, expected)

        result = pd.interval_range(start=start.asm8, end=end.asm8)
        tm.assert_index_equal(result, expected)

        # equivalent freq with timestamp
        equiv_freq = [
            "D",
            Day(),
            pd.Timedelta(days=1),
            timedelta(days=1),
            pd.DateOffset(days=1),
        ]
        for freq in equiv_freq:
            result = pd.interval_range(start=start, end=end, freq=freq)
            tm.assert_index_equal(result, expected)

        # equivalent timedelta-like start/end
        start, end = (
            pd.Timedelta(days=1).as_unit("us"),
            pd.Timedelta(days=10).as_unit("us"),
        )
        expected = pd.interval_range(start=start, end=end)

        result = pd.interval_range(
            start=start.to_pytimedelta(), end=end.to_pytimedelta()
        )
        tm.assert_index_equal(result, expected)

        result = pd.interval_range(start=start.asm8, end=end.asm8)
        tm.assert_index_equal(result, expected)

        # equivalent freq with timedelta
        equiv_freq = ["D", Day(), pd.Timedelta(days=1), timedelta(days=1)]
        for freq in equiv_freq:
            result = pd.interval_range(start=start, end=end, freq=freq)
            tm.assert_index_equal(result, expected)

    def test_errors(self):
        # not enough params
        msg = (
            "Of the four parameters: start, end, periods, and freq, "
            "exactly three must be specified"
        )

        with pytest.raises(ValueError, match=msg):
            pd.interval_range(start=0)

        with pytest.raises(ValueError, match=msg):
            pd.interval_range(end=5)

        with pytest.raises(ValueError, match=msg):
            pd.interval_range(periods=2)

        with pytest.raises(ValueError, match=msg):
            pd.interval_range()

        # too many params
        with pytest.raises(ValueError, match=msg):
            pd.interval_range(start=0, end=5, periods=6, freq=1.5)

        # mixed units
        msg = "start, end, freq need to be type compatible"
        with pytest.raises(TypeError, match=msg):
            pd.interval_range(start=0, end=pd.Timestamp("20130101"), freq=2)

        with pytest.raises(TypeError, match=msg):
            pd.interval_range(start=0, end=pd.Timedelta("1 day"), freq=2)

        with pytest.raises(TypeError, match=msg):
            pd.interval_range(start=0, end=10, freq="D")

        with pytest.raises(TypeError, match=msg):
            pd.interval_range(start=pd.Timestamp("20130101"), end=10, freq="D")

        with pytest.raises(TypeError, match=msg):
            pd.interval_range(
                start=pd.Timestamp("20130101"), end=pd.Timedelta("1 day"), freq="D"
            )

        with pytest.raises(TypeError, match=msg):
            pd.interval_range(
                start=pd.Timestamp("20130101"), end=pd.Timestamp("20130110"), freq=2
            )

        with pytest.raises(TypeError, match=msg):
            pd.interval_range(start=pd.Timedelta("1 day"), end=10, freq="D")

        with pytest.raises(TypeError, match=msg):
            pd.interval_range(
                start=pd.Timedelta("1 day"), end=pd.Timestamp("20130110"), freq="D"
            )

        with pytest.raises(TypeError, match=msg):
            pd.interval_range(
                start=pd.Timedelta("1 day"), end=pd.Timedelta("10 days"), freq=2
            )

        # invalid periods
        msg = "periods must be an integer, got foo"
        with pytest.raises(TypeError, match=msg):
            pd.interval_range(start=0, periods="foo")

        # invalid start
        msg = "start must be numeric or datetime-like, got foo"
        with pytest.raises(ValueError, match=msg):
            pd.interval_range(start="foo", periods=10)

        # invalid end
        msg = r"end must be numeric or datetime-like, got \(0, 1\]"
        with pytest.raises(ValueError, match=msg):
            pd.interval_range(end=pd.Interval(0, 1), periods=10)

        # invalid freq for datetime-like
        msg = "freq must be numeric or convertible to DateOffset, got foo"
        with pytest.raises(ValueError, match=msg):
            pd.interval_range(start=0, end=10, freq="foo")

        with pytest.raises(ValueError, match=msg):
            pd.interval_range(start=pd.Timestamp("20130101"), periods=10, freq="foo")

        with pytest.raises(ValueError, match=msg):
            pd.interval_range(end=pd.Timedelta("1 day"), periods=10, freq="foo")

        # mixed tz
        start = pd.Timestamp("2017-01-01", tz="US/Eastern")
        end = pd.Timestamp("2017-01-07", tz="US/Pacific")
        msg = "Start and end cannot both be tz-aware with different timezones"
        with pytest.raises(TypeError, match=msg):
            pd.interval_range(start=start, end=end)

    def test_float_freq(self):
        # GH 54477
        result = pd.interval_range(0, 1, freq=0.1)
        expected = pd.IntervalIndex.from_breaks([0 + 0.1 * n for n in range(11)])
        tm.assert_index_equal(result, expected)

        result = pd.interval_range(0, 1, freq=0.6)
        expected = pd.IntervalIndex.from_breaks([0, 0.6])
        tm.assert_index_equal(result, expected)

    def test_interval_range_float32_start_int_freq(self):
        # GH 58964
        result = pd.interval_range(start=np.float32(0), end=2, freq=1)
        expected = pd.IntervalIndex.from_tuples(
            [(0.0, 1.0), (1.0, 2.0)], dtype="interval[float64, right]"
        )
        tm.assert_index_equal(result, expected)
