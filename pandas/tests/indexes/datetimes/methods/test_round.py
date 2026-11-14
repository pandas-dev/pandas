import numpy as np
import pytest

from pandas._libs.tslibs import to_offset
from pandas._libs.tslibs.offsets import INVALID_FREQ_ERR_MSG
from pandas._libs.tslibs.timezones import tz_compare

import pandas as pd
import pandas._testing as tm


class TestDatetimeIndexRound:
    def test_round_daily(self):
        dti = pd.date_range("20130101 09:10:11", periods=5)
        result = dti.round("D")
        expected = pd.date_range("20130101", periods=5)
        tm.assert_index_equal(result, expected, check_freq=False)

        dti = dti.tz_localize("UTC").tz_convert("US/Eastern")
        result = dti.round("D")
        expected = pd.date_range("20130101", periods=5).tz_localize("US/Eastern")
        tm.assert_index_equal(result, expected, check_freq=False)

        result = dti.round("s")
        tm.assert_index_equal(result, dti)

    @pytest.mark.parametrize(
        "freq, error_msg",
        [
            ("YE", "<YearEnd: month=12> is a non-fixed frequency"),
            ("ME", "<MonthEnd> is a non-fixed frequency"),
            ("foobar", "Invalid frequency: foobar"),
        ],
    )
    def test_round_invalid(self, freq, error_msg):
        dti = pd.date_range("20130101 09:10:11", periods=5)
        dti = dti.tz_localize("UTC").tz_convert("US/Eastern")
        with pytest.raises(ValueError, match=error_msg):
            dti.round(freq)

    def test_round(self, tz_naive_fixture, unit):
        tz = tz_naive_fixture
        rng = pd.date_range(
            start="2016-01-01", periods=5, freq="30Min", tz=tz, unit=unit
        )
        elt = rng[1]

        expected_rng = pd.DatetimeIndex(
            [
                pd.Timestamp("2016-01-01 00:00:00", tz=tz),
                pd.Timestamp("2016-01-01 00:00:00", tz=tz),
                pd.Timestamp("2016-01-01 01:00:00", tz=tz),
                pd.Timestamp("2016-01-01 02:00:00", tz=tz),
                pd.Timestamp("2016-01-01 02:00:00", tz=tz),
            ]
        ).as_unit(unit)
        expected_elt = expected_rng[1]

        result = rng.round(freq="h")
        tm.assert_index_equal(result, expected_rng)
        assert elt.round(freq="h") == expected_elt

        msg = INVALID_FREQ_ERR_MSG
        with pytest.raises(ValueError, match=msg):
            rng.round(freq="foo")
        with pytest.raises(ValueError, match=msg):
            elt.round(freq="foo")

        msg = "<MonthEnd> is a non-fixed frequency"
        with pytest.raises(ValueError, match=msg):
            rng.round(freq="ME")
        with pytest.raises(ValueError, match=msg):
            elt.round(freq="ME")

    def test_round2(self, tz_naive_fixture):
        tz = tz_naive_fixture
        # GH#14440 & GH#15578
        index = pd.DatetimeIndex(["2016-10-17 12:00:00.0015"], tz=tz).as_unit("ns")
        result = index.round("ms")
        expected = pd.DatetimeIndex(["2016-10-17 12:00:00.002000"], tz=tz).as_unit("ns")
        tm.assert_index_equal(result, expected)

        for freq in ["us", "ns"]:
            tm.assert_index_equal(index, index.round(freq))

    def test_round3(self, tz_naive_fixture):
        tz = tz_naive_fixture
        index = pd.DatetimeIndex(["2016-10-17 12:00:00.00149"], tz=tz).as_unit("ns")
        result = index.round("ms")
        expected = pd.DatetimeIndex(["2016-10-17 12:00:00.001000"], tz=tz).as_unit("ns")
        tm.assert_index_equal(result, expected)

    def test_round4(self, tz_naive_fixture):
        index = pd.DatetimeIndex(["2016-10-17 12:00:00.001501031"], dtype="M8[ns]")
        result = index.round("10ns")
        expected = pd.DatetimeIndex(["2016-10-17 12:00:00.001501030"], dtype="M8[ns]")
        tm.assert_index_equal(result, expected)

        ts = "2016-10-17 12:00:00.001501031"
        dti = pd.DatetimeIndex([ts], dtype="M8[ns]")
        with tm.assert_produces_warning(False):
            dti.round("1010ns")

    @pytest.mark.parametrize("method", ["round", "floor", "ceil"])
    @pytest.mark.parametrize(
        "unit, freq, freq_unit",
        [
            ("s", "700ms", "ms"),
            ("s", "2500ms", "ms"),
            ("ms", "3us", "us"),
            ("us", "300ns", "ns"),
        ],
    )
    def test_round_freq_not_multiple_of_resolution(self, method, unit, freq, freq_unit):
        # GH#67978 freq is neither a multiple nor a divisor of one unit of the
        #  index resolution, so the result would not be representable.  The
        #  value is deliberately off the freq grid so the second half below is
        #  not a no-op for any of the parametrizations.
        dti = pd.DatetimeIndex(["2020-01-01 00:00:06.500000001"])

        msg = rf"freq=.* is incompatible with unit={unit}"
        with pytest.raises(ValueError, match=msg):
            getattr(dti.as_unit(unit), method)(freq)

        finer = dti.as_unit(freq_unit)
        result = getattr(finer, method)(freq)
        assert not result.equals(finer)
        tm.assert_index_equal(
            result, getattr(finer.as_unit("ns"), method)(freq).as_unit(freq_unit)
        )

    @pytest.mark.parametrize("method", ["round", "floor", "ceil"])
    @pytest.mark.parametrize("unit, freq", [("s", "250ms"), ("us", "250ns")])
    def test_round_freq_divides_resolution(self, method, unit, freq):
        # GH#67978 freq divides one unit of the index resolution evenly, so
        #  every representable value is already a multiple of freq
        dti = pd.DatetimeIndex(["2020-01-01 00:00:06"]).as_unit(unit)
        tm.assert_index_equal(getattr(dti, method)(freq), dti)

    def test_no_rounding_occurs(self, tz_naive_fixture):
        # GH 21262
        tz = tz_naive_fixture
        rng = pd.date_range(
            start="2016-01-01", periods=5, freq="2Min", tz=tz, unit="ns"
        )

        expected_rng = pd.DatetimeIndex(
            [
                pd.Timestamp("2016-01-01 00:00:00", tz=tz),
                pd.Timestamp("2016-01-01 00:02:00", tz=tz),
                pd.Timestamp("2016-01-01 00:04:00", tz=tz),
                pd.Timestamp("2016-01-01 00:06:00", tz=tz),
                pd.Timestamp("2016-01-01 00:08:00", tz=tz),
            ]
        ).as_unit("ns")

        result = rng.round(freq="2min")
        tm.assert_index_equal(result, expected_rng)

    @pytest.mark.parametrize(
        "test_input, rounder, freq, expected",
        [
            (["2117-01-01 00:00:45"], "floor", "15s", ["2117-01-01 00:00:45"]),
            (["2117-01-01 00:00:45"], "ceil", "15s", ["2117-01-01 00:00:45"]),
            (
                ["2117-01-01 00:00:45.000000012"],
                "floor",
                "10ns",
                ["2117-01-01 00:00:45.000000010"],
            ),
            (
                ["1823-01-01 00:00:01.000000012"],
                "ceil",
                "10ns",
                ["1823-01-01 00:00:01.000000020"],
            ),
            (["1823-01-01 00:00:01"], "floor", "1s", ["1823-01-01 00:00:01"]),
            (["1823-01-01 00:00:01"], "ceil", "1s", ["1823-01-01 00:00:01"]),
            (["2018-01-01 00:15:00"], "ceil", "15min", ["2018-01-01 00:15:00"]),
            (["2018-01-01 00:15:00"], "floor", "15min", ["2018-01-01 00:15:00"]),
            (["1823-01-01 03:00:00"], "ceil", "3h", ["1823-01-01 03:00:00"]),
            (["1823-01-01 03:00:00"], "floor", "3h", ["1823-01-01 03:00:00"]),
            (
                ("NaT", "1823-01-01 00:00:01"),
                "floor",
                "1s",
                ("NaT", "1823-01-01 00:00:01"),
            ),
            (
                ("NaT", "1823-01-01 00:00:01"),
                "ceil",
                "1s",
                ("NaT", "1823-01-01 00:00:01"),
            ),
        ],
    )
    def test_ceil_floor_edge(self, test_input, rounder, freq, expected):
        dt = pd.DatetimeIndex(list(test_input))
        func = getattr(dt, rounder)
        result = func(freq)
        expected = pd.DatetimeIndex(list(expected))
        assert expected.equals(result)

    @pytest.mark.parametrize(
        "start, index_freq, periods",
        [("2018-01-01", "12h", 25), ("2018-01-01 0:0:0.124999", "1ns", 1000)],
    )
    @pytest.mark.parametrize(
        "round_freq",
        [
            "2ns",
            "3ns",
            "4ns",
            "5ns",
            "6ns",
            "7ns",
            "250ns",
            "500ns",
            "750ns",
            "1us",
            "19us",
            "250us",
            "500us",
            "750us",
            "1s",
            "2s",
            "3s",
            "12h",
            "1D",
        ],
    )
    def test_round_int64(self, start, index_freq, periods, round_freq):
        dt = pd.date_range(start=start, freq=index_freq, periods=periods, unit="ns")
        unit = to_offset(round_freq).nanos

        # test floor
        result = dt.floor(round_freq)
        diff = dt.asi8 - result.asi8
        mod = result.asi8 % unit
        assert (mod == 0).all(), f"floor not a {round_freq} multiple"
        assert (0 <= diff).all() and (diff < unit).all(), "floor error"

        # test ceil
        result = dt.ceil(round_freq)
        diff = result.asi8 - dt.asi8
        mod = result.asi8 % unit
        assert (mod == 0).all(), f"ceil not a {round_freq} multiple"
        assert (0 <= diff).all() and (diff < unit).all(), "ceil error"

        # test round
        result = dt.round(round_freq)
        diff = abs(result.asi8 - dt.asi8)
        mod = result.asi8 % unit
        assert (mod == 0).all(), f"round not a {round_freq} multiple"
        assert (diff <= unit // 2).all(), "round error"
        if unit % 2 == 0:
            assert (result.asi8[diff == unit // 2] % 2 == 0).all(), (
                "round half to even error"
            )


@pytest.mark.parametrize(
    "freq, values, floors, ceils",
    [
        (
            "ME",
            ["2017-01-23 12:34:56", "2017-01-31", "2017-01-31 00:00:01"],
            ["2016-12-31", "2017-01-31", "2017-01-31"],
            ["2017-01-31", "2017-01-31", "2017-02-28"],
        ),
        (
            "MS",
            ["2017-01-23 12:34:56", "2017-01-01", "2016-12-31 23:59:59"],
            ["2017-01-01", "2017-01-01", "2016-12-01"],
            ["2017-02-01", "2017-01-01", "2017-01-01"],
        ),
        (
            "QE-NOV",
            ["2017-01-23 12:34:56", "2017-02-28 00:00:01"],
            ["2016-11-30", "2017-02-28"],
            ["2017-02-28", "2017-05-31"],
        ),
        (
            "W-SUN",
            ["2017-01-23 12:34:56", "2017-01-22"],
            ["2017-01-22", "2017-01-22"],
            ["2017-01-29", "2017-01-22"],
        ),
        (
            "B",
            ["2017-01-23 12:34:56", "2017-01-28 09:00"],
            ["2017-01-23", "2017-01-27"],
            ["2017-01-24", "2017-01-30"],
        ),
        ("YS", ["2017-01-23 12:34:56"], ["2017-01-01"], ["2018-01-01"]),
        ("SMS", ["2017-01-23 12:34:56"], ["2017-01-15"], ["2017-02-01"]),
        ("BME", ["2017-01-23 12:34:56"], ["2016-12-30"], ["2017-01-31"]),
        # the offsets below have no vectorized arithmetic
        ("WOM-1WED", ["2017-01-23 12:34:56"], ["2017-01-04"], ["2017-02-01"]),
        ("LWOM-MON", ["2017-01-23 12:34:56"], ["2016-12-26"], ["2017-01-30"]),
        ("CBMS", ["2017-01-23 12:34:56"], ["2017-01-02"], ["2017-02-01"]),
        (
            pd.offsets.CustomBusinessDay(holidays=["2017-01-23"]),
            ["2017-01-23 12:34:56"],
            ["2017-01-20"],
            ["2017-01-24"],
        ),
        (
            pd.offsets.FY5253(weekday=4, startingMonth=12, variation="last"),
            ["2017-01-23 12:34:56"],
            ["2016-12-30"],
            ["2017-12-29"],
        ),
        (pd.offsets.Easter(), ["2017-01-23 12:34:56"], ["2016-03-27"], ["2017-04-16"]),
    ],
)
def test_floor_ceil_non_fixed_freq(freq, values, floors, ceils, tz_naive_fixture, unit):
    # GH#15303
    tz = tz_naive_fixture
    dti = pd.DatetimeIndex([*values, pd.NaT], tz=tz).as_unit(unit)

    for method, expected_values in [("floor", floors), ("ceil", ceils)]:
        expected = pd.DatetimeIndex([*expected_values, pd.NaT], tz=tz).as_unit(unit)
        result = getattr(dti, method)(freq)
        tm.assert_index_equal(result, expected)
        result = getattr(pd.Series(dti).dt, method)(freq)
        tm.assert_series_equal(result, pd.Series(expected))
        for ts, exp in zip(dti[:-1], expected[:-1], strict=True):
            res = getattr(ts, method)(freq)
            assert res == exp
            assert res.unit == exp.unit
            assert tz_compare(res.tz, exp.tz)


@pytest.mark.parametrize(
    "method, value", [("floor", "2017-08-13 12:00"), ("ceil", "2017-08-12 12:00")]
)
def test_floor_ceil_non_fixed_freq_nonexistent(method, value):
    # GH#15303 both snap to Sunday 2017-08-13 00:00, which does not exist here
    dti = pd.DatetimeIndex([value], tz="America/Santiago")
    with pytest.raises(ValueError, match="nonexistent"):
        getattr(dti, method)("W-SUN")
    with pytest.raises(ValueError, match="nonexistent"):
        getattr(dti[0], method)("W-SUN")

    expected = pd.DatetimeIndex(["2017-08-13 01:00"], tz="America/Santiago")
    result = getattr(dti, method)("W-SUN", nonexistent="shift_forward")
    tm.assert_index_equal(result, expected)
    assert getattr(dti[0], method)("W-SUN", nonexistent="shift_forward") == expected[0]


def test_floor_non_fixed_freq_last_ns_day():
    # GH#15303 the next midnight after 2262-04-11 is out of nanosecond bounds
    dti = pd.DatetimeIndex(["2262-04-11 12:00", "2262-04-10 12:00", pd.NaT]).as_unit(
        "ns"
    )
    expected = pd.DatetimeIndex(["2262-04-01", "2262-04-01", pd.NaT]).as_unit("ns")
    tm.assert_index_equal(dti.floor("MS"), expected)
    assert dti[0].floor("MS") == expected[0]


@pytest.mark.parametrize(
    "freq",
    [
        "2ME",
        "-1ME",
        "0MS",
        pd.offsets.Week(),
        pd.offsets.Week(n=-1, weekday=2),
        pd.offsets.DateOffset(months=1),
        pd.offsets.BusinessHour(),
    ],
)
@pytest.mark.parametrize("method", ["floor", "ceil"])
def test_floor_ceil_unsupported_non_fixed_freq_raises(freq, method):
    # GH#15303
    dti = pd.DatetimeIndex(["2017-01-23 12:34:56"])
    msg = "floor and ceil support only a single step of a date-anchored frequency"
    with pytest.raises(ValueError, match=msg):
        getattr(dti, method)(freq)
    with pytest.raises(ValueError, match=msg):
        getattr(dti[0], method)(freq)


@pytest.mark.parametrize(
    "freq",
    [
        pd.offsets.BusinessDay(offset=pd.Timedelta(hours=1)),
        pd.offsets.CustomBusinessDay(offset=pd.Timedelta(days=1)),
        pd.offsets.CustomBusinessMonthEnd(offset=pd.Timedelta(hours=1)),
    ],
)
@pytest.mark.parametrize("method", ["floor", "ceil"])
def test_floor_ceil_business_offset_with_offset_raises(freq, method):
    # GH#15303
    dti = pd.DatetimeIndex(["2017-01-23 12:34:56"])
    msg = "floor and ceil do not support .* with a nonzero 'offset'"
    with pytest.raises(ValueError, match=msg):
        getattr(dti, method)(freq)
    with pytest.raises(ValueError, match=msg):
        getattr(dti[0], method)(freq)


@pytest.mark.parametrize("is_dst, utc_offset", [(True, "-04:00"), (False, "-05:00")])
def test_floor_non_fixed_freq_ambiguous(is_dst, utc_offset):
    # GH#15303 Sunday 2017-11-05 00:00 is ambiguous here
    tz = "America/Havana"
    dti = pd.DatetimeIndex(["2017-11-05 12:00"], tz=tz)
    msg = "Cannot infer dst time"
    with pytest.raises(ValueError, match=msg):
        dti.floor("W-SUN")
    with pytest.raises(ValueError, match=msg):
        dti[0].floor("W-SUN")

    expected = pd.Timestamp(f"2017-11-05 00:00{utc_offset}").tz_convert(
        "America/Havana"
    )
    result = dti.floor("W-SUN", ambiguous=np.array([is_dst]))
    tm.assert_index_equal(result, pd.DatetimeIndex([expected]))
    assert dti[0].floor("W-SUN", ambiguous=is_dst) == expected
