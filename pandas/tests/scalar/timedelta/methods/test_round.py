import numpy as np
import pytest

from pandas._libs import lib
from pandas._libs.tslibs import (
    NaT,
    iNaT,
)
from pandas.errors import OutOfBoundsTimedelta

import pandas as pd


class TestTimedeltaRound:
    @pytest.mark.parametrize(
        "freq,s1,s2",
        [
            # This first case has s1, s2 being the same as t1,t2 below
            (
                "ns",
                "1 days 02:34:56.789123456",
                "-1 days 02:34:56.789123456",
            ),
            (
                "us",
                "1 days 02:34:56.789123000",
                "-1 days 02:34:56.789123000",
            ),
            (
                "ms",
                "1 days 02:34:56.789000000",
                "-1 days 02:34:56.789000000",
            ),
            ("s", "1 days 02:34:57", "-1 days 02:34:57"),
            ("2s", "1 days 02:34:56", "-1 days 02:34:56"),
            ("5s", "1 days 02:34:55", "-1 days 02:34:55"),
            ("min", "1 days 02:35:00", "-1 days 02:35:00"),
            ("12min", "1 days 02:36:00", "-1 days 02:36:00"),
            ("h", "1 days 03:00:00", "-1 days 03:00:00"),
            ("D", "1 days", "-1 days"),
        ],
    )
    def test_round(self, freq, s1, s2):
        s1 = pd.Timedelta(s1)
        s2 = pd.Timedelta(s2)
        t1 = pd.Timedelta("1 days 02:34:56.789123456")
        t2 = pd.Timedelta("-1 days 02:34:56.789123456")

        r1 = t1.round(freq)
        assert r1 == s1
        r2 = t2.round(freq)
        assert r2 == s2

    def test_round_invalid(self):
        t1 = pd.Timedelta("1 days 02:34:56.789123456")

        for freq, msg in [
            ("YE", "<YearEnd: month=12> is a non-fixed frequency"),
            ("ME", "<MonthEnd> is a non-fixed frequency"),
            ("foobar", "Invalid frequency: foobar"),
        ]:
            with pytest.raises(ValueError, match=msg):
                t1.round(freq)

    def test_round_implementation_bounds(self):
        # See also: analogous test for Timestamp
        # GH#38964
        result = pd.Timedelta.min.ceil("s")
        expected = pd.Timedelta.min + pd.Timedelta(seconds=1) - pd.Timedelta(145224193)
        assert result == expected

        result = pd.Timedelta.max.floor("s")
        expected = pd.Timedelta.max - pd.Timedelta(854775807)
        assert result == expected

        msg = (
            r"Cannot round -106752 days \+00:12:43.145224193 to freq=s without overflow"
        )
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            pd.Timedelta.min.floor("s")
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            pd.Timedelta.min.round("s")

        msg = "Cannot round 106751 days 23:47:16.854775807 to freq=s without overflow"
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            pd.Timedelta.max.ceil("s")
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            pd.Timedelta.max.round("s")

    @pytest.mark.parametrize(
        "val",
        [
            iNaT + 1,
            -1,
            0,
            1,
            lib.i8max,
            10**9 - 1,
            10**9,
            10**9 + 1,
            60 * 10**9 - 1,
            24 * 3600 * 10**9 - 1,
            -(10**9 + 1),
            -(60 * 10**9 - 1),
            -(24 * 3600 * 10**9 - 1),
        ],
    )
    @pytest.mark.parametrize(
        "method", [pd.Timedelta.round, pd.Timedelta.floor, pd.Timedelta.ceil]
    )
    def test_round_sanity(self, val, method):
        cls = pd.Timedelta
        err_cls = OutOfBoundsTimedelta

        val = np.int64(val)
        td = cls(val)

        def checker(ts, nanos, unit):
            # First check that we do raise in cases where we should
            if nanos == 1:
                pass
            else:
                div, mod = divmod(ts._value, nanos)
                diff = int(nanos - mod)
                lb = ts._value - mod
                assert lb <= ts._value  # i.e. no overflows with python ints
                ub = ts._value + diff
                assert ub > ts._value  # i.e. no overflows with python ints

                msg = "without overflow"
                if mod == 0:
                    # We should never be raising in this
                    pass
                elif method is cls.ceil:
                    if ub > cls.max._value:
                        with pytest.raises(err_cls, match=msg):
                            method(ts, unit)
                        return
                elif method is cls.floor:
                    if lb < cls.min._value:
                        with pytest.raises(err_cls, match=msg):
                            method(ts, unit)
                        return
                elif mod >= diff:
                    if ub > cls.max._value:
                        with pytest.raises(err_cls, match=msg):
                            method(ts, unit)
                        return
                elif lb < cls.min._value:
                    with pytest.raises(err_cls, match=msg):
                        method(ts, unit)
                    return

            res = method(ts, unit)

            td = res - ts
            diff = abs(td._value)
            assert diff < nanos
            assert res._value % nanos == 0

            if method is cls.round:
                assert diff <= nanos / 2
            elif method is cls.floor:
                assert res <= ts
            elif method is cls.ceil:
                assert res >= ts

        nanos = 1
        checker(td, nanos, "ns")

        nanos = 1000
        checker(td, nanos, "us")

        nanos = 1_000_000
        checker(td, nanos, "ms")

        nanos = 1_000_000_000
        checker(td, nanos, "s")

        nanos = 60 * 1_000_000_000
        checker(td, nanos, "min")

        nanos = 60 * 60 * 1_000_000_000
        checker(td, nanos, "h")

        nanos = 24 * 60 * 60 * 1_000_000_000
        checker(td, nanos, "D")

    def test_round_non_nano(self, unit):
        td = pd.Timedelta("1 days 02:34:57").as_unit(unit)

        res = td.round("min")
        assert res == pd.Timedelta("1 days 02:35:00")
        assert res._creso == td._creso

        res = td.floor("min")
        assert res == pd.Timedelta("1 days 02:34:00")
        assert res._creso == td._creso

        res = td.ceil("min")
        assert res == pd.Timedelta("1 days 02:35:00")
        assert res._creso == td._creso

    @pytest.mark.parametrize(
        "timedelta,frequency,expected_ceil,expected_round,expected_floor",
        [
            (
                pd.Timedelta("1001ms"),
                pd.Timedelta("1s"),
                pd.Timedelta("2s"),
                pd.Timedelta("1s"),
                pd.Timedelta("1s"),
            ),
            (
                pd.Timedelta("1001ms"),
                pd.Timedelta("1ms"),
                pd.Timedelta("1001ms"),
                pd.Timedelta("1001ms"),
                pd.Timedelta("1001ms"),
            ),
            (
                pd.Timedelta("1 days 2 min 3 us 42 ns"),
                pd.Timedelta("1s"),
                pd.Timedelta("1 days 2 min 1s"),
                pd.Timedelta("1 days 2 min"),
                pd.Timedelta("1 days 2 min"),
            ),
            (
                pd.Timedelta("5 hours 9 minutes 15.13 seconds"),
                pd.Timedelta("1 hour"),
                pd.Timedelta("6 hours"),
                pd.Timedelta("5 hours"),
                pd.Timedelta("5 hours"),
            ),
            (
                pd.Timedelta("5 hours 9 minutes 15.13 seconds"),
                pd.Timedelta("1 hour 30 min"),
                pd.Timedelta("6 hours"),
                pd.Timedelta("4 hours 30 minutes"),
                pd.Timedelta("4 hours 30 minutes"),
            ),
            # Edge cases derived from TestTimestampRound.test_ceil_floor_edge
            (
                pd.Timedelta("1 days 45 seconds"),
                pd.Timedelta("15s"),
                pd.Timedelta("1 days 45 seconds"),
                pd.Timedelta("1 days 45 seconds"),
                pd.Timedelta("1 days 45 seconds"),
            ),
            (
                pd.Timedelta("1 days 45.000000012 seconds"),
                pd.Timedelta("10ns"),
                pd.Timedelta("1 days 45.000000020 seconds"),
                pd.Timedelta("1 days 45.000000010 seconds"),
                pd.Timedelta("1 days 45.000000010 seconds"),
            ),
            (
                pd.Timedelta("1 days 1.000000012 seconds"),
                pd.Timedelta("10ns"),
                pd.Timedelta("1 days 1.000000020 seconds"),
                pd.Timedelta("1 days 1.000000010 seconds"),
                pd.Timedelta("1 days 1.000000010 seconds"),
            ),
        ],
    )
    def test_rounding_with_timedelta_freq(
        self, timedelta, frequency, expected_ceil, expected_round, expected_floor
    ):
        # GH#63687 - Timedelta rounding methods accept Timedelta arguments
        assert timedelta.ceil(frequency) == expected_ceil
        assert timedelta.round(frequency) == expected_round
        assert timedelta.floor(frequency) == expected_floor

    def test_rounding_nat_frequency(self):
        td = pd.Timedelta("1001ms")

        with pytest.raises(TypeError, match="Argument 'freq' has incorrect type"):
            td.ceil(NaT)
        with pytest.raises(TypeError, match="Argument 'freq' has incorrect type"):
            td.floor(NaT)
        with pytest.raises(TypeError, match="Argument 'freq' has incorrect type"):
            td.round(NaT)

    def test_rounding_nat_timedelta(self):
        freq = pd.Timedelta("1s")

        assert NaT.ceil(freq) is NaT
        assert NaT.floor(freq) is NaT
        assert NaT.round(freq) is NaT

    def test_round_freq_finer_than_resolution(self):
        # GH#64828
        td = pd.Timedelta(1.0, unit="days").as_unit("s")
        assert td.unit == "s"
        assert td.round("100ms") == pd.Timedelta("1 days 00:00:00")
        assert td.floor("100ms") == pd.Timedelta("1 days 00:00:00")
        assert td.ceil("100ms") == pd.Timedelta("1 days 00:00:00")
