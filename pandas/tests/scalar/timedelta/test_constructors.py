from datetime import timedelta
from itertools import (
    chain,
    zip_longest,
)
import re

import numpy as np
import pytest

from pandas._libs.tslibs import OutOfBoundsTimedelta

from pandas import (
    NA,
    NaT,
    Timedelta,
    offsets,
    to_timedelta,
)

TD_KWARGS_UNITS = {
    "weeks": ("w",),
    "days": ("d", "day", "days"),
    "hours": ("h", "hr", "hour", "hours"),
    "minutes": ("m", "t", "min", "minute", "minutes"),
    "seconds": ("s", "sec", "second", "seconds"),
    "milliseconds": ("l", "ms", "milli", "millis", "millisecond", "milliseconds"),
    "microseconds": ("u", "us", "µs", "micro", "micros", "microsecond", "microseconds"),
    "nanoseconds": ("n", "ns", "nano", "nanos", "nanosecond", "nanoseconds"),
}
TD_MAX_PER_KWARG = {
    "nanoseconds": Timedelta.max.value,
    "microseconds": Timedelta.max.value // 1_000,
    "milliseconds": Timedelta.max.value // 1_000_000,
    "seconds": Timedelta.max.value // 1_000_000_000,
    "minutes": Timedelta.max.value // (1_000_000_000 * 60),
    "hours": Timedelta.max.value // (1_000_000_000 * 60 * 60),
    "days": Timedelta.max.value // (1_000_000_000 * 60 * 60 * 24),
    "weeks": Timedelta.max.value // (1_000_000_000 * 60 * 60 * 24 * 7),
}
TD_MIN_PER_KWARG = {
    "nanoseconds": Timedelta.min.value,
    "microseconds": Timedelta.min.value // 1_000,
    "milliseconds": Timedelta.min.value // 1_000_000,
    "seconds": Timedelta.min.value // 1_000_000_000,
    "minutes": Timedelta.min.value // (1_000_000_000 * 60),
    "hours": Timedelta.min.value // (1_000_000_000 * 60 * 60),
    "days": Timedelta.min.value // (1_000_000_000 * 60 * 60 * 24),
    "weeks": Timedelta.min.value // (1_000_000_000 * 60 * 60 * 24 * 7),
}
TD_MAX_PER_UNIT = dict(
    chain.from_iterable(
        zip_longest(units, (TD_MAX_PER_KWARG[k],), fillvalue=TD_MAX_PER_KWARG[k])
        for k, units in TD_KWARGS_UNITS.items()
    )
)
TD_MIN_PER_UNIT = dict(
    chain.from_iterable(
        zip_longest(units, (TD_MIN_PER_KWARG[k],), fillvalue=TD_MIN_PER_KWARG[k])
        for k, units in TD_KWARGS_UNITS.items()
    )
)
TD_KWARGS_NP_TD64_UNITS = dict(
    zip(TD_MAX_PER_KWARG, ("ns", "us", "ms", "s", "m", "h", "D", "W"))
)
NP_TD64_MAX_PER_UNIT = dict(
    zip(("ns", "us", "ms", "s", "m", "h", "D", "W"), TD_MAX_PER_KWARG.values())
)
NP_TD64_MIN_PER_UNIT = dict(
    zip(("ns", "us", "ms", "s", "m", "h", "D", "W"), TD_MIN_PER_KWARG.values())
)


skip_ns = lambda d: {k: v for k, v in d.items() if not k.startswith("n")}


def test_construction():
    expected = np.timedelta64(10, "D").astype("m8[ns]").view("i8")
    assert Timedelta(10, unit="d").value == expected
    assert Timedelta(10.0, unit="d").value == expected
    assert Timedelta("10 days").value == expected
    assert Timedelta(days=10).value == expected
    assert Timedelta(days=10.0).value == expected

    expected += np.timedelta64(10, "s").astype("m8[ns]").view("i8")
    assert Timedelta("10 days 00:00:10").value == expected
    assert Timedelta(days=10, seconds=10).value == expected
    assert Timedelta(days=10, milliseconds=10 * 1000).value == expected
    assert Timedelta(days=10, microseconds=10 * 1000 * 1000).value == expected

    # rounding cases
    assert Timedelta(82739999850000).value == 82739999850000
    assert "0 days 22:58:59.999850" in str(Timedelta(82739999850000))
    assert Timedelta(123072001000000).value == 123072001000000
    assert "1 days 10:11:12.001" in str(Timedelta(123072001000000))

    # string conversion with/without leading zero
    # GH#9570
    assert Timedelta("0:00:00") == timedelta(hours=0)
    assert Timedelta("00:00:00") == timedelta(hours=0)
    assert Timedelta("-1:00:00") == -timedelta(hours=1)
    assert Timedelta("-01:00:00") == -timedelta(hours=1)

    # more strings & abbrevs
    # GH#8190
    assert Timedelta("1 h") == timedelta(hours=1)
    assert Timedelta("1 hour") == timedelta(hours=1)
    assert Timedelta("1 hr") == timedelta(hours=1)
    assert Timedelta("1 hours") == timedelta(hours=1)
    assert Timedelta("-1 hours") == -timedelta(hours=1)
    assert Timedelta("1 m") == timedelta(minutes=1)
    assert Timedelta("1.5 m") == timedelta(seconds=90)
    assert Timedelta("1 minute") == timedelta(minutes=1)
    assert Timedelta("1 minutes") == timedelta(minutes=1)
    assert Timedelta("1 s") == timedelta(seconds=1)
    assert Timedelta("1 second") == timedelta(seconds=1)
    assert Timedelta("1 seconds") == timedelta(seconds=1)
    assert Timedelta("1 ms") == timedelta(milliseconds=1)
    assert Timedelta("1 milli") == timedelta(milliseconds=1)
    assert Timedelta("1 millisecond") == timedelta(milliseconds=1)
    assert Timedelta("1 us") == timedelta(microseconds=1)
    assert Timedelta("1 µs") == timedelta(microseconds=1)
    assert Timedelta("1 micros") == timedelta(microseconds=1)
    assert Timedelta("1 microsecond") == timedelta(microseconds=1)
    assert Timedelta("1.5 microsecond") == Timedelta("00:00:00.000001500")
    assert Timedelta("1 ns") == Timedelta("00:00:00.000000001")
    assert Timedelta("1 nano") == Timedelta("00:00:00.000000001")
    assert Timedelta("1 nanosecond") == Timedelta("00:00:00.000000001")

    # combos
    assert Timedelta("10 days 1 hour") == timedelta(days=10, hours=1)
    assert Timedelta("10 days 1 h") == timedelta(days=10, hours=1)
    assert Timedelta("10 days 1 h 1m 1s") == timedelta(
        days=10, hours=1, minutes=1, seconds=1
    )
    assert Timedelta("-10 days 1 h 1m 1s") == -timedelta(
        days=10, hours=1, minutes=1, seconds=1
    )
    assert Timedelta("-10 days 1 h 1m 1s") == -timedelta(
        days=10, hours=1, minutes=1, seconds=1
    )
    assert Timedelta("-10 days 1 h 1m 1s 3us") == -timedelta(
        days=10, hours=1, minutes=1, seconds=1, microseconds=3
    )
    assert Timedelta("-10 days 1 h 1.5m 1s 3us") == -timedelta(
        days=10, hours=1, minutes=1, seconds=31, microseconds=3
    )

    # floats
    expected = np.timedelta64(10, "s").astype("m8[ns]").view("i8") + np.timedelta64(
        500, "ms"
    ).astype("m8[ns]").view("i8")
    assert Timedelta(10.5, unit="s").value == expected

    # offset
    assert to_timedelta(offsets.Hour(2)) == Timedelta(hours=2)
    assert Timedelta(offsets.Hour(2)) == Timedelta(hours=2)
    assert Timedelta(offsets.Second(2)) == Timedelta(seconds=2)

    # GH#11995: unicode
    expected = Timedelta("1H")
    result = Timedelta("1H")
    assert result == expected
    assert to_timedelta(offsets.Hour(2)) == Timedelta("0 days, 02:00:00")


@pytest.mark.parametrize("unit", ("ps", "ns"))
def test_from_np_td64_ignores_unit(unit: str):
    """
    Ignore the unit, as it may cause silently overflows leading to incorrect results,
    and in non-overflow cases is irrelevant GH#46827.
    """
    td64 = np.timedelta64(NP_TD64_MAX_PER_UNIT["h"], "h")
    msg = re.escape(
        "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
    )

    assert Timedelta(td64, unit=unit) == Timedelta(td64)

    with pytest.raises(OutOfBoundsTimedelta, match=msg):
        Timedelta(td64 * 2, unit=unit)


@pytest.mark.parametrize(("td_kwarg", "np_unit"), TD_KWARGS_NP_TD64_UNITS.items())
@pytest.mark.parametrize(
    "np_dtype",
    (np.int64, np.int32, np.int16, np.float64, np.float32, np.float16),
)
def test_td_construction_with_np_dtypes(np_dtype: type, td_kwarg: str, np_unit: str):
    # GH#8757: test construction with np dtypes
    expected_ns = np.timedelta64(1, np_unit).astype("m8[ns]").view("i8")
    assert Timedelta(**{td_kwarg: np_dtype(1)}).value == expected_ns


@pytest.mark.parametrize(
    "val",
    [
        "1s",
        "-1s",
        "1us",
        "-1us",
        "1 day",
        "-1 day",
        "-23:59:59.999999",
        "-1 days +23:59:59.999999",
        "-1ns",
        "1ns",
        "-23:59:59.999999999",
    ],
)
def test_td_from_repr_roundtrip(val):
    # round-trip both for string and value
    td = Timedelta(val)
    assert Timedelta(td.value) == td

    assert Timedelta(str(td)) == td
    assert Timedelta(td._repr_base(format="all")) == td
    assert Timedelta(td._repr_base()) == td


@pytest.mark.parametrize(
    "fmt,exp",
    [
        (
            "P6DT0H50M3.010010012S",
            Timedelta(
                days=6,
                minutes=50,
                seconds=3,
                milliseconds=10,
                microseconds=10,
                nanoseconds=12,
            ),
        ),
        (
            "P-6DT0H50M3.010010012S",
            Timedelta(
                days=-6,
                minutes=50,
                seconds=3,
                milliseconds=10,
                microseconds=10,
                nanoseconds=12,
            ),
        ),
        ("P4DT12H30M5S", Timedelta(days=4, hours=12, minutes=30, seconds=5)),
        ("P0DT0H0M0.000000123S", Timedelta(nanoseconds=123)),
        ("P0DT0H0M0.00001S", Timedelta(microseconds=10)),
        ("P0DT0H0M0.001S", Timedelta(milliseconds=1)),
        ("P0DT0H1M0S", Timedelta(minutes=1)),
        ("P1DT25H61M61S", Timedelta(days=1, hours=25, minutes=61, seconds=61)),
        ("PT1S", Timedelta(seconds=1)),
        ("PT0S", Timedelta(seconds=0)),
        ("P1WT0S", Timedelta(days=7, seconds=0)),
        ("P1D", Timedelta(days=1)),
        ("P1DT1H", Timedelta(days=1, hours=1)),
        ("P1W", Timedelta(days=7)),
        ("PT300S", Timedelta(seconds=300)),
        ("P1DT0H0M00000000000S", Timedelta(days=1)),
        ("PT-6H3M", Timedelta(hours=-6, minutes=3)),
        ("-PT6H3M", Timedelta(hours=-6, minutes=-3)),
        ("-PT-6H+3M", Timedelta(hours=6, minutes=-3)),
    ],
)
def test_iso_constructor(fmt, exp):
    assert Timedelta(fmt) == exp


@pytest.mark.parametrize(
    "constructed_td, conversion",
    [
        (Timedelta(nanoseconds=100), "100ns"),
        (
            Timedelta(
                days=1,
                hours=1,
                minutes=1,
                weeks=1,
                seconds=1,
                milliseconds=1,
                microseconds=1,
                nanoseconds=1,
            ),
            694861001001001,
        ),
        (Timedelta(microseconds=1) + Timedelta(nanoseconds=1), "1us1ns"),
        (Timedelta(microseconds=1) - Timedelta(nanoseconds=1), "999ns"),
        (Timedelta(microseconds=1) + 5 * Timedelta(nanoseconds=-2), "990ns"),
    ],
)
def test_td_constructor_on_nanoseconds(constructed_td, conversion):
    # GH#9273
    assert constructed_td == Timedelta(conversion)


@pytest.mark.parametrize(
    ("args", "kwargs"),
    [
        ((), {}),
        (("ps",), {}),
        (("ns",), {}),
        (("ms",), {}),
        ((), {"seconds": 3}),
        (("ns",), {"minutes": 2}),
    ],
)
def test_other_args_ignored_if_timedelta_value_passed(args: tuple, kwargs: dict):
    original = Timedelta(1)
    new = Timedelta(original, *args, **kwargs)

    assert new == original
    if not any((args, kwargs)):
        assert new is original


@pytest.mark.parametrize(
    "value",
    (
        None,
        np.nan,
        NaT,
        pytest.param(
            NA,
            marks=pytest.mark.xfail(
                reason="constructor fails",
                raises=ValueError,
                strict=True,
            ),
        ),
    ),
    ids=("None", "np.nan", "pd.NaT", "pd.NA"),
)
def test_returns_nat_for_most_na_values(value):
    assert Timedelta(value) is NaT


class TestInvalidArgCombosFormats:
    def test_raises_if_no_args_passed(self):
        msg = re.escape(
            "cannot construct a Timedelta without a value/unit or descriptive keywords "
            "(days,seconds....)"
        )

        with pytest.raises(ValueError, match=msg):
            Timedelta()

    @pytest.mark.parametrize("unit", ("years", "months", "day", "ps"))
    def test_raises_for_invalid_kwarg(self, unit: str):
        msg = re.escape(
            "cannot construct a Timedelta from the passed arguments, allowed keywords "
            "are ('weeks', 'days', 'hours', 'minutes', 'seconds', 'milliseconds', "
            "'microseconds', 'nanoseconds')"
        )

        with pytest.raises(ValueError, match=msg):
            Timedelta(**{unit: 1})  # type: ignore[arg-type]

    def test_raises_if_kwarg_has_str_value(self):
        msg = "Invalid type <class 'str'>. Must be int or float."

        with pytest.raises(TypeError, match=msg):
            Timedelta(nanoseconds="1")

    @pytest.mark.parametrize(
        ("constructor", "value", "unit", "msg"),
        (
            (Timedelta, "10s", "ms", "the value is a str"),
            (to_timedelta, "10s", "ms", "the input is/contains a str"),
            (to_timedelta, ["1", "2", "3"], "s", "the input contains a str"),
        ),
        ids=("Timedelta", "to_timedelta-scalar", "to_timedelta-sequence"),
    )
    def test_raises_if_both_str_value_and_unit_passed(
        self,
        constructor,
        value,
        unit,
        msg,
    ):
        msg = "unit must not be specified if " + msg

        with pytest.raises(ValueError, match=msg):
            constructor(value, unit=unit)

    @pytest.mark.parametrize(
        "value",
        [
            "PPPPPPPPPPPP",
            "PDTHMS",
            "P0DT999H999M999S",
            "P1DT0H0M0.0000000000000S",
            "P1DT0H0M0.S",
            "P",
            "-P",
        ],
    )
    def test_raises_for_invalid_iso_like_str_value(self, value):
        msg = f"Invalid ISO 8601 Duration format - {value}"

        with pytest.raises(ValueError, match=msg):
            Timedelta(value)

    def test_raises_if_str_value_contains_no_units(self):
        msg = "no units specified"

        with pytest.raises(ValueError, match=msg):
            Timedelta("3.1415")

    @pytest.mark.parametrize(
        ("value", "msg"),
        (
            ("us", "unit abbreviation w/o a number"),
            ("seconds", "unit abbreviation w/o a number"),
            ("garbage", "unit abbreviation w/o a number"),
            # GH39710 Timedelta input string with only symbols and no digits raises
            ("+", "symbols w/o a number"),
            ("-", "symbols w/o a number"),
        ),
    )
    def test_raises_if_str_value_contains_no_numeric_component(
        self,
        value: str,
        msg: str,
    ):
        with pytest.raises(ValueError, match=msg):
            Timedelta(value)

    @pytest.mark.parametrize(
        "value",
        (
            "--",
            # Currently invalid as it has a - on the hh:mm:dd part
            # (only allowed on the days)
            "-10 days -1 h 1.5m 1s 3us",
            "10 days -1 h 1.5m 1s 3us",
        ),
    )
    def test_raises_for_str_value_with_minus_sign(self, value: str):
        msg = "only leading negative signs are allowed"
        with pytest.raises(ValueError, match=msg):
            Timedelta(value)

    @pytest.mark.parametrize("unit", ["Y", "y", "M"])
    def test_raises_if_ambiguous_units_passed(self, unit: str):
        msg = (
            "Units 'M', 'Y', and 'y' are no longer supported, as they do not "
            "represent unambiguous timedelta values durations."
        )

        with pytest.raises(ValueError, match=msg):
            Timedelta(1, unit)


class TestOverflow:
    @pytest.mark.parametrize(("unit", "max_val"), TD_MAX_PER_UNIT.items())
    def test_int_plus_units_too_big(self, unit: str, max_val: int, request):
        if unit == "w":
            mark = pytest.mark.xfail(
                reason="does not raise",
                raises=pytest.fail.Exception,
                strict=True,
            )
            request.node.add_marker(mark)

        too_big = max_val + 1
        msg = re.escape(
            "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
        )

        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            Timedelta(too_big, unit=unit)

    @pytest.mark.parametrize(("unit", "min_val"), skip_ns(TD_MIN_PER_UNIT).items())
    def test_int_plus_units_too_small(self, unit: str, min_val: int, request):
        if unit == "w":
            mark = pytest.mark.xfail(
                reason="does not raise",
                raises=pytest.fail.Exception,
                strict=True,
            )
            request.node.add_marker(mark)

        too_small = min_val - 1
        msg = re.escape(
            "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
        )

        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            Timedelta(too_small, unit=unit)

    @pytest.mark.parametrize(("kwarg", "max_val"), TD_MAX_PER_KWARG.items())
    def test_kwarg_too_big(self, kwarg: str, max_val: int):
        too_big = max_val + 1
        msg = re.escape(
            "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
        )

        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            assert Timedelta(**{kwarg: too_big})  # type: ignore[arg-type]

    @pytest.mark.parametrize(("kwarg", "min_val"), skip_ns(TD_MIN_PER_KWARG).items())
    def test_kwarg_too_small(self, kwarg: str, min_val: int):
        too_small = min_val - 1
        msg = re.escape(
            "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
        )

        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            Timedelta(**{kwarg: too_small})  # type: ignore[arg-type]

    @pytest.mark.parametrize(("kwarg", "max_val"), skip_ns(TD_MAX_PER_KWARG).items())
    def test_from_timedelta_too_big(self, kwarg: str, max_val: int):
        too_big = timedelta(**{kwarg: max_val + 1})
        msg = re.escape(
            "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
        )

        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            Timedelta(too_big)

    @pytest.mark.parametrize(("kwarg", "min_val"), skip_ns(TD_MIN_PER_KWARG).items())
    def test_from_timedelta_too_small(self, kwarg: str, min_val: int):
        too_small = timedelta(**{kwarg: min_val - 1})
        msg = re.escape(
            "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
        )

        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            Timedelta(too_small)

    @pytest.mark.parametrize(("unit", "max_val"), skip_ns(NP_TD64_MAX_PER_UNIT).items())
    def test_from_np_td64_too_big(self, unit: str, max_val: int):
        too_big = np.timedelta64(max_val + 1, unit)
        msg = re.escape(
            "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
        )

        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            Timedelta(too_big)

    @pytest.mark.parametrize(("unit", "min_val"), skip_ns(NP_TD64_MIN_PER_UNIT).items())
    def test_from_np_td64_too_small(self, unit: str, min_val: int):
        too_small = np.timedelta64(min_val - 1, unit)
        msg = re.escape(
            "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
        )

        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            Timedelta(too_small)

    def test_too_small_by_1ns_returns_nat(self):
        too_small = Timedelta.min.value - 1
        too_small_np_td = np.timedelta64(too_small)

        assert isinstance(too_small, int)
        assert isinstance(too_small_np_td, np.timedelta64)

        assert Timedelta(too_small, "ns") is NaT
        assert Timedelta(nanoseconds=too_small) is NaT
        assert Timedelta(too_small_np_td) is NaT
