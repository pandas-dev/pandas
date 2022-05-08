"""
Most Timedelta scalar tests; See test_arithmetic for tests of binary operations with a
Timedelta scalar.
"""

from __future__ import annotations

from datetime import timedelta
from itertools import (
    chain,
    product,
    zip_longest,
)
import operator
import re

from hypothesis import (
    given,
    strategies as st,
)
import numpy as np
import pytest

from pandas._libs import lib
from pandas._libs.tslibs import (
    OutOfBoundsTimedelta,
    iNaT,
)

import pandas as pd
from pandas import (
    NA,
    NaT,
    Timedelta,
    TimedeltaIndex,
    offsets,
    to_timedelta,
)
import pandas._testing as tm

TD_UNITS = (
    ("n", "ns", "nano", "nanos", "nanosecond", "nanoseconds"),
    ("u", "us", "µs", "micro", "micros", "microsecond", "microseconds"),
    ("l", "ms", "milli", "millis", "millisecond", "milliseconds"),
    ("s", "sec", "second", "seconds"),
    ("m", "t", "min", "minute", "minutes"),
    ("h", "hr", "hour", "hours"),
    ("d", "day", "days"),
    ("w",),
)
TD_UNITS_UNIQUE = tuple(map(operator.itemgetter(0), TD_UNITS))
TD_KWARGS = (
    "nanoseconds",
    "microseconds",
    "milliseconds",
    "seconds",
    "minutes",
    "hours",
    "days",
    "weeks",
)
TD_COMPONENTS = tuple(reversed(TD_KWARGS[:-1]))
TD64_UNITS = ("ns", "us", "ms", "s", "m", "h", "D", "W")

TD_KWARGS_TD_UNITS = dict(zip(TD_KWARGS, TD_UNITS))
TD_UNITS_TD_KWARGS = dict(
    chain.from_iterable(
        zip_longest(units, (kwarg,), fillvalue=kwarg)
        for kwarg, units in TD_KWARGS_TD_UNITS.items()
    )
)
TD_KWARGS_TD64_UNITS = dict(zip(TD_KWARGS, TD64_UNITS))
TD_UNITS_TD64_UNITS = dict(
    chain.from_iterable(
        zip_longest(td_units, (TD64_UNITS[i],), fillvalue=TD64_UNITS[i])
        for i, td_units in enumerate(TD_UNITS)
    )
)

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
    "microseconds": Timedelta.min.value // 1_000 + 1,
    "milliseconds": Timedelta.min.value // 1_000_000 + 1,
    "seconds": Timedelta.min.value // 1_000_000_000 + 1,
    "minutes": Timedelta.min.value // (1_000_000_000 * 60) + 1,
    "hours": Timedelta.min.value // (1_000_000_000 * 60 * 60) + 1,
    "days": Timedelta.min.value // (1_000_000_000 * 60 * 60 * 24) + 1,
    "weeks": Timedelta.min.value // (1_000_000_000 * 60 * 60 * 24 * 7) + 1,
}
# simplified to include only one key corresponding to each unit
TD_MAX_PER_UNIT = dict(zip(TD_UNITS_UNIQUE, TD_MAX_PER_KWARG.values()))
TD_MIN_PER_UNIT = dict(zip(TD_UNITS_UNIQUE, TD_MIN_PER_KWARG.values()))
TD64_MAX_PER_UNIT = dict(zip(TD64_UNITS, TD_MAX_PER_KWARG.values()))
TD64_MIN_PER_UNIT = dict(zip(TD64_UNITS, TD_MIN_PER_KWARG.values()))

xfail_does_not_raise = pytest.mark.xfail(
    reason="does not raise",
    raises=pytest.fail.Exception,
    strict=True,
)
skip_ns = lambda s: (u for u in s if not u.startswith("n"))


@pytest.fixture(name="timedelta_kwarg", params=skip_ns(TD_KWARGS))
def fixture_timedelta_kwarg(request) -> str:
    return request.param


@pytest.fixture(name="td_max_per_unit", params=TD_MAX_PER_UNIT)
def fixture_td_max_per_unit(request) -> tuple:
    unit = request.param
    if request.cls is TestOverflow and unit == "w":
        request.applymarker(xfail_does_not_raise)

    return unit, TD_MAX_PER_UNIT[unit]


@pytest.fixture(name="td_min_per_unit", params=TD_MIN_PER_UNIT)
def fixture_td_min_per_unit(request) -> tuple:
    unit = request.param
    if request.cls is TestOverflow and unit == "w":
        request.applymarker(xfail_does_not_raise)

    return unit, TD_MIN_PER_UNIT[unit]


@pytest.fixture(name="td_max_per_kwarg", params=TD_MAX_PER_KWARG)
def fixture_td_max_per_kwarg(request) -> tuple:
    kwarg = request.param
    return kwarg, TD_MAX_PER_KWARG[kwarg]


@pytest.fixture(name="td_min_per_kwarg", params=TD_MIN_PER_KWARG)
def fixture_td_min_per_kwarg(request) -> tuple:
    kwarg = request.param
    return kwarg, TD_MIN_PER_KWARG[kwarg]


@pytest.fixture(name="td64_max_per_unit", params=skip_ns(TD64_MAX_PER_UNIT))
def fixture_td64_max_per_unit(request) -> tuple:
    unit = request.param
    return unit, TD64_MAX_PER_UNIT[unit]


@pytest.fixture(name="td64_min_per_unit", params=skip_ns(TD64_MIN_PER_UNIT))
def fixture_td64_min_per_unit(request) -> tuple:
    unit = request.param
    return unit, TD64_MIN_PER_UNIT[unit]


@pytest.fixture(name="td_overflow_msg")
def fixture_td_overflow_msg() -> str:
    return re.escape(
        "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
    )


@pytest.fixture(name="non_nano_reso", params=(7, 8, 9))
def fixture_non_nano_reso(request):
    """7, 8, 9 correspond to second, millisecond, and microsecond, respectively"""
    return request.param


@pytest.fixture(name="non_nano_td")
def fixture_non_nano_td(non_nano_reso: int) -> Timedelta:
    # microsecond that would be just out of bounds for nano
    us = np.int64((TD_MAX_PER_KWARG["days"] + 1) * 86_400 * 1_000_000)
    values = {
        9: us,
        8: us // 1000,
        7: us // 1_000_000,
    }

    return Timedelta._from_value_and_reso(values[non_nano_reso], non_nano_reso)


class TestConstruction:
    """
    Tests of the public constructor, Timedelta.__new__().
    """

    def test_type(self):
        td = Timedelta(1)

        assert isinstance(td, Timedelta)
        assert isinstance(td, timedelta)

    @pytest.mark.parametrize("td_unit, td64_unit", TD_UNITS_TD64_UNITS.items())
    def test_from_value_and_unit(
        self,
        td_unit: str,
        td64_unit: str,
        any_real_numpy_dtype: str,
    ):
        """GH#8757: test construction with np dtypes"""
        expected_ns = np.timedelta64(1, td64_unit).astype("m8[ns]").view("i8")
        one = np.dtype(any_real_numpy_dtype).type(1)
        td = Timedelta(one, td_unit)

        assert td.value == expected_ns

    @pytest.mark.parametrize("subset", map(slice, range(1, len(TD_UNITS_UNIQUE))))
    def test_from_str(self, subset: slice):
        """GH#8190"""
        td64s = tuple(np.timedelta64(1, u) for u in TD64_UNITS[subset])
        str_value = " ".join(tuple(f"1 {u}" for u in TD_UNITS_UNIQUE[subset]))
        expected_ns = np.sum(td64s).astype("m8[ns]").view("i8")
        td = Timedelta(str_value)
        neg_td = Timedelta("-" + str_value)

        assert td.value == expected_ns
        assert neg_td.value == -1 * expected_ns

    @pytest.mark.parametrize(
        "value, expected_hours",
        (
            ("0:00:00", 0),
            ("1:00:00", 1),
        ),
    )
    def test_from_str_with_without_leading_zero(self, value: str, expected_hours: int):
        """GH#9570"""
        expected_ns = np.timedelta64(expected_hours, "h").astype("m8[ns]").view("i8")
        td0 = Timedelta(value)
        td1 = Timedelta("0" + value)

        assert td0.value == expected_ns
        assert td1.value == expected_ns

    @pytest.mark.parametrize(
        ("value", "expected"),
        (
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
        ),
    )
    def test_from_isoformat_str(self, value: str, expected: Timedelta):
        assert Timedelta(value) == expected

    @pytest.mark.parametrize("subset", map(slice, range(1, len(TD_KWARGS))))
    def test_from_kwargs(self, subset: slice, any_real_numpy_dtype: str):
        td64s = tuple(np.timedelta64(1, u) for u in TD64_UNITS[subset])
        kwargs = {u: np.dtype(any_real_numpy_dtype).type(1) for u in TD_KWARGS[subset]}
        expected_ns = np.sum(td64s).astype("m8[ns]").view("i8")
        td = Timedelta(**kwargs)

        assert td.value == expected_ns

    @pytest.mark.parametrize("td_unit, td_kwarg", TD_UNITS_TD_KWARGS.items())
    def test_kwarg_unit_equivalence(self, request, td_unit: str, td_kwarg: str):
        if td_kwarg == "weeks":
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="this one isn't valid",
                    raises=ValueError,
                    strict=True,
                )
            )

        from_unit = Timedelta(1, td_unit)
        from_kwarg = Timedelta(**{td_kwarg: 1})  # type: ignore[arg-type]
        from_str_unit = Timedelta(f"1 {td_unit}")
        from_str_kwarg = Timedelta(f"1 {td_kwarg}")

        assert from_unit == from_kwarg == from_str_unit == from_str_kwarg

    @pytest.mark.parametrize(
        "value, td_unit, expected_ns",
        (
            (9.123, "us", 9123),
            (9.123456, "ms", 9123456),
            (9.123456789, "s", 9123456789),
        ),
    )
    def test_float_values_not_rounded(
        self,
        value: float,
        td_unit: str,
        expected_ns: int,
    ):
        """GH#12690"""
        td_kwarg = TD_UNITS_TD_KWARGS[td_unit]
        from_float = Timedelta(value, td_unit)
        from_str = Timedelta(f"{value} {td_unit}")
        from_kwarg = Timedelta(**{td_kwarg: value})  # type: ignore[arg-type]

        assert from_float.value == expected_ns
        assert from_str.value == expected_ns
        assert from_kwarg.value == expected_ns

    def test_from_offset(self, tick_classes):
        offset = tick_classes(1)
        assert Timedelta(offset).value == offset.nanos

    @pytest.mark.parametrize("td_unit", TD_UNITS)
    def test_from_td64_ignores_unit(self, td_unit: str, td_overflow_msg: str):
        """
        Ignore the unit, as it may cause silently overflows leading to incorrect
        results, and in non-overflow cases is irrelevant GH#46827.
        """
        td64 = np.timedelta64(TD64_MAX_PER_UNIT["h"], "h")

        assert Timedelta(td64, td_unit) == Timedelta(td64)
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta(td64 * 2, td_unit)

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
    def test_from_td_ignores_other_args(self, args: tuple, kwargs: dict):
        original = Timedelta(1)
        new = Timedelta(original, *args, **kwargs)

        assert new == original
        if not any((args, kwargs)):
            assert new is original

    def test_from_timedelta(self, timedelta_kwarg: str):
        kwargs = {timedelta_kwarg: 1}
        assert Timedelta(**kwargs) == timedelta(**kwargs)  # type: ignore[arg-type]

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
    def test_from_na_value_returns_nat(self, value):
        assert Timedelta(value) is NaT

    def test_raises_if_no_args_passed(self):
        msg = (
            "cannot construct a Timedelta without a value/unit or descriptive keywords"
        )

        with pytest.raises(ValueError, match=msg):
            Timedelta()

    @pytest.mark.parametrize("unit", ("years", "months", "day", "ps", "reso", "_reso"))
    def test_raises_for_invalid_kwarg(self, unit: str):
        msg = "cannot construct a Timedelta from the passed arguments"

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
    def test_raises_for_invalid_isolike_str_value(self, value):
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
    def test_raises_if_str_value_has_no_numeric_component(self, value: str, msg: str):
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
    def test_raises_for_str_value_with_second_minus_sign(self, value: str):
        msg = "only leading negative signs are allowed"
        with pytest.raises(ValueError, match=msg):
            Timedelta(value)

    @pytest.mark.parametrize(
        ("unit", "func"),
        product(("Y", "y", "M"), (Timedelta, to_timedelta)),
    )
    def test_warns_or_raises_if_ambiguous_unit_passed(self, unit: str, func):
        msg = "Units 'M', 'Y', and 'y' are no longer supported"

        with pytest.raises(ValueError, match=msg):
            func(1, unit)

    def test_reso_invariant_if_td_created_via_public_api(self, td_max_per_unit: tuple):
        unit, max_value = td_max_per_unit
        td_small = Timedelta(1, unit)
        td_max = Timedelta(max_value, unit)
        msg = "attribute '_reso' of 'pandas._libs.tslibs.timedeltas._Timedelta'"

        assert getattr(td_small, "_reso") == 10
        assert getattr(td_max, "_reso") == 10
        with pytest.raises(AttributeError, match=msg):
            setattr(td_max, "_reso", 9)

    def test_reso_configurable_via_private_api(self, non_nano_reso: int):
        td = Timedelta._from_value_and_reso(np.int64(1), non_nano_reso)
        assert td.value == 1
        assert getattr(td, "_reso") == non_nano_reso


class TestOverflow:
    def test_value_unit_too_big(self, td_max_per_unit: tuple, td_overflow_msg: str):
        unit, value = td_max_per_unit

        assert Timedelta(value, unit) <= Timedelta.max
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta(value + 1, unit)

    def test_value_unit_too_small(self, td_min_per_unit: tuple, td_overflow_msg: str):
        unit, value = td_min_per_unit
        too_small = value - 1

        assert Timedelta(value, unit) >= Timedelta.min
        if unit == "n":
            result = Timedelta(too_small, unit)
            assert result is NaT  # type: ignore[comparison-overlap]
            too_small -= 1
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta(too_small, unit)

    def test_kwarg_too_big(self, td_max_per_kwarg: tuple, td_overflow_msg: str):
        kwarg, value = td_max_per_kwarg

        assert Timedelta(**{kwarg: value}) <= Timedelta.max
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            assert Timedelta(**{kwarg: value + 1})

    def test_kwarg_too_small(self, td_min_per_kwarg: tuple, td_overflow_msg: str):
        kwarg, value = td_min_per_kwarg
        too_small = value - 1

        assert Timedelta(**{kwarg: value}) >= Timedelta.min
        if kwarg == "nanoseconds":
            result = Timedelta(**{kwarg: too_small})
            assert result is NaT  # type: ignore[comparison-overlap]
            too_small -= 1
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta(**{kwarg: too_small})

    def test_from_timedelta_too_big(self, timedelta_kwarg: str, td_overflow_msg: str):
        max_val = TD_MAX_PER_KWARG[timedelta_kwarg]

        assert Timedelta(timedelta(**{timedelta_kwarg: max_val})) <= Timedelta.max
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta(timedelta(**{timedelta_kwarg: max_val + 1}))

    def test_from_timedelta_too_small(self, timedelta_kwarg: str, td_overflow_msg: str):
        min_val = TD_MIN_PER_KWARG[timedelta_kwarg]

        assert Timedelta(timedelta(**{timedelta_kwarg: min_val})) >= Timedelta.min
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta(timedelta(**{timedelta_kwarg: min_val - 1}))

    def test_from_td64_too_big(self, td64_max_per_unit: tuple, td_overflow_msg: str):
        unit, value = td64_max_per_unit

        assert Timedelta(np.timedelta64(value, unit)) <= Timedelta.max
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta(np.timedelta64(value + 1, unit))

    def test_from_td64_too_small(self, td64_min_per_unit: tuple, td_overflow_msg: str):
        unit, value = td64_min_per_unit

        assert Timedelta(np.timedelta64(value, unit)) >= Timedelta.min
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta(np.timedelta64(value - 1, unit))


class TestNonNano:
    """
    WIP.
    """

    def test_unary_non_nano(self, non_nano_td, non_nano_reso):
        assert abs(non_nano_td)._reso == non_nano_reso
        assert (-non_nano_td)._reso == non_nano_reso
        assert (+non_nano_td)._reso == non_nano_reso

    def test_sub_preserves_reso(self, non_nano_td, non_nano_reso):
        res = non_nano_td - non_nano_td
        expected = Timedelta._from_value_and_reso(0, non_nano_reso)
        assert res == expected
        assert res._reso == non_nano_reso

    def test_mul_preserves_reso(self, non_nano_td, non_nano_reso):
        # The non_nano_td fixture should always be far from the implementation
        #  bound, so doubling does not risk overflow.
        res = non_nano_td * 2
        assert res.value == non_nano_td.value * 2
        assert res._reso == non_nano_reso

    def test_cmp_cross_reso(self, non_nano_td):
        # numpy gets this wrong because of silent overflow
        assert Timedelta.max < non_nano_td
        assert non_nano_td > Timedelta.max
        assert not Timedelta.max == non_nano_td
        assert non_nano_td != Timedelta.max

    def test_to_pytimedelta(self, non_nano_td):
        res = non_nano_td.to_pytimedelta()
        expected = timedelta(days=106752)
        assert type(res) is timedelta
        assert res == expected

    @pytest.mark.parametrize(
        "converter",
        (
            operator.methodcaller("to_timedelta64"),
            operator.methodcaller("to_numpy"),
            operator.attrgetter("asm8"),
        ),
    )
    def test_to_timedelta64(self, non_nano_td, converter):
        td64 = converter(non_nano_td)
        reso_dtype = {7: "m8[s]", 8: "m8[ms]", 9: "m8[us]"}

        assert isinstance(td64, np.timedelta64)
        assert td64.view("i8") == non_nano_td.value
        assert td64.dtype == reso_dtype[non_nano_td._reso]


class TestUnaryOps:
    def test_invert(self):
        td = Timedelta(10, unit="d")

        msg = "bad operand type for unary ~"
        with pytest.raises(TypeError, match=msg):
            ~td

        # check this matches pytimedelta and timedelta64
        with pytest.raises(TypeError, match=msg):
            ~(td.to_pytimedelta())

        umsg = "ufunc 'invert' not supported for the input types"
        with pytest.raises(TypeError, match=umsg):
            ~(td.to_timedelta64())

    def test_unary_ops(self):
        td = Timedelta(10, unit="d")

        # __neg__, __pos__
        assert -td == Timedelta(-10, unit="d")
        assert -td == Timedelta("-10d")
        assert +td == Timedelta(10, unit="d")

        # __abs__, __abs__(__neg__)
        assert abs(td) == td
        assert abs(-td) == td
        assert abs(-td) == Timedelta("10d")


class TestAttributes:
    def test_min_max_correspond_to_int64_boundaries(self):
        """GH#12727"""
        assert Timedelta.min.value == iNaT + 1
        assert Timedelta.max.value == lib.i8max

    def test_fields(self):
        """GH#10050: compat with datetime.timedelta; GH#31354"""
        fields = ("days", "seconds", "microseconds", "nanoseconds")
        td = Timedelta("1 days, 10:11:12")

        assert td.days == 1
        assert td.seconds == 10 * 3600 + 11 * 60 + 12
        assert td.microseconds == 0
        assert td.nanoseconds == 0
        assert all(isinstance(v, int) for v in operator.attrgetter(*fields)(td))
        assert td.days * 24 * 3600 * int(1e9) + td.seconds * int(1e9) == td.value

    @pytest.mark.parametrize("field", ("hours", "minutes", "milliseconds"))
    def test_fields_not_exposed(self, field: str):
        msg = f"'Timedelta' object has no attribute '{field}'"

        with pytest.raises(AttributeError, match=msg):
            getattr(Timedelta.max, field)

    @pytest.mark.parametrize(
        "td, expected_values",
        (
            (Timedelta("-1 us"), (-1, 23, 59, 59, 999, 999, 0)),
            (Timedelta("-1 days 1 us"), (-2, 23, 59, 59, 999, 999, 0)),
        ),
    )
    def test_components(self, td, expected_values: tuple[int]):
        values = operator.attrgetter(*TD_COMPONENTS)(td.components)

        assert values == expected_values
        assert all(isinstance(v, int) for v in values)

    def test_resolution_string(self):
        assert Timedelta(days=1).resolution_string == "D"
        assert Timedelta(hours=1).resolution_string == "H"
        assert Timedelta(minutes=1).resolution_string == "T"
        assert Timedelta(seconds=1).resolution_string == "S"
        assert Timedelta(milliseconds=1).resolution_string == "L"
        assert Timedelta(microseconds=1).resolution_string == "U"
        assert Timedelta(nanoseconds=1).resolution_string == "N"

    @pytest.mark.parametrize("td_units", TD_UNITS)
    def test_resolution_is_class_attr(self, td_units: str):
        """GH#21344; mirrors datetime.timedelta"""
        td = Timedelta(1, td_units[0])

        assert td.resolution is Timedelta.resolution
        assert Timedelta.resolution == Timedelta(1, "ns")

    def test_asm8_is_alias_for_to_timedelta64(self):
        result = Timedelta.max.asm8

        assert result == Timedelta.max.to_timedelta64()
        assert isinstance(result, np.timedelta64)

    @pytest.mark.parametrize(
        "attr, expected_value",
        (("delta", 1), ("freq", None), ("is_populated", False)),
    )
    def test_deprecated_attrs(self, attr: str, expected_value):
        """GH#46430, GH#46476"""
        td = Timedelta(1, "ns")
        msg = f"Timedelta.{attr}"
        with tm.assert_produces_warning(FutureWarning, match=msg):
            getattr(td, attr) == expected_value

        with pytest.raises(AttributeError, match="is not writable"):
            setattr(td, attr, "coconut")


class TestMethods:
    @pytest.mark.parametrize(
        "value, expected",
        (
            (
                "1 days, 10:11:12.123456789",
                1 * 86400 + 10 * 3600 + 11 * 60 + 12.123456,
            ),
            ("30S", 30.0),
            ("0", 0.0),
            ("-2S", -2.0),
            ("5.324S", 5.324),
        ),
    )
    def test_total_seconds(self, value: str, expected: float):
        # see gh-10939
        td = Timedelta(value)
        assert td.total_seconds() == expected

    def test_to_pytimedelta(self):
        td = Timedelta("1 days, 10:11:12.012345")
        py_td = td.to_pytimedelta()

        assert py_td == td
        assert Timedelta(py_td) == td
        assert isinstance(py_td, timedelta)
        assert not isinstance(py_td, Timedelta)

    @pytest.mark.parametrize(
        "td, expected",
        (
            (Timedelta(500, "ns"), timedelta(0)),
            (Timedelta(501, "ns"), timedelta(microseconds=1)),
        ),
    )
    def test_to_pytimedelta_rounds_ns(self, td: Timedelta, expected: timedelta):
        assert td.to_pytimedelta() == expected

    def test_to_timedelta64(self):
        td64 = Timedelta.max.to_timedelta64()

        assert td64 == Timedelta.max
        assert Timedelta(td64) == Timedelta.max
        assert isinstance(td64, np.timedelta64)

    def test_to_numpy(self):
        """GH#24653: alias .to_numpy() for scalars"""
        assert Timedelta.max.to_timedelta64() == Timedelta.max.to_numpy()

    @pytest.mark.parametrize(
        "args, kwargs",
        (
            (("m8[ns]",), {}),
            ((), {"copy": True}),
            (("m8[ns]",), {"copy": True}),
        ),
    )
    def test_to_numpy_raises_if_args_passed(self, args, kwargs):
        # GH#44460
        msg = "dtype and copy arguments are ignored"
        with pytest.raises(ValueError, match=msg):
            Timedelta.max.to_numpy(*args, **kwargs)

    @pytest.mark.parametrize(
        "freq,s1,s2",
        [
            # This first case has s1, s2 being the same as t1,t2 below
            (
                "N",
                Timedelta("1 days 02:34:56.789123456"),
                Timedelta("-1 days 02:34:56.789123456"),
            ),
            (
                "U",
                Timedelta("1 days 02:34:56.789123000"),
                Timedelta("-1 days 02:34:56.789123000"),
            ),
            (
                "L",
                Timedelta("1 days 02:34:56.789000000"),
                Timedelta("-1 days 02:34:56.789000000"),
            ),
            ("S", Timedelta("1 days 02:34:57"), Timedelta("-1 days 02:34:57")),
            ("2S", Timedelta("1 days 02:34:56"), Timedelta("-1 days 02:34:56")),
            ("5S", Timedelta("1 days 02:34:55"), Timedelta("-1 days 02:34:55")),
            ("T", Timedelta("1 days 02:35:00"), Timedelta("-1 days 02:35:00")),
            ("12T", Timedelta("1 days 02:36:00"), Timedelta("-1 days 02:36:00")),
            ("H", Timedelta("1 days 03:00:00"), Timedelta("-1 days 03:00:00")),
            ("d", Timedelta("1 days"), Timedelta("-1 days")),
        ],
    )
    def test_round(self, freq, s1, s2):

        t1 = Timedelta("1 days 02:34:56.789123456")
        t2 = Timedelta("-1 days 02:34:56.789123456")

        r1 = t1.round(freq)
        assert r1 == s1
        r2 = t2.round(freq)
        assert r2 == s2

    def test_round_invalid(self):
        t1 = Timedelta("1 days 02:34:56.789123456")

        for freq, msg in [
            ("Y", "<YearEnd: month=12> is a non-fixed frequency"),
            ("M", "<MonthEnd> is a non-fixed frequency"),
            ("foobar", "Invalid frequency: foobar"),
        ]:
            with pytest.raises(ValueError, match=msg):
                t1.round(freq)

    def test_round_implementation_bounds(self):
        # See also: analogous test for Timestamp
        # GH#38964
        result = Timedelta.min.ceil("s")
        expected = Timedelta.min + Timedelta(seconds=1) - Timedelta(145224193)
        assert result == expected

        result = Timedelta.max.floor("s")
        expected = Timedelta.max - Timedelta(854775807)
        assert result == expected

        with pytest.raises(OverflowError, match="value too large"):
            Timedelta.min.floor("s")

        # the second message here shows up in windows builds
        msg = "|".join(
            ["Python int too large to convert to C long", "int too big to convert"]
        )
        with pytest.raises(OverflowError, match=msg):
            Timedelta.max.ceil("s")

    @given(val=st.integers(min_value=iNaT + 1, max_value=lib.i8max))
    @pytest.mark.parametrize(
        "method", [Timedelta.round, Timedelta.floor, Timedelta.ceil]
    )
    def test_round_sanity(self, val, method):
        val = np.int64(val)
        td = Timedelta(val)

        assert method(td, "ns") == td

        res = method(td, "us")
        nanos = 1000
        assert np.abs((res - td).value) < nanos
        assert res.value % nanos == 0

        res = method(td, "ms")
        nanos = 1_000_000
        assert np.abs((res - td).value) < nanos
        assert res.value % nanos == 0

        res = method(td, "s")
        nanos = 1_000_000_000
        assert np.abs((res - td).value) < nanos
        assert res.value % nanos == 0

        res = method(td, "min")
        nanos = 60 * 1_000_000_000
        assert np.abs((res - td).value) < nanos
        assert res.value % nanos == 0

        res = method(td, "h")
        nanos = 60 * 60 * 1_000_000_000
        assert np.abs((res - td).value) < nanos
        assert res.value % nanos == 0

        res = method(td, "D")
        nanos = 24 * 60 * 60 * 1_000_000_000
        assert np.abs((res - td).value) < nanos
        assert res.value % nanos == 0

    def test_pickle(self):
        assert Timedelta.max == tm.round_trip_pickle(Timedelta.max)

    @pytest.mark.parametrize("num_days", range(20))
    def test_hash_equals_timedelta_hash(self, num_days: int):
        """GH#11129"""
        kwargs = {"days": num_days, "seconds": 1}
        td = Timedelta(**kwargs)  # type: ignore[arg-type]

        assert hash(td) == hash(timedelta(**kwargs))

    @pytest.mark.parametrize("ns", (1, 500))
    def test_hash_differs_from_timedelta_hash_if_ns_lost(self, ns: int):
        td = Timedelta(ns, "ns")
        assert hash(td) != hash(td.to_pytimedelta())

    @pytest.mark.parametrize("td_kwarg", TD_KWARGS)
    def test_only_zero_value_falsy(self, td_kwarg):
        """GH#21484"""
        assert bool(Timedelta(**{td_kwarg: 0})) is False
        assert bool(Timedelta(**{td_kwarg: 1})) is True
        assert bool(Timedelta(**{td_kwarg: -1})) is True

    @pytest.mark.parametrize(
        "td, expected_iso",
        [
            (
                Timedelta(days=6, milliseconds=123, nanoseconds=45),
                "P6DT0H0M0.123000045S",
            ),
            (Timedelta(days=4, hours=12, minutes=30, seconds=5), "P4DT12H30M5S"),
            (Timedelta(nanoseconds=123), "P0DT0H0M0.000000123S"),
            # trim nano
            (Timedelta(microseconds=10), "P0DT0H0M0.00001S"),
            # trim micro
            (Timedelta(milliseconds=1), "P0DT0H0M0.001S"),
            # don't strip every 0
            (Timedelta(minutes=1), "P0DT0H1M0S"),
        ],
    )
    def test_isoformat(self, td, expected_iso):
        assert td.isoformat() == expected_iso

    @pytest.mark.parametrize(
        ("value, expected"),
        (
            ("1 W", "7 days 00:00:00"),
            ("-1 W", "-7 days +00:00:00"),
            ("1 D", "1 days 00:00:00"),
            ("-1 D", "-1 days +00:00:00"),
            ("1 H", "0 days 01:00:00"),
            ("-1 H", "-1 days +23:00:00"),
            ("1 m", "0 days 00:01:00"),
            ("-1 m", "-1 days +23:59:00"),
            ("1 m", "0 days 00:01:00"),
            ("-1 m", "-1 days +23:59:00"),
            ("1 s", "0 days 00:00:01"),
            ("-1 s", "-1 days +23:59:59"),
            ("1 ms", "0 days 00:00:00.001000"),
            ("-1 ms", "-1 days +23:59:59.999000"),
            ("1 us", "0 days 00:00:00.000001"),
            ("-1 us", "-1 days +23:59:59.999999"),
            ("1 ns", "0 days 00:00:00.000000001"),
            ("-1 ns", "-1 days +23:59:59.999999999"),
        ),
    )
    def test_str_and_repr(self, value: str, expected: str):
        expected_repr = f"Timedelta('{expected}')"
        td = Timedelta(value)

        assert str(td) == expected
        assert repr(td) == expected_repr
        assert Timedelta(expected) == td


class TestToTimedelta:
    """Move elsewhere"""

    def test_iso_conversion(self):
        # GH #21877
        expected = Timedelta(1, unit="s")
        assert to_timedelta("P0DT0H0M1S") == expected

    def test_nat_converters(self):
        result = to_timedelta("nat").to_numpy()
        assert result.dtype.kind == "M"
        assert result.astype("int64") == iNaT

        result = to_timedelta("nan").to_numpy()
        assert result.dtype.kind == "M"
        assert result.astype("int64") == iNaT

    def test_contains(self):
        # Checking for any NaT-like objects
        # GH 13603
        td = to_timedelta(range(5), unit="d") + offsets.Hour(1)
        for v in [NaT, None, float("nan"), np.nan]:
            assert not (v in td)

        td = to_timedelta([NaT])
        for v in [NaT, None, float("nan"), np.nan]:
            assert v in td

        # invalid
        msg = "have leftover units"
        with pytest.raises(ValueError, match=msg):
            Timedelta("- 1days, 00")

    @pytest.mark.parametrize("unit, np_unit", TD_UNITS_TD64_UNITS.items())
    @pytest.mark.parametrize("wrapper", [np.array, list, pd.Index])
    def test_unit_parser(self, unit, np_unit, wrapper):
        # validate all units, GH 6855, GH 21762
        # array-likes
        expected = TimedeltaIndex(
            [np.timedelta64(i, np_unit) for i in np.arange(5).tolist()]
        )
        result = to_timedelta(wrapper(range(5)), unit=unit)
        tm.assert_index_equal(result, expected)
        result = TimedeltaIndex(wrapper(range(5)), unit=unit)
        tm.assert_index_equal(result, expected)

        str_repr = [f"{x}{unit}" for x in np.arange(5)]
        result = to_timedelta(wrapper(str_repr))
        tm.assert_index_equal(result, expected)
        result = to_timedelta(wrapper(str_repr))
        tm.assert_index_equal(result, expected)


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


def test_nan_total_seconds():
    # put elsewhere? a test of NaT, not Timedelta, behavior
    rng = Timedelta(np.nan)
    assert np.isnan(rng.total_seconds())
