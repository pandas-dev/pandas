from datetime import (
    datetime,
    timedelta,
)

from dateutil.tz.tz import tzlocal
import numpy as np
import pytest

from pandas._libs.tslibs import (
    OutOfBoundsDatetime,
    Timestamp,
)
from pandas.compat import (
    IS64,
    WASM,
    is_platform_windows,
)

from pandas import (
    DatetimeIndex,
    Period,
    Timedelta,
)
import pandas._testing as tm

from pandas.tseries.offsets import (
    FY5253,
    BDay,
    BMonthBegin,
    BMonthEnd,
    BQuarterBegin,
    BQuarterEnd,
    BusinessHour,
    BYearBegin,
    BYearEnd,
    CBMonthBegin,
    CBMonthEnd,
    CDay,
    CustomBusinessHour,
    DateOffset,
    FY5253Quarter,
    HalfYearEnd,
    LastWeekOfMonth,
    MonthBegin,
    MonthEnd,
    QuarterEnd,
    SemiMonthBegin,
    SemiMonthEnd,
    Week,
    WeekOfMonth,
    YearBegin,
    YearEnd,
)


def _get_offset(klass, value=1, normalize=False):
    # create instance from offset class
    if klass is FY5253:
        klass = klass(
            n=value,
            startingMonth=1,
            weekday=1,
            variation="last",
            normalize=normalize,
        )
    elif klass is FY5253Quarter:
        klass = klass(
            n=value,
            startingMonth=1,
            weekday=1,
            qtr_with_extra_week=1,
            variation="last",
            normalize=normalize,
        )
    elif klass is LastWeekOfMonth:
        klass = klass(n=value, weekday=5, normalize=normalize)
    elif klass is WeekOfMonth:
        klass = klass(n=value, week=1, weekday=5, normalize=normalize)
    elif klass is Week:
        klass = klass(n=value, weekday=5, normalize=normalize)
    elif klass is DateOffset:
        klass = klass(days=value, normalize=normalize)
    else:
        klass = klass(value, normalize=normalize)
    return klass


@pytest.fixture(
    params=[
        BDay,
        BusinessHour,
        BMonthEnd,
        BMonthBegin,
        BQuarterEnd,
        BQuarterBegin,
        BYearEnd,
        BYearBegin,
        CDay,
        CustomBusinessHour,
        CBMonthEnd,
        CBMonthBegin,
        MonthEnd,
        MonthBegin,
        SemiMonthBegin,
        SemiMonthEnd,
        QuarterEnd,
        LastWeekOfMonth,
        WeekOfMonth,
        Week,
        YearBegin,
        YearEnd,
        FY5253,
        FY5253Quarter,
        DateOffset,
    ]
)
def _offset(request):
    return request.param


@pytest.mark.skipif(WASM, reason="OverflowError received on WASM")
def test_apply_out_of_range(request, tz_naive_fixture, _offset):
    tz = tz_naive_fixture

    # try to create an out-of-bounds result timestamp; if we can't create
    # the offset skip
    try:
        if _offset in (BusinessHour, CustomBusinessHour):
            # Using 10000 in BusinessHour fails in tz check because of DST
            # difference
            offset = _get_offset(_offset, value=100000)
        else:
            offset = _get_offset(_offset, value=10000)

        result = Timestamp("20080101") + offset
        assert isinstance(result, datetime)
        assert result.tzinfo is None

        # Check tz is preserved
        t = Timestamp("20080101", tz=tz)
        result = t + offset
        assert isinstance(result, datetime)
        if tz is not None:
            assert t.tzinfo is not None

        if (
            isinstance(tz, tzlocal)
            and ((not IS64) or WASM)
            and _offset is not DateOffset
        ):
            # If we hit OutOfBoundsDatetime on non-64 bit machines
            # we'll drop out of the try clause before the next test
            request.applymarker(
                pytest.mark.xfail(reason="OverflowError inside tzlocal past 2038")
            )
        elif (
            isinstance(tz, tzlocal)
            and is_platform_windows()
            and _offset in (QuarterEnd, BQuarterBegin, BQuarterEnd, FY5253Quarter)
        ):
            request.applymarker(
                pytest.mark.xfail(reason="After GH#49737 t.tzinfo is None on CI")
            )
        assert str(t.tzinfo) == str(result.tzinfo), (t.tzinfo, result.tzinfo)

    except OutOfBoundsDatetime:
        pass
    except (ValueError, KeyError, NotImplementedError, OverflowError, OSError):
        # we are creating an invalid offset
        # so ignore
        # - NotImplementedError is raised for tz-aware timestamps outside Python's
        #   range that are created by adding the offset
        # - OverflowError is raised in the same case on 32-bit systems with tzlocal
        # - OSError is raised in the same case on Windows with tzlocal
        pass


@pytest.mark.parametrize(
    "offset",
    [
        DateOffset(months=1),
        MonthBegin(),
        MonthEnd(),
        QuarterEnd(),
        YearEnd(),
        SemiMonthBegin(),
        SemiMonthEnd(),
    ],
)
def test_apply_array_out_of_bounds_raises(offset):
    # GH#66549 the vectorized month/quarter shifts raised a bare OverflowError
    #  naming a C source line, where the scalar path raises OutOfBoundsDatetime
    dti = DatetimeIndex([Timestamp.max])

    with pytest.raises(OutOfBoundsDatetime, match="Out of bounds nanosecond"):
        dti + offset

    with pytest.raises(OutOfBoundsDatetime):
        Timestamp.max + offset


@pytest.mark.parametrize(
    "offset",
    [
        QuarterEnd(2**62),
        YearEnd(2**62),
        YearBegin(-(2**63)),
        BYearEnd(2**63 - 1),
        HalfYearEnd(-(2**63)),
    ],
)
def test_apply_array_quarters_large_n_raises(offset):
    # GH#66549 the quarter shift computed `modby * n` in wrapping int64, so
    #  e.g. YearEnd(2**62) shifted by 12 * 2**62 == 0 months instead of raising
    dti = DatetimeIndex(["1970-06-15"])

    with pytest.raises(OutOfBoundsDatetime, match="Out of bounds"):
        dti + offset

    with pytest.raises(OutOfBoundsDatetime, match="Out of bounds"):
        dti[0] + offset


def test_apply_array_quarters_large_n_names_same_year_as_scalar():
    # GH#66549 the wrapped shift amount also has to be redone in Python ints
    #  for the message to name the year the scalar path names
    dti = DatetimeIndex(["1970-06-15"])
    offset = YearEnd(2**62)

    with pytest.raises(OutOfBoundsDatetime, match="year 4611686018427389873"):
        dti + offset

    with pytest.raises(OutOfBoundsDatetime, match="4611686018427389873-12-31"):
        dti[0] + offset


def test_shift_month_year_overflow_names_shifted_year():
    # GH#66549 the scalar path formed `dts.year + dy` in int64, so an n this
    #  large wrapped the year negative in the error message
    ts = Timestamp("1970-06-15")

    with pytest.raises(OutOfBoundsDatetime, match="year 9223372036854777776"):
        ts + BYearEnd(2**63 - 1)


def test_apply_array_quarters_large_n_keeps_nat():
    # GH#66549 a representable shift is unaffected by the overflow check
    dti = DatetimeIndex(["1970-06-15", "NaT"])

    result = dti + YearEnd(2)

    expected = DatetimeIndex(["1971-12-31", "NaT"]).as_unit(result.unit)
    tm.assert_index_equal(result, expected)


def test_apply_array_quarters_large_n_all_nat():
    # GH#66549 an all-NaT input has nothing to shift, so even an n that no
    #  element could survive comes back all-NaT rather than raising
    dti = DatetimeIndex(["NaT", "NaT"])

    result = dti + YearEnd(2**62)

    tm.assert_index_equal(result, dti.as_unit(result.unit))


def test_apply_array_out_of_bounds_raises_non_nano():
    # GH#66549 the message names the resolution the result overflowed
    dti = DatetimeIndex(np.array([np.iinfo(np.int64).max], dtype="M8[s]"))

    with pytest.raises(OutOfBoundsDatetime, match="Out of bounds second"):
        dti + MonthEnd()


@pytest.mark.parametrize(
    "offset, expected",
    [
        (DateOffset(months=1, days=-31), "2262-04-10 23:47:16.854775807"),
        (DateOffset(years=1, days=-366), "2262-04-10 23:47:16.854775807"),
        (DateOffset(months=2, days=-62), "2262-04-10 23:47:16.854775807"),
        (DateOffset(months=1, hours=-745), "2262-04-10 22:47:16.854775807"),
    ],
)
def test_apply_array_composed_out_of_bounds_intermediate(offset, expected):
    # GH#66549 the month shift used to be materialized before the timedelta
    #  component was applied, so a composed offset whose final result is
    #  representable was rejected on the array path
    dti = DatetimeIndex([Timestamp.max])

    result = dti + offset
    expected_dti = DatetimeIndex([Timestamp(expected)])
    tm.assert_index_equal(result, expected_dti)
    assert Timestamp.max + offset == expected_dti[0]


def test_apply_array_composed_out_of_bounds_intermediate_below_min():
    # GH#66549 same going the other way off the bottom of the range
    dti = DatetimeIndex([Timestamp.min])
    offset = DateOffset(months=-1, days=32)

    result = dti + offset
    expected = DatetimeIndex([Timestamp("1677-09-22 00:12:43.145224193")])
    tm.assert_index_equal(result, expected)
    assert Timestamp.min + offset == expected[0]


def test_apply_array_composed_still_out_of_bounds():
    # GH#66549 a composed offset whose result is genuinely unrepresentable
    #  still raises, and names the resolution rather than a C source line
    dti = DatetimeIndex([Timestamp.max])

    with pytest.raises(OutOfBoundsDatetime, match="Out of bounds nanosecond"):
        dti + DateOffset(months=1, days=-1)


def test_apply_array_composed_out_of_bounds_intermediate_non_nano():
    # GH#66549 the intermediate is out of the second-resolution range
    dti = DatetimeIndex(np.array([np.iinfo(np.int64).max], dtype="M8[s]"))

    result = dti + DateOffset(months=1, days=-31)

    tm.assert_index_equal(result, dti)


def test_apply_array_composed_keeps_nat():
    # GH#66549 the fused path leaves missing values missing
    dti = DatetimeIndex(["2000-01-31", "NaT"])

    result = dti + DateOffset(months=1, days=1)

    expected = DatetimeIndex(["2000-03-01", "NaT"]).as_unit(result.unit)
    tm.assert_index_equal(result, expected)


def test_offsets_compare_equal(_offset):
    # root cause of GH#456: __ne__ was not implemented
    offset1 = _offset()
    offset2 = _offset()
    assert not offset1 != offset2
    assert offset1 == offset2


@pytest.mark.parametrize(
    "date, offset2",
    [
        [Timestamp(2008, 1, 1), BDay(2)],
        [Timestamp(2014, 7, 1, 10, 00), BusinessHour(n=3)],
        [
            Timestamp(2014, 7, 1, 10),
            CustomBusinessHour(
                holidays=["2014-06-27", Timestamp(2014, 6, 30), Timestamp("2014-07-02")]
            ),
        ],
        [Timestamp(2008, 1, 2), SemiMonthEnd(2)],
        [Timestamp(2008, 1, 2), SemiMonthBegin(2)],
        [Timestamp(2008, 1, 2), Week(2)],
        [Timestamp(2008, 1, 2), WeekOfMonth(2)],
        [Timestamp(2008, 1, 2), LastWeekOfMonth(2)],
    ],
)
def test_rsub(date, offset2):
    assert date - offset2 == (-offset2) + date


@pytest.mark.parametrize(
    "date, offset2",
    [
        [Timestamp(2008, 1, 1), BDay(2)],
        [Timestamp(2014, 7, 1, 10, 00), BusinessHour(n=3)],
        [
            Timestamp(2014, 7, 1, 10),
            CustomBusinessHour(
                holidays=["2014-06-27", Timestamp(2014, 6, 30), Timestamp("2014-07-02")]
            ),
        ],
        [Timestamp(2008, 1, 2), SemiMonthEnd(2)],
        [Timestamp(2008, 1, 2), SemiMonthBegin(2)],
        [Timestamp(2008, 1, 2), Week(2)],
        [Timestamp(2008, 1, 2), WeekOfMonth(2)],
        [Timestamp(2008, 1, 2), LastWeekOfMonth(2)],
    ],
)
def test_radd(date, offset2):
    assert date + offset2 == offset2 + date


@pytest.mark.parametrize(
    "date, offset_box, offset2",
    [
        [Timestamp(2008, 1, 1), BDay, BDay(2)],
        [Timestamp(2008, 1, 2), SemiMonthEnd, SemiMonthEnd(2)],
        [Timestamp(2008, 1, 2), SemiMonthBegin, SemiMonthBegin(2)],
        [Timestamp(2008, 1, 2), Week, Week(2)],
        [Timestamp(2008, 1, 2), WeekOfMonth, WeekOfMonth(2)],
        [Timestamp(2008, 1, 2), LastWeekOfMonth, LastWeekOfMonth(2)],
    ],
)
def test_sub(date, offset_box, offset2):
    off = offset2
    msg = "Cannot subtract datetime from offset"
    with pytest.raises(TypeError, match=msg):
        off - date

    assert 2 * off - off == off
    assert date - offset2 == date + offset_box(-2)
    assert date - offset2 == date - (2 * off - off)


@pytest.mark.parametrize(
    "offset_box, offset1",
    [
        [BDay, BDay()],
        [LastWeekOfMonth, LastWeekOfMonth()],
        [WeekOfMonth, WeekOfMonth()],
        [Week, Week()],
        [SemiMonthBegin, SemiMonthBegin()],
        [SemiMonthEnd, SemiMonthEnd()],
        [CustomBusinessHour, CustomBusinessHour(weekmask="Tue Wed Thu Fri")],
        [BusinessHour, BusinessHour()],
    ],
)
def test_Mult1(offset_box, offset1):
    dt = Timestamp(2008, 1, 2)
    assert dt + 10 * offset1 == dt + offset_box(10)
    assert dt + 5 * offset1 == dt + offset_box(5)


@pytest.mark.parametrize("other", ["foo", 3, 3.5, None])
def test_add_unsupported_type_raises(_offset, other):
    # GH#36590 the standard "unsupported operand type(s)" TypeError, not a
    #  cython signature error leaking out of the offset's _add_datetime
    off = _get_offset(_offset)
    msg = "unsupported operand type"
    with pytest.raises(TypeError, match=msg):
        off + other


def test_add_period_defers_to_period():
    # declining the operation lets python fall back to Period.__radd__
    assert MonthEnd() + Period("2022-01", freq="M") == Period("2022-02", freq="M")


def test_add_offset_raises(_offset):
    # GH#36590 offset composition is not supported (GH#10902); the error should
    #  be an ordinary one, not a cython signature error naming datetime.datetime
    off = _get_offset(_offset)
    msg = "|".join(["unsupported operand type", "Cannot add"])
    with pytest.raises(TypeError, match=msg):
        off + BMonthEnd()


@pytest.mark.parametrize(
    "other", [timedelta(days=1), np.timedelta64(1, "D"), Timedelta(days=1)]
)
@pytest.mark.parametrize("offset", [MonthEnd(), YearEnd(), DateOffset(months=1)])
def test_add_timedelta_to_anchored_offset_raises(offset, other):
    # GH#36590 anchored offsets have no timedelta behaviour to fall back on
    # older numpy versions raise a UFuncTypeError (a TypeError subclass) for
    #  np.timedelta64 instead of deferring to the python message
    msg = "|".join(["unsupported operand type", "cannot use operands with types"])
    with pytest.raises(TypeError, match=msg):
        offset + other


@pytest.mark.parametrize(
    "other", [timedelta(hours=3), np.timedelta64(3, "h"), Timedelta(hours=3)]
)
def test_add_timedelta_to_business_day(other):
    # a timedelta folds into the business day offset's `offset` attribute
    expected = BDay(offset=timedelta(hours=3))
    assert BDay() + other == expected
    assert other + BDay() == expected


def test_compare_str(_offset):
    # GH#23524
    # comparing to strings that cannot be cast to DateOffsets should
    #  not raise for __eq__ or __ne__
    off = _get_offset(_offset)

    assert not off == "infer"
    assert off != "foo"
    # Note: inequalities are only implemented for Tick subclasses;
    #  tests for this are in test_ticks
