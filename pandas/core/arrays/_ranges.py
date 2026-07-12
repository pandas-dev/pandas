"""
Helper functions to generate range-like data for DatetimeArray
(and possibly TimedeltaArray/PeriodArray)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pandas._libs.lib import i8max
from pandas._libs.tslibs import (
    BaseOffset,
    Day,
    OutOfBoundsDatetime,
    Timedelta,
    Timestamp,
    iNaT,
)

from pandas.core.construction import range_to_ndarray

if TYPE_CHECKING:
    from pandas._typing import (
        TimeUnit,
        npt,
    )


def generate_regular_range(
    start: Timestamp | Timedelta | None,
    end: Timestamp | Timedelta | None,
    periods: int | None,
    freq: BaseOffset,
    unit: TimeUnit = "ns",
) -> npt.NDArray[np.intp]:
    """
    Generate a range of dates or timestamps with the spans between dates
    described by the given `freq` DateOffset.

    Parameters
    ----------
    start : Timedelta, Timestamp or None
        First point of produced date range.
    end : Timedelta, Timestamp or None
        Last point of produced date range.
    periods : int or None
        Number of periods in produced date range.
    freq : Tick
        Describes space between dates in produced date range.
    unit : {'s', 'ms', 'us', 'ns'}, default "ns"
        The resolution the output is meant to represent.

    Returns
    -------
    ndarray[np.int64]
        Representing the given resolution.
    """
    istart = start._value if start is not None else None
    iend = end._value if end is not None else None
    if isinstance(freq, Day):
        # In contexts without a timezone, a Day offset is unambiguously
        #  interpretable as Timedelta-like.
        td = Timedelta(days=freq.n)
    else:
        freq.nanos  # raises if non-fixed frequency
        td = Timedelta(freq)
    b: int
    e: int
    try:
        td = td.as_unit(unit, round_ok=False)
    except ValueError as err:
        raise ValueError(
            f"freq={freq} is incompatible with unit={unit}. "
            "Use a lower freq or a higher unit instead."
        ) from err
    stride = int(td._value)

    if periods is None and istart is not None and iend is not None:
        b = istart
        # cannot just use e = Timestamp(end) + 1 because arange breaks when
        # stride is too large, see GH10887
        e = b + (iend - b) // stride * stride + stride // 2 + 1
    elif istart is not None and periods is not None:
        b = istart
        e = _generate_range_overflow_safe(b, periods, stride, side="start")
    elif iend is not None and periods is not None:
        e = iend + stride
        b = _generate_range_overflow_safe(e, periods, stride, side="end")
    else:
        raise ValueError(
            "at least 'start' or 'end' should be specified if a 'period' is given."
        )

    return range_to_ndarray(range(b, e, stride))


def generate_daily_offset_range(
    start: Timestamp | None,
    end: Timestamp | None,
    periods: int | None,
    freq: BaseOffset,
    unit: TimeUnit = "ns",
) -> npt.NDArray[np.int64]:
    """
    Generate a range for offsets whose on-offset dates are a subset of a
    daily grid, by generating a daily range and filtering.

    This is a performance optimization (GH#16463) for offsets like BusinessDay
    that implement ``_get_daily_offset_mask``. Instead of iteratively applying
    the offset, we generate a daily range and filter it vectorized.

    Parameters
    ----------
    start : Timestamp or None
    end : Timestamp or None
    periods : int or None
    freq : BaseOffset
        Must have ``_supports_daily_offset_mask`` and ``n >= 1``
        (caller ensures both).
    unit : str, default "ns"

    Returns
    -------
    ndarray[int64]
        The filtered i8 values.
    """
    abs_n = freq.n  # caller ensures n >= 1

    # Resolve the endpoints with offset arithmetic so the on-offset count is
    # exact regardless of holidays/weekmasks. The old ``needed * 7 + 6``-day
    # buffer assumed >=1 on-offset day per 7 calendar days, which fails for e.g.
    # a CustomBusinessDay whose holidays blank out whole weeks -- the result was
    # then silently short (GH#64648 post-merge). Arithmetic runs on normalized
    # dates with the time-of-day re-attached, so ``normalize=True`` offsets
    # don't strip it from the output (GH#44025).
    if periods is not None and end is not None:
        # end + periods: anchor at the last on-offset date <= end so that
        # ``periods`` on-offset dates are returned (GH#64834).
        tod: Timedelta = end - end.normalize()
        # pyright can't see through the NaT-returning Timestamp constructor.
        anchor: Timestamp = Timestamp(  # pyright: ignore[reportAssignmentType]
            freq.rollback(end.normalize())
        )
        start = (anchor - (periods - 1) * freq + tod).as_unit(unit)
        end = (anchor + tod).as_unit(unit)
    else:
        # start is guaranteed non-None here by the caller (exactly two of
        # start/end/periods are given, and this branch is not end+periods).
        assert start is not None
        # Roll an off-offset start forward to the first on-offset date (n >= 1).
        tod = start - start.normalize()
        anchor = Timestamp(  # pyright: ignore[reportAssignmentType]
            freq.rollforward(start.normalize())
        )
        start = (anchor + tod).as_unit(unit)
        if periods is not None:
            # start + periods.
            end = (anchor + (periods - 1) * freq + tod).as_unit(unit)
        elif (
            # start + end. GH#64790: align end's time-of-day to start's so the
            # last on-offset date isn't dropped when end's is earlier. Mask
            # offsets preserve time-of-day unless ``freq.normalize``; testing
            # that attribute rather than probing ``freq._apply(start)`` avoids
            # raising near Timestamp.max (GH#64648).
            tod and not freq.normalize and end >= start  # type: ignore[operator]
        ):
            try:
                end = (end.normalize() + tod).as_unit(unit)  # type: ignore[union-attr]
            except OutOfBoundsDatetime:
                # the boundary element this alignment would admit is beyond
                # Timestamp.max, so keeping the raw end excludes it correctly
                pass

    i8values = generate_regular_range(start, end, None, Day(), unit=unit)
    dt64 = i8values.view(f"datetime64[{unit}]")
    i8values = i8values[freq._get_daily_offset_mask(dt64)]  # type: ignore[attr-defined]

    if abs_n > 1:
        # start is the first on-offset date and end the last, so a forward
        # stride anchored at the front lands on end (GH#64648, GH#65604).
        i8values = i8values[::abs_n]

    return i8values


def _generate_range_overflow_safe(
    endpoint: int, periods: int, stride: int, side: str = "start"
) -> int:
    """
    Calculate the second endpoint for passing to np.arange, checking
    to avoid an integer overflow.  Catch OverflowError and re-raise
    as OutOfBoundsDatetime.

    Parameters
    ----------
    endpoint : int
        nanosecond timestamp of the known endpoint of the desired range
    periods : int
        number of periods in the desired range
    stride : int
        nanoseconds between periods in the desired range
    side : {'start', 'end'}
        which end of the range `endpoint` refers to

    Returns
    -------
    other_end : int

    Raises
    ------
    OutOfBoundsDatetime
    """
    # GH#14187 raise instead of incorrectly wrapping around
    assert side in ["start", "end"]

    i64max = np.uint64(i8max)
    msg = f"Cannot generate range with {side}={endpoint} and periods={periods}"

    with np.errstate(over="raise"):
        # if periods * strides cannot be multiplied within the *uint64* bounds,
        #  we cannot salvage the operation by recursing, so raise
        try:
            addend = np.uint64(periods) * np.uint64(np.abs(stride))
        except FloatingPointError as err:
            raise OutOfBoundsDatetime(msg) from err

    if np.abs(addend) <= i64max:
        # relatively easy case without casting concerns
        return _generate_range_overflow_safe_signed(endpoint, periods, stride, side)

    elif (endpoint > 0 and side == "start" and stride > 0) or (
        endpoint < 0 < stride and side == "end"
    ):
        # no chance of not-overflowing
        raise OutOfBoundsDatetime(msg)

    elif side == "end" and endpoint - stride <= i64max < endpoint:
        # in _generate_regular_range we added `stride` thereby overflowing
        #  the bounds.  Adjust to fix this.
        return _generate_range_overflow_safe(
            endpoint - stride, periods - 1, stride, side
        )

    # split into smaller pieces
    mid_periods = periods // 2
    remaining = periods - mid_periods
    assert 0 < remaining < periods, (remaining, periods, endpoint, stride)

    midpoint = int(_generate_range_overflow_safe(endpoint, mid_periods, stride, side))
    return _generate_range_overflow_safe(midpoint, remaining, stride, side)


def _generate_range_overflow_safe_signed(
    endpoint: int, periods: int, stride: int, side: str
) -> int:
    """
    A special case for _generate_range_overflow_safe where `periods * stride`
    can be calculated without overflowing int64 bounds.
    """
    assert side in ["start", "end"]
    if side == "end":
        stride *= -1

    with np.errstate(over="raise"):
        addend = np.int64(periods) * np.int64(stride)
        try:
            # easy case with no overflows
            result = np.int64(endpoint) + addend
            if result == iNaT:
                # Putting this into a DatetimeArray/TimedeltaArray
                #  would incorrectly be interpreted as NaT
                raise OverflowError
            return int(result)
        except (FloatingPointError, OverflowError):
            # with endpoint negative and addend positive we risk
            #  FloatingPointError; with reversed signed we risk OverflowError
            pass

        # if stride and endpoint had opposite signs, then endpoint + addend
        #  should never overflow.  so they must have the same signs
        assert (stride > 0 and endpoint >= 0) or (stride < 0 and endpoint <= 0)

        if stride > 0:
            # watch out for very special case in which we just slightly
            #  exceed implementation bounds, but when passing the result to
            #  np.arange will get a result slightly within the bounds

            uresult = np.uint64(endpoint) + np.uint64(addend)
            i64max = np.uint64(i8max)
            assert uresult > i64max
            if uresult <= i64max + np.uint64(stride):
                return int(uresult)

    raise OutOfBoundsDatetime(
        f"Cannot generate range with {side}={endpoint} and periods={periods}"
    )
