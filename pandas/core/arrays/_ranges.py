"""
Helper functions to generate range-like data for DatetimeArray
(and possibly TimedeltaArray/PeriodArray)
"""

from typing import Tuple

import numpy as np

from pandas._libs.tslibs import OutOfBoundsDatetime, Timestamp

from pandas.tseries.offsets import DateOffset, Tick, generate_range


def generate_regular_range(
    start: Timestamp, end: Timestamp, periods: int, freq: DateOffset
) -> Tuple[np.ndarray, str]:
    """
    Generate a range of dates with the spans between dates described by
    the given `freq` DateOffset.

    Parameters
    ----------
    start : Timestamp or None
        first point of produced date range
    end : Timestamp or None
        last point of produced date range
    periods : int
        number of periods in produced date range
    freq : DateOffset
        describes space between dates in produced date range

    Returns
    -------
    ndarray[np.int64] representing nanosecond unix timestamps
    """
    if isinstance(freq, Tick):
        stride = freq.nanos
        if periods is None:
            b = Timestamp(start).value
            # cannot just use e = Timestamp(end) + 1 because arange breaks when
            # stride is too large, see GH10887
            e = b + (Timestamp(end).value - b) // stride * stride + stride // 2 + 1
            # end.tz == start.tz by this point due to _generate implementation
            tz = start.tz
        elif start is not None:
            b = Timestamp(start).value
            e = _generate_range_overflow_safe(b, periods, stride, side="start")
            tz = start.tz
        elif end is not None:
            e = Timestamp(end).value + stride
            b = _generate_range_overflow_safe(e, periods, stride, side="end")
            tz = end.tz
        else:
            raise ValueError(
                "at least 'start' or 'end' should be specified "
                "if a 'period' is given."
            )

        with np.errstate(over="raise"):
            # If the range is sufficiently large, np.arange may overflow
            #  and incorrectly return an empty array if not caught.
            try:
                values = np.arange(b, e, stride, dtype=np.int64)
            except FloatingPointError:
                xdr = [b]
                while xdr[-1] != e:
                    xdr.append(xdr[-1] + stride)
                values = np.array(xdr[:-1], dtype=np.int64)

    else:
        tz = None
        # start and end should have the same timezone by this point
        if start is not None:
            tz = start.tz
        elif end is not None:
            tz = end.tz

        xdr = generate_range(start=start, end=end, periods=periods, offset=freq)

        values = np.array([x.value for x in xdr], dtype=np.int64)

    return values, tz


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

    i64max = np.uint64(np.iinfo(np.int64).max)
    msg = f"Cannot generate range with {side}={endpoint} and periods={periods}"

    with np.errstate(over="raise"):
        # if periods * strides cannot be multiplied within the *uint64* bounds,
        #  we cannot salvage the operation by recursing, so raise
        try:
            addend = np.uint64(periods) * np.uint64(np.abs(stride))
        except FloatingPointError:
            raise OutOfBoundsDatetime(msg)

    if np.abs(addend) <= i64max:
        # relatively easy case without casting concerns
        return _generate_range_overflow_safe_signed(endpoint, periods, stride, side)

    elif (endpoint > 0 and side == "start" and stride > 0) or (
        endpoint < 0 and side == "end" and stride > 0
    ):
        # no chance of not-overflowing
        raise OutOfBoundsDatetime(msg)

    elif side == "end" and endpoint > i64max and endpoint - stride <= i64max:
        # in _generate_regular_range we added `stride` thereby overflowing
        #  the bounds.  Adjust to fix this.
        return _generate_range_overflow_safe(
            endpoint - stride, periods - 1, stride, side
        )

    # split into smaller pieces
    mid_periods = periods // 2
    remaining = periods - mid_periods
    assert 0 < remaining < periods, (remaining, periods, endpoint, stride)

    midpoint = _generate_range_overflow_safe(endpoint, mid_periods, stride, side)
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
            return np.int64(endpoint) + addend
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
            result = np.uint64(endpoint) + np.uint64(addend)
            i64max = np.uint64(np.iinfo(np.int64).max)
            assert result > i64max
            if result <= i64max + np.uint64(stride):
                return result

    raise OutOfBoundsDatetime(
        f"Cannot generate range with {side}={endpoint} and periods={periods}"
    )
