# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

import datetime as dt
import operator as op
import sys
import warnings
import zoneinfo
from functools import cache, partial
from importlib import resources
from pathlib import Path
from typing import TYPE_CHECKING, Annotated, overload

from hypothesis.errors import InvalidArgument
from hypothesis.internal.validation import check_type, check_valid_interval
from hypothesis.strategies._internal.core import sampled_from
from hypothesis.strategies._internal.lazy import unwrap_strategies
from hypothesis.strategies._internal.misc import just, none, nothing
from hypothesis.strategies._internal.strategies import (
    FilteredStrategy,
    OneOfStrategy,
    SampledFromStrategy,
    SearchStrategy,
)
from hypothesis.strategies._internal.utils import defines_strategy

if TYPE_CHECKING:
    from annotated_types import Timezone

    NaiveDatetime = Annotated[dt.datetime, Timezone(None)]
    AwareDatetime = Annotated[dt.datetime, Timezone(...)]
elif at := sys.modules.get("annotated_types"):
    NaiveDatetime = Annotated[dt.datetime, at.Timezone(None)]
    AwareDatetime = Annotated[dt.datetime, at.Timezone(...)]
else:
    NaiveDatetime = AwareDatetime = dt.datetime

DATENAMES = ("year", "month", "day")
TIMENAMES = ("hour", "minute", "second", "microsecond")

_MICROSECOND = dt.timedelta(microseconds=1)


def _comparator_bound(condition):
    """Return ``(op, bound)`` for filter conditions like ``partial(op, bound)``,
    with a single positional argument to one of the five comparison operators."""
    if (
        isinstance(condition, partial)
        and len(condition.args) == 1
        and not condition.keywords
        and condition.func in (op.lt, op.le, op.eq, op.ge, op.gt)
    ):
        return condition.func, condition.args[0]
    return None


def _narrowed_bounds(func, arg, min_value, max_value, shift):
    """Narrow [min_value, max_value] to satisfy the condition ``func(arg, x)``.

    ``shift(value, steps)`` moves value by that many of the smallest representable
    steps, raising OverflowError if the result would be unrepresentable.  Returns
    the narrowed (min_value, max_value), or None if no values can satisfy the
    condition.
    """
    if func in (op.lt, op.gt):
        try:
            arg = shift(arg, 1 if func is op.lt else -1)
        except OverflowError:  # gt the maximum value, or lt the minimum
            return None
    lo, hi = {
        # We're talking about op(arg, x) - the reverse of our usual intuition!
        op.lt: (arg, max_value),  # lambda x: arg < x
        op.le: (arg, max_value),  # lambda x: arg <= x
        op.eq: (arg, arg),  #       lambda x: arg == x
        op.ge: (min_value, arg),  # lambda x: arg >= x
        op.gt: (min_value, arg),  # lambda x: arg > x
    }[func]
    lo = max(lo, min_value)
    hi = min(hi, max_value)
    if hi < lo:
        return None
    return lo, hi


def _timezones_kind(strat):
    """Classify the values a timezones= strategy can generate: "none" if only
    None, "aware" if only tzinfo instances, or "unknown" if we can't tell."""
    strat = unwrap_strategies(strat)
    if isinstance(strat, SampledFromStrategy) and all(
        name == "filter" for name, _ in strat._transformations
    ):
        kinds = {
            "none" if e is None else "aware" if isinstance(e, dt.tzinfo) else "unknown"
            for e in strat.elements
        }
        return kinds.pop() if len(kinds) == 1 else "unknown"
    if isinstance(strat, OneOfStrategy):
        kinds = {_timezones_kind(s) for s in strat.original_strategies}
        return kinds.pop() if len(kinds) == 1 else "unknown"
    return "unknown"


def is_pytz_timezone(tz):
    if not isinstance(tz, dt.tzinfo):
        return False
    module = type(tz).__module__
    return module == "pytz" or module.startswith("pytz.")


def replace_tzinfo(value, timezone):
    if is_pytz_timezone(timezone):
        # Pytz timezones are a little complicated, and using the .replace method
        # can cause some weird issues, so we use their special "localize" instead.
        #
        # We use the fold attribute as a convenient boolean for is_dst, even though
        # they're semantically distinct.  For ambiguous or imaginary hours, fold says
        # whether you should use the offset that applies before the gap (fold=0) or
        # the offset that applies after the gap (fold=1). is_dst says whether you
        # should choose the side that is "DST" or "STD" (STD->STD or DST->DST
        # transitions are unclear as you might expect).
        #
        # WARNING: this is INCORRECT for timezones with negative DST offsets such as
        #       "Europe/Dublin", but it's unclear what we could do instead beyond
        #       documenting the problem and recommending use of `dateutil` instead.
        return timezone.localize(value, is_dst=not value.fold)
    return value.replace(tzinfo=timezone)


def _instant(value):
    """A sort key ordering aware datetimes by the moment they refer to.

    Unlike comparison of datetimes which share a tzinfo - which falls back to
    ignoring both the timezone and the fold attribute - this respects the fold,
    and unlike .astimezone() it cannot overflow near datetime.min/max.
    """
    return value.replace(tzinfo=None) - dt.datetime.min - value.utcoffset()


def _ambiguous(value, tz):
    # Whether the naive value is inside a DST fold, i.e. is a wall time which
    # occurs twice in tz, so that its utcoffset depends on the fold attribute.
    return (
        replace_tzinfo(value.replace(fold=0), tz).utcoffset()
        != replace_tzinfo(value.replace(fold=1), tz).utcoffset()
    )


def datetime_does_not_exist(value):
    """This function tests whether the given datetime can be round-tripped to and
    from UTC.  It is an exact inverse of (and very similar to) the dateutil method
    https://dateutil.readthedocs.io/en/stable/tz.html#dateutil.tz.datetime_exists
    """
    # Naive datetimes cannot be imaginary, but we need this special case because
    # chaining .astimezone() ends with *the system local timezone*, not None.
    # See bug report in https://github.com/HypothesisWorks/hypothesis/issues/2662
    if value.tzinfo is None:
        return False
    try:
        # Does the naive portion of the datetime change when round-tripped to
        # UTC?  If so, or if this overflows, we say that it does not exist.
        roundtrip = value.astimezone(dt.timezone.utc).astimezone(value.tzinfo)
    except OverflowError:
        # Overflows at datetime.min or datetime.max boundary condition.
        # Rejecting these is acceptable, because timezones are close to
        # meaningless before ~1900 and subject to a lot of change by
        # 9999, so it should be a very small fraction of possible values.
        return True

    if (
        value.tzinfo is not roundtrip.tzinfo
        and value.utcoffset() != roundtrip.utcoffset()
    ):
        # This only ever occurs during imaginary (i.e. nonexistent) datetimes,
        # and only for pytz timezones which do not follow PEP-495 semantics.
        # (may exclude a few other edge cases, but you should use zoneinfo anyway)
        return True

    assert value.tzinfo is roundtrip.tzinfo, "so only the naive portions are compared"
    return value != roundtrip


def _num_days_in_month(year, month):
    """Branchless equivalent of ``monthrange(year, month)[1]`` for valid inputs.

    Written using only arithmetic and (in)equality, with no branching or indexing.
    This avoids concretizing the input or adding more path constraints than necessary.
    """
    leap = (year % 4 == 0) * (1 - (year % 100 == 0) * (year % 400 != 0))
    is_feb = month == 2
    is_30_day = 1 - (month != 4) * (month != 6) * (month != 9) * (month != 11)
    return 31 - is_30_day - is_feb * (3 - leap)


def draw_capped_multipart(
    data, min_value, max_value, duration_names=DATENAMES + TIMENAMES
):
    assert isinstance(min_value, (dt.date, dt.time, dt.datetime))
    assert type(min_value) == type(max_value)
    assert min_value <= max_value

    # cap_{low, high} records whether every field drawn so far has equalled
    # ``min_value``'s / ``max_value``'s, i.e. whether that bound is still "active" and
    # constrains the next field.
    #
    # cap_{low, high} are conceptually booleans. We define them as integers and interpret
    # boolean operations on them as multiplication, so that we don't concretize or
    # branch under symbolic backends. See
    # https://github.com/HypothesisWorks/hypothesis/issues/4759.
    cap_low = 1
    cap_high = 1
    result = {}
    for name in duration_names:
        natural_low = getattr(dt.datetime.min, name)
        if name == "day":
            natural_high = _num_days_in_month(result["year"], result["month"])
        else:
            natural_high = getattr(dt.datetime.max, name)
        # equivalent to:
        #   low  = min_value.<name> if cap_low  else natural_low
        #   high = max_value.<name> if cap_high else natural_high
        low = natural_low + cap_low * (getattr(min_value, name) - natural_low)
        high = natural_high + cap_high * (getattr(max_value, name) - natural_high)
        if name == "year":
            val = data.draw_integer(low, high, shrink_towards=2000)
        else:
            val = data.draw_integer(low, high)
        result[name] = val
        cap_low = cap_low * (val == low)
        cap_high = cap_high * (val == high)
    if hasattr(min_value, "fold"):
        # The `fold` attribute is ignored in comparison of naive datetimes.
        # In tz-aware datetimes it would require *very* invasive changes to
        # the logic above, and be very sensitive to the specific timezone
        # (at the cost of efficient shrinking and mutation), so at least for
        # now we stick with the status quo and generate it independently.
        result["fold"] = data.draw_integer(0, 1)
    return result


def _shift_datetime(value, steps):
    return value + steps * _MICROSECOND


class DatetimeStrategy(SearchStrategy):
    def __init__(self, min_value, max_value, timezones_strat, allow_imaginary):
        super().__init__()
        assert isinstance(timezones_strat, SearchStrategy)
        assert isinstance(allow_imaginary, bool)
        self.aware = (min_value is not None and min_value.tzinfo is not None) or (
            max_value is not None and max_value.tzinfo is not None
        )
        if self.aware:
            for value in (min_value, max_value):
                assert value is None or (
                    isinstance(value, dt.datetime) and value.tzinfo is not None
                )
            # The instants bounding this strategy, as _instant() sort keys.
            # UTC offsets are less than a day, so a None bound is replaced by
            # a key which lies outside the representable range.
            self.min_instant = (
                dt.timedelta(days=-2) if min_value is None else _instant(min_value)
            )
            self.max_instant = (
                dt.datetime.max - dt.datetime.min + dt.timedelta(days=2)
                if max_value is None
                else _instant(max_value)
            )
            assert self.min_instant <= self.max_instant
        else:
            for value in (min_value, max_value):
                assert isinstance(value, dt.datetime)
                assert value.tzinfo is None
            assert min_value <= max_value
        self.min_value = min_value
        self.max_value = max_value
        self.tz_strat = timezones_strat
        self.allow_imaginary = allow_imaginary

    def do_draw(self, data):
        # We start by drawing a timezone, and an initial datetime.
        tz = data.draw(self.tz_strat)
        if self.aware:
            if not isinstance(tz, dt.tzinfo):
                raise InvalidArgument(
                    f"Drew {tz!r} from the timezones strategy {self.tz_strat!r}, "
                    "but with aware min_value/max_value bounds the timezones "
                    "strategy must only generate tzinfo objects (not None)"
                )
            result = self.draw_aware_datetime(data, tz)
        else:
            result = self.draw_naive_datetime_and_combine(data, tz)

        # TODO: with some probability, systematically search for one of
        #   - an imaginary time (if allowed),
        #   - a time within 24hrs of a leap second (if there any are within bounds),
        #   - other subtle, little-known, or nasty issues as described in
        #     https://github.com/HypothesisWorks/hypothesis/issues/69

        # If we happened to end up with a disallowed imaginary time, reject it.
        if (not self.allow_imaginary) and datetime_does_not_exist(result):
            data.mark_invalid(f"{result} does not exist (usually a DST transition)")
        return result

    def in_bounds(self, value):
        return self.min_instant <= _instant(value) <= self.max_instant

    def draw_aware_datetime(self, data, tz):
        def wall_clock(bound, extreme):
            if bound is None:
                return extreme
            try:
                return bound.astimezone(tz).replace(tzinfo=None)
            except OverflowError:
                # UTC offsets are less than a day, so an overflowing bound
                # must be within a day of datetime.min/max, converting to a
                # moment beyond them.  If every wall time representable in tz
                # is on the in-bounds side, the bound is simply vacuous here;
                # otherwise nothing in tz is in bounds.
                near_min = bound.replace(tzinfo=None) - dt.datetime.min < dt.timedelta(
                    days=2
                )
                if near_min == (extreme is dt.datetime.min):
                    return extreme
                data.mark_invalid(f"{bound!r} is not representable in {tz!r}")

        min_local = wall_clock(self.min_value, dt.datetime.min)
        max_local = wall_clock(self.max_value, dt.datetime.max)
        if min_local > max_local or (
            max_local - min_local <= dt.timedelta(days=1)
            and (_ambiguous(min_local, tz) or _ambiguous(max_local, tz))
        ):
            # A large fraction of the wall times between bounds inside or close
            # to a DST fold would risk rejection below - and bounds inside the
            # same fold may even be in inverted wall-clock order, like
            # 01:59 EDT < 01:01 EST - so we recurse to draw in UTC, where wall
            # times are unambiguous and ordered, and convert.  This is the
            # standard draw with the standard shrink order, except that
            # simplicity is judged on the UTC wall time rather than the local.
            value = self.draw_aware_datetime(data, dt.timezone.utc)
            try:
                return value.astimezone(tz)
            except OverflowError:
                data.mark_invalid(f"{value!r} is not representable in {tz!r}")
        result = draw_capped_multipart(data, min_local, max_local)
        value = replace_tzinfo(dt.datetime(**result), timezone=tz)
        if not self.in_bounds(value):
            # An ambiguous wall time next to a bound, with the out-of-bounds fold.
            data.mark_invalid(f"{value!r} is outside the bounds")
        return value

    def draw_naive_datetime_and_combine(self, data, tz):
        result = draw_capped_multipart(data, self.min_value, self.max_value)
        try:
            return replace_tzinfo(dt.datetime(**result), timezone=tz)
        except (ValueError, OverflowError):
            data.mark_invalid(
                f"Failed to draw a datetime between {self.min_value!r} and "
                f"{self.max_value!r} with timezone from {self.tz_strat!r}."
            )

    def filter(self, condition):
        if (parsed := _comparator_bound(condition)) is not None and isinstance(
            arg := parsed[1], dt.datetime
        ):
            func = parsed[0]
            try:
                bound_aware = arg.utcoffset() is not None
            except Exception:
                # A tzinfo whose utcoffset() raises; comparing against this
                # bound will raise the same error at draw time.
                return super().filter(condition)
            if not bound_aware:
                # The bound compares as naive (either no tzinfo, or a tzinfo
                # without a UTC offset), so we can only rewrite it into the
                # naive wall-clock bounds if it really is naive and every
                # generated value is too.
                if (
                    arg.tzinfo is None
                    and not self.aware
                    and _timezones_kind(self.tz_strat) == "none"
                ):
                    bounds = _narrowed_bounds(
                        func, arg, self.min_value, self.max_value, _shift_datetime
                    )
                    if bounds is None:
                        return nothing()
                    if bounds == (self.min_value, self.max_value):
                        return self
                    return datetimes(
                        *bounds,
                        timezones=self.tz_strat,
                        allow_imaginary=self.allow_imaginary,
                    )
            else:
                # An aware bound constrains the instant of generated values,
                # so we narrow our aware bounds to the closed interval of
                # satisfying instants - retaining strict predicates below,
                # which then reject at most the boundary instant per timezone.
                # We compare bounds by their _instant() key, since comparison
                # of datetimes which share a tzinfo would fall back to
                # wall-clock order, ignoring the fold.
                if self.aware:
                    min_value, max_value = self.min_value, self.max_value
                elif (self.min_value, self.max_value) == (
                    dt.datetime.min,
                    dt.datetime.max,
                ) and _timezones_kind(self.tz_strat) != "none":
                    # An unbounded naive-mode strategy whose values are all
                    # aware: promote to aware mode, bounded by the filter.
                    min_value = max_value = None
                else:
                    return super().filter(condition)
                key = _instant(arg)
                if func in (op.lt, op.le, op.eq) and (
                    min_value is None or _instant(min_value) < key
                ):
                    min_value = arg
                if func in (op.gt, op.ge, op.eq) and (
                    max_value is None or key < _instant(max_value)
                ):
                    max_value = arg
                if min_value is not None and max_value is not None:
                    lo, hi = _instant(min_value), _instant(max_value)
                    if hi < lo or (func in (op.lt, op.gt) and lo == hi == key):
                        # Only aware-mode strategies can reach this, and they
                        # generate only aware values (or raise for a bad
                        # timezones strategy), so this is provably empty.
                        return nothing()
                if min_value is self.min_value and max_value is self.max_value:
                    result = self
                else:
                    result = DatetimeStrategy(
                        min_value, max_value, self.tz_strat, self.allow_imaginary
                    )
                if func in (op.lt, op.gt):
                    return FilteredStrategy(result, (condition,))
                return result
        return super().filter(condition)


@overload
def datetimes(
    min_value: NaiveDatetime | None = None,
    max_value: NaiveDatetime | None = None,
    *,
    timezones: SearchStrategy[None] | None = None,
) -> SearchStrategy[NaiveDatetime]:  # pragma: no cover
    ...


@overload
def datetimes(
    min_value: dt.datetime | None = None,
    max_value: dt.datetime | None = None,
    *,
    timezones: SearchStrategy[dt.tzinfo],
    allow_imaginary: bool = True,
) -> SearchStrategy[AwareDatetime]:  # pragma: no cover
    ...


@overload
def datetimes(
    min_value: None = None,
    max_value: None = None,
    *,
    timezones: SearchStrategy[dt.tzinfo | None],
    allow_imaginary: bool = True,
) -> SearchStrategy[dt.datetime]:  # pragma: no cover
    ...


@defines_strategy(force_reusable_values=True)
def datetimes(
    min_value: dt.datetime | None = None,
    max_value: dt.datetime | None = None,
    *,
    timezones: SearchStrategy[dt.tzinfo | None] | None = None,
    allow_imaginary: bool = True,
) -> SearchStrategy[dt.datetime]:
    """datetimes(min_value=None, max_value=None, *, timezones=None, allow_imaginary=True)

    A strategy for generating datetimes, which may be timezone-aware.

    If ``min_value`` and ``max_value`` are naive datetimes, or omitted, this
    strategy works by drawing a naive datetime between them - defaulting to
    ``datetime.min`` and ``datetime.max`` respectively - and then attaching
    a timezone drawn from ``timezones``, which defaults to
    :func:`~hypothesis.strategies.none`.

    If instead both bounds are timezone-aware, they are treated as moments in
    time, and ``timezones`` defaults to :func:`~hypothesis.strategies.timezones`.
    Each generated datetime is aware, in a timezone drawn from ``timezones`` -
    which must not generate ``None`` - and lies between the two moments.
    Passing one aware and one naive bound is an error.

    ``timezones`` must be a strategy that generates either ``None``, for naive
    datetimes, or :class:`~python:datetime.tzinfo` objects for 'aware' datetimes.
    You can construct your own, though we recommend using one of these built-in
    strategies:

    * with the standard library: :func:`hypothesis.strategies.timezones`;
    * with :pypi:`dateutil <python-dateutil>`:
      :func:`hypothesis.extra.dateutil.timezones`; or
    * with :pypi:`pytz`: :func:`hypothesis.extra.pytz.timezones`.

    You may pass ``allow_imaginary=False`` to filter out "imaginary" datetimes
    which did not (or will not) occur due to daylight savings, leap seconds,
    timezone and calendar adjustments, etc.  Imaginary datetimes are allowed
    by default, because malformed timestamps are a common source of bugs.

    .. note::

        Arithmetic and comparisons on timezone-aware datetimes can be very
        surprising around daylight-savings changes.  See `this CPython issue
        <https://github.com/python/cpython/issues/116035>`__ for details
        and discussion.

    Examples from this strategy shrink towards midnight on January 1st 2000,
    local time.
    """
    check_type(bool, allow_imaginary, "allow_imaginary")
    if min_value is not None:
        check_type(dt.datetime, min_value, "min_value")
    if max_value is not None:
        check_type(dt.datetime, max_value, "max_value")
    if timezones is not None and not isinstance(timezones, SearchStrategy):
        raise InvalidArgument(
            f"{timezones=} must be a SearchStrategy that can "
            "provide tzinfo for datetimes (either None or dt.tzinfo objects)"
        )
    if (min_value is None or min_value.tzinfo is None) and (
        max_value is None or max_value.tzinfo is None
    ):
        min_value = dt.datetime.min if min_value is None else min_value
        max_value = dt.datetime.max if max_value is None else max_value
        if timezones is None:
            timezones = none()
        check_valid_interval(min_value, max_value, "min_value", "max_value")
    else:
        # Aware bounds describe moments in time; we check both are aware here,
        # and then at draw time convert them to the drawn timezone and proceed
        # as in the naive case.
        for name, value in [("min_value", min_value), ("max_value", max_value)]:
            if value is not None and value.tzinfo is None:
                raise InvalidArgument(
                    f"{name}={value!r} is naive, but the other bound is "
                    "timezone-aware; the bounds must be both naive or both aware"
                )
        if timezones is None:
            timezones = _timezones()
        # Compare explicitly as moments in time: comparison of datetimes which
        # share a tzinfo falls back to wall-clock order, ignoring the fold.
        if (
            min_value is not None
            and max_value is not None
            and _instant(max_value) < _instant(min_value)
        ):
            raise InvalidArgument(
                f"Cannot have {max_value=} < {min_value=}, comparing as "
                "moments in time"
            )
    return DatetimeStrategy(min_value, max_value, timezones, allow_imaginary)


_ARBITRARY_DATE = dt.date(2000, 1, 1)


def _shift_time(value, steps):
    # dt.time supports no arithmetic, so we go via a datetime on a fixed day
    # and treat crossing midnight as overflowing the representable range.
    shifted = dt.datetime.combine(_ARBITRARY_DATE, value) + steps * _MICROSECOND
    if shifted.date() != _ARBITRARY_DATE:
        raise OverflowError
    return shifted.time()


class TimeStrategy(SearchStrategy):
    def __init__(self, min_value, max_value, timezones_strat):
        super().__init__()
        self.min_value = min_value
        self.max_value = max_value
        self.tz_strat = timezones_strat

    def do_draw(self, data):
        result = draw_capped_multipart(data, self.min_value, self.max_value, TIMENAMES)
        tz = data.draw(self.tz_strat)
        return dt.time(**result, tzinfo=tz)

    def filter(self, condition):
        # We only rewrite naive times: ordering aware times works in terms of
        # utcoffset(), which is None for e.g. ZoneInfo tzinfos on a time - so
        # such values compare as naive anyway, and rewriting fixed-offset aware
        # times isn't worth the extra complexity.
        if (
            (parsed := _comparator_bound(condition)) is not None
            and isinstance(arg := parsed[1], dt.time)
            and arg.tzinfo is None
            and _timezones_kind(self.tz_strat) == "none"
        ):
            bounds = _narrowed_bounds(
                parsed[0], arg, self.min_value, self.max_value, _shift_time
            )
            if bounds is None:
                return nothing()
            if bounds == (self.min_value, self.max_value):
                return self
            return times(*bounds, timezones=self.tz_strat)
        return super().filter(condition)


@defines_strategy(force_reusable_values=True)
def times(
    min_value: dt.time = dt.time.min,
    max_value: dt.time = dt.time.max,
    *,
    timezones: SearchStrategy[dt.tzinfo | None] = none(),
) -> SearchStrategy[dt.time]:
    """times(min_value=datetime.time.min, max_value=datetime.time.max, *, timezones=none())

    A strategy for times between ``min_value`` and ``max_value``.

    The ``timezones`` argument is handled as for :py:func:`datetimes`.

    Examples from this strategy shrink towards midnight, with the timezone
    component shrinking as for the strategy that provided it.
    """
    check_type(dt.time, min_value, "min_value")
    check_type(dt.time, max_value, "max_value")
    if min_value.tzinfo is not None:
        raise InvalidArgument(f"{min_value=} must not have tzinfo")
    if max_value.tzinfo is not None:
        raise InvalidArgument(f"{max_value=} must not have tzinfo")
    check_valid_interval(min_value, max_value, "min_value", "max_value")
    return TimeStrategy(min_value, max_value, timezones)


def _shift_date(value, steps):
    return value + steps * dt.timedelta(days=1)


class DateStrategy(SearchStrategy):
    def __init__(self, min_value, max_value):
        super().__init__()
        assert isinstance(min_value, dt.date)
        assert isinstance(max_value, dt.date)
        assert min_value < max_value
        self.min_value = min_value
        self.max_value = max_value

    def do_draw(self, data):
        return dt.date(
            **draw_capped_multipart(data, self.min_value, self.max_value, DATENAMES)
        )

    def filter(self, condition):
        if (
            (parsed := _comparator_bound(condition)) is not None
            # datetime is a date subclass, but not comparable with dates
            and isinstance(arg := parsed[1], dt.date)
            and not isinstance(arg, dt.datetime)
        ):
            bounds = _narrowed_bounds(
                parsed[0], arg, self.min_value, self.max_value, _shift_date
            )
            if bounds is None:
                return nothing()
            if bounds == (self.min_value, self.max_value):
                return self
            return dates(*bounds)

        return super().filter(condition)


@defines_strategy(force_reusable_values=True)
def dates(
    min_value: dt.date = dt.date.min, max_value: dt.date = dt.date.max
) -> SearchStrategy[dt.date]:
    """dates(min_value=datetime.date.min, max_value=datetime.date.max)

    A strategy for dates between ``min_value`` and ``max_value``.

    Examples from this strategy shrink towards January 1st 2000.
    """
    check_type(dt.date, min_value, "min_value")
    check_type(dt.date, max_value, "max_value")
    # datetime is a subclass of date, so check_type() accepts it - but a datetime
    # bound is almost certainly a mistake, and breaks our drawing logic downstream.
    if isinstance(min_value, dt.datetime):
        raise InvalidArgument(f"{min_value=} is a datetime, but expected a date")
    if isinstance(max_value, dt.datetime):
        raise InvalidArgument(f"{max_value=} is a datetime, but expected a date")
    check_valid_interval(min_value, max_value, "min_value", "max_value")
    if min_value == max_value:
        return just(min_value)
    return DateStrategy(min_value, max_value)


class TimedeltaStrategy(SearchStrategy):
    def __init__(self, min_value, max_value):
        super().__init__()
        assert isinstance(min_value, dt.timedelta)
        assert isinstance(max_value, dt.timedelta)
        assert min_value < max_value
        self.min_value = min_value
        self.max_value = max_value

    def do_draw(self, data):
        result = {}
        low_bound = True
        high_bound = True
        for name in ("days", "seconds", "microseconds"):
            low = getattr(self.min_value if low_bound else dt.timedelta.min, name)
            high = getattr(self.max_value if high_bound else dt.timedelta.max, name)
            val = data.draw_integer(low, high)
            result[name] = val
            low_bound = low_bound and val == low
            high_bound = high_bound and val == high
        return dt.timedelta(**result)


@defines_strategy(force_reusable_values=True)
def timedeltas(
    min_value: dt.timedelta = dt.timedelta.min,
    max_value: dt.timedelta = dt.timedelta.max,
) -> SearchStrategy[dt.timedelta]:
    """timedeltas(min_value=datetime.timedelta.min, max_value=datetime.timedelta.max)

    A strategy for timedeltas between ``min_value`` and ``max_value``.

    Examples from this strategy shrink towards zero.
    """
    check_type(dt.timedelta, min_value, "min_value")
    check_type(dt.timedelta, max_value, "max_value")
    check_valid_interval(min_value, max_value, "min_value", "max_value")
    if min_value == max_value:
        return just(min_value)
    return TimedeltaStrategy(min_value=min_value, max_value=max_value)


@cache
def _valid_key_cacheable(tzpath, key):
    assert isinstance(tzpath, tuple)  # zoneinfo changed, better update this function!
    for root in tzpath:
        if Path(root).joinpath(key).exists():  # pragma: no branch
            # No branch because most systems only have one TZPATH component.
            return True
    else:  # pragma: no cover
        # This branch is only taken for names which are known to zoneinfo
        # but not present on the filesystem, i.e. on Windows with tzdata,
        # and so is never executed by our coverage tests.
        *package_loc, resource_name = key.split("/")
        package = "tzdata.zoneinfo." + ".".join(package_loc)
        try:
            return (resources.files(package) / resource_name).exists()
        except ModuleNotFoundError:
            return False


@defines_strategy(force_reusable_values=True)
def timezone_keys(
    *,
    # allow_alias: bool = True,
    # allow_deprecated: bool = True,
    allow_prefix: bool = True,
) -> SearchStrategy[str]:
    """A strategy for :wikipedia:`IANA timezone names <List_of_tz_database_time_zones>`.

    As well as timezone names like ``"UTC"``, ``"Australia/Sydney"``, or
    ``"America/New_York"``, this strategy can generate:

    - Aliases such as ``"Antarctica/McMurdo"``, which links to ``"Pacific/Auckland"``.
    - Deprecated names such as ``"Antarctica/South_Pole"``, which *also* links to
      ``"Pacific/Auckland"``.  Note that most but
      not all deprecated timezone names are also aliases.
    - Timezone names with the ``"posix/"`` or ``"right/"`` prefixes, unless
      ``allow_prefix=False``.

    These strings are provided separately from Tzinfo objects - such as ZoneInfo
    instances from the timezones() strategy - to facilitate testing of timezone
    logic without needing workarounds to access non-canonical names.

    .. note::

        `The tzdata package is required on Windows
        <https://docs.python.org/3/library/zoneinfo.html#data-sources>`__.
        ``pip install hypothesis[zoneinfo]`` installs it, if and only if needed.

    On Windows, you may need to access IANA timezone data via the :pypi:`tzdata`
    package.  For non-IANA timezones, such as Windows-native names or GNU TZ
    strings, we recommend using :func:`~hypothesis.strategies.sampled_from` with
    the :pypi:`dateutil <python-dateutil>` package, e.g.
    :meth:`dateutil:dateutil.tz.tzwin.list`.
    """
    # check_type(bool, allow_alias, "allow_alias")
    # check_type(bool, allow_deprecated, "allow_deprecated")
    check_type(bool, allow_prefix, "allow_prefix")

    with warnings.catch_warnings():
        try:
            warnings.simplefilter("ignore", EncodingWarning)
        except NameError:  # pragma: no cover
            pass
        # On Python 3.12 (and others?), `available_timezones()` opens files
        # without specifying an encoding - which our selftests make an error.
        available_timezones = ("UTC", *sorted(zoneinfo.available_timezones()))

    # TODO: filter out alias and deprecated names if disallowed

    # When prefixes are allowed, we first choose a key and then flatmap to get our
    # choice with one of the available prefixes.  That in turn means that we need
    # some logic to determine which prefixes are available for a given key:

    def valid_key(key):
        return key == "UTC" or _valid_key_cacheable(zoneinfo.TZPATH, key)

    # TODO: work out how to place a higher priority on "weird" timezones
    # For details see https://github.com/HypothesisWorks/hypothesis/issues/2414
    strategy = sampled_from([key for key in available_timezones if valid_key(key)])

    if not allow_prefix:
        return strategy

    def sample_with_prefixes(zone):
        keys_with_prefixes = (zone, f"posix/{zone}", f"right/{zone}")
        return sampled_from([key for key in keys_with_prefixes if valid_key(key)])

    return strategy.flatmap(sample_with_prefixes)


@defines_strategy(force_reusable_values=True)
def timezones(*, no_cache: bool = False) -> SearchStrategy["zoneinfo.ZoneInfo"]:
    """A strategy for :class:`python:zoneinfo.ZoneInfo` objects.

    If ``no_cache=True``, the generated instances are constructed using
    :meth:`ZoneInfo.no_cache <python:zoneinfo.ZoneInfo.no_cache>` instead
    of the usual constructor.  This may change the semantics of your datetimes
    in surprising ways, so only use it if you know that you need to!

    .. note::

        `The tzdata package is required on Windows
        <https://docs.python.org/3/library/zoneinfo.html#data-sources>`__.
        ``pip install hypothesis[zoneinfo]`` installs it, if and only if needed.
    """
    check_type(bool, no_cache, "no_cache")
    return timezone_keys().map(
        zoneinfo.ZoneInfo.no_cache if no_cache else zoneinfo.ZoneInfo
    )


# In datetimes() above, the ``timezones`` argument shadows this module's
# timezones() strategy, so we refer to it by this alias instead.
_timezones = timezones
