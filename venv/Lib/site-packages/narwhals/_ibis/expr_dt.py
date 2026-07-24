from __future__ import annotations

from typing import TYPE_CHECKING, Any

import ibis

from narwhals._duration import Interval
from narwhals._ibis.utils import UNITS_DICT_BUCKET, UNITS_DICT_TRUNCATE
from narwhals._sql.expr_dt import SQLExprDateTimeNamesSpace
from narwhals._utils import not_implemented

if TYPE_CHECKING:
    from collections.abc import Callable

    import ibis.expr.types as ir

    from narwhals._ibis.expr import IbisExpr
    from narwhals._ibis.utils import BucketUnit, TruncateUnit


class IbisExprDateTimeNamespace(SQLExprDateTimeNamesSpace["IbisExpr"]):
    def millisecond(self) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.millisecond())

    def microsecond(self) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.microsecond())

    def to_string(self, format: str) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.strftime(format))

    def weekday(self) -> IbisExpr:
        # Ibis uses 0-6 for Monday-Sunday. Add 1 to match polars.
        return self.compliant._with_callable(lambda expr: expr.day_of_week.index() + 1)

    def _bucket(self, kwds: dict[BucketUnit, Any], /) -> Callable[..., ir.TimestampValue]:
        def fn(expr: ir.TimestampValue) -> ir.TimestampValue:
            return expr.bucket(**kwds)

        return fn

    def _truncate(self, unit: TruncateUnit, /) -> Callable[..., ir.TimestampValue]:
        def fn(expr: ir.TimestampValue) -> ir.TimestampValue:
            return expr.truncate(unit)

        return fn

    def truncate(self, every: str) -> IbisExpr:
        interval = Interval.parse(every)
        multiple, unit = interval.multiple, interval.unit
        if unit == "q":
            multiple, unit = 3 * multiple, "mo"
        if multiple != 1:
            if self.compliant._backend_version < (7, 1):  # pragma: no cover
                msg = "Truncating datetimes with multiples of the unit is only supported in Ibis >= 7.1."
                raise NotImplementedError(msg)
            fn = self._bucket({UNITS_DICT_BUCKET[unit]: multiple})
        else:
            fn = self._truncate(UNITS_DICT_TRUNCATE[unit])
        return self.compliant._with_callable(fn)

    def offset_by(self, by: str) -> IbisExpr:
        interval = Interval.parse_no_constraints(by)
        multiple, unit = interval.multiple, interval.unit
        if unit == "ns":
            msg = "Offsetting by nanoseconds is not yet supported for ibis."
            raise NotImplementedError(msg)
        # `ibis.interval` does not accept `quarters` on all supported ibis
        # versions (`truncate` avoids passing it through for the same reason),
        # so express a quarter offset as the equivalent number of months.
        offset_multiple, offset_unit = multiple, unit
        if unit == "q":
            offset_multiple, offset_unit = 3 * multiple, "mo"
        kwds: dict[BucketUnit, Any] = {UNITS_DICT_BUCKET[offset_unit]: offset_multiple}
        offset = ibis.interval(**kwds)

        def fn(expr: ir.TimestampValue) -> ir.TimestampValue:
            # Ibis stores timezone-aware data as UTC, so calendar offsets are
            # applied to the underlying instant rather than wall-clock time.
            # The result would differ from other backends across a DST
            # transition, so refuse it rather than return a surprising value.
            tz = getattr(expr.type(), "timezone", None)
            if unit in {"y", "q", "mo", "d"} and tz is not None:
                msg = (
                    f"Offsetting timezone-aware data by {UNITS_DICT_BUCKET[unit]} is "
                    "not supported for ibis, as the result would be incorrect across "
                    "daylight saving time transitions."
                )
                raise NotImplementedError(msg)
            return expr.add(offset)

        return self.compliant._with_callable(fn)

    def replace_time_zone(self, time_zone: str | None) -> IbisExpr:
        if time_zone is None:
            return self.compliant._with_callable(lambda expr: expr.cast("timestamp"))
        msg = "`replace_time_zone` with non-null `time_zone` not yet implemented for Ibis"  # pragma: no cover
        raise NotImplementedError(msg)

    nanosecond = not_implemented()
    total_minutes = not_implemented()
    total_seconds = not_implemented()
    total_milliseconds = not_implemented()
    total_microseconds = not_implemented()
    total_nanoseconds = not_implemented()
    convert_time_zone = not_implemented()
    timestamp = not_implemented()
