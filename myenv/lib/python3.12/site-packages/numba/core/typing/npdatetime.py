"""
Typing declarations for np.timedelta64.
"""


from itertools import product
import operator

from numba.core import types
from numba.core.typing.templates import (AttributeTemplate, ConcreteTemplate,
                                         AbstractTemplate, infer_global, infer,
                                         infer_getattr, signature)
from numba.np import npdatetime_helpers
from numba.np.numpy_support import numpy_version


# timedelta64-only operations

class TimedeltaUnaryOp(AbstractTemplate):

    def generic(self, args, kws):
        if len(args) == 2:
            # Guard against binary + and -
            return
        op, = args
        if not isinstance(op, types.NPTimedelta):
            return
        return signature(op, op)


class TimedeltaBinOp(AbstractTemplate):

    def generic(self, args, kws):
        if len(args) == 1:
            # Guard against unary + and -
            return
        left, right = args
        if not all(isinstance(tp, types.NPTimedelta) for tp in args):
            return
        if npdatetime_helpers.can_cast_timedelta_units(left.unit, right.unit):
            return signature(right, left, right)
        elif npdatetime_helpers.can_cast_timedelta_units(right.unit, left.unit):
            return signature(left, left, right)


class TimedeltaCmpOp(AbstractTemplate):

    def generic(self, args, kws):
        # For equality comparisons, all units are inter-comparable
        left, right = args
        if not all(isinstance(tp, types.NPTimedelta) for tp in args):
            return
        return signature(types.boolean, left, right)


class TimedeltaOrderedCmpOp(AbstractTemplate):

    def generic(self, args, kws):
        # For ordered comparisons, units must be compatible
        left, right = args
        if not all(isinstance(tp, types.NPTimedelta) for tp in args):
            return
        if (npdatetime_helpers.can_cast_timedelta_units(left.unit, right.unit) or
            npdatetime_helpers.can_cast_timedelta_units(right.unit, left.unit)):
            return signature(types.boolean, left, right)


class TimedeltaMixOp(AbstractTemplate):

    def generic(self, args, kws):
        """
        (timedelta64, {int, float}) -> timedelta64
        ({int, float}, timedelta64) -> timedelta64
        """
        left, right = args
        if isinstance(right, types.NPTimedelta):
            td, other = right, left
            sig_factory = lambda other: signature(td, other, td)
        elif isinstance(left, types.NPTimedelta):
            td, other = left, right
            sig_factory = lambda other: signature(td, td, other)
        else:
            return
        if not isinstance(other, (types.Float, types.Integer)):
            return
        # Force integer types to convert to signed because it matches
        # timedelta64 semantics better.
        if isinstance(other, types.Integer):
            other = types.int64
        return sig_factory(other)


class TimedeltaDivOp(AbstractTemplate):

    def generic(self, args, kws):
        """
        (timedelta64, {int, float}) -> timedelta64
        (timedelta64, timedelta64) -> float
        """
        left, right = args
        if not isinstance(left, types.NPTimedelta):
            return
        if isinstance(right, types.NPTimedelta):
            if (npdatetime_helpers.can_cast_timedelta_units(left.unit, right.unit)
                or npdatetime_helpers.can_cast_timedelta_units(right.unit, left.unit)):
                return signature(types.float64, left, right)
        elif isinstance(right, (types.Float)):
            return signature(left, left, right)
        elif isinstance(right, (types.Integer)):
            # Force integer types to convert to signed because it matches
            # timedelta64 semantics better.
            return signature(left, left, types.int64)


@infer_global(operator.pos)
class TimedeltaUnaryPos(TimedeltaUnaryOp):
    key = operator.pos

@infer_global(operator.neg)
class TimedeltaUnaryNeg(TimedeltaUnaryOp):
    key = operator.neg

@infer_global(operator.add)
@infer_global(operator.iadd)
class TimedeltaBinAdd(TimedeltaBinOp):
    key = operator.add

@infer_global(operator.sub)
@infer_global(operator.isub)
class TimedeltaBinSub(TimedeltaBinOp):
    key = operator.sub

@infer_global(operator.mul)
@infer_global(operator.imul)
class TimedeltaBinMult(TimedeltaMixOp):
    key = operator.mul

@infer_global(operator.truediv)
@infer_global(operator.itruediv)
class TimedeltaTrueDiv(TimedeltaDivOp):
    key = operator.truediv

@infer_global(operator.floordiv)
@infer_global(operator.ifloordiv)
class TimedeltaFloorDiv(TimedeltaDivOp):
    key = operator.floordiv

if numpy_version >= (1, 25):
    @infer_global(operator.eq)
    class TimedeltaCmpEq(TimedeltaOrderedCmpOp):
        key = operator.eq

    @infer_global(operator.ne)
    class TimedeltaCmpNe(TimedeltaOrderedCmpOp):
        key = operator.ne
else:
    @infer_global(operator.eq)
    class TimedeltaCmpEq(TimedeltaCmpOp):
        key = operator.eq

    @infer_global(operator.ne)
    class TimedeltaCmpNe(TimedeltaCmpOp):
        key = operator.ne

@infer_global(operator.lt)
class TimedeltaCmpLt(TimedeltaOrderedCmpOp):
    key = operator.lt

@infer_global(operator.le)
class TimedeltaCmpLE(TimedeltaOrderedCmpOp):
    key = operator.le

@infer_global(operator.gt)
class TimedeltaCmpGt(TimedeltaOrderedCmpOp):
    key = operator.gt

@infer_global(operator.ge)
class TimedeltaCmpGE(TimedeltaOrderedCmpOp):
    key = operator.ge


@infer_global(abs)
class TimedeltaAbs(TimedeltaUnaryOp):
    pass


# datetime64 operations

@infer_global(operator.add)
@infer_global(operator.iadd)
class DatetimePlusTimedelta(AbstractTemplate):
    key = operator.add

    def generic(self, args, kws):
        if len(args) == 1:
            # Guard against unary +
            return
        left, right = args
        if isinstance(right, types.NPTimedelta):
            dt = left
            td = right
        elif isinstance(left, types.NPTimedelta):
            dt = right
            td = left
        else:
            return
        if isinstance(dt, types.NPDatetime):
            unit = npdatetime_helpers.combine_datetime_timedelta_units(dt.unit,
                                                                       td.unit)
            if unit is not None:
                return signature(types.NPDatetime(unit), left, right)

@infer_global(operator.sub)
@infer_global(operator.isub)
class DatetimeMinusTimedelta(AbstractTemplate):
    key = operator.sub

    def generic(self, args, kws):
        if len(args) == 1:
            # Guard against unary -
            return
        dt, td = args
        if isinstance(dt, types.NPDatetime) and isinstance(td,
                                                           types.NPTimedelta):
            unit = npdatetime_helpers.combine_datetime_timedelta_units(dt.unit,
                                                                       td.unit)
            if unit is not None:
                return signature(types.NPDatetime(unit), dt, td)

@infer_global(operator.sub)
class DatetimeMinusDatetime(AbstractTemplate):
    key = operator.sub

    def generic(self, args, kws):
        if len(args) == 1:
            # Guard against unary -
            return
        left, right = args
        if isinstance(left, types.NPDatetime) and isinstance(right,
                                                             types.NPDatetime):
            unit = npdatetime_helpers.get_best_unit(left.unit, right.unit)
            return signature(types.NPTimedelta(unit), left, right)


class DatetimeCmpOp(AbstractTemplate):

    def generic(self, args, kws):
        # For datetime64 comparisons, all units are inter-comparable
        left, right = args
        if not all(isinstance(tp, types.NPDatetime) for tp in args):
            return
        return signature(types.boolean, left, right)


@infer_global(operator.eq)
class DatetimeCmpEq(DatetimeCmpOp):
    key = operator.eq

@infer_global(operator.ne)
class DatetimeCmpNe(DatetimeCmpOp):
    key = operator.ne

@infer_global(operator.lt)
class DatetimeCmpLt(DatetimeCmpOp):
    key = operator.lt

@infer_global(operator.le)
class DatetimeCmpLE(DatetimeCmpOp):
    key = operator.le

@infer_global(operator.gt)
class DatetimeCmpGt(DatetimeCmpOp):
    key = operator.gt

@infer_global(operator.ge)
class DatetimeCmpGE(DatetimeCmpOp):
    key = operator.ge


@infer_global(npdatetime_helpers.datetime_minimum)
@infer_global(npdatetime_helpers.datetime_maximum)
class DatetimeMinMax(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        assert len(args) == 2
        error_msg = "DatetimeMinMax requires both arguments to be NPDatetime type or both arguments to be NPTimedelta types"
        assert isinstance(args[0], (types.NPDatetime, types.NPTimedelta)), error_msg
        if isinstance(args[0], types.NPDatetime):
            if not isinstance(args[1], types.NPDatetime):
                raise TypeError(error_msg)
        else:
            if not isinstance(args[1], types.NPTimedelta):
                raise TypeError(error_msg)
        return signature(args[0], *args)
