"""
Implementation of operations on numpy timedelta64.
"""

import numpy as np
import operator

import llvmlite.ir
from llvmlite.ir import Constant

from numba.core import types, cgutils
from numba.core.cgutils import create_constant_array
from numba.core.imputils import (lower_builtin, lower_constant,
                                 impl_ret_untracked, lower_cast)
from numba.np import npdatetime_helpers, numpy_support, npyfuncs
from numba.extending import overload_method
from numba.core.config import IS_32BITS
from numba.core.errors import LoweringError

# datetime64 and timedelta64 use the same internal representation
DATETIME64 = TIMEDELTA64 = llvmlite.ir.IntType(64)
NAT = Constant(TIMEDELTA64, npdatetime_helpers.NAT)

TIMEDELTA_BINOP_SIG = (types.NPTimedelta,) * 2


def scale_by_constant(builder, val, factor):
    """
    Multiply *val* by the constant *factor*.
    """
    return builder.mul(val, Constant(TIMEDELTA64, factor))


def unscale_by_constant(builder, val, factor):
    """
    Divide *val* by the constant *factor*.
    """
    return builder.sdiv(val, Constant(TIMEDELTA64, factor))


def add_constant(builder, val, const):
    """
    Add constant *const* to *val*.
    """
    return builder.add(val, Constant(TIMEDELTA64, const))


def scale_timedelta(context, builder, val, srcty, destty):
    """
    Scale the timedelta64 *val* from *srcty* to *destty*
    (both numba.types.NPTimedelta instances)
    """
    factor = npdatetime_helpers.get_timedelta_conversion_factor(
        srcty.unit, destty.unit)
    if factor is None:
        # This can happen when using explicit output in a ufunc.
        msg = f"cannot convert timedelta64 from {srcty.unit} to {destty.unit}"
        raise LoweringError(msg)
    return scale_by_constant(builder, val, factor)


def normalize_timedeltas(context, builder, left, right, leftty, rightty):
    """
    Scale either *left* or *right* to the other's unit, in order to have
    homogeneous units.
    """
    factor = npdatetime_helpers.get_timedelta_conversion_factor(
        leftty.unit, rightty.unit)
    if factor is not None:
        return scale_by_constant(builder, left, factor), right
    factor = npdatetime_helpers.get_timedelta_conversion_factor(
        rightty.unit, leftty.unit)
    if factor is not None:
        return left, scale_by_constant(builder, right, factor)
    # Typing should not let this happen, except on == and != operators
    raise RuntimeError("cannot normalize %r and %r" % (leftty, rightty))


def alloc_timedelta_result(builder, name='ret'):
    """
    Allocate a NaT-initialized datetime64 (or timedelta64) result slot.
    """
    ret = cgutils.alloca_once(builder, TIMEDELTA64, name=name)
    builder.store(NAT, ret)
    return ret


def alloc_boolean_result(builder, name='ret'):
    """
    Allocate an uninitialized boolean result slot.
    """
    ret = cgutils.alloca_once(builder, llvmlite.ir.IntType(1), name=name)
    return ret


def is_not_nat(builder, val):
    """
    Return a predicate which is true if *val* is not NaT.
    """
    return builder.icmp_unsigned('!=', val, NAT)


def are_not_nat(builder, vals):
    """
    Return a predicate which is true if all of *vals* are not NaT.
    """
    assert len(vals) >= 1
    pred = is_not_nat(builder, vals[0])
    for val in vals[1:]:
        pred = builder.and_(pred, is_not_nat(builder, val))
    return pred


normal_year_months = create_constant_array(
    TIMEDELTA64,
    [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
leap_year_months = create_constant_array(
    TIMEDELTA64,
    [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
normal_year_months_acc = create_constant_array(
    TIMEDELTA64,
    [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334])
leap_year_months_acc = create_constant_array(
    TIMEDELTA64,
    [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335])


@lower_constant(types.NPDatetime)
@lower_constant(types.NPTimedelta)
def datetime_constant(context, builder, ty, pyval):
    return DATETIME64(pyval.astype(np.int64))


# Arithmetic operators on timedelta64

@lower_builtin(operator.pos, types.NPTimedelta)
def timedelta_pos_impl(context, builder, sig, args):
    res = args[0]
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(operator.neg, types.NPTimedelta)
def timedelta_neg_impl(context, builder, sig, args):
    res = builder.neg(args[0])
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(abs, types.NPTimedelta)
def timedelta_abs_impl(context, builder, sig, args):
    val, = args
    ret = alloc_timedelta_result(builder)
    with builder.if_else(cgutils.is_scalar_neg(builder, val)) as (then, otherwise):
        with then:
            builder.store(builder.neg(val), ret)
        with otherwise:
            builder.store(val, ret)
    res = builder.load(ret)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def timedelta_sign_impl(context, builder, sig, args):
    """
    np.sign(timedelta64)
    """
    val, = args
    ret = alloc_timedelta_result(builder)
    zero = Constant(TIMEDELTA64, 0)
    with builder.if_else(builder.icmp_signed('>', val, zero)
                         ) as (gt_zero, le_zero):
        with gt_zero:
            builder.store(Constant(TIMEDELTA64, 1), ret)
        with le_zero:
            with builder.if_else(builder.icmp_unsigned('==', val, zero)
                                 ) as (eq_zero, lt_zero):
                with eq_zero:
                    builder.store(Constant(TIMEDELTA64, 0), ret)
                with lt_zero:
                    builder.store(Constant(TIMEDELTA64, -1), ret)
    res = builder.load(ret)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(operator.add, *TIMEDELTA_BINOP_SIG)
@lower_builtin(operator.iadd, *TIMEDELTA_BINOP_SIG)
def timedelta_add_impl(context, builder, sig, args):
    [va, vb] = args
    [ta, tb] = sig.args
    ret = alloc_timedelta_result(builder)
    with cgutils.if_likely(builder, are_not_nat(builder, [va, vb])):
        va = scale_timedelta(context, builder, va, ta, sig.return_type)
        vb = scale_timedelta(context, builder, vb, tb, sig.return_type)
        builder.store(builder.add(va, vb), ret)
    res = builder.load(ret)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(operator.sub, *TIMEDELTA_BINOP_SIG)
@lower_builtin(operator.isub, *TIMEDELTA_BINOP_SIG)
def timedelta_sub_impl(context, builder, sig, args):
    [va, vb] = args
    [ta, tb] = sig.args
    ret = alloc_timedelta_result(builder)
    with cgutils.if_likely(builder, are_not_nat(builder, [va, vb])):
        va = scale_timedelta(context, builder, va, ta, sig.return_type)
        vb = scale_timedelta(context, builder, vb, tb, sig.return_type)
        builder.store(builder.sub(va, vb), ret)
    res = builder.load(ret)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def _timedelta_times_number(context, builder, td_arg, td_type,
                            number_arg, number_type, return_type):
    ret = alloc_timedelta_result(builder)
    with cgutils.if_likely(builder, is_not_nat(builder, td_arg)):
        if isinstance(number_type, types.Float):
            val = builder.sitofp(td_arg, number_arg.type)
            val = builder.fmul(val, number_arg)
            val = _cast_to_timedelta(context, builder, val)
        else:
            val = builder.mul(td_arg, number_arg)
        # The scaling is required for ufunc np.multiply() with an explicit
        # output in a different unit.
        val = scale_timedelta(context, builder, val, td_type, return_type)
        builder.store(val, ret)
    return builder.load(ret)


@lower_builtin(operator.mul, types.NPTimedelta, types.Integer)
@lower_builtin(operator.imul, types.NPTimedelta, types.Integer)
@lower_builtin(operator.mul, types.NPTimedelta, types.Float)
@lower_builtin(operator.imul, types.NPTimedelta, types.Float)
def timedelta_times_number(context, builder, sig, args):
    res = _timedelta_times_number(context, builder,
                                  args[0], sig.args[0], args[1], sig.args[1],
                                  sig.return_type)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(operator.mul, types.Integer, types.NPTimedelta)
@lower_builtin(operator.imul, types.Integer, types.NPTimedelta)
@lower_builtin(operator.mul, types.Float, types.NPTimedelta)
@lower_builtin(operator.imul, types.Float, types.NPTimedelta)
def number_times_timedelta(context, builder, sig, args):
    res = _timedelta_times_number(context, builder,
                                  args[1], sig.args[1], args[0], sig.args[0],
                                  sig.return_type)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(operator.truediv, types.NPTimedelta, types.Integer)
@lower_builtin(operator.itruediv, types.NPTimedelta, types.Integer)
@lower_builtin(operator.floordiv, types.NPTimedelta, types.Integer)
@lower_builtin(operator.ifloordiv, types.NPTimedelta, types.Integer)
@lower_builtin(operator.truediv, types.NPTimedelta, types.Float)
@lower_builtin(operator.itruediv, types.NPTimedelta, types.Float)
@lower_builtin(operator.floordiv, types.NPTimedelta, types.Float)
@lower_builtin(operator.ifloordiv, types.NPTimedelta, types.Float)
def timedelta_over_number(context, builder, sig, args):
    td_arg, number_arg = args
    number_type = sig.args[1]
    ret = alloc_timedelta_result(builder)
    ok = builder.and_(is_not_nat(builder, td_arg),
                      builder.not_(cgutils.is_scalar_zero_or_nan(builder, number_arg)))
    with cgutils.if_likely(builder, ok):
        # Denominator is non-zero, non-NaN
        if isinstance(number_type, types.Float):
            val = builder.sitofp(td_arg, number_arg.type)
            val = builder.fdiv(val, number_arg)
            val = _cast_to_timedelta(context, builder, val)
        else:
            val = builder.sdiv(td_arg, number_arg)
        # The scaling is required for ufuncs np.*divide() with an explicit
        # output in a different unit.
        val = scale_timedelta(context, builder, val,
                              sig.args[0], sig.return_type)
        builder.store(val, ret)
    res = builder.load(ret)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(operator.truediv, *TIMEDELTA_BINOP_SIG)
@lower_builtin(operator.itruediv, *TIMEDELTA_BINOP_SIG)
def timedelta_over_timedelta(context, builder, sig, args):
    [va, vb] = args
    [ta, tb] = sig.args
    not_nan = are_not_nat(builder, [va, vb])
    ll_ret_type = context.get_value_type(sig.return_type)
    ret = cgutils.alloca_once(builder, ll_ret_type, name='ret')
    builder.store(Constant(ll_ret_type, float('nan')), ret)
    with cgutils.if_likely(builder, not_nan):
        va, vb = normalize_timedeltas(context, builder, va, vb, ta, tb)
        va = builder.sitofp(va, ll_ret_type)
        vb = builder.sitofp(vb, ll_ret_type)
        builder.store(builder.fdiv(va, vb), ret)
    res = builder.load(ret)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(operator.floordiv, *TIMEDELTA_BINOP_SIG)
def timedelta_floor_div_timedelta(context, builder, sig, args):
    [va, vb] = args
    [ta, tb] = sig.args
    ll_ret_type = context.get_value_type(sig.return_type)
    not_nan = are_not_nat(builder, [va, vb])
    ret = cgutils.alloca_once(builder, ll_ret_type, name='ret')
    zero = Constant(ll_ret_type, 0)
    one = Constant(ll_ret_type, 1)
    builder.store(zero, ret)
    with cgutils.if_likely(builder, not_nan):
        va, vb = normalize_timedeltas(context, builder, va, vb, ta, tb)
        # is the denominator zero or NaT?
        denom_ok = builder.not_(builder.icmp_signed('==', vb, zero))
        with cgutils.if_likely(builder, denom_ok):
            # is either arg negative?
            vaneg = builder.icmp_signed('<', va, zero)
            neg = builder.or_(vaneg, builder.icmp_signed('<', vb, zero))
            with builder.if_else(neg) as (then, otherwise):
                with then:  # one or more value negative
                    with builder.if_else(vaneg) as (negthen, negotherwise):
                        with negthen:
                            top = builder.sub(va, one)
                            div = builder.sdiv(top, vb)
                            builder.store(div, ret)
                        with negotherwise:
                            top = builder.add(va, one)
                            div = builder.sdiv(top, vb)
                            builder.store(div, ret)
                with otherwise:
                    div = builder.sdiv(va, vb)
                    builder.store(div, ret)
    res = builder.load(ret)
    return impl_ret_untracked(context, builder, sig.return_type, res)

def timedelta_mod_timedelta(context, builder, sig, args):
    # inspired by https://github.com/numpy/numpy/blob/fe8072a12d65e43bd2e0b0f9ad67ab0108cc54b3/numpy/core/src/umath/loops.c.src#L1424
    # alg is basically as `a % b`:
    # if a or b is NaT return NaT
    # elseif b is 0 return NaT
    # else pretend a and b are int and do pythonic int modulus

    [va, vb] = args
    [ta, tb] = sig.args
    not_nan = are_not_nat(builder, [va, vb])
    ll_ret_type = context.get_value_type(sig.return_type)
    ret = alloc_timedelta_result(builder)
    builder.store(NAT, ret)
    zero = Constant(ll_ret_type, 0)
    with cgutils.if_likely(builder, not_nan):
        va, vb = normalize_timedeltas(context, builder, va, vb, ta, tb)
        # is the denominator zero or NaT?
        denom_ok = builder.not_(builder.icmp_signed('==', vb, zero))
        with cgutils.if_likely(builder, denom_ok):
            # is either arg negative?
            vapos = builder.icmp_signed('>', va, zero)
            vbpos = builder.icmp_signed('>', vb, zero)
            rem = builder.srem(va, vb)
            cond = builder.or_(builder.and_(vapos, vbpos),
                               builder.icmp_signed('==', rem, zero))
            with builder.if_else(cond) as (then, otherwise):
                with then:
                    builder.store(rem, ret)
                with otherwise:
                    builder.store(builder.add(rem, vb), ret)

    res = builder.load(ret)
    return impl_ret_untracked(context, builder, sig.return_type, res)

# Comparison operators on timedelta64


def _create_timedelta_comparison_impl(ll_op, default_value):
    def impl(context, builder, sig, args):
        [va, vb] = args
        [ta, tb] = sig.args
        ret = alloc_boolean_result(builder)
        with builder.if_else(are_not_nat(builder, [va, vb])) as (then, otherwise):
            with then:
                try:
                    norm_a, norm_b = normalize_timedeltas(
                        context, builder, va, vb, ta, tb)
                except RuntimeError:
                    # Cannot normalize units => the values are unequal (except if NaT)
                    builder.store(default_value, ret)
                else:
                    builder.store(builder.icmp_unsigned(ll_op, norm_a, norm_b), ret)
            with otherwise:
                # NaT ==/>=/>/</<= NaT is False
                # NaT != <anything, including NaT> is True
                if ll_op == '!=':
                    builder.store(cgutils.true_bit, ret)
                else:
                    builder.store(cgutils.false_bit, ret)
        res = builder.load(ret)
        return impl_ret_untracked(context, builder, sig.return_type, res)

    return impl


def _create_timedelta_ordering_impl(ll_op):
    def impl(context, builder, sig, args):
        [va, vb] = args
        [ta, tb] = sig.args
        ret = alloc_boolean_result(builder)
        with builder.if_else(are_not_nat(builder, [va, vb])) as (then, otherwise):
            with then:
                norm_a, norm_b = normalize_timedeltas(
                    context, builder, va, vb, ta, tb)
                builder.store(builder.icmp_signed(ll_op, norm_a, norm_b), ret)
            with otherwise:
                # NaT >=/>/</<= NaT is False
                builder.store(cgutils.false_bit, ret)
        res = builder.load(ret)
        return impl_ret_untracked(context, builder, sig.return_type, res)

    return impl


timedelta_eq_timedelta_impl = _create_timedelta_comparison_impl(
    '==', cgutils.false_bit)
timedelta_ne_timedelta_impl = _create_timedelta_comparison_impl(
    '!=', cgutils.true_bit)
timedelta_lt_timedelta_impl = _create_timedelta_ordering_impl('<')
timedelta_le_timedelta_impl = _create_timedelta_ordering_impl('<=')
timedelta_gt_timedelta_impl = _create_timedelta_ordering_impl('>')
timedelta_ge_timedelta_impl = _create_timedelta_ordering_impl('>=')

for op_, func in [(operator.eq, timedelta_eq_timedelta_impl),
                  (operator.ne, timedelta_ne_timedelta_impl),
                  (operator.lt,  timedelta_lt_timedelta_impl),
                  (operator.le, timedelta_le_timedelta_impl),
                  (operator.gt,  timedelta_gt_timedelta_impl),
                  (operator.ge, timedelta_ge_timedelta_impl)]:
    lower_builtin(op_, *TIMEDELTA_BINOP_SIG)(func)


# Arithmetic on datetime64

def is_leap_year(builder, year_val):
    """
    Return a predicate indicating whether *year_val* (offset by 1970) is a
    leap year.
    """
    actual_year = builder.add(year_val, Constant(DATETIME64, 1970))
    multiple_of_4 = cgutils.is_null(
        builder, builder.and_(actual_year, Constant(DATETIME64, 3)))
    not_multiple_of_100 = cgutils.is_not_null(
        builder, builder.srem(actual_year, Constant(DATETIME64, 100)))
    multiple_of_400 = cgutils.is_null(
        builder, builder.srem(actual_year, Constant(DATETIME64, 400)))
    return builder.and_(multiple_of_4,
                        builder.or_(not_multiple_of_100, multiple_of_400))


def year_to_days(builder, year_val):
    """
    Given a year *year_val* (offset to 1970), return the number of days
    since the 1970 epoch.
    """
    # The algorithm below is copied from Numpy's get_datetimestruct_days()
    # (src/multiarray/datetime.c)
    ret = cgutils.alloca_once(builder, TIMEDELTA64)
    # First approximation
    days = scale_by_constant(builder, year_val, 365)
    # Adjust for leap years
    with builder.if_else(cgutils.is_neg_int(builder, year_val)) \
            as (if_neg, if_pos):
        with if_pos:
            # At or after 1970:
            # 1968 is the closest leap year before 1970.
            # Exclude the current year, so add 1.
            from_1968 = add_constant(builder, year_val, 1)
            # Add one day for each 4 years
            p_days = builder.add(days,
                                 unscale_by_constant(builder, from_1968, 4))
            # 1900 is the closest previous year divisible by 100
            from_1900 = add_constant(builder, from_1968, 68)
            # Subtract one day for each 100 years
            p_days = builder.sub(p_days,
                                 unscale_by_constant(builder, from_1900, 100))
            # 1600 is the closest previous year divisible by 400
            from_1600 = add_constant(builder, from_1900, 300)
            # Add one day for each 400 years
            p_days = builder.add(p_days,
                                 unscale_by_constant(builder, from_1600, 400))
            builder.store(p_days, ret)
        with if_neg:
            # Before 1970:
            # NOTE `year_val` is negative, and so will be `from_1972` and `from_2000`.
            # 1972 is the closest later year after 1970.
            # Include the current year, so subtract 2.
            from_1972 = add_constant(builder, year_val, -2)
            # Subtract one day for each 4 years (`from_1972` is negative)
            n_days = builder.add(days,
                                 unscale_by_constant(builder, from_1972, 4))
            # 2000 is the closest later year divisible by 100
            from_2000 = add_constant(builder, from_1972, -28)
            # Add one day for each 100 years
            n_days = builder.sub(n_days,
                                 unscale_by_constant(builder, from_2000, 100))
            # 2000 is also the closest later year divisible by 400
            # Subtract one day for each 400 years
            n_days = builder.add(n_days,
                                 unscale_by_constant(builder, from_2000, 400))
            builder.store(n_days, ret)
    return builder.load(ret)


def reduce_datetime_for_unit(builder, dt_val, src_unit, dest_unit):
    dest_unit_code = npdatetime_helpers.DATETIME_UNITS[dest_unit]
    src_unit_code = npdatetime_helpers.DATETIME_UNITS[src_unit]
    if dest_unit_code < 2 or src_unit_code >= 2:
        return dt_val, src_unit
    # Need to compute the day ordinal for *dt_val*
    if src_unit_code == 0:
        # Years to days
        year_val = dt_val
        days_val = year_to_days(builder, year_val)

    else:
        # Months to days
        leap_array = cgutils.global_constant(builder, "leap_year_months_acc",
                                             leap_year_months_acc)
        normal_array = cgutils.global_constant(builder, "normal_year_months_acc",
                                               normal_year_months_acc)

        days = cgutils.alloca_once(builder, TIMEDELTA64)

        # First compute year number and month number
        year, month = cgutils.divmod_by_constant(builder, dt_val, 12)

        # Then deduce the number of days
        with builder.if_else(is_leap_year(builder, year)) as (then, otherwise):
            with then:
                addend = builder.load(cgutils.gep(builder, leap_array,
                                                  0, month, inbounds=True))
                builder.store(addend, days)
            with otherwise:
                addend = builder.load(cgutils.gep(builder, normal_array,
                                                  0, month, inbounds=True))
                builder.store(addend, days)

        days_val = year_to_days(builder, year)
        days_val = builder.add(days_val, builder.load(days))

    if dest_unit_code == 2:
        # Need to scale back to weeks
        weeks, _ = cgutils.divmod_by_constant(builder, days_val, 7)
        return weeks, 'W'
    else:
        return days_val, 'D'


def convert_datetime_for_arith(builder, dt_val, src_unit, dest_unit):
    """
    Convert datetime *dt_val* from *src_unit* to *dest_unit*.
    """
    # First partial conversion to days or weeks, if necessary.
    dt_val, dt_unit = reduce_datetime_for_unit(
        builder, dt_val, src_unit, dest_unit)
    # Then multiply by the remaining constant factor.
    dt_factor = npdatetime_helpers.get_timedelta_conversion_factor(dt_unit, dest_unit)
    if dt_factor is None:
        # This can happen when using explicit output in a ufunc.
        raise LoweringError("cannot convert datetime64 from %r to %r"
                                       % (src_unit, dest_unit))
    return scale_by_constant(builder, dt_val, dt_factor)


def _datetime_timedelta_arith(ll_op_name):
    def impl(context, builder, dt_arg, dt_unit,
             td_arg, td_unit, ret_unit):
        ret = alloc_timedelta_result(builder)
        with cgutils.if_likely(builder, are_not_nat(builder, [dt_arg, td_arg])):
            dt_arg = convert_datetime_for_arith(builder, dt_arg,
                                                dt_unit, ret_unit)
            td_factor = npdatetime_helpers.get_timedelta_conversion_factor(
                td_unit, ret_unit)
            td_arg = scale_by_constant(builder, td_arg, td_factor)
            ret_val = getattr(builder, ll_op_name)(dt_arg, td_arg)
            builder.store(ret_val, ret)
        return builder.load(ret)
    return impl


_datetime_plus_timedelta = _datetime_timedelta_arith('add')
_datetime_minus_timedelta = _datetime_timedelta_arith('sub')

# datetime64 + timedelta64


@lower_builtin(operator.add, types.NPDatetime, types.NPTimedelta)
@lower_builtin(operator.iadd, types.NPDatetime, types.NPTimedelta)
def datetime_plus_timedelta(context, builder, sig, args):
    dt_arg, td_arg = args
    dt_type, td_type = sig.args
    res = _datetime_plus_timedelta(context, builder,
                                   dt_arg, dt_type.unit,
                                   td_arg, td_type.unit,
                                   sig.return_type.unit)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(operator.add, types.NPTimedelta, types.NPDatetime)
@lower_builtin(operator.iadd, types.NPTimedelta, types.NPDatetime)
def timedelta_plus_datetime(context, builder, sig, args):
    td_arg, dt_arg = args
    td_type, dt_type = sig.args
    res = _datetime_plus_timedelta(context, builder,
                                   dt_arg, dt_type.unit,
                                   td_arg, td_type.unit,
                                   sig.return_type.unit)
    return impl_ret_untracked(context, builder, sig.return_type, res)

# datetime64 - timedelta64


@lower_builtin(operator.sub, types.NPDatetime, types.NPTimedelta)
@lower_builtin(operator.isub, types.NPDatetime, types.NPTimedelta)
def datetime_minus_timedelta(context, builder, sig, args):
    dt_arg, td_arg = args
    dt_type, td_type = sig.args
    res = _datetime_minus_timedelta(context, builder,
                                    dt_arg, dt_type.unit,
                                    td_arg, td_type.unit,
                                    sig.return_type.unit)
    return impl_ret_untracked(context, builder, sig.return_type, res)

# datetime64 - datetime64


@lower_builtin(operator.sub, types.NPDatetime, types.NPDatetime)
def datetime_minus_datetime(context, builder, sig, args):
    va, vb = args
    ta, tb = sig.args
    unit_a = ta.unit
    unit_b = tb.unit
    ret_unit = sig.return_type.unit
    ret = alloc_timedelta_result(builder)
    with cgutils.if_likely(builder, are_not_nat(builder, [va, vb])):
        va = convert_datetime_for_arith(builder, va, unit_a, ret_unit)
        vb = convert_datetime_for_arith(builder, vb, unit_b, ret_unit)
        ret_val = builder.sub(va, vb)
        builder.store(ret_val, ret)
    res = builder.load(ret)
    return impl_ret_untracked(context, builder, sig.return_type, res)

# datetime64 comparisons


def _create_datetime_comparison_impl(ll_op):
    def impl(context, builder, sig, args):
        va, vb = args
        ta, tb = sig.args
        unit_a = ta.unit
        unit_b = tb.unit
        ret_unit = npdatetime_helpers.get_best_unit(unit_a, unit_b)
        ret = alloc_boolean_result(builder)
        with builder.if_else(are_not_nat(builder, [va, vb])) as (then, otherwise):
            with then:
                norm_a = convert_datetime_for_arith(
                    builder, va, unit_a, ret_unit)
                norm_b = convert_datetime_for_arith(
                    builder, vb, unit_b, ret_unit)
                ret_val = builder.icmp_signed(ll_op, norm_a, norm_b)
                builder.store(ret_val, ret)
            with otherwise:
                if ll_op == '!=':
                    ret_val = cgutils.true_bit
                else:
                    ret_val = cgutils.false_bit
                builder.store(ret_val, ret)
        res = builder.load(ret)
        return impl_ret_untracked(context, builder, sig.return_type, res)

    return impl


datetime_eq_datetime_impl = _create_datetime_comparison_impl('==')
datetime_ne_datetime_impl = _create_datetime_comparison_impl('!=')
datetime_lt_datetime_impl = _create_datetime_comparison_impl('<')
datetime_le_datetime_impl = _create_datetime_comparison_impl('<=')
datetime_gt_datetime_impl = _create_datetime_comparison_impl('>')
datetime_ge_datetime_impl = _create_datetime_comparison_impl('>=')

for op, func in [(operator.eq, datetime_eq_datetime_impl),
                 (operator.ne, datetime_ne_datetime_impl),
                 (operator.lt, datetime_lt_datetime_impl),
                 (operator.le, datetime_le_datetime_impl),
                 (operator.gt, datetime_gt_datetime_impl),
                 (operator.ge, datetime_ge_datetime_impl)]:
    lower_builtin(op, *[types.NPDatetime]*2)(func)


########################################################################
# datetime/timedelta fmax/fmin maximum/minimum support

def _gen_datetime_max_impl(NAT_DOMINATES):
    def datetime_max_impl(context, builder, sig, args):
        # note this could be optimizing relying on the actual value of NAT
        # but as NumPy doesn't rely on this, this seems more resilient
        in1, in2 = args
        in1_not_nat = is_not_nat(builder, in1)
        in2_not_nat = is_not_nat(builder, in2)
        in1_ge_in2 = builder.icmp_signed('>=', in1, in2)
        res = builder.select(in1_ge_in2, in1, in2)
        if NAT_DOMINATES:
            # NaT now dominates, like NaN
            in1, in2 = in2, in1
        res = builder.select(in1_not_nat, res, in2)
        res = builder.select(in2_not_nat, res, in1)

        return impl_ret_untracked(context, builder, sig.return_type, res)
    return datetime_max_impl

datetime_maximum_impl = _gen_datetime_max_impl(True)
datetime_fmax_impl = _gen_datetime_max_impl(False)

def _gen_datetime_min_impl(NAT_DOMINATES):
    def datetime_min_impl(context, builder, sig, args):
        # note this could be optimizing relying on the actual value of NAT
        # but as NumPy doesn't rely on this, this seems more resilient
        in1, in2 = args
        in1_not_nat = is_not_nat(builder, in1)
        in2_not_nat = is_not_nat(builder, in2)
        in1_le_in2 = builder.icmp_signed('<=', in1, in2)
        res = builder.select(in1_le_in2, in1, in2)
        if NAT_DOMINATES:
            # NaT now dominates, like NaN
            in1, in2 = in2, in1
        res = builder.select(in1_not_nat, res, in2)
        res = builder.select(in2_not_nat, res, in1)

        return impl_ret_untracked(context, builder, sig.return_type, res)
    return datetime_min_impl

datetime_minimum_impl = _gen_datetime_min_impl(True)
datetime_fmin_impl = _gen_datetime_min_impl(False)

def _gen_timedelta_max_impl(NAT_DOMINATES):
    def timedelta_max_impl(context, builder, sig, args):
        # note this could be optimizing relying on the actual value of NAT
        # but as NumPy doesn't rely on this, this seems more resilient
        in1, in2 = args
        in1_not_nat = is_not_nat(builder, in1)
        in2_not_nat = is_not_nat(builder, in2)
        in1_ge_in2 = builder.icmp_signed('>=', in1, in2)
        res = builder.select(in1_ge_in2, in1, in2)
        if NAT_DOMINATES:
            # NaT now dominates, like NaN
            in1, in2 = in2, in1
        res = builder.select(in1_not_nat, res, in2)
        res = builder.select(in2_not_nat, res, in1)

        return impl_ret_untracked(context, builder, sig.return_type, res)
    return timedelta_max_impl

timedelta_maximum_impl = _gen_timedelta_max_impl(True)
timedelta_fmax_impl = _gen_timedelta_max_impl(False)

def _gen_timedelta_min_impl(NAT_DOMINATES):
    def timedelta_min_impl(context, builder, sig, args):
        # note this could be optimizing relying on the actual value of NAT
        # but as NumPy doesn't rely on this, this seems more resilient
        in1, in2 = args
        in1_not_nat = is_not_nat(builder, in1)
        in2_not_nat = is_not_nat(builder, in2)
        in1_le_in2 = builder.icmp_signed('<=', in1, in2)
        res = builder.select(in1_le_in2, in1, in2)
        if NAT_DOMINATES:
            # NaT now dominates, like NaN
            in1, in2 = in2, in1
        res = builder.select(in1_not_nat, res, in2)
        res = builder.select(in2_not_nat, res, in1)

        return impl_ret_untracked(context, builder, sig.return_type, res)
    return timedelta_min_impl

timedelta_minimum_impl = _gen_timedelta_min_impl(True)
timedelta_fmin_impl = _gen_timedelta_min_impl(False)

def _cast_to_timedelta(context, builder, val):
    temp = builder.alloca(TIMEDELTA64)
    val_is_nan = builder.fcmp_unordered('uno', val, val)
    with builder.if_else(val_is_nan) as (
            then, els):
        with then:
            # NaN does not guarantee to cast to NAT.
            # We should store NAT explicitly.
            builder.store(NAT, temp)
        with els:
            builder.store(builder.fptosi(val, TIMEDELTA64), temp)
    return builder.load(temp)


@lower_builtin(np.isnat, types.NPDatetime)
@lower_builtin(np.isnat, types.NPTimedelta)
def _np_isnat_impl(context, builder, sig, args):
    return npyfuncs.np_datetime_isnat_impl(context, builder, sig, args)


@lower_cast(types.NPDatetime, types.Integer)
@lower_cast(types.NPTimedelta, types.Integer)
def _cast_npdatetime_int64(context, builder, fromty, toty, val):
    if toty.bitwidth != 64: # all date time types are 64 bit
        msg = f"Cannot cast {fromty} to {toty} as {toty} is not 64 bits wide."
        raise ValueError(msg)
    return val


@overload_method(types.NPTimedelta, '__hash__')
@overload_method(types.NPDatetime, '__hash__')
def ol_hash_npdatetime(x):
    if numpy_support.numpy_version >= (2, 2)\
          and isinstance(x, types.NPTimedelta) and not x.unit:
        raise ValueError("Can't hash generic timedelta64")

    if IS_32BITS:
        def impl(x):
            x = np.int64(x)
            if x < 2**31 - 1:  # x < LONG_MAX
                y = np.int32(x)
            else:
                hi = (np.int64(x) & 0xffffffff00000000) >> 32
                lo = (np.int64(x) & 0x00000000ffffffff)
                y = np.int32(lo + (1000003) * hi)
            if y == -1:
                y = np.int32(-2)
            return y
    else:
        def impl(x):
            if np.int64(x) == -1:
                return np.int64(-2)
            return np.int64(x)
    return impl


lower_builtin(npdatetime_helpers.datetime_minimum, types.NPDatetime, types.NPDatetime)(datetime_minimum_impl)
lower_builtin(npdatetime_helpers.datetime_minimum, types.NPTimedelta, types.NPTimedelta)(datetime_minimum_impl)
lower_builtin(npdatetime_helpers.datetime_maximum, types.NPDatetime, types.NPDatetime)(datetime_maximum_impl)
lower_builtin(npdatetime_helpers.datetime_maximum, types.NPTimedelta, types.NPTimedelta)(datetime_maximum_impl)
