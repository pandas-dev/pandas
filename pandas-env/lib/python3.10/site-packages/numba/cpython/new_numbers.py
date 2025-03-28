import math
import numbers
import itertools

import numpy as np
import operator

from llvmlite import ir
from llvmlite.ir import Constant

from numba.core.imputils import (lower_builtin, lower_getattr,
                                 lower_cast, lower_constant,
                                 impl_ret_untracked)
from numba.core import typing, types, errors, cgutils
from numba.core.extending import overload_method
from numba.cpython.unsafe.numbers import viewer

def _int_arith_flags(rettype):
    """
    Return the modifier flags for integer arithmetic.
    """
    if rettype.signed:
        # Ignore the effects of signed overflow.  This is important for
        # optimization of some indexing operations.  For example
        # array[i+1] could see `i+1` trigger a signed overflow and
        # give a negative number.  With Python's indexing, a negative
        # index is treated differently: its resolution has a runtime cost.
        # Telling LLVM to ignore signed overflows allows it to optimize
        # away the check for a negative `i+1` if it knows `i` is positive.
        return ['nsw']
    else:
        return []


def int_add_impl(context, builder, sig, args):
    [va, vb] = args
    [ta, tb] = sig.args
    a = context.cast(builder, va, ta, sig.return_type)
    b = context.cast(builder, vb, tb, sig.return_type)
    res = builder.add(a, b, flags=_int_arith_flags(sig.return_type))
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_sub_impl(context, builder, sig, args):
    [va, vb] = args
    [ta, tb] = sig.args
    a = context.cast(builder, va, ta, sig.return_type)
    b = context.cast(builder, vb, tb, sig.return_type)
    res = builder.sub(a, b, flags=_int_arith_flags(sig.return_type))
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_mul_impl(context, builder, sig, args):
    [va, vb] = args
    [ta, tb] = sig.args
    a = context.cast(builder, va, ta, sig.return_type)
    b = context.cast(builder, vb, tb, sig.return_type)
    res = builder.mul(a, b, flags=_int_arith_flags(sig.return_type))
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_divmod_signed(context, builder, ty, x, y):
    """
    Reference Objects/intobject.c
    xdivy = x / y;
    xmody = (long)(x - (unsigned long)xdivy * y);
    /* If the signs of x and y differ, and the remainder is non-0,
     * C89 doesn't define whether xdivy is now the floor or the
     * ceiling of the infinitely precise quotient.  We want the floor,
     * and we have it iff the remainder's sign matches y's.
     */
    if (xmody && ((y ^ xmody) < 0) /* i.e. and signs differ */) {
        xmody += y;
        --xdivy;
        assert(xmody && ((y ^ xmody) >= 0));
    }
    *p_xdivy = xdivy;
    *p_xmody = xmody;
    """
    assert x.type == y.type

    ZERO = y.type(0)
    ONE = y.type(1)

    # NOTE: On x86 at least, dividing the lowest representable integer
    # (e.g. 0x80000000 for int32) by -1 causes a SIFGPE (division overflow),
    # causing the process to crash.
    # We return 0, 0 instead (more or less like Numpy).

    resdiv = cgutils.alloca_once_value(builder, ZERO)
    resmod = cgutils.alloca_once_value(builder, ZERO)

    is_overflow = builder.and_(
        builder.icmp_signed('==', x, x.type(ty.minval)),
        builder.icmp_signed('==', y, y.type(-1)))

    with builder.if_then(builder.not_(is_overflow), likely=True):
        # Note LLVM will optimize this to a single divmod instruction,
        # if available on the target CPU (e.g. x86).
        xdivy = builder.sdiv(x, y)
        xmody = builder.srem(x, y)

        y_xor_xmody_ltz = builder.icmp_signed('<', builder.xor(y, xmody), ZERO)
        xmody_istrue = builder.icmp_signed('!=', xmody, ZERO)
        cond = builder.and_(xmody_istrue, y_xor_xmody_ltz)

        with builder.if_else(cond) as (if_different_signs, if_same_signs):
            with if_same_signs:
                builder.store(xdivy, resdiv)
                builder.store(xmody, resmod)

            with if_different_signs:
                builder.store(builder.sub(xdivy, ONE), resdiv)
                builder.store(builder.add(xmody, y), resmod)

    return builder.load(resdiv), builder.load(resmod)


def int_divmod(context, builder, ty, x, y):
    """
    Integer divmod(x, y).  The caller must ensure that y != 0.
    """
    if ty.signed:
        return int_divmod_signed(context, builder, ty, x, y)
    else:
        return builder.udiv(x, y), builder.urem(x, y)


def _int_divmod_impl(context, builder, sig, args, zerodiv_message):
    va, vb = args
    ta, tb = sig.args

    ty = sig.return_type
    if isinstance(ty, types.UniTuple):
        ty = ty.dtype
    a = context.cast(builder, va, ta, ty)
    b = context.cast(builder, vb, tb, ty)
    quot = cgutils.alloca_once(builder, a.type, name="quot")
    rem = cgutils.alloca_once(builder, a.type, name="rem")

    with builder.if_else(cgutils.is_scalar_zero(builder, b), likely=False
                         ) as (if_zero, if_non_zero):
        with if_zero:
            if not context.error_model.fp_zero_division(
                builder, (zerodiv_message,)):
                # No exception raised => return 0
                # XXX We should also set the FPU exception status, but
                # there's no easy way to do that from LLVM.
                builder.store(b, quot)
                builder.store(b, rem)
        with if_non_zero:
            q, r = int_divmod(context, builder, ty, a, b)
            builder.store(q, quot)
            builder.store(r, rem)

    return quot, rem


@lower_builtin(divmod, types.Integer, types.Integer)
def int_divmod_impl(context, builder, sig, args):
    quot, rem = _int_divmod_impl(context, builder, sig, args,
                                 "integer divmod by zero")

    return cgutils.pack_array(builder,
                              (builder.load(quot), builder.load(rem)))


@lower_builtin(operator.floordiv, types.Integer, types.Integer)
@lower_builtin(operator.ifloordiv, types.Integer, types.Integer)
def int_floordiv_impl(context, builder, sig, args):
    quot, rem = _int_divmod_impl(context, builder, sig, args,
                                 "integer division by zero")
    return builder.load(quot)


@lower_builtin(operator.truediv, types.Integer, types.Integer)
@lower_builtin(operator.itruediv, types.Integer, types.Integer)
def int_truediv_impl(context, builder, sig, args):
    [va, vb] = args
    [ta, tb] = sig.args
    a = context.cast(builder, va, ta, sig.return_type)
    b = context.cast(builder, vb, tb, sig.return_type)
    with cgutils.if_zero(builder, b):
        context.error_model.fp_zero_division(builder, ("division by zero",))
    res = builder.fdiv(a, b)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(operator.mod, types.Integer, types.Integer)
@lower_builtin(operator.imod, types.Integer, types.Integer)
def int_rem_impl(context, builder, sig, args):
    quot, rem = _int_divmod_impl(context, builder, sig, args,
                                 "integer modulo by zero")
    return builder.load(rem)


def _get_power_zerodiv_return(context, return_type):
    if (isinstance(return_type, types.Integer)
        and not context.error_model.raise_on_fp_zero_division):
        # If not raising, return 0x8000... when computing 0 ** <negative number>
        return -1 << (return_type.bitwidth - 1)
    else:
        return False


def int_power_impl(context, builder, sig, args):
    """
    a ^ b, where a is an integer or real, and b an integer
    """
    is_integer = isinstance(sig.args[0], types.Integer)
    tp = sig.return_type
    zerodiv_return = _get_power_zerodiv_return(context, tp)

    def int_power(a, b):
        # Ensure computations are done with a large enough width
        r = tp(1)
        a = tp(a)
        if b < 0:
            invert = True
            exp = -b
            if exp < 0:
                raise OverflowError
            if is_integer:
                if a == 0:
                    if zerodiv_return:
                        return zerodiv_return
                    else:
                        raise ZeroDivisionError("0 cannot be raised to a negative power")
                if a != 1 and a != -1:
                    return 0
        else:
            invert = False
            exp = b
        if exp > 0x10000:
            # Optimization cutoff: fallback on the generic algorithm
            return math.pow(a, float(b))
        while exp != 0:
            if exp & 1:
                r *= a
            exp >>= 1
            a *= a

        return 1.0 / r if invert else r

    res = context.compile_internal(builder, int_power, sig, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(operator.pow, types.Integer, types.IntegerLiteral)
@lower_builtin(operator.ipow, types.Integer, types.IntegerLiteral)
@lower_builtin(operator.pow, types.Float, types.IntegerLiteral)
@lower_builtin(operator.ipow, types.Float, types.IntegerLiteral)
def static_power_impl(context, builder, sig, args):
    """
    a ^ b, where a is an integer or real, and b a constant integer
    """
    exp = sig.args[1].value
    if not isinstance(exp, numbers.Integral):
        raise NotImplementedError
    if abs(exp) > 0x10000:
        # Optimization cutoff: fallback on the generic algorithm above
        raise NotImplementedError
    invert = exp < 0
    exp = abs(exp)

    tp = sig.return_type
    is_integer = isinstance(tp, types.Integer)
    zerodiv_return = _get_power_zerodiv_return(context, tp)

    val = context.cast(builder, args[0], sig.args[0], tp)
    lty = val.type

    def mul(a, b):
        if is_integer:
            return builder.mul(a, b)
        else:
            return builder.fmul(a, b)

    # Unroll the exponentiation loop
    res = lty(1)
    a = val
    while exp != 0:
        if exp & 1:
            res = mul(res, val)
        exp >>= 1
        val = mul(val, val)

    if invert:
        # If the exponent was negative, fix the result by inverting it
        if is_integer:
            # Integer inversion
            def invert_impl(a):
                if a == 0:
                    if zerodiv_return:
                        return zerodiv_return
                    else:
                        raise ZeroDivisionError("0 cannot be raised to a negative power")
                if a != 1 and a != -1:
                    return 0
                else:
                    return a

        else:
            # Real inversion
            def invert_impl(a):
                return 1.0 / a

        res = context.compile_internal(builder, invert_impl,
                                       typing.signature(tp, tp), (res,))

    return res


def int_slt_impl(context, builder, sig, args):
    res = builder.icmp_signed('<', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_sle_impl(context, builder, sig, args):
    res = builder.icmp_signed('<=', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_sgt_impl(context, builder, sig, args):
    res = builder.icmp_signed('>', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_sge_impl(context, builder, sig, args):
    res = builder.icmp_signed('>=', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_ult_impl(context, builder, sig, args):
    res = builder.icmp_unsigned('<', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_ule_impl(context, builder, sig, args):
    res = builder.icmp_unsigned('<=', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_ugt_impl(context, builder, sig, args):
    res = builder.icmp_unsigned('>', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_uge_impl(context, builder, sig, args):
    res = builder.icmp_unsigned('>=', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_eq_impl(context, builder, sig, args):
    res = builder.icmp_unsigned('==', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_ne_impl(context, builder, sig, args):
    res = builder.icmp_unsigned('!=', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_signed_unsigned_cmp(op):
    def impl(context, builder, sig, args):
        (left, right) = args
        # This code is translated from the NumPy source.
        # What we're going to do is divide the range of a signed value at zero.
        # If the signed value is less than zero, then we can treat zero as the
        # unsigned value since the unsigned value is necessarily zero or larger
        # and any signed comparison between a negative value and zero/infinity
        # will yield the same result. If the signed value is greater than or
        # equal to zero, then we can safely cast it to an unsigned value and do
        # the expected unsigned-unsigned comparison operation.
        # Original: https://github.com/numpy/numpy/pull/23713
        cmp_zero = builder.icmp_signed('<', left, Constant(left.type, 0))
        lt_zero = builder.icmp_signed(op, left, Constant(left.type, 0))
        ge_zero = builder.icmp_unsigned(op, left, right)
        res = builder.select(cmp_zero, lt_zero, ge_zero)
        return impl_ret_untracked(context, builder, sig.return_type, res)
    return impl


def int_unsigned_signed_cmp(op):
    def impl(context, builder, sig, args):
        (left, right) = args
        # See the function `int_signed_unsigned_cmp` for implementation notes.
        cmp_zero = builder.icmp_signed('<', right, Constant(right.type, 0))
        lt_zero = builder.icmp_signed(op, Constant(right.type, 0), right)
        ge_zero = builder.icmp_unsigned(op, left, right)
        res = builder.select(cmp_zero, lt_zero, ge_zero)
        return impl_ret_untracked(context, builder, sig.return_type, res)
    return impl


def int_abs_impl(context, builder, sig, args):
    [x] = args
    ZERO = Constant(x.type, None)
    ltz = builder.icmp_signed('<', x, ZERO)
    negated = builder.neg(x)
    res = builder.select(ltz, negated, x)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def uint_abs_impl(context, builder, sig, args):
    [x] = args
    return impl_ret_untracked(context, builder, sig.return_type, x)


def int_shl_impl(context, builder, sig, args):
    [valty, amtty] = sig.args
    [val, amt] = args
    val = context.cast(builder, val, valty, sig.return_type)
    amt = context.cast(builder, amt, amtty, sig.return_type)
    res = builder.shl(val, amt)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_shr_impl(context, builder, sig, args):
    [valty, amtty] = sig.args
    [val, amt] = args
    val = context.cast(builder, val, valty, sig.return_type)
    amt = context.cast(builder, amt, amtty, sig.return_type)
    if sig.return_type.signed:
        res = builder.ashr(val, amt)
    else:
        res = builder.lshr(val, amt)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_and_impl(context, builder, sig, args):
    [at, bt] = sig.args
    [av, bv] = args
    cav = context.cast(builder, av, at, sig.return_type)
    cbc = context.cast(builder, bv, bt, sig.return_type)
    res = builder.and_(cav, cbc)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_or_impl(context, builder, sig, args):
    [at, bt] = sig.args
    [av, bv] = args
    cav = context.cast(builder, av, at, sig.return_type)
    cbc = context.cast(builder, bv, bt, sig.return_type)
    res = builder.or_(cav, cbc)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_xor_impl(context, builder, sig, args):
    [at, bt] = sig.args
    [av, bv] = args
    cav = context.cast(builder, av, at, sig.return_type)
    cbc = context.cast(builder, bv, bt, sig.return_type)
    res = builder.xor(cav, cbc)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_negate_impl(context, builder, sig, args):
    [typ] = sig.args
    [val] = args
    # Negate before upcasting, for unsigned numbers
    res = builder.neg(val)
    res = context.cast(builder, res, typ, sig.return_type)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_positive_impl(context, builder, sig, args):
    [typ] = sig.args
    [val] = args
    res = context.cast(builder, val, typ, sig.return_type)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_invert_impl(context, builder, sig, args):
    [typ] = sig.args
    [val] = args
    # Invert before upcasting, for unsigned numbers
    res = builder.xor(val, Constant(val.type, int('1' * val.type.width, 2)))
    res = context.cast(builder, res, typ, sig.return_type)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def int_sign_impl(context, builder, sig, args):
    """
    np.sign(int)
    """
    [x] = args
    POS = Constant(x.type, 1)
    NEG = Constant(x.type, -1)
    ZERO = Constant(x.type, 0)

    cmp_zero = builder.icmp_unsigned('==', x, ZERO)
    cmp_pos = builder.icmp_signed('>', x, ZERO)

    presult = cgutils.alloca_once(builder, x.type)

    bb_zero = builder.append_basic_block(".zero")
    bb_postest = builder.append_basic_block(".postest")
    bb_pos = builder.append_basic_block(".pos")
    bb_neg = builder.append_basic_block(".neg")
    bb_exit = builder.append_basic_block(".exit")

    builder.cbranch(cmp_zero, bb_zero, bb_postest)

    with builder.goto_block(bb_zero):
        builder.store(ZERO, presult)
        builder.branch(bb_exit)

    with builder.goto_block(bb_postest):
        builder.cbranch(cmp_pos, bb_pos, bb_neg)

    with builder.goto_block(bb_pos):
        builder.store(POS, presult)
        builder.branch(bb_exit)

    with builder.goto_block(bb_neg):
        builder.store(NEG, presult)
        builder.branch(bb_exit)

    builder.position_at_end(bb_exit)
    res = builder.load(presult)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def bool_negate_impl(context, builder, sig, args):
    [typ] = sig.args
    [val] = args
    res = context.cast(builder, val, typ, sig.return_type)
    res = builder.neg(res)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def bool_unary_positive_impl(context, builder, sig, args):
    [typ] = sig.args
    [val] = args
    res = context.cast(builder, val, typ, sig.return_type)
    return impl_ret_untracked(context, builder, sig.return_type, res)


lower_builtin(operator.eq, types.py_bool, types.py_bool)(int_eq_impl)
lower_builtin(operator.ne, types.py_bool, types.py_bool)(int_ne_impl)
lower_builtin(operator.lt, types.py_bool, types.py_bool)(int_ult_impl)
lower_builtin(operator.le, types.py_bool, types.py_bool)(int_ule_impl)
lower_builtin(operator.gt, types.py_bool, types.py_bool)(int_ugt_impl)
lower_builtin(operator.ge, types.py_bool, types.py_bool)(int_uge_impl)
lower_builtin(operator.neg, types.py_bool)(bool_negate_impl)
lower_builtin(operator.pos, types.py_bool)(bool_unary_positive_impl)


def _implement_integer_operators():
    ty = types.Integer

    lower_builtin(operator.add, ty, ty)(int_add_impl)
    lower_builtin(operator.iadd, ty, ty)(int_add_impl)
    lower_builtin(operator.sub, ty, ty)(int_sub_impl)
    lower_builtin(operator.isub, ty, ty)(int_sub_impl)
    lower_builtin(operator.mul, ty, ty)(int_mul_impl)
    lower_builtin(operator.imul, ty, ty)(int_mul_impl)
    lower_builtin(operator.eq, ty, ty)(int_eq_impl)
    lower_builtin(operator.ne, ty, ty)(int_ne_impl)

    lower_builtin(operator.lshift, ty, ty)(int_shl_impl)
    lower_builtin(operator.ilshift, ty, ty)(int_shl_impl)
    lower_builtin(operator.rshift, ty, ty)(int_shr_impl)
    lower_builtin(operator.irshift, ty, ty)(int_shr_impl)

    lower_builtin(operator.neg, ty)(int_negate_impl)
    lower_builtin(operator.pos, ty)(int_positive_impl)

    lower_builtin(operator.pow, ty, ty)(int_power_impl)
    lower_builtin(operator.ipow, ty, ty)(int_power_impl)
    lower_builtin(pow, ty, ty)(int_power_impl)

    for ty in types.np_unsigned_domain:
        lower_builtin(operator.lt, ty, ty)(int_ult_impl)
        lower_builtin(operator.le, ty, ty)(int_ule_impl)
        lower_builtin(operator.gt, ty, ty)(int_ugt_impl)
        lower_builtin(operator.ge, ty, ty)(int_uge_impl)
        lower_builtin(operator.pow, types.Float, ty)(int_power_impl)
        lower_builtin(operator.ipow, types.Float, ty)(int_power_impl)
        lower_builtin(pow, types.Float, ty)(int_power_impl)
        lower_builtin(abs, ty)(uint_abs_impl)

    lower_builtin(operator.lt, types.IntegerLiteral, types.IntegerLiteral)(int_slt_impl)
    lower_builtin(operator.gt, types.IntegerLiteral, types.IntegerLiteral)(int_slt_impl)
    lower_builtin(operator.le, types.IntegerLiteral, types.IntegerLiteral)(int_slt_impl)
    lower_builtin(operator.ge, types.IntegerLiteral, types.IntegerLiteral)(int_slt_impl)
    for ty in itertools.chain(types.py_signed_domain,
                              types.np_signed_domain):
        lower_builtin(operator.lt, ty, ty)(int_slt_impl)
        lower_builtin(operator.le, ty, ty)(int_sle_impl)
        lower_builtin(operator.gt, ty, ty)(int_sgt_impl)
        lower_builtin(operator.ge, ty, ty)(int_sge_impl)
        lower_builtin(operator.pow, types.Float, ty)(int_power_impl)
        lower_builtin(operator.ipow, types.Float, ty)(int_power_impl)
        lower_builtin(pow, types.Float, ty)(int_power_impl)
        lower_builtin(abs, ty)(int_abs_impl)

def _implement_bitwise_operators():
    for ty in (types.Boolean, types.Integer):
        lower_builtin(operator.and_, ty, ty)(int_and_impl)
        lower_builtin(operator.iand, ty, ty)(int_and_impl)
        lower_builtin(operator.or_, ty, ty)(int_or_impl)
        lower_builtin(operator.ior, ty, ty)(int_or_impl)
        lower_builtin(operator.xor, ty, ty)(int_xor_impl)
        lower_builtin(operator.ixor, ty, ty)(int_xor_impl)

        lower_builtin(operator.invert, ty)(int_invert_impl)

_implement_integer_operators()

_implement_bitwise_operators()


def real_add_impl(context, builder, sig, args):
    res = builder.fadd(*args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_sub_impl(context, builder, sig, args):
    res = builder.fsub(*args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_mul_impl(context, builder, sig, args):
    res = builder.fmul(*args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_div_impl(context, builder, sig, args):
    with cgutils.if_zero(builder, args[1]):
        context.error_model.fp_zero_division(builder, ("division by zero",))
    res = builder.fdiv(*args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_divmod(context, builder, x, y):
    assert x.type == y.type
    floatty = x.type

    module = builder.module
    fname = context.mangler(".numba.python.rem", [x.type])
    fnty = ir.FunctionType(floatty, (floatty, floatty, ir.PointerType(floatty)))
    fn = cgutils.get_or_insert_function(module, fnty, fname)

    if fn.is_declaration:
        fn.linkage = 'linkonce_odr'
        fnbuilder = ir.IRBuilder(fn.append_basic_block('entry'))
        fx, fy, pmod = fn.args
        div, mod = real_divmod_func_body(context, fnbuilder, fx, fy)
        fnbuilder.store(mod, pmod)
        fnbuilder.ret(div)

    pmod = cgutils.alloca_once(builder, floatty)
    quotient = builder.call(fn, (x, y, pmod))
    return quotient, builder.load(pmod)


def real_divmod_func_body(context, builder, vx, wx):
    # Reference Objects/floatobject.c
    #
    # float_divmod(PyObject *v, PyObject *w)
    # {
    #     double vx, wx;
    #     double div, mod, floordiv;
    #     CONVERT_TO_DOUBLE(v, vx);
    #     CONVERT_TO_DOUBLE(w, wx);
    #     mod = fmod(vx, wx);
    #     /* fmod is typically exact, so vx-mod is *mathematically* an
    #        exact multiple of wx.  But this is fp arithmetic, and fp
    #        vx - mod is an approximation; the result is that div may
    #        not be an exact integral value after the division, although
    #        it will always be very close to one.
    #     */
    #     div = (vx - mod) / wx;
    #     if (mod) {
    #         /* ensure the remainder has the same sign as the denominator */
    #         if ((wx < 0) != (mod < 0)) {
    #             mod += wx;
    #             div -= 1.0;
    #         }
    #     }
    #     else {
    #         /* the remainder is zero, and in the presence of signed zeroes
    #            fmod returns different results across platforms; ensure
    #            it has the same sign as the denominator; we'd like to do
    #            "mod = wx * 0.0", but that may get optimized away */
    #         mod *= mod;  /* hide "mod = +0" from optimizer */
    #         if (wx < 0.0)
    #             mod = -mod;
    #     }
    #     /* snap quotient to nearest integral value */
    #     if (div) {
    #         floordiv = floor(div);
    #         if (div - floordiv > 0.5)
    #             floordiv += 1.0;
    #     }
    #     else {
    #         /* div is zero - get the same sign as the true quotient */
    #         div *= div;             /* hide "div = +0" from optimizers */
    #         floordiv = div * vx / wx; /* zero w/ sign of vx/wx */
    #     }
    #     return Py_BuildValue("(dd)", floordiv, mod);
    # }
    pmod = cgutils.alloca_once(builder, vx.type)
    pdiv = cgutils.alloca_once(builder, vx.type)
    pfloordiv = cgutils.alloca_once(builder, vx.type)

    mod = builder.frem(vx, wx)
    div = builder.fdiv(builder.fsub(vx, mod), wx)

    builder.store(mod, pmod)
    builder.store(div, pdiv)

    # Note the use of negative zero for proper negating with `ZERO - x`
    ZERO = vx.type(0.0)
    NZERO = vx.type(-0.0)
    ONE = vx.type(1.0)
    mod_istrue = builder.fcmp_unordered('!=', mod, ZERO)
    wx_ltz = builder.fcmp_ordered('<', wx, ZERO)
    mod_ltz = builder.fcmp_ordered('<', mod, ZERO)

    with builder.if_else(mod_istrue, likely=True) as (if_nonzero_mod, if_zero_mod):
        with if_nonzero_mod:
            # `mod` is non-zero or NaN
            # Ensure the remainder has the same sign as the denominator
            wx_ltz_ne_mod_ltz = builder.icmp_unsigned('!=', wx_ltz, mod_ltz)

            with builder.if_then(wx_ltz_ne_mod_ltz):
                builder.store(builder.fsub(div, ONE), pdiv)
                builder.store(builder.fadd(mod, wx), pmod)

        with if_zero_mod:
            # `mod` is zero, select the proper sign depending on
            # the denominator's sign
            mod = builder.select(wx_ltz, NZERO, ZERO)
            builder.store(mod, pmod)

    del mod, div

    div = builder.load(pdiv)
    div_istrue = builder.fcmp_ordered('!=', div, ZERO)

    with builder.if_then(div_istrue):
        realtypemap = {'float': types.float32,
                       'double': types.float64}
        realtype = realtypemap[str(wx.type)]
        floorfn = context.get_function(math.floor,
                                       typing.signature(realtype, realtype))
        floordiv = floorfn(builder, [div])
        floordivdiff = builder.fsub(div, floordiv)
        floordivincr = builder.fadd(floordiv, ONE)
        HALF = Constant(wx.type, 0.5)
        pred = builder.fcmp_ordered('>', floordivdiff, HALF)
        floordiv = builder.select(pred, floordivincr, floordiv)
        builder.store(floordiv, pfloordiv)

    with cgutils.ifnot(builder, div_istrue):
        div = builder.fmul(div, div)
        builder.store(div, pdiv)
        floordiv = builder.fdiv(builder.fmul(div, vx), wx)
        builder.store(floordiv, pfloordiv)

    return builder.load(pfloordiv), builder.load(pmod)


@lower_builtin(divmod, types.Float, types.Float)
def real_divmod_impl(context, builder, sig, args, loc=None):
    x, y = args
    quot = cgutils.alloca_once(builder, x.type, name="quot")
    rem = cgutils.alloca_once(builder, x.type, name="rem")

    with builder.if_else(cgutils.is_scalar_zero(builder, y), likely=False
                         ) as (if_zero, if_non_zero):
        with if_zero:
            if not context.error_model.fp_zero_division(
                builder, ("modulo by zero",), loc):
                # No exception raised => compute the nan result,
                # and set the FP exception word for Numpy warnings.
                q = builder.fdiv(x, y)
                r = builder.frem(x, y)
                builder.store(q, quot)
                builder.store(r, rem)
        with if_non_zero:
            q, r = real_divmod(context, builder, x, y)
            builder.store(q, quot)
            builder.store(r, rem)

    return cgutils.pack_array(builder,
                              (builder.load(quot), builder.load(rem)))


def real_mod_impl(context, builder, sig, args, loc=None):
    x, y = args
    res = cgutils.alloca_once(builder, x.type)
    with builder.if_else(cgutils.is_scalar_zero(builder, y), likely=False
                         ) as (if_zero, if_non_zero):
        with if_zero:
            if not context.error_model.fp_zero_division(
                builder, ("modulo by zero",), loc):
                # No exception raised => compute the nan result,
                # and set the FP exception word for Numpy warnings.
                rem = builder.frem(x, y)
                builder.store(rem, res)
        with if_non_zero:
            _, rem = real_divmod(context, builder, x, y)
            builder.store(rem, res)
    return impl_ret_untracked(context, builder, sig.return_type,
                              builder.load(res))


def real_floordiv_impl(context, builder, sig, args, loc=None):
    x, y = args
    res = cgutils.alloca_once(builder, x.type)
    with builder.if_else(cgutils.is_scalar_zero(builder, y), likely=False
                         ) as (if_zero, if_non_zero):
        with if_zero:
            if not context.error_model.fp_zero_division(
                builder, ("division by zero",), loc):
                # No exception raised => compute the +/-inf or nan result,
                # and set the FP exception word for Numpy warnings.
                quot = builder.fdiv(x, y)
                builder.store(quot, res)
        with if_non_zero:
            quot, _ = real_divmod(context, builder, x, y)
            builder.store(quot, res)
    return impl_ret_untracked(context, builder, sig.return_type,
                              builder.load(res))


def real_power_impl(context, builder, sig, args):
    x, y = args
    module = builder.module
    if context.implement_powi_as_math_call:
        imp = context.get_function(math.pow, sig)
        res = imp(builder, args)
    else:
        fn = module.declare_intrinsic('llvm.pow', [y.type])
        res = builder.call(fn, (x, y))
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_lt_impl(context, builder, sig, args):
    res = builder.fcmp_ordered('<', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_le_impl(context, builder, sig, args):
    res = builder.fcmp_ordered('<=', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_gt_impl(context, builder, sig, args):
    res = builder.fcmp_ordered('>', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_ge_impl(context, builder, sig, args):
    res = builder.fcmp_ordered('>=', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_eq_impl(context, builder, sig, args):
    res = builder.fcmp_ordered('==', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_ne_impl(context, builder, sig, args):
    res = builder.fcmp_unordered('!=', *args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_abs_impl(context, builder, sig, args):
    [ty] = sig.args
    sig = typing.signature(ty, ty)
    impl = context.get_function(math.fabs, sig)
    return impl(builder, args)


def real_negate_impl(context, builder, sig, args):
    from numba.cpython import mathimpl
    res = mathimpl.negate_real(builder, args[0])
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_positive_impl(context, builder, sig, args):
    [typ] = sig.args
    [val] = args
    res = context.cast(builder, val, typ, sig.return_type)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def real_sign_impl(context, builder, sig, args):
    """
    np.sign(float)
    """
    [x] = args
    POS = Constant(x.type, 1)
    NEG = Constant(x.type, -1)
    ZERO = Constant(x.type, 0)

    presult = cgutils.alloca_once(builder, x.type)

    is_pos = builder.fcmp_ordered('>', x, ZERO)
    is_neg = builder.fcmp_ordered('<', x, ZERO)

    with builder.if_else(is_pos) as (gt_zero, not_gt_zero):
        with gt_zero:
            builder.store(POS, presult)
        with not_gt_zero:
            with builder.if_else(is_neg) as (lt_zero, not_lt_zero):
                with lt_zero:
                    builder.store(NEG, presult)
                with not_lt_zero:
                    # For both NaN and 0, the result of sign() is simply
                    # the input value.
                    builder.store(x, presult)

    res = builder.load(presult)
    return impl_ret_untracked(context, builder, sig.return_type, res)


ty = types.Float

lower_builtin(operator.add, ty, ty)(real_add_impl)
lower_builtin(operator.iadd, ty, ty)(real_add_impl)
lower_builtin(operator.sub, ty, ty)(real_sub_impl)
lower_builtin(operator.isub, ty, ty)(real_sub_impl)
lower_builtin(operator.mul, ty, ty)(real_mul_impl)
lower_builtin(operator.imul, ty, ty)(real_mul_impl)
lower_builtin(operator.floordiv, ty, ty)(real_floordiv_impl)
lower_builtin(operator.ifloordiv, ty, ty)(real_floordiv_impl)
lower_builtin(operator.truediv, ty, ty)(real_div_impl)
lower_builtin(operator.itruediv, ty, ty)(real_div_impl)
lower_builtin(operator.mod, ty, ty)(real_mod_impl)
lower_builtin(operator.imod, ty, ty)(real_mod_impl)
lower_builtin(operator.pow, ty, ty)(real_power_impl)
lower_builtin(operator.ipow, ty, ty)(real_power_impl)
lower_builtin(pow, ty, ty)(real_power_impl)

lower_builtin(operator.eq, ty, ty)(real_eq_impl)
lower_builtin(operator.ne, ty, ty)(real_ne_impl)
lower_builtin(operator.lt, ty, ty)(real_lt_impl)
lower_builtin(operator.le, ty, ty)(real_le_impl)
lower_builtin(operator.gt, ty, ty)(real_gt_impl)
lower_builtin(operator.ge, ty, ty)(real_ge_impl)

lower_builtin(abs, ty)(real_abs_impl)

lower_builtin(operator.neg, ty)(real_negate_impl)
lower_builtin(operator.pos, ty)(real_positive_impl)

del ty


@lower_getattr(types.Complex, "real")
def complex_real_impl(context, builder, typ, value):
    cplx = context.make_complex(builder, typ, value=value)
    res = cplx.real
    return impl_ret_untracked(context, builder, typ, res)

@lower_getattr(types.Complex, "imag")
def complex_imag_impl(context, builder, typ, value):
    cplx = context.make_complex(builder, typ, value=value)
    res = cplx.imag
    return impl_ret_untracked(context, builder, typ, res)

@lower_builtin("complex.conjugate", types.Complex)
def complex_conjugate_impl(context, builder, sig, args):
    from numba.cpython import mathimpl
    z = context.make_complex(builder, sig.args[0], args[0])
    z.imag = mathimpl.negate_real(builder, z.imag)
    res = z._getvalue()
    return impl_ret_untracked(context, builder, sig.return_type, res)

def real_real_impl(context, builder, typ, value):
    return impl_ret_untracked(context, builder, typ, value)

def real_imag_impl(context, builder, typ, value):
    res = cgutils.get_null_value(value.type)
    return impl_ret_untracked(context, builder, typ, res)

def real_conjugate_impl(context, builder, sig, args):
    return impl_ret_untracked(context, builder, sig.return_type, args[0])

for cls in (types.Float, types.Integer):
    lower_getattr(cls, "real")(real_real_impl)
    lower_getattr(cls, "imag")(real_imag_impl)
    lower_builtin("complex.conjugate", cls)(real_conjugate_impl)


@lower_builtin(operator.pow, types.Complex, types.Complex)
@lower_builtin(operator.ipow, types.Complex, types.Complex)
@lower_builtin(pow, types.Complex, types.Complex)
def complex_power_impl(context, builder, sig, args):
    [ca, cb] = args
    ty = sig.args[0]
    fty = ty.underlying_float
    a = context.make_helper(builder, ty, value=ca)
    b = context.make_helper(builder, ty, value=cb)
    c = context.make_helper(builder, ty)
    module = builder.module
    pa = a._getpointer()
    pb = b._getpointer()
    pc = c._getpointer()

    # Optimize for square because cpow loses a lot of precision
    TWO = context.get_constant(fty, 2)
    ZERO = context.get_constant(fty, 0)

    b_real_is_two = builder.fcmp_ordered('==', b.real, TWO)
    b_imag_is_zero = builder.fcmp_ordered('==', b.imag, ZERO)
    b_is_two = builder.and_(b_real_is_two, b_imag_is_zero)

    with builder.if_else(b_is_two) as (then, otherwise):
        with then:
            # Lower as multiplication
            res = complex_mul_impl(context, builder, sig, (ca, ca))
            cres = context.make_helper(builder, ty, value=res)
            c.real = cres.real
            c.imag = cres.imag

        with otherwise:
            # Lower with call to external function
            func_name = {
                types.complex64: "numba_cpowf",
                types.complex128: "numba_cpow",
                }[ty]
            fnty = ir.FunctionType(ir.VoidType(), [pa.type] * 3)
            cpow = cgutils.get_or_insert_function(module, fnty, func_name)
            builder.call(cpow, (pa, pb, pc))

    res = builder.load(pc)
    return impl_ret_untracked(context, builder, sig.return_type, res)

def complex_add_impl(context, builder, sig, args):
    [cx, cy] = args
    ty = sig.args[0]
    x = context.make_complex(builder, ty, value=cx)
    y = context.make_complex(builder, ty, value=cy)
    z = context.make_complex(builder, ty)
    a = x.real
    b = x.imag
    c = y.real
    d = y.imag
    z.real = builder.fadd(a, c)
    z.imag = builder.fadd(b, d)
    res = z._getvalue()
    return impl_ret_untracked(context, builder, sig.return_type, res)


def complex_sub_impl(context, builder, sig, args):
    [cx, cy] = args
    ty = sig.args[0]
    x = context.make_complex(builder, ty, value=cx)
    y = context.make_complex(builder, ty, value=cy)
    z = context.make_complex(builder, ty)
    a = x.real
    b = x.imag
    c = y.real
    d = y.imag
    z.real = builder.fsub(a, c)
    z.imag = builder.fsub(b, d)
    res = z._getvalue()
    return impl_ret_untracked(context, builder, sig.return_type, res)


def complex_mul_impl(context, builder, sig, args):
    """
    (a+bi)(c+di)=(ac-bd)+i(ad+bc)
    """
    [cx, cy] = args
    ty = sig.args[0]
    x = context.make_complex(builder, ty, value=cx)
    y = context.make_complex(builder, ty, value=cy)
    z = context.make_complex(builder, ty)
    a = x.real
    b = x.imag
    c = y.real
    d = y.imag
    ac = builder.fmul(a, c)
    bd = builder.fmul(b, d)
    ad = builder.fmul(a, d)
    bc = builder.fmul(b, c)
    z.real = builder.fsub(ac, bd)
    z.imag = builder.fadd(ad, bc)
    res = z._getvalue()
    return impl_ret_untracked(context, builder, sig.return_type, res)


NAN = float('nan')

def complex_div_impl(context, builder, sig, args):
    def complex_div(a, b):
        # This is CPython's algorithm (in _Py_c_quot()).
        areal = a.real
        aimag = a.imag
        breal = b.real
        bimag = b.imag
        if not breal and not bimag:
            raise ZeroDivisionError("complex division by zero")
        if abs(breal) >= abs(bimag):
            # Divide tops and bottom by b.real
            if not breal:
                return complex(NAN, NAN)
            ratio = bimag / breal
            denom = breal + bimag * ratio
            return complex(
                (areal + aimag * ratio) / denom,
                (aimag - areal * ratio) / denom)
        else:
            # Divide tops and bottom by b.imag
            if not bimag:
                return complex(NAN, NAN)
            ratio = breal / bimag
            denom = breal * ratio + bimag
            return complex(
                (a.real * ratio + a.imag) / denom,
                (a.imag * ratio - a.real) / denom)

    res = context.compile_internal(builder, complex_div, sig, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def complex_negate_impl(context, builder, sig, args):
    from numba.cpython import mathimpl
    [typ] = sig.args
    [val] = args
    cmplx = context.make_complex(builder, typ, value=val)
    res = context.make_complex(builder, typ)
    res.real = mathimpl.negate_real(builder, cmplx.real)
    res.imag = mathimpl.negate_real(builder, cmplx.imag)
    res = res._getvalue()
    return impl_ret_untracked(context, builder, sig.return_type, res)


def complex_positive_impl(context, builder, sig, args):
    [val] = args
    return impl_ret_untracked(context, builder, sig.return_type, val)


def complex_eq_impl(context, builder, sig, args):
    [cx, cy] = args
    typ = sig.args[0]
    x = context.make_complex(builder, typ, value=cx)
    y = context.make_complex(builder, typ, value=cy)

    reals_are_eq = builder.fcmp_ordered('==', x.real, y.real)
    imags_are_eq = builder.fcmp_ordered('==', x.imag, y.imag)
    res = builder.and_(reals_are_eq, imags_are_eq)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def complex_ne_impl(context, builder, sig, args):
    [cx, cy] = args
    typ = sig.args[0]
    x = context.make_complex(builder, typ, value=cx)
    y = context.make_complex(builder, typ, value=cy)

    reals_are_ne = builder.fcmp_unordered('!=', x.real, y.real)
    imags_are_ne = builder.fcmp_unordered('!=', x.imag, y.imag)
    res = builder.or_(reals_are_ne, imags_are_ne)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def complex_abs_impl(context, builder, sig, args):
    """
    abs(z) := hypot(z.real, z.imag)
    """
    def complex_abs(z):
        return math.hypot(z.real, z.imag)

    res = context.compile_internal(builder, complex_abs, sig, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


ty = types.Complex

lower_builtin(operator.add, ty, ty)(complex_add_impl)
lower_builtin(operator.iadd, ty, ty)(complex_add_impl)
lower_builtin(operator.sub, ty, ty)(complex_sub_impl)
lower_builtin(operator.isub, ty, ty)(complex_sub_impl)
lower_builtin(operator.mul, ty, ty)(complex_mul_impl)
lower_builtin(operator.imul, ty, ty)(complex_mul_impl)
lower_builtin(operator.truediv, ty, ty)(complex_div_impl)
lower_builtin(operator.itruediv, ty, ty)(complex_div_impl)
lower_builtin(operator.neg, ty)(complex_negate_impl)
lower_builtin(operator.pos, ty)(complex_positive_impl)
# Complex modulo is deprecated in python3

lower_builtin(operator.eq, ty, ty)(complex_eq_impl)
lower_builtin(operator.ne, ty, ty)(complex_ne_impl)

lower_builtin(abs, ty)(complex_abs_impl)

del ty


@lower_builtin("number.item", types.Boolean)
@lower_builtin("number.item", types.Number)
def number_item_impl(context, builder, sig, args):
    """
    The no-op .item() method on booleans and numbers.
    """
    return args[0]


#------------------------------------------------------------------------------


def number_not_impl(context, builder, sig, args):
    [typ] = sig.args
    [val] = args
    istrue = context.cast(builder, val, typ, sig.return_type)
    res = builder.not_(istrue)
    return impl_ret_untracked(context, builder, sig.return_type, res)

@lower_builtin(bool, types.Boolean)
def bool_as_bool(context, builder, sig, args):
    [val] = args
    return val

@lower_builtin(bool, types.Integer)
def int_as_bool(context, builder, sig, args):
    [val] = args
    return builder.icmp_unsigned('!=', val, Constant(val.type, 0))

@lower_builtin(bool, types.Float)
def float_as_bool(context, builder, sig, args):
    [val] = args
    return builder.fcmp_unordered('!=', val, Constant(val.type, 0.0))

@lower_builtin(bool, types.Complex)
def complex_as_bool(context, builder, sig, args):
    [typ] = sig.args
    [val] = args
    cmplx = context.make_complex(builder, typ, val)
    real, imag = cmplx.real, cmplx.imag
    zero = Constant(real.type, 0.0)
    real_istrue = builder.fcmp_unordered('!=', real, zero)
    imag_istrue = builder.fcmp_unordered('!=', imag, zero)
    return builder.or_(real_istrue, imag_istrue)


for ty in (types.Integer, types.Float, types.Complex):
    lower_builtin(operator.not_, ty)(number_not_impl)


lower_builtin(operator.not_, types.py_bool)(number_not_impl)


#------------------------------------------------------------------------------
# Hashing numbers, see hashing.py

#-------------------------------------------------------------------------------
# Implicit casts between numerics

@lower_cast(types.IntegerLiteral, types.Integer)
@lower_cast(types.IntegerLiteral, types.Float)
@lower_cast(types.IntegerLiteral, types.Complex)
def literal_int_to_number(context, builder, fromty, toty, val):
    lit = context.get_constant_generic(
        builder,
        fromty.literal_type,
        fromty.literal_value,
        )
    return context.cast(builder, lit, fromty.literal_type, toty)


@lower_cast(types.Integer, types.Integer)
def integer_to_integer(context, builder, fromty, toty, val):
    if toty.bitwidth == fromty.bitwidth:
        # Just a change of signedness
        return val
    elif toty.bitwidth < fromty.bitwidth:
        # Downcast
        return builder.trunc(val, context.get_value_type(toty))
    elif fromty.signed:
        # Signed upcast
        return builder.sext(val, context.get_value_type(toty))
    else:
        # Unsigned upcast
        return builder.zext(val, context.get_value_type(toty))

@lower_cast(types.Integer, types.voidptr)
def integer_to_voidptr(context, builder, fromty, toty, val):
    return builder.inttoptr(val, context.get_value_type(toty))

@lower_cast(types.Float, types.Float)
def float_to_float(context, builder, fromty, toty, val):
    lty = context.get_value_type(toty)
    if fromty.bitwidth < toty.bitwidth:
        return builder.fpext(val, lty)
    else:
        return builder.fptrunc(val, lty)

@lower_cast(types.Integer, types.Float)
def integer_to_float(context, builder, fromty, toty, val):
    lty = context.get_value_type(toty)
    if fromty.signed:
        return builder.sitofp(val, lty)
    else:
        return builder.uitofp(val, lty)

@lower_cast(types.Float, types.Integer)
def float_to_integer(context, builder, fromty, toty, val):
    lty = context.get_value_type(toty)
    if toty.signed:
        return builder.fptosi(val, lty)
    else:
        return builder.fptoui(val, lty)

@lower_cast(types.Float, types.Complex)
@lower_cast(types.Integer, types.Complex)
def non_complex_to_complex(context, builder, fromty, toty, val):
    real = context.cast(builder, val, fromty, toty.underlying_float)
    imag = context.get_constant(toty.underlying_float, 0)

    cmplx = context.make_complex(builder, toty)
    cmplx.real = real
    cmplx.imag = imag
    return cmplx._getvalue()

@lower_cast(types.Complex, types.Complex)
def complex_to_complex(context, builder, fromty, toty, val):
    srcty = fromty.underlying_float
    dstty = toty.underlying_float

    src = context.make_complex(builder, fromty, value=val)
    dst = context.make_complex(builder, toty)
    dst.real = context.cast(builder, src.real, srcty, dstty)
    dst.imag = context.cast(builder, src.imag, srcty, dstty)
    return dst._getvalue()

@lower_cast(types.Any, types.Boolean)
def any_to_boolean(context, builder, fromty, toty, val):
    return context.is_true(builder, fromty, val)

@lower_cast(types.Boolean, types.Number)
def boolean_to_any(context, builder, fromty, toty, val):
    # Casting from boolean to anything first casts to int32
    asint = builder.zext(val, ir.IntType(32))
    return context.cast(builder, asint, types.int32, toty)

@lower_cast(types.IntegerLiteral, types.Boolean)
@lower_cast(types.BooleanLiteral, types.Boolean)
def literal_int_to_boolean(context, builder, fromty, toty, val):
    lit = context.get_constant_generic(
        builder,
        fromty.literal_type,
        fromty.literal_value,
        )
    return context.is_true(builder, fromty.literal_type, lit)

#-------------------------------------------------------------------------------
# Constants

@lower_constant(types.Complex)
def constant_complex(context, builder, ty, pyval):
    fty = ty.underlying_float
    real = context.get_constant_generic(builder, fty, pyval.real)
    imag = context.get_constant_generic(builder, fty, pyval.imag)
    return Constant.literal_struct((real, imag))

@lower_constant(types.Integer)
@lower_constant(types.Float)
@lower_constant(types.Boolean)
def constant_integer(context, builder, ty, pyval):
    # See https://github.com/numba/numba/issues/6979
    # llvmlite ir.IntType specialises the formatting of the constant for a
    # cpython bool. A NumPy np.bool_ is not a cpython bool so force it to be one
    # so that the constant renders correctly!
    if isinstance(pyval, np.bool_):
        pyval = bool(pyval)
    lty = context.get_value_type(ty)
    return lty(pyval)


#-------------------------------------------------------------------------------
# View

def scalar_view(scalar, viewty):
    """ Typing for the np scalar 'view' method. """
    if (isinstance(scalar, (types.Float, types.Integer))
            and isinstance(viewty, types.abstract.DTypeSpec)):
        if scalar.bitwidth != viewty.dtype.bitwidth:
            raise errors.TypingError(
                "Changing the dtype of a 0d array is only supported if the "
                "itemsize is unchanged")

        def impl(scalar, viewty):
            return viewer(scalar, viewty)
        return impl


overload_method(types.Float, 'view')(scalar_view)
overload_method(types.Integer, 'view')(scalar_view)
