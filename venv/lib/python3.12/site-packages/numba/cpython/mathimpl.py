"""
Provide math calls that uses intrinsics or libc math functions.
"""

import math
import operator
import sys
import numpy as np

import llvmlite.ir
from llvmlite.ir import Constant

from numba.core.imputils import Registry, impl_ret_untracked
from numba import typeof
from numba.core import types, utils, config, cgutils
from numba.core.extending import overload
from numba.core.typing import signature
from numba.cpython.unsafe.numbers import trailing_zeros


registry = Registry('mathimpl')
lower = registry.lower


# Helpers, shared with cmathimpl.
_NP_FLT_FINFO = np.finfo(np.dtype('float32'))
FLT_MAX = _NP_FLT_FINFO.max
FLT_MIN = _NP_FLT_FINFO.tiny

_NP_DBL_FINFO = np.finfo(np.dtype('float64'))
DBL_MAX = _NP_DBL_FINFO.max
DBL_MIN = _NP_DBL_FINFO.tiny

FLOAT_ABS_MASK = 0x7fffffff
FLOAT_SIGN_MASK = 0x80000000
DOUBLE_ABS_MASK = 0x7fffffffffffffff
DOUBLE_SIGN_MASK = 0x8000000000000000


def is_nan(builder, val):
    """
    Return a condition testing whether *val* is a NaN.
    """
    return builder.fcmp_unordered('uno', val, val)

def is_inf(builder, val):
    """
    Return a condition testing whether *val* is an infinite.
    """
    pos_inf = Constant(val.type, float("+inf"))
    neg_inf = Constant(val.type, float("-inf"))
    isposinf = builder.fcmp_ordered('==', val, pos_inf)
    isneginf = builder.fcmp_ordered('==', val, neg_inf)
    return builder.or_(isposinf, isneginf)

def is_finite(builder, val):
    """
    Return a condition testing whether *val* is a finite.
    """
    # is_finite(x)  <=>  x - x != NaN
    val_minus_val = builder.fsub(val, val)
    return builder.fcmp_ordered('ord', val_minus_val, val_minus_val)

def f64_as_int64(builder, val):
    """
    Bitcast a double into a 64-bit integer.
    """
    assert val.type == llvmlite.ir.DoubleType()
    return builder.bitcast(val, llvmlite.ir.IntType(64))

def int64_as_f64(builder, val):
    """
    Bitcast a 64-bit integer into a double.
    """
    assert val.type == llvmlite.ir.IntType(64)
    return builder.bitcast(val, llvmlite.ir.DoubleType())

def f32_as_int32(builder, val):
    """
    Bitcast a float into a 32-bit integer.
    """
    assert val.type == llvmlite.ir.FloatType()
    return builder.bitcast(val, llvmlite.ir.IntType(32))

def int32_as_f32(builder, val):
    """
    Bitcast a 32-bit integer into a float.
    """
    assert val.type == llvmlite.ir.IntType(32)
    return builder.bitcast(val, llvmlite.ir.FloatType())

def negate_real(builder, val):
    """
    Negate real number *val*, with proper handling of zeros.
    """
    # The negative zero forces LLVM to handle signed zeros properly.
    return builder.fsub(Constant(val.type, -0.0), val)

def call_fp_intrinsic(builder, name, args):
    """
    Call a LLVM intrinsic floating-point operation.
    """
    mod = builder.module
    intr = mod.declare_intrinsic(name, [a.type for a in args])
    return builder.call(intr, args)


def _unary_int_input_wrapper_impl(wrapped_impl):
    """
    Return an implementation factory to convert the single integral input
    argument to a float64, then defer to the *wrapped_impl*.
    """
    def implementer(context, builder, sig, args):
        val, = args
        input_type = sig.args[0]
        fpval = context.cast(builder, val, input_type, types.float64)
        inner_sig = signature(types.float64, types.float64)
        res = wrapped_impl(context, builder, inner_sig, (fpval,))
        return context.cast(builder, res, types.float64, sig.return_type)

    return implementer

def unary_math_int_impl(fn, float_impl):
    impl = _unary_int_input_wrapper_impl(float_impl)
    lower(fn, types.Integer)(impl)

def unary_math_intr(fn, intrcode):
    """
    Implement the math function *fn* using the LLVM intrinsic *intrcode*.
    """
    @lower(fn, types.Float)
    def float_impl(context, builder, sig, args):
        res = call_fp_intrinsic(builder, intrcode, args)
        return impl_ret_untracked(context, builder, sig.return_type, res)

    unary_math_int_impl(fn, float_impl)
    return float_impl

def unary_math_extern(fn, f32extern, f64extern, int_restype=False):
    """
    Register implementations of Python function *fn* using the
    external function named *f32extern* and *f64extern* (for float32
    and float64 inputs, respectively).
    If *int_restype* is true, then the function's return value should be
    integral, otherwise floating-point.
    """
    f_restype = types.int64 if int_restype else None

    def float_impl(context, builder, sig, args):
        """
        Implement *fn* for a types.Float input.
        """
        [val] = args
        mod = builder.module
        input_type = sig.args[0]
        lty = context.get_value_type(input_type)
        func_name = {
            types.float32: f32extern,
            types.float64: f64extern,
            }[input_type]
        fnty = llvmlite.ir.FunctionType(lty, [lty])
        fn = cgutils.insert_pure_function(builder.module, fnty, name=func_name)
        res = builder.call(fn, (val,))
        res = context.cast(builder, res, input_type, sig.return_type)
        return impl_ret_untracked(context, builder, sig.return_type, res)

    lower(fn, types.Float)(float_impl)

    # Implement wrapper for integer inputs
    unary_math_int_impl(fn, float_impl)

    return float_impl


unary_math_intr(math.fabs, 'llvm.fabs')
exp_impl = unary_math_intr(math.exp, 'llvm.exp')
log_impl = unary_math_intr(math.log, 'llvm.log')
log10_impl = unary_math_intr(math.log10, 'llvm.log10')
log2_impl = unary_math_intr(math.log2, 'llvm.log2')
sin_impl = unary_math_intr(math.sin, 'llvm.sin')
cos_impl = unary_math_intr(math.cos, 'llvm.cos')

log1p_impl = unary_math_extern(math.log1p, "log1pf", "log1p")
expm1_impl = unary_math_extern(math.expm1, "expm1f", "expm1")
erf_impl = unary_math_extern(math.erf, "erff", "erf")
erfc_impl = unary_math_extern(math.erfc, "erfcf", "erfc")

tan_impl = unary_math_extern(math.tan, "tanf", "tan")
asin_impl = unary_math_extern(math.asin, "asinf", "asin")
acos_impl = unary_math_extern(math.acos, "acosf", "acos")
atan_impl = unary_math_extern(math.atan, "atanf", "atan")

asinh_impl = unary_math_extern(math.asinh, "asinhf", "asinh")
acosh_impl = unary_math_extern(math.acosh, "acoshf", "acosh")
atanh_impl = unary_math_extern(math.atanh, "atanhf", "atanh")
sinh_impl = unary_math_extern(math.sinh, "sinhf", "sinh")
cosh_impl = unary_math_extern(math.cosh, "coshf", "cosh")
tanh_impl = unary_math_extern(math.tanh, "tanhf", "tanh")

log2_impl = unary_math_extern(math.log2, "log2f", "log2")
ceil_impl = unary_math_extern(math.ceil, "ceilf", "ceil", True)
floor_impl = unary_math_extern(math.floor, "floorf", "floor", True)

gamma_impl = unary_math_extern(math.gamma, "numba_gammaf", "numba_gamma") # work-around
sqrt_impl = unary_math_extern(math.sqrt, "sqrtf", "sqrt")
trunc_impl = unary_math_extern(math.trunc, "truncf", "trunc", True)
lgamma_impl = unary_math_extern(math.lgamma, "lgammaf", "lgamma")


@lower(math.isnan, types.Float)
def isnan_float_impl(context, builder, sig, args):
    [val] = args
    res = is_nan(builder, val)
    return impl_ret_untracked(context, builder, sig.return_type, res)

@lower(math.isnan, types.Integer)
def isnan_int_impl(context, builder, sig, args):
    res = cgutils.false_bit
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower(math.isinf, types.Float)
def isinf_float_impl(context, builder, sig, args):
    [val] = args
    res = is_inf(builder, val)
    return impl_ret_untracked(context, builder, sig.return_type, res)

@lower(math.isinf, types.Integer)
def isinf_int_impl(context, builder, sig, args):
    res = cgutils.false_bit
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower(math.isfinite, types.Float)
def isfinite_float_impl(context, builder, sig, args):
    [val] = args
    res = is_finite(builder, val)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower(math.isfinite, types.Integer)
def isfinite_int_impl(context, builder, sig, args):
    res = cgutils.true_bit
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower(math.copysign, types.Float, types.Float)
def copysign_float_impl(context, builder, sig, args):
    lty = args[0].type
    mod = builder.module
    fn = cgutils.get_or_insert_function(mod, llvmlite.ir.FunctionType(lty, (lty, lty)),
                                        'llvm.copysign.%s' % lty.intrinsic_name)
    res = builder.call(fn, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


# -----------------------------------------------------------------------------


@lower(math.frexp, types.Float)
def frexp_impl(context, builder, sig, args):
    val, = args
    fltty = context.get_data_type(sig.args[0])
    intty = context.get_data_type(sig.return_type[1])
    expptr = cgutils.alloca_once(builder, intty, name='exp')
    fnty = llvmlite.ir.FunctionType(fltty, (fltty, llvmlite.ir.PointerType(intty)))
    fname = {
        "float": "numba_frexpf",
        "double": "numba_frexp",
        }[str(fltty)]
    fn = cgutils.get_or_insert_function(builder.module, fnty, fname)
    res = builder.call(fn, (val, expptr))
    res = cgutils.make_anonymous_struct(builder, (res, builder.load(expptr)))
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower(math.ldexp, types.Float, types.intc)
def ldexp_impl(context, builder, sig, args):
    val, exp = args
    fltty, intty = map(context.get_data_type, sig.args)
    fnty = llvmlite.ir.FunctionType(fltty, (fltty, intty))
    fname = {
        "float": "numba_ldexpf",
        "double": "numba_ldexp",
        }[str(fltty)]
    fn = cgutils.insert_pure_function(builder.module, fnty, name=fname)
    res = builder.call(fn, (val, exp))
    return impl_ret_untracked(context, builder, sig.return_type, res)


# -----------------------------------------------------------------------------


@lower(math.atan2, types.int64, types.int64)
def atan2_s64_impl(context, builder, sig, args):
    [y, x] = args
    y = builder.sitofp(y, llvmlite.ir.DoubleType())
    x = builder.sitofp(x, llvmlite.ir.DoubleType())
    fsig = signature(types.float64, types.float64, types.float64)
    return atan2_float_impl(context, builder, fsig, (y, x))

@lower(math.atan2, types.uint64, types.uint64)
def atan2_u64_impl(context, builder, sig, args):
    [y, x] = args
    y = builder.uitofp(y, llvmlite.ir.DoubleType())
    x = builder.uitofp(x, llvmlite.ir.DoubleType())
    fsig = signature(types.float64, types.float64, types.float64)
    return atan2_float_impl(context, builder, fsig, (y, x))

@lower(math.atan2, types.Float, types.Float)
def atan2_float_impl(context, builder, sig, args):
    assert len(args) == 2
    mod = builder.module
    ty = sig.args[0]
    lty = context.get_value_type(ty)
    func_name = {
        types.float32: "atan2f",
        types.float64: "atan2"
        }[ty]
    fnty = llvmlite.ir.FunctionType(lty, (lty, lty))
    fn = cgutils.insert_pure_function(builder.module, fnty, name=func_name)
    res = builder.call(fn, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


# -----------------------------------------------------------------------------


@lower(math.hypot, types.int64, types.int64)
def hypot_s64_impl(context, builder, sig, args):
    [x, y] = args
    y = builder.sitofp(y, llvmlite.ir.DoubleType())
    x = builder.sitofp(x, llvmlite.ir.DoubleType())
    fsig = signature(types.float64, types.float64, types.float64)
    res = hypot_float_impl(context, builder, fsig, (x, y))
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower(math.hypot, types.uint64, types.uint64)
def hypot_u64_impl(context, builder, sig, args):
    [x, y] = args
    y = builder.sitofp(y, llvmlite.ir.DoubleType())
    x = builder.sitofp(x, llvmlite.ir.DoubleType())
    fsig = signature(types.float64, types.float64, types.float64)
    res = hypot_float_impl(context, builder, fsig, (x, y))
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower(math.hypot, types.Float, types.Float)
def hypot_float_impl(context, builder, sig, args):
    xty, yty = sig.args
    assert xty == yty == sig.return_type
    x, y = args

    # Windows has alternate names for hypot/hypotf, see
    # https://msdn.microsoft.com/fr-fr/library/a9yb3dbt%28v=vs.80%29.aspx
    fname = {
        types.float32: "_hypotf" if sys.platform == 'win32' else "hypotf",
        types.float64: "_hypot" if sys.platform == 'win32' else "hypot",
    }[xty]
    plat_hypot = types.ExternalFunction(fname, sig)

    if sys.platform == 'win32' and config.MACHINE_BITS == 32:
        inf = xty(float('inf'))

        def hypot_impl(x, y):
            if math.isinf(x) or math.isinf(y):
                return inf
            return plat_hypot(x, y)
    else:
        def hypot_impl(x, y):
            return plat_hypot(x, y)

    res = context.compile_internal(builder, hypot_impl, sig, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


# -----------------------------------------------------------------------------

@lower(math.radians, types.Float)
def radians_float_impl(context, builder, sig, args):
    [x] = args
    coef = context.get_constant(sig.return_type, math.pi / 180)
    res = builder.fmul(x, coef)
    return impl_ret_untracked(context, builder, sig.return_type, res)

unary_math_int_impl(math.radians, radians_float_impl)

# -----------------------------------------------------------------------------

@lower(math.degrees, types.Float)
def degrees_float_impl(context, builder, sig, args):
    [x] = args
    coef = context.get_constant(sig.return_type, 180 / math.pi)
    res = builder.fmul(x, coef)
    return impl_ret_untracked(context, builder, sig.return_type, res)

unary_math_int_impl(math.degrees, degrees_float_impl)

# -----------------------------------------------------------------------------

@lower(math.pow, types.Float, types.Float)
@lower(math.pow, types.Float, types.Integer)
def pow_impl(context, builder, sig, args):
    impl = context.get_function(operator.pow, sig)
    return impl(builder, args)

# -----------------------------------------------------------------------------

@lower(math.nextafter, types.Float, types.Float)
def nextafter_impl(context, builder, sig, args):
    assert len(args) == 2
    ty = sig.args[0]
    lty = context.get_value_type(ty)
    func_name = {
        types.float32: "nextafterf",
        types.float64: "nextafter"
        }[ty]
    fnty = llvmlite.ir.FunctionType(lty, (lty, lty))
    fn = cgutils.insert_pure_function(builder.module, fnty, name=func_name)
    res = builder.call(fn, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)

# -----------------------------------------------------------------------------

def _unsigned(T):
    """Convert integer to unsigned integer of equivalent width."""
    pass

@overload(_unsigned)
def _unsigned_impl(T):
    if T in types.unsigned_domain:
        return lambda T: T
    elif T in types.signed_domain:
        newT = getattr(types, 'uint{}'.format(T.bitwidth))
        return lambda T: newT(T)


def gcd_impl(context, builder, sig, args):
    xty, yty = sig.args
    assert xty == yty == sig.return_type
    x, y = args

    def gcd(a, b):
        """
        Stein's algorithm, heavily cribbed from Julia implementation.
        """
        T = type(a)
        if a == 0: return abs(b)
        if b == 0: return abs(a)
        za = trailing_zeros(a)
        zb = trailing_zeros(b)
        k = min(za, zb)
        # Uses np.*_shift instead of operators due to return types
        u = _unsigned(abs(np.right_shift(a, za)))
        v = _unsigned(abs(np.right_shift(b, zb)))
        while u != v:
            if u > v:
                u, v = v, u
            v -= u
            v = np.right_shift(v, trailing_zeros(v))
        r = np.left_shift(T(u), k)
        return r

    res = context.compile_internal(builder, gcd, sig, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


lower(math.gcd, types.Integer, types.Integer)(gcd_impl)
