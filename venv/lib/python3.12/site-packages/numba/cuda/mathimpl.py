import math
import operator
from llvmlite import ir
from numba.core import types, typing, cgutils, targetconfig
from numba.core.imputils import Registry
from numba.types import float32, float64, int64, uint64
from numba.cuda import libdevice
from numba import cuda

registry = Registry()
lower = registry.lower


booleans = []
booleans += [('isnand', 'isnanf', math.isnan)]
booleans += [('isinfd', 'isinff', math.isinf)]
booleans += [('isfinited', 'finitef', math.isfinite)]

unarys = []
unarys += [('ceil', 'ceilf', math.ceil)]
unarys += [('floor', 'floorf', math.floor)]
unarys += [('fabs', 'fabsf', math.fabs)]
unarys += [('exp', 'expf', math.exp)]
unarys += [('expm1', 'expm1f', math.expm1)]
unarys += [('erf', 'erff', math.erf)]
unarys += [('erfc', 'erfcf', math.erfc)]
unarys += [('tgamma', 'tgammaf', math.gamma)]
unarys += [('lgamma', 'lgammaf', math.lgamma)]
unarys += [('sqrt', 'sqrtf', math.sqrt)]
unarys += [('log', 'logf', math.log)]
unarys += [('log2', 'log2f', math.log2)]
unarys += [('log10', 'log10f', math.log10)]
unarys += [('log1p', 'log1pf', math.log1p)]
unarys += [('acosh', 'acoshf', math.acosh)]
unarys += [('acos', 'acosf', math.acos)]
unarys += [('cos', 'cosf', math.cos)]
unarys += [('cosh', 'coshf', math.cosh)]
unarys += [('asinh', 'asinhf', math.asinh)]
unarys += [('asin', 'asinf', math.asin)]
unarys += [('sin', 'sinf', math.sin)]
unarys += [('sinh', 'sinhf', math.sinh)]
unarys += [('atan', 'atanf', math.atan)]
unarys += [('atanh', 'atanhf', math.atanh)]
unarys += [('tan', 'tanf', math.tan)]
unarys += [('trunc', 'truncf', math.trunc)]

unarys_fastmath = {}
unarys_fastmath['cosf'] = 'fast_cosf'
unarys_fastmath['sinf'] = 'fast_sinf'
unarys_fastmath['tanf'] = 'fast_tanf'
unarys_fastmath['expf'] = 'fast_expf'
unarys_fastmath['log2f'] = 'fast_log2f'
unarys_fastmath['log10f'] = 'fast_log10f'
unarys_fastmath['logf'] = 'fast_logf'

binarys = []
binarys += [('copysign', 'copysignf', math.copysign)]
binarys += [('atan2', 'atan2f', math.atan2)]
binarys += [('pow', 'powf', math.pow)]
binarys += [('fmod', 'fmodf', math.fmod)]
binarys += [('hypot', 'hypotf', math.hypot)]
binarys += [('remainder', 'remainderf', math.remainder)]

binarys_fastmath = {}
binarys_fastmath['powf'] = 'fast_powf'


@lower(math.isinf, types.Integer)
@lower(math.isnan, types.Integer)
def math_isinf_isnan_int(context, builder, sig, args):
    return context.get_constant(types.boolean, 0)


@lower(operator.truediv, types.float32, types.float32)
def maybe_fast_truediv(context, builder, sig, args):
    if context.fastmath:
        sig = typing.signature(float32, float32, float32)
        impl = context.get_function(libdevice.fast_fdividef, sig)
        return impl(builder, args)
    else:
        with cgutils.if_zero(builder, args[1]):
            context.error_model.fp_zero_division(builder, ("division by zero",))
        res = builder.fdiv(*args)
        return res


@lower(math.isfinite, types.Integer)
def math_isfinite_int(context, builder, sig, args):
    return context.get_constant(types.boolean, 1)


@lower(math.sin, types.float16)
def fp16_sin_impl(context, builder, sig, args):
    def fp16_sin(x):
        return cuda.fp16.hsin(x)

    return context.compile_internal(builder, fp16_sin, sig, args)


@lower(math.cos, types.float16)
def fp16_cos_impl(context, builder, sig, args):
    def fp16_cos(x):
        return cuda.fp16.hcos(x)

    return context.compile_internal(builder, fp16_cos, sig, args)


@lower(math.log, types.float16)
def fp16_log_impl(context, builder, sig, args):
    def fp16_log(x):
        return cuda.fp16.hlog(x)

    return context.compile_internal(builder, fp16_log, sig, args)


@lower(math.log10, types.float16)
def fp16_log10_impl(context, builder, sig, args):
    def fp16_log10(x):
        return cuda.fp16.hlog10(x)

    return context.compile_internal(builder, fp16_log10, sig, args)


@lower(math.log2, types.float16)
def fp16_log2_impl(context, builder, sig, args):
    def fp16_log2(x):
        return cuda.fp16.hlog2(x)

    return context.compile_internal(builder, fp16_log2, sig, args)


@lower(math.exp, types.float16)
def fp16_exp_impl(context, builder, sig, args):
    def fp16_exp(x):
        return cuda.fp16.hexp(x)

    return context.compile_internal(builder, fp16_exp, sig, args)


@lower(math.floor, types.float16)
def fp16_floor_impl(context, builder, sig, args):
    def fp16_floor(x):
        return cuda.fp16.hfloor(x)

    return context.compile_internal(builder, fp16_floor, sig, args)


@lower(math.ceil, types.float16)
def fp16_ceil_impl(context, builder, sig, args):
    def fp16_ceil(x):
        return cuda.fp16.hceil(x)

    return context.compile_internal(builder, fp16_ceil, sig, args)


@lower(math.sqrt, types.float16)
def fp16_sqrt_impl(context, builder, sig, args):
    def fp16_sqrt(x):
        return cuda.fp16.hsqrt(x)

    return context.compile_internal(builder, fp16_sqrt, sig, args)


@lower(math.fabs, types.float16)
def fp16_fabs_impl(context, builder, sig, args):
    def fp16_fabs(x):
        return cuda.fp16.habs(x)

    return context.compile_internal(builder, fp16_fabs, sig, args)


@lower(math.trunc, types.float16)
def fp16_trunc_impl(context, builder, sig, args):
    def fp16_trunc(x):
        return cuda.fp16.htrunc(x)

    return context.compile_internal(builder, fp16_trunc, sig, args)


def impl_boolean(key, ty, libfunc):
    def lower_boolean_impl(context, builder, sig, args):
        libfunc_impl = context.get_function(libfunc,
                                            typing.signature(types.int32, ty))
        result = libfunc_impl(builder, args)
        return context.cast(builder, result, types.int32, types.boolean)

    lower(key, ty)(lower_boolean_impl)


def get_lower_unary_impl(key, ty, libfunc):
    def lower_unary_impl(context, builder, sig, args):
        actual_libfunc = libfunc
        fast_replacement = None
        if ty == float32 and context.fastmath:
            fast_replacement = unarys_fastmath.get(libfunc.__name__)

        if fast_replacement is not None:
            actual_libfunc = getattr(libdevice, fast_replacement)

        libfunc_impl = context.get_function(actual_libfunc,
                                            typing.signature(ty, ty))
        return libfunc_impl(builder, args)
    return lower_unary_impl


def get_unary_impl_for_fn_and_ty(fn, ty):
    # tanh is a special case - because it is not registered like the other
    # unary implementations, it does not appear in the unarys list. However,
    # its implementation can be looked up by key like the other
    # implementations, so we add it to the list we search here.
    tanh_impls = ('tanh', 'tanhf', math.tanh)
    for fname64, fname32, key in unarys + [tanh_impls]:
        if fn == key:
            if ty == float32:
                impl = getattr(libdevice, fname32)
            elif ty == float64:
                impl = getattr(libdevice, fname64)

            return get_lower_unary_impl(key, ty, impl)

    raise RuntimeError(f"Implementation of {fn} for {ty} not found")


def impl_unary(key, ty, libfunc):
    lower_unary_impl = get_lower_unary_impl(key, ty, libfunc)
    lower(key, ty)(lower_unary_impl)


def impl_unary_int(key, ty, libfunc):
    def lower_unary_int_impl(context, builder, sig, args):
        if sig.args[0] == int64:
            convert = builder.sitofp
        elif sig.args[0] == uint64:
            convert = builder.uitofp
        else:
            m = 'Only 64-bit integers are supported for generic unary int ops'
            raise TypeError(m)

        arg = convert(args[0], ir.DoubleType())
        sig = typing.signature(float64, float64)
        libfunc_impl = context.get_function(libfunc, sig)
        return libfunc_impl(builder, [arg])

    lower(key, ty)(lower_unary_int_impl)


def get_lower_binary_impl(key, ty, libfunc):
    def lower_binary_impl(context, builder, sig, args):
        actual_libfunc = libfunc
        fast_replacement = None
        if ty == float32 and context.fastmath:
            fast_replacement = binarys_fastmath.get(libfunc.__name__)

        if fast_replacement is not None:
            actual_libfunc = getattr(libdevice, fast_replacement)

        libfunc_impl = context.get_function(actual_libfunc,
                                            typing.signature(ty, ty, ty))
        return libfunc_impl(builder, args)
    return lower_binary_impl


def get_binary_impl_for_fn_and_ty(fn, ty):
    for fname64, fname32, key in binarys:
        if fn == key:
            if ty == float32:
                impl = getattr(libdevice, fname32)
            elif ty == float64:
                impl = getattr(libdevice, fname64)

            return get_lower_binary_impl(key, ty, impl)

    raise RuntimeError(f"Implementation of {fn} for {ty} not found")


def impl_binary(key, ty, libfunc):
    lower_binary_impl = get_lower_binary_impl(key, ty, libfunc)
    lower(key, ty, ty)(lower_binary_impl)


def impl_binary_int(key, ty, libfunc):
    def lower_binary_int_impl(context, builder, sig, args):
        if sig.args[0] == int64:
            convert = builder.sitofp
        elif sig.args[0] == uint64:
            convert = builder.uitofp
        else:
            m = 'Only 64-bit integers are supported for generic binary int ops'
            raise TypeError(m)

        args = [convert(arg, ir.DoubleType()) for arg in args]
        sig = typing.signature(float64, float64, float64)
        libfunc_impl = context.get_function(libfunc, sig)
        return libfunc_impl(builder, args)

    lower(key, ty, ty)(lower_binary_int_impl)


for fname64, fname32, key in booleans:
    impl32 = getattr(libdevice, fname32)
    impl64 = getattr(libdevice, fname64)
    impl_boolean(key, float32, impl32)
    impl_boolean(key, float64, impl64)


for fname64, fname32, key in unarys:
    impl32 = getattr(libdevice, fname32)
    impl64 = getattr(libdevice, fname64)
    impl_unary(key, float32, impl32)
    impl_unary(key, float64, impl64)
    impl_unary_int(key, int64, impl64)
    impl_unary_int(key, uint64, impl64)


for fname64, fname32, key in binarys:
    impl32 = getattr(libdevice, fname32)
    impl64 = getattr(libdevice, fname64)
    impl_binary(key, float32, impl32)
    impl_binary(key, float64, impl64)
    impl_binary_int(key, int64, impl64)
    impl_binary_int(key, uint64, impl64)


def impl_pow_int(ty, libfunc):
    def lower_pow_impl_int(context, builder, sig, args):
        powi_sig = typing.signature(ty, ty, types.int32)
        libfunc_impl = context.get_function(libfunc, powi_sig)
        return libfunc_impl(builder, args)

    lower(math.pow, ty, types.int32)(lower_pow_impl_int)


impl_pow_int(types.float32, libdevice.powif)
impl_pow_int(types.float64, libdevice.powi)


def impl_modf(ty, libfunc):
    retty = types.UniTuple(ty, 2)

    def lower_modf_impl(context, builder, sig, args):
        modf_sig = typing.signature(retty, ty)
        libfunc_impl = context.get_function(libfunc, modf_sig)
        return libfunc_impl(builder, args)

    lower(math.modf, ty)(lower_modf_impl)


impl_modf(types.float32, libdevice.modff)
impl_modf(types.float64, libdevice.modf)


def impl_frexp(ty, libfunc):
    retty = types.Tuple((ty, types.int32))

    def lower_frexp_impl(context, builder, sig, args):
        frexp_sig = typing.signature(retty, ty)
        libfunc_impl = context.get_function(libfunc, frexp_sig)
        return libfunc_impl(builder, args)

    lower(math.frexp, ty)(lower_frexp_impl)


impl_frexp(types.float32, libdevice.frexpf)
impl_frexp(types.float64, libdevice.frexp)


def impl_ldexp(ty, libfunc):
    def lower_ldexp_impl(context, builder, sig, args):
        ldexp_sig = typing.signature(ty, ty, types.int32)
        libfunc_impl = context.get_function(libfunc, ldexp_sig)
        return libfunc_impl(builder, args)

    lower(math.ldexp, ty, types.int32)(lower_ldexp_impl)


impl_ldexp(types.float32, libdevice.ldexpf)
impl_ldexp(types.float64, libdevice.ldexp)


def impl_tanh(ty, libfunc):
    def lower_tanh_impl(context, builder, sig, args):
        def get_compute_capability():
            flags = targetconfig.ConfigStack().top()
            return flags.compute_capability

        def tanh_impl_libdevice():
            tanh_sig = typing.signature(ty, ty)
            libfunc_impl = context.get_function(libfunc, tanh_sig)
            return libfunc_impl(builder, args)

        def tanhf_impl_fastmath():
            fnty = ir.FunctionType(ir.FloatType(), [ir.FloatType()])
            asm = ir.InlineAsm(fnty, 'tanh.approx.f32 $0, $1;', '=f,f')
            return builder.call(asm, args)

        if ty == float32 and context.fastmath:
            cc = get_compute_capability()
            if cc >= (7,5):
                return tanhf_impl_fastmath()

        return tanh_impl_libdevice()

    lower(math.tanh, ty)(lower_tanh_impl)


impl_tanh(types.float32, libdevice.tanhf)
impl_tanh(types.float64, libdevice.tanh)

impl_unary_int(math.tanh, int64, libdevice.tanh)
impl_unary_int(math.tanh, uint64, libdevice.tanh)

# Complex power implementations - translations of _Py_c_pow from CPython
# https://github.com/python/cpython/blob/a755410e054e1e2390de5830befc08fe80706c66/Objects/complexobject.c#L123-L151
#
# The complex64 variant casts all constants and some variables to ensure that
# as much computation is done in single precision as possible. A small number
# of operations are still done in 64-bit, but these come from libdevice code.


def cpow_implement(fty, cty):
    def core(context, builder, sig, args):
        def cpow_internal(a, b):

            if b.real == fty(0.0) and b.imag == fty(0.0):
                return cty(1.0) + cty(0.0j)
            elif a.real == fty(0.0) and b.real == fty(0.0):
                return cty(0.0) + cty(0.0j)

            vabs = math.hypot(a.real, a.imag)
            len = math.pow(vabs, b.real)
            at = math.atan2(a.imag, a.real)
            phase = at * b.real
            if b.imag != fty(0.0):
                len /= math.exp(at * b.imag)
                phase += b.imag * math.log(vabs)

            return len * (cty(math.cos(phase)) +
                          cty(math.sin(phase) * cty(1.0j)))

        return context.compile_internal(builder, cpow_internal, sig, args)

    lower(operator.pow, cty, cty)(core)
    lower(operator.ipow, cty, cty)(core)
    lower(pow, cty, cty)(core)


cpow_implement(types.float32, types.complex64)
cpow_implement(types.float64, types.complex128)
