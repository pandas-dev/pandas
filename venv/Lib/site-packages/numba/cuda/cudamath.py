import math
from numba.core import types
from numba.core.typing.templates import ConcreteTemplate, signature, Registry


registry = Registry()
infer_global = registry.register_global


@infer_global(math.acos)
@infer_global(math.acosh)
@infer_global(math.asin)
@infer_global(math.asinh)
@infer_global(math.atan)
@infer_global(math.atanh)
@infer_global(math.cosh)
@infer_global(math.degrees)
@infer_global(math.erf)
@infer_global(math.erfc)
@infer_global(math.expm1)
@infer_global(math.gamma)
@infer_global(math.lgamma)
@infer_global(math.log1p)
@infer_global(math.radians)
@infer_global(math.sinh)
@infer_global(math.tanh)
@infer_global(math.tan)
class Math_unary(ConcreteTemplate):
    cases = [
        signature(types.float64, types.int64),
        signature(types.float64, types.uint64),
        signature(types.float32, types.float32),
        signature(types.float64, types.float64),
    ]


@infer_global(math.sin)
@infer_global(math.cos)
@infer_global(math.ceil)
@infer_global(math.floor)
@infer_global(math.sqrt)
@infer_global(math.log)
@infer_global(math.log2)
@infer_global(math.log10)
@infer_global(math.exp)
@infer_global(math.fabs)
@infer_global(math.trunc)
class Math_unary_with_fp16(ConcreteTemplate):
    cases = [
        signature(types.float64, types.int64),
        signature(types.float64, types.uint64),
        signature(types.float32, types.float32),
        signature(types.float64, types.float64),
        signature(types.float16, types.float16),
    ]


@infer_global(math.atan2)
class Math_atan2(ConcreteTemplate):
    key = math.atan2
    cases = [
        signature(types.float64, types.int64, types.int64),
        signature(types.float64, types.uint64, types.uint64),
        signature(types.float32, types.float32, types.float32),
        signature(types.float64, types.float64, types.float64),
    ]


@infer_global(math.hypot)
class Math_hypot(ConcreteTemplate):
    key = math.hypot
    cases = [
        signature(types.float64, types.int64, types.int64),
        signature(types.float64, types.uint64, types.uint64),
        signature(types.float32, types.float32, types.float32),
        signature(types.float64, types.float64, types.float64),
    ]


@infer_global(math.copysign)
@infer_global(math.fmod)
class Math_binary(ConcreteTemplate):
    cases = [
        signature(types.float32, types.float32, types.float32),
        signature(types.float64, types.float64, types.float64),
    ]


@infer_global(math.remainder)
class Math_remainder(ConcreteTemplate):
    cases = [
        signature(types.float32, types.float32, types.float32),
        signature(types.float64, types.float64, types.float64),
    ]


@infer_global(math.pow)
class Math_pow(ConcreteTemplate):
    cases = [
        signature(types.float32, types.float32, types.float32),
        signature(types.float64, types.float64, types.float64),
        signature(types.float32, types.float32, types.int32),
        signature(types.float64, types.float64, types.int32),
    ]


@infer_global(math.frexp)
class Math_frexp(ConcreteTemplate):
    cases = [
        signature(types.Tuple([types.float32, types.int32]), types.float32),
        signature(types.Tuple([types.float64, types.int32]), types.float64),
    ]


@infer_global(math.ldexp)
class Math_ldexp(ConcreteTemplate):
    cases = [
        signature(types.float32, types.float32, types.int32),
        signature(types.float64, types.float64, types.int32),
    ]


@infer_global(math.isinf)
@infer_global(math.isnan)
@infer_global(math.isfinite)
class Math_isnan(ConcreteTemplate):
    cases = [
        signature(types.boolean, types.int64),
        signature(types.boolean, types.uint64),
        signature(types.boolean, types.float32),
        signature(types.boolean, types.float64),
    ]


@infer_global(math.modf)
class Math_modf(ConcreteTemplate):
    cases = [
        signature(types.UniTuple(types.float64, 2), types.float64),
        signature(types.UniTuple(types.float32, 2), types.float32)
    ]
