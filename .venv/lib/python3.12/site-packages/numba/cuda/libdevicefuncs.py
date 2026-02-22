from collections import namedtuple
from textwrap import indent

from numba.types import float32, float64, int16, int32, int64, void, Tuple
from numba.core.typing.templates import signature

arg = namedtuple("arg", ("name", "ty", "is_ptr"))

functions = {
    "__nv_abs": (int32, [arg(name="x", ty=int32, is_ptr=False)]),
    "__nv_acos": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_acosf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_acosh": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_acoshf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_asin": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_asinf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_asinh": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_asinhf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_atan": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_atan2": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_atan2f": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_atanf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_atanh": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_atanhf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_brev": (int32, [arg(name="x", ty=int32, is_ptr=False)]),
    "__nv_brevll": (int64, [arg(name="x", ty=int64, is_ptr=False)]),
    "__nv_byte_perm": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
            arg(name="z", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_cbrt": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_cbrtf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_ceil": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_ceilf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_clz": (int32, [arg(name="x", ty=int32, is_ptr=False)]),
    "__nv_clzll": (int32, [arg(name="x", ty=int64, is_ptr=False)]),
    "__nv_copysign": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_copysignf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_cos": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_cosf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_cosh": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_coshf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_cospi": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_cospif": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_dadd_rd": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_dadd_rn": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_dadd_ru": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_dadd_rz": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_ddiv_rd": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_ddiv_rn": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_ddiv_ru": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_ddiv_rz": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_dmul_rd": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_dmul_rn": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_dmul_ru": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_dmul_rz": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_double2float_rd": (
        float32,
        [arg(name="d", ty=float64, is_ptr=False)],
    ),
    "__nv_double2float_rn": (
        float32,
        [arg(name="d", ty=float64, is_ptr=False)],
    ),
    "__nv_double2float_ru": (
        float32,
        [arg(name="d", ty=float64, is_ptr=False)],
    ),
    "__nv_double2float_rz": (
        float32,
        [arg(name="d", ty=float64, is_ptr=False)],
    ),
    "__nv_double2hiint": (int32, [arg(name="d", ty=float64, is_ptr=False)]),
    "__nv_double2int_rd": (int32, [arg(name="d", ty=float64, is_ptr=False)]),
    "__nv_double2int_rn": (int32, [arg(name="d", ty=float64, is_ptr=False)]),
    "__nv_double2int_ru": (int32, [arg(name="d", ty=float64, is_ptr=False)]),
    "__nv_double2int_rz": (int32, [arg(name="d", ty=float64, is_ptr=False)]),
    "__nv_double2ll_rd": (int64, [arg(name="f", ty=float64, is_ptr=False)]),
    "__nv_double2ll_rn": (int64, [arg(name="f", ty=float64, is_ptr=False)]),
    "__nv_double2ll_ru": (int64, [arg(name="f", ty=float64, is_ptr=False)]),
    "__nv_double2ll_rz": (int64, [arg(name="f", ty=float64, is_ptr=False)]),
    "__nv_double2loint": (int32, [arg(name="d", ty=float64, is_ptr=False)]),
    "__nv_double2uint_rd": (int32, [arg(name="d", ty=float64, is_ptr=False)]),
    "__nv_double2uint_rn": (int32, [arg(name="d", ty=float64, is_ptr=False)]),
    "__nv_double2uint_ru": (int32, [arg(name="d", ty=float64, is_ptr=False)]),
    "__nv_double2uint_rz": (int32, [arg(name="d", ty=float64, is_ptr=False)]),
    "__nv_double2ull_rd": (int64, [arg(name="f", ty=float64, is_ptr=False)]),
    "__nv_double2ull_rn": (int64, [arg(name="f", ty=float64, is_ptr=False)]),
    "__nv_double2ull_ru": (int64, [arg(name="f", ty=float64, is_ptr=False)]),
    "__nv_double2ull_rz": (int64, [arg(name="f", ty=float64, is_ptr=False)]),
    "__nv_double_as_longlong": (
        int64,
        [arg(name="x", ty=float64, is_ptr=False)],
    ),
    "__nv_drcp_rd": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_drcp_rn": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_drcp_ru": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_drcp_rz": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_dsqrt_rd": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_dsqrt_rn": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_dsqrt_ru": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_dsqrt_rz": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_erf": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_erfc": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_erfcf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_erfcinv": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_erfcinvf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_erfcx": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_erfcxf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_erff": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_erfinv": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_erfinvf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_exp": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_exp10": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_exp10f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_exp2": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_exp2f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_expf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_expm1": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_expm1f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fabs": (float64, [arg(name="f", ty=float64, is_ptr=False)]),
    "__nv_fabsf": (float32, [arg(name="f", ty=float32, is_ptr=False)]),
    "__nv_fadd_rd": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fadd_rn": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fadd_ru": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fadd_rz": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fast_cosf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fast_exp10f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fast_expf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fast_fdividef": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fast_log10f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fast_log2f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fast_logf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fast_powf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fast_sincosf": (
        void,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="sptr", ty=float32, is_ptr=True),
            arg(name="cptr", ty=float32, is_ptr=True),
        ],
    ),
    "__nv_fast_sinf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fast_tanf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fdim": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_fdimf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fdiv_rd": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fdiv_rn": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fdiv_ru": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fdiv_rz": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_ffs": (int32, [arg(name="x", ty=int32, is_ptr=False)]),
    "__nv_ffsll": (int32, [arg(name="x", ty=int64, is_ptr=False)]),
    "__nv_finitef": (int32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_float2half_rn": (int16, [arg(name="f", ty=float32, is_ptr=False)]),
    "__nv_float2int_rd": (int32, [arg(name="in", ty=float32, is_ptr=False)]),
    "__nv_float2int_rn": (int32, [arg(name="in", ty=float32, is_ptr=False)]),
    "__nv_float2int_ru": (int32, [arg(name="in", ty=float32, is_ptr=False)]),
    "__nv_float2int_rz": (int32, [arg(name="in", ty=float32, is_ptr=False)]),
    "__nv_float2ll_rd": (int64, [arg(name="f", ty=float32, is_ptr=False)]),
    "__nv_float2ll_rn": (int64, [arg(name="f", ty=float32, is_ptr=False)]),
    "__nv_float2ll_ru": (int64, [arg(name="f", ty=float32, is_ptr=False)]),
    "__nv_float2ll_rz": (int64, [arg(name="f", ty=float32, is_ptr=False)]),
    "__nv_float2uint_rd": (int32, [arg(name="in", ty=float32, is_ptr=False)]),
    "__nv_float2uint_rn": (int32, [arg(name="in", ty=float32, is_ptr=False)]),
    "__nv_float2uint_ru": (int32, [arg(name="in", ty=float32, is_ptr=False)]),
    "__nv_float2uint_rz": (int32, [arg(name="in", ty=float32, is_ptr=False)]),
    "__nv_float2ull_rd": (int64, [arg(name="f", ty=float32, is_ptr=False)]),
    "__nv_float2ull_rn": (int64, [arg(name="f", ty=float32, is_ptr=False)]),
    "__nv_float2ull_ru": (int64, [arg(name="f", ty=float32, is_ptr=False)]),
    "__nv_float2ull_rz": (int64, [arg(name="f", ty=float32, is_ptr=False)]),
    "__nv_float_as_int": (int32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_floor": (float64, [arg(name="f", ty=float64, is_ptr=False)]),
    "__nv_floorf": (float32, [arg(name="f", ty=float32, is_ptr=False)]),
    "__nv_fma": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
            arg(name="z", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_fma_rd": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
            arg(name="z", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_fma_rn": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
            arg(name="z", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_fma_ru": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
            arg(name="z", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_fma_rz": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
            arg(name="z", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_fmaf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
            arg(name="z", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fmaf_rd": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
            arg(name="z", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fmaf_rn": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
            arg(name="z", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fmaf_ru": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
            arg(name="z", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fmaf_rz": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
            arg(name="z", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fmax": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_fmaxf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fmin": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_fminf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fmod": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_fmodf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fmul_rd": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fmul_rn": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fmul_ru": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fmul_rz": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_frcp_rd": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_frcp_rn": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_frcp_ru": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_frcp_rz": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_frexp": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="b", ty=int32, is_ptr=True),
        ],
    ),
    "__nv_frexpf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="b", ty=int32, is_ptr=True),
        ],
    ),
    "__nv_frsqrt_rn": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fsqrt_rd": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fsqrt_rn": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fsqrt_ru": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fsqrt_rz": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_fsub_rd": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fsub_rn": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fsub_ru": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_fsub_rz": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_hadd": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_half2float": (float32, [arg(name="h", ty=int16, is_ptr=False)]),
    "__nv_hiloint2double": (
        float64,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_hypot": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_hypotf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_ilogb": (int32, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_ilogbf": (int32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_int2double_rn": (float64, [arg(name="i", ty=int32, is_ptr=False)]),
    "__nv_int2float_rd": (float32, [arg(name="in", ty=int32, is_ptr=False)]),
    "__nv_int2float_rn": (float32, [arg(name="in", ty=int32, is_ptr=False)]),
    "__nv_int2float_ru": (float32, [arg(name="in", ty=int32, is_ptr=False)]),
    "__nv_int2float_rz": (float32, [arg(name="in", ty=int32, is_ptr=False)]),
    "__nv_int_as_float": (float32, [arg(name="x", ty=int32, is_ptr=False)]),
    "__nv_isfinited": (int32, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_isinfd": (int32, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_isinff": (int32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_isnand": (int32, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_isnanf": (int32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_j0": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_j0f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_j1": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_j1f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_jn": (
        float64,
        [
            arg(name="n", ty=int32, is_ptr=False),
            arg(name="x", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_jnf": (
        float32,
        [
            arg(name="n", ty=int32, is_ptr=False),
            arg(name="x", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_ldexp": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_ldexpf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_lgamma": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_lgammaf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_ll2double_rd": (float64, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ll2double_rn": (float64, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ll2double_ru": (float64, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ll2double_rz": (float64, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ll2float_rd": (float32, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ll2float_rn": (float32, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ll2float_ru": (float32, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ll2float_rz": (float32, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_llabs": (int64, [arg(name="x", ty=int64, is_ptr=False)]),
    "__nv_llmax": (
        int64,
        [
            arg(name="x", ty=int64, is_ptr=False),
            arg(name="y", ty=int64, is_ptr=False),
        ],
    ),
    "__nv_llmin": (
        int64,
        [
            arg(name="x", ty=int64, is_ptr=False),
            arg(name="y", ty=int64, is_ptr=False),
        ],
    ),
    "__nv_llrint": (int64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_llrintf": (int64, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_llround": (int64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_llroundf": (int64, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_log": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_log10": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_log10f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_log1p": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_log1pf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_log2": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_log2f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_logb": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_logbf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_logf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_longlong_as_double": (
        float64,
        [arg(name="x", ty=int64, is_ptr=False)],
    ),
    "__nv_max": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_min": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_modf": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="b", ty=float64, is_ptr=True),
        ],
    ),
    "__nv_modff": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="b", ty=float32, is_ptr=True),
        ],
    ),
    "__nv_mul24": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_mul64hi": (
        int64,
        [
            arg(name="x", ty=int64, is_ptr=False),
            arg(name="y", ty=int64, is_ptr=False),
        ],
    ),
    "__nv_mulhi": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    # __nv_nan and __nv_nanf are excluded - they return a representation of a
    # quiet NaN, but the argument they take seems to be undocumented, and
    # follows a strange form - it is not an output like every other pointer
    # argument. If a NaN is required, one can be obtained in CUDA Python by
    # other means, e.g. `math.nan`. They are left in this list for completeness
    # / reference.
    # "__nv_nan": (float64, [arg(name="tagp", ty=int8, is_ptr=True)]),
    # "__nv_nanf": (float32, [arg(name="tagp", ty=int8, is_ptr=True)]),
    "__nv_nearbyint": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_nearbyintf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_nextafter": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_nextafterf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_normcdf": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_normcdff": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_normcdfinv": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_normcdfinvf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_popc": (int32, [arg(name="x", ty=int32, is_ptr=False)]),
    "__nv_popcll": (int32, [arg(name="x", ty=int64, is_ptr=False)]),
    "__nv_pow": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_powf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_powi": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_powif": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_rcbrt": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_rcbrtf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_remainder": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_remainderf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
        ],
    ),
    "__nv_remquo": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=float64, is_ptr=False),
            arg(name="c", ty=int32, is_ptr=True),
        ],
    ),
    "__nv_remquof": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=float32, is_ptr=False),
            arg(name="quo", ty=int32, is_ptr=True),
        ],
    ),
    "__nv_rhadd": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_rint": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_rintf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_round": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_roundf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_rsqrt": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_rsqrtf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_sad": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
            arg(name="z", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_saturatef": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_scalbn": (
        float64,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_scalbnf": (
        float32,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_signbitd": (int32, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_signbitf": (int32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_sin": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_sincos": (
        void,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="sptr", ty=float64, is_ptr=True),
            arg(name="cptr", ty=float64, is_ptr=True),
        ],
    ),
    "__nv_sincosf": (
        void,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="sptr", ty=float32, is_ptr=True),
            arg(name="cptr", ty=float32, is_ptr=True),
        ],
    ),
    "__nv_sincospi": (
        void,
        [
            arg(name="x", ty=float64, is_ptr=False),
            arg(name="sptr", ty=float64, is_ptr=True),
            arg(name="cptr", ty=float64, is_ptr=True),
        ],
    ),
    "__nv_sincospif": (
        void,
        [
            arg(name="x", ty=float32, is_ptr=False),
            arg(name="sptr", ty=float32, is_ptr=True),
            arg(name="cptr", ty=float32, is_ptr=True),
        ],
    ),
    "__nv_sinf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_sinh": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_sinhf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_sinpi": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_sinpif": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_sqrt": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_sqrtf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_tan": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_tanf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_tanh": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_tanhf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_tgamma": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_tgammaf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_trunc": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_truncf": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_uhadd": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_uint2double_rn": (float64, [arg(name="i", ty=int32, is_ptr=False)]),
    "__nv_uint2float_rd": (float32, [arg(name="in", ty=int32, is_ptr=False)]),
    "__nv_uint2float_rn": (float32, [arg(name="in", ty=int32, is_ptr=False)]),
    "__nv_uint2float_ru": (float32, [arg(name="in", ty=int32, is_ptr=False)]),
    "__nv_uint2float_rz": (float32, [arg(name="in", ty=int32, is_ptr=False)]),
    "__nv_ull2double_rd": (float64, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ull2double_rn": (float64, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ull2double_ru": (float64, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ull2double_rz": (float64, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ull2float_rd": (float32, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ull2float_rn": (float32, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ull2float_ru": (float32, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ull2float_rz": (float32, [arg(name="l", ty=int64, is_ptr=False)]),
    "__nv_ullmax": (
        int64,
        [
            arg(name="x", ty=int64, is_ptr=False),
            arg(name="y", ty=int64, is_ptr=False),
        ],
    ),
    "__nv_ullmin": (
        int64,
        [
            arg(name="x", ty=int64, is_ptr=False),
            arg(name="y", ty=int64, is_ptr=False),
        ],
    ),
    "__nv_umax": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_umin": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_umul24": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_umul64hi": (
        int64,
        [
            arg(name="x", ty=int64, is_ptr=False),
            arg(name="y", ty=int64, is_ptr=False),
        ],
    ),
    "__nv_umulhi": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_urhadd": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_usad": (
        int32,
        [
            arg(name="x", ty=int32, is_ptr=False),
            arg(name="y", ty=int32, is_ptr=False),
            arg(name="z", ty=int32, is_ptr=False),
        ],
    ),
    "__nv_y0": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_y0f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_y1": (float64, [arg(name="x", ty=float64, is_ptr=False)]),
    "__nv_y1f": (float32, [arg(name="x", ty=float32, is_ptr=False)]),
    "__nv_yn": (
        float64,
        [
            arg(name="n", ty=int32, is_ptr=False),
            arg(name="x", ty=float64, is_ptr=False),
        ],
    ),
    "__nv_ynf": (
        float32,
        [
            arg(name="n", ty=int32, is_ptr=False),
            arg(name="x", ty=float32, is_ptr=False),
        ],
    ),
}


def create_signature(retty, args):
    """
    Given the return type and arguments for a libdevice function, return the
    signature of the stub function used to call it from CUDA Python.
    """

    # Any pointer arguments should be part of the return type.
    return_types = [arg.ty for arg in args if arg.is_ptr]
    # If the return type is void, there is no point adding it to the list of
    # return types.
    if retty != void:
        return_types.insert(0, retty)

    if len(return_types) > 1:
        retty = Tuple(return_types)
    else:
        retty = return_types[0]

    argtypes = [arg.ty for arg in args if not arg.is_ptr]

    return signature(retty, *argtypes)


# The following code generates the stubs for libdevice functions.
#
# Stubs can be regenerated (e.g. if the functions dict above is modified) with:
#
# python -c "from numba.cuda.libdevicefuncs import generate_stubs; \
#            generate_stubs()" > numba/cuda/libdevice.py

docstring_template = """
See https://docs.nvidia.com/cuda/libdevice-users-guide/{func}.html

{param_types}
:rtype: {retty}
"""

param_template = """\
:param {a.name}: Argument.
:type {a.name}: {a.ty}"""


def generate_stubs():
    for name, (retty, args) in functions.items():
        # Some libdevice functions have arguments called `in`, which causes a
        # syntax error in Python, so we rename these to `x`.
        def argname(arg):
            if arg.name == "in":
                return "x"
            else:
                return arg.name

        argnames = [argname(a) for a in args if not a.is_ptr]
        argstr = ", ".join(argnames)
        signature = create_signature(retty, args)

        param_types = "\n".join(
            [param_template.format(a=a) for a in args if not a.is_ptr]
        )
        docstring = docstring_template.format(
            param_types=param_types, retty=signature.return_type, func=name
        )
        docstring = indent(docstring, "    ")
        print(f'def {name[5:]}({argstr}):\n    """{docstring}"""\n\n')
