"""Contains information on how to translate different ufuncs for the CUDA
target. It is a database of different ufuncs and how each of its loops maps to
a function that implements the inner kernel of that ufunc (the inner kernel
being the per-element function).

Use get_ufunc_info() to get the information related to a ufunc.
"""

import math
import numpy as np
from functools import lru_cache
from numba.core import typing


def get_ufunc_info(ufunc_key):
    return ufunc_db()[ufunc_key]


@lru_cache
def ufunc_db():
    # Imports here are at function scope to avoid circular imports
    from numba.cpython import cmathimpl, mathimpl, numbers
    from numba.np import npyfuncs
    from numba.np.numpy_support import numpy_version
    from numba.cuda.mathimpl import (get_unary_impl_for_fn_and_ty,
                                     get_binary_impl_for_fn_and_ty)

    def np_unary_impl(fn, context, builder, sig, args):
        npyfuncs._check_arity_and_homogeneity(sig, args, 1)
        impl = get_unary_impl_for_fn_and_ty(fn, sig.args[0])
        return impl(context, builder, sig, args)

    def np_binary_impl(fn, context, builder, sig, args):
        npyfuncs._check_arity_and_homogeneity(sig, args, 2)
        impl = get_binary_impl_for_fn_and_ty(fn, sig.args[0])
        return impl(context, builder, sig, args)

    def np_real_log_impl(context, builder, sig, args):
        return np_unary_impl(math.log, context, builder, sig, args)

    def np_real_log2_impl(context, builder, sig, args):
        return np_unary_impl(math.log2, context, builder, sig, args)

    def np_real_log10_impl(context, builder, sig, args):
        return np_unary_impl(math.log10, context, builder, sig, args)

    def np_real_sin_impl(context, builder, sig, args):
        return np_unary_impl(math.sin, context, builder, sig, args)

    def np_real_cos_impl(context, builder, sig, args):
        return np_unary_impl(math.cos, context, builder, sig, args)

    def np_real_tan_impl(context, builder, sig, args):
        return np_unary_impl(math.tan, context, builder, sig, args)

    def np_real_asin_impl(context, builder, sig, args):
        return np_unary_impl(math.asin, context, builder, sig, args)

    def np_real_acos_impl(context, builder, sig, args):
        return np_unary_impl(math.acos, context, builder, sig, args)

    def np_real_atan_impl(context, builder, sig, args):
        return np_unary_impl(math.atan, context, builder, sig, args)

    def np_real_atan2_impl(context, builder, sig, args):
        return np_binary_impl(math.atan2, context, builder, sig, args)

    def np_real_hypot_impl(context, builder, sig, args):
        return np_binary_impl(math.hypot, context, builder, sig, args)

    def np_real_sinh_impl(context, builder, sig, args):
        return np_unary_impl(math.sinh, context, builder, sig, args)

    def np_complex_sinh_impl(context, builder, sig, args):
        # npymath does not provide a complex sinh. The code in funcs.inc.src
        # is translated here...
        npyfuncs._check_arity_and_homogeneity(sig, args, 1)

        ty = sig.args[0]
        fty = ty.underlying_float
        fsig1 = typing.signature(*[fty] * 2)
        x = context.make_complex(builder, ty, args[0])
        out = context.make_complex(builder, ty)
        xr = x.real
        xi = x.imag

        sxi = np_real_sin_impl(context, builder, fsig1, [xi])
        shxr = np_real_sinh_impl(context, builder, fsig1, [xr])
        cxi = np_real_cos_impl(context, builder, fsig1, [xi])
        chxr = np_real_cosh_impl(context, builder, fsig1, [xr])

        out.real = builder.fmul(cxi, shxr)
        out.imag = builder.fmul(sxi, chxr)

        return out._getvalue()

    def np_real_cosh_impl(context, builder, sig, args):
        return np_unary_impl(math.cosh, context, builder, sig, args)

    def np_complex_cosh_impl(context, builder, sig, args):
        # npymath does not provide a complex cosh. The code in funcs.inc.src
        # is translated here...
        npyfuncs._check_arity_and_homogeneity(sig, args, 1)

        ty = sig.args[0]
        fty = ty.underlying_float
        fsig1 = typing.signature(*[fty] * 2)
        x = context.make_complex(builder, ty, args[0])
        out = context.make_complex(builder, ty)
        xr = x.real
        xi = x.imag

        cxi = np_real_cos_impl(context, builder, fsig1, [xi])
        chxr = np_real_cosh_impl(context, builder, fsig1, [xr])
        sxi = np_real_sin_impl(context, builder, fsig1, [xi])
        shxr = np_real_sinh_impl(context, builder, fsig1, [xr])

        out.real = builder.fmul(cxi, chxr)
        out.imag = builder.fmul(sxi, shxr)

        return out._getvalue()

    def np_real_tanh_impl(context, builder, sig, args):
        return np_unary_impl(math.tanh, context, builder, sig, args)

    def np_complex_tanh_impl(context, builder, sig, args):
        # npymath does not provide complex tan functions. The code
        # in funcs.inc.src for tanh is translated here...
        npyfuncs._check_arity_and_homogeneity(sig, args, 1)

        ty = sig.args[0]
        fty = ty.underlying_float
        fsig1 = typing.signature(*[fty] * 2)
        ONE = context.get_constant(fty, 1.0)
        x = context.make_complex(builder, ty, args[0])
        out = context.make_complex(builder, ty)

        xr = x.real
        xi = x.imag
        si = np_real_sin_impl(context, builder, fsig1, [xi])
        ci = np_real_cos_impl(context, builder, fsig1, [xi])
        shr = np_real_sinh_impl(context, builder, fsig1, [xr])
        chr_ = np_real_cosh_impl(context, builder, fsig1, [xr])
        rs = builder.fmul(ci, shr)
        is_ = builder.fmul(si, chr_)
        rc = builder.fmul(ci, chr_)
        # Note: opposite sign for `ic` from code in funcs.inc.src
        ic = builder.fmul(si, shr)
        sqr_rc = builder.fmul(rc, rc)
        sqr_ic = builder.fmul(ic, ic)
        d = builder.fadd(sqr_rc, sqr_ic)
        inv_d = builder.fdiv(ONE, d)
        rs_rc = builder.fmul(rs, rc)
        is_ic = builder.fmul(is_, ic)
        is_rc = builder.fmul(is_, rc)
        rs_ic = builder.fmul(rs, ic)
        numr = builder.fadd(rs_rc, is_ic)
        numi = builder.fsub(is_rc, rs_ic)
        out.real = builder.fmul(numr, inv_d)
        out.imag = builder.fmul(numi, inv_d)

        return out._getvalue()

    def np_real_asinh_impl(context, builder, sig, args):
        return np_unary_impl(math.asinh, context, builder, sig, args)

    def np_real_acosh_impl(context, builder, sig, args):
        return np_unary_impl(math.acosh, context, builder, sig, args)

    def np_real_atanh_impl(context, builder, sig, args):
        return np_unary_impl(math.atanh, context, builder, sig, args)

    db = {}

    db[np.sin] = {
        'f->f': np_real_sin_impl,
        'd->d': np_real_sin_impl,
        'F->F': npyfuncs.np_complex_sin_impl,
        'D->D': npyfuncs.np_complex_sin_impl,
    }

    db[np.cos] = {
        'f->f': np_real_cos_impl,
        'd->d': np_real_cos_impl,
        'F->F': npyfuncs.np_complex_cos_impl,
        'D->D': npyfuncs.np_complex_cos_impl,
    }

    db[np.tan] = {
        'f->f': np_real_tan_impl,
        'd->d': np_real_tan_impl,
        'F->F': cmathimpl.tan_impl,
        'D->D': cmathimpl.tan_impl,
    }

    db[np.arcsin] = {
        'f->f': np_real_asin_impl,
        'd->d': np_real_asin_impl,
        'F->F': cmathimpl.asin_impl,
        'D->D': cmathimpl.asin_impl,
    }

    db[np.arccos] = {
        'f->f': np_real_acos_impl,
        'd->d': np_real_acos_impl,
        'F->F': cmathimpl.acos_impl,
        'D->D': cmathimpl.acos_impl,
    }

    db[np.arctan] = {
        'f->f': np_real_atan_impl,
        'd->d': np_real_atan_impl,
        'F->F': cmathimpl.atan_impl,
        'D->D': cmathimpl.atan_impl,
    }

    db[np.arctan2] = {
        'ff->f': np_real_atan2_impl,
        'dd->d': np_real_atan2_impl,
    }

    db[np.hypot] = {
        'ff->f': np_real_hypot_impl,
        'dd->d': np_real_hypot_impl,
    }

    db[np.sinh] = {
        'f->f': np_real_sinh_impl,
        'd->d': np_real_sinh_impl,
        'F->F': np_complex_sinh_impl,
        'D->D': np_complex_sinh_impl,
    }

    db[np.cosh] = {
        'f->f': np_real_cosh_impl,
        'd->d': np_real_cosh_impl,
        'F->F': np_complex_cosh_impl,
        'D->D': np_complex_cosh_impl,
    }

    db[np.tanh] = {
        'f->f': np_real_tanh_impl,
        'd->d': np_real_tanh_impl,
        'F->F': np_complex_tanh_impl,
        'D->D': np_complex_tanh_impl,
    }

    db[np.arcsinh] = {
        'f->f': np_real_asinh_impl,
        'd->d': np_real_asinh_impl,
        'F->F': cmathimpl.asinh_impl,
        'D->D': cmathimpl.asinh_impl,
    }

    db[np.arccosh] = {
        'f->f': np_real_acosh_impl,
        'd->d': np_real_acosh_impl,
        'F->F': npyfuncs.np_complex_acosh_impl,
        'D->D': npyfuncs.np_complex_acosh_impl,
    }

    db[np.arctanh] = {
        'f->f': np_real_atanh_impl,
        'd->d': np_real_atanh_impl,
        'F->F': cmathimpl.atanh_impl,
        'D->D': cmathimpl.atanh_impl,
    }

    db[np.deg2rad] = {
        'f->f': mathimpl.radians_float_impl,
        'd->d': mathimpl.radians_float_impl,
    }

    db[np.radians] = db[np.deg2rad]

    db[np.rad2deg] = {
        'f->f': mathimpl.degrees_float_impl,
        'd->d': mathimpl.degrees_float_impl,
    }

    db[np.degrees] = db[np.rad2deg]

    db[np.greater] = {
        '??->?': numbers.int_ugt_impl,
        'bb->?': numbers.int_sgt_impl,
        'BB->?': numbers.int_ugt_impl,
        'hh->?': numbers.int_sgt_impl,
        'HH->?': numbers.int_ugt_impl,
        'ii->?': numbers.int_sgt_impl,
        'II->?': numbers.int_ugt_impl,
        'll->?': numbers.int_sgt_impl,
        'LL->?': numbers.int_ugt_impl,
        'qq->?': numbers.int_sgt_impl,
        'QQ->?': numbers.int_ugt_impl,
        'ff->?': numbers.real_gt_impl,
        'dd->?': numbers.real_gt_impl,
        'FF->?': npyfuncs.np_complex_gt_impl,
        'DD->?': npyfuncs.np_complex_gt_impl,
    }
    if numpy_version >= (1, 25):
        db[np.greater].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('>'),
            'Qq->?': numbers.int_unsigned_signed_cmp('>')})

    db[np.greater_equal] = {
        '??->?': numbers.int_uge_impl,
        'bb->?': numbers.int_sge_impl,
        'BB->?': numbers.int_uge_impl,
        'hh->?': numbers.int_sge_impl,
        'HH->?': numbers.int_uge_impl,
        'ii->?': numbers.int_sge_impl,
        'II->?': numbers.int_uge_impl,
        'll->?': numbers.int_sge_impl,
        'LL->?': numbers.int_uge_impl,
        'qq->?': numbers.int_sge_impl,
        'QQ->?': numbers.int_uge_impl,
        'ff->?': numbers.real_ge_impl,
        'dd->?': numbers.real_ge_impl,
        'FF->?': npyfuncs.np_complex_ge_impl,
        'DD->?': npyfuncs.np_complex_ge_impl,
    }
    if numpy_version >= (1, 25):
        db[np.greater_equal].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('>='),
            'Qq->?': numbers.int_unsigned_signed_cmp('>=')})

    db[np.less] = {
        '??->?': numbers.int_ult_impl,
        'bb->?': numbers.int_slt_impl,
        'BB->?': numbers.int_ult_impl,
        'hh->?': numbers.int_slt_impl,
        'HH->?': numbers.int_ult_impl,
        'ii->?': numbers.int_slt_impl,
        'II->?': numbers.int_ult_impl,
        'll->?': numbers.int_slt_impl,
        'LL->?': numbers.int_ult_impl,
        'qq->?': numbers.int_slt_impl,
        'QQ->?': numbers.int_ult_impl,
        'ff->?': numbers.real_lt_impl,
        'dd->?': numbers.real_lt_impl,
        'FF->?': npyfuncs.np_complex_lt_impl,
        'DD->?': npyfuncs.np_complex_lt_impl,
    }
    if numpy_version >= (1, 25):
        db[np.less].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('<'),
            'Qq->?': numbers.int_unsigned_signed_cmp('<')})

    db[np.less_equal] = {
        '??->?': numbers.int_ule_impl,
        'bb->?': numbers.int_sle_impl,
        'BB->?': numbers.int_ule_impl,
        'hh->?': numbers.int_sle_impl,
        'HH->?': numbers.int_ule_impl,
        'ii->?': numbers.int_sle_impl,
        'II->?': numbers.int_ule_impl,
        'll->?': numbers.int_sle_impl,
        'LL->?': numbers.int_ule_impl,
        'qq->?': numbers.int_sle_impl,
        'QQ->?': numbers.int_ule_impl,
        'ff->?': numbers.real_le_impl,
        'dd->?': numbers.real_le_impl,
        'FF->?': npyfuncs.np_complex_le_impl,
        'DD->?': npyfuncs.np_complex_le_impl,
    }
    if numpy_version >= (1, 25):
        db[np.less_equal].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('<='),
            'Qq->?': numbers.int_unsigned_signed_cmp('<=')})

    db[np.not_equal] = {
        '??->?': numbers.int_ne_impl,
        'bb->?': numbers.int_ne_impl,
        'BB->?': numbers.int_ne_impl,
        'hh->?': numbers.int_ne_impl,
        'HH->?': numbers.int_ne_impl,
        'ii->?': numbers.int_ne_impl,
        'II->?': numbers.int_ne_impl,
        'll->?': numbers.int_ne_impl,
        'LL->?': numbers.int_ne_impl,
        'qq->?': numbers.int_ne_impl,
        'QQ->?': numbers.int_ne_impl,
        'ff->?': numbers.real_ne_impl,
        'dd->?': numbers.real_ne_impl,
        'FF->?': npyfuncs.np_complex_ne_impl,
        'DD->?': npyfuncs.np_complex_ne_impl,
    }
    if numpy_version >= (1, 25):
        db[np.not_equal].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('!='),
            'Qq->?': numbers.int_unsigned_signed_cmp('!=')})

    db[np.equal] = {
        '??->?': numbers.int_eq_impl,
        'bb->?': numbers.int_eq_impl,
        'BB->?': numbers.int_eq_impl,
        'hh->?': numbers.int_eq_impl,
        'HH->?': numbers.int_eq_impl,
        'ii->?': numbers.int_eq_impl,
        'II->?': numbers.int_eq_impl,
        'll->?': numbers.int_eq_impl,
        'LL->?': numbers.int_eq_impl,
        'qq->?': numbers.int_eq_impl,
        'QQ->?': numbers.int_eq_impl,
        'ff->?': numbers.real_eq_impl,
        'dd->?': numbers.real_eq_impl,
        'FF->?': npyfuncs.np_complex_eq_impl,
        'DD->?': npyfuncs.np_complex_eq_impl,
    }
    if numpy_version >= (1, 25):
        db[np.equal].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('=='),
            'Qq->?': numbers.int_unsigned_signed_cmp('==')})

    db[np.logical_and] = {
        '??->?': npyfuncs.np_logical_and_impl,
        'bb->?': npyfuncs.np_logical_and_impl,
        'BB->?': npyfuncs.np_logical_and_impl,
        'hh->?': npyfuncs.np_logical_and_impl,
        'HH->?': npyfuncs.np_logical_and_impl,
        'ii->?': npyfuncs.np_logical_and_impl,
        'II->?': npyfuncs.np_logical_and_impl,
        'll->?': npyfuncs.np_logical_and_impl,
        'LL->?': npyfuncs.np_logical_and_impl,
        'qq->?': npyfuncs.np_logical_and_impl,
        'QQ->?': npyfuncs.np_logical_and_impl,
        'ff->?': npyfuncs.np_logical_and_impl,
        'dd->?': npyfuncs.np_logical_and_impl,
        'FF->?': npyfuncs.np_complex_logical_and_impl,
        'DD->?': npyfuncs.np_complex_logical_and_impl,
    }

    db[np.logical_or] = {
        '??->?': npyfuncs.np_logical_or_impl,
        'bb->?': npyfuncs.np_logical_or_impl,
        'BB->?': npyfuncs.np_logical_or_impl,
        'hh->?': npyfuncs.np_logical_or_impl,
        'HH->?': npyfuncs.np_logical_or_impl,
        'ii->?': npyfuncs.np_logical_or_impl,
        'II->?': npyfuncs.np_logical_or_impl,
        'll->?': npyfuncs.np_logical_or_impl,
        'LL->?': npyfuncs.np_logical_or_impl,
        'qq->?': npyfuncs.np_logical_or_impl,
        'QQ->?': npyfuncs.np_logical_or_impl,
        'ff->?': npyfuncs.np_logical_or_impl,
        'dd->?': npyfuncs.np_logical_or_impl,
        'FF->?': npyfuncs.np_complex_logical_or_impl,
        'DD->?': npyfuncs.np_complex_logical_or_impl,
    }

    db[np.logical_xor] = {
        '??->?': npyfuncs.np_logical_xor_impl,
        'bb->?': npyfuncs.np_logical_xor_impl,
        'BB->?': npyfuncs.np_logical_xor_impl,
        'hh->?': npyfuncs.np_logical_xor_impl,
        'HH->?': npyfuncs.np_logical_xor_impl,
        'ii->?': npyfuncs.np_logical_xor_impl,
        'II->?': npyfuncs.np_logical_xor_impl,
        'll->?': npyfuncs.np_logical_xor_impl,
        'LL->?': npyfuncs.np_logical_xor_impl,
        'qq->?': npyfuncs.np_logical_xor_impl,
        'QQ->?': npyfuncs.np_logical_xor_impl,
        'ff->?': npyfuncs.np_logical_xor_impl,
        'dd->?': npyfuncs.np_logical_xor_impl,
        'FF->?': npyfuncs.np_complex_logical_xor_impl,
        'DD->?': npyfuncs.np_complex_logical_xor_impl,
    }

    db[np.logical_not] = {
        '?->?': npyfuncs.np_logical_not_impl,
        'b->?': npyfuncs.np_logical_not_impl,
        'B->?': npyfuncs.np_logical_not_impl,
        'h->?': npyfuncs.np_logical_not_impl,
        'H->?': npyfuncs.np_logical_not_impl,
        'i->?': npyfuncs.np_logical_not_impl,
        'I->?': npyfuncs.np_logical_not_impl,
        'l->?': npyfuncs.np_logical_not_impl,
        'L->?': npyfuncs.np_logical_not_impl,
        'q->?': npyfuncs.np_logical_not_impl,
        'Q->?': npyfuncs.np_logical_not_impl,
        'f->?': npyfuncs.np_logical_not_impl,
        'd->?': npyfuncs.np_logical_not_impl,
        'F->?': npyfuncs.np_complex_logical_not_impl,
        'D->?': npyfuncs.np_complex_logical_not_impl,
    }

    db[np.maximum] = {
        '??->?': npyfuncs.np_logical_or_impl,
        'bb->b': npyfuncs.np_int_smax_impl,
        'BB->B': npyfuncs.np_int_umax_impl,
        'hh->h': npyfuncs.np_int_smax_impl,
        'HH->H': npyfuncs.np_int_umax_impl,
        'ii->i': npyfuncs.np_int_smax_impl,
        'II->I': npyfuncs.np_int_umax_impl,
        'll->l': npyfuncs.np_int_smax_impl,
        'LL->L': npyfuncs.np_int_umax_impl,
        'qq->q': npyfuncs.np_int_smax_impl,
        'QQ->Q': npyfuncs.np_int_umax_impl,
        'ff->f': npyfuncs.np_real_maximum_impl,
        'dd->d': npyfuncs.np_real_maximum_impl,
        'FF->F': npyfuncs.np_complex_maximum_impl,
        'DD->D': npyfuncs.np_complex_maximum_impl,
    }

    db[np.minimum] = {
        '??->?': npyfuncs.np_logical_and_impl,
        'bb->b': npyfuncs.np_int_smin_impl,
        'BB->B': npyfuncs.np_int_umin_impl,
        'hh->h': npyfuncs.np_int_smin_impl,
        'HH->H': npyfuncs.np_int_umin_impl,
        'ii->i': npyfuncs.np_int_smin_impl,
        'II->I': npyfuncs.np_int_umin_impl,
        'll->l': npyfuncs.np_int_smin_impl,
        'LL->L': npyfuncs.np_int_umin_impl,
        'qq->q': npyfuncs.np_int_smin_impl,
        'QQ->Q': npyfuncs.np_int_umin_impl,
        'ff->f': npyfuncs.np_real_minimum_impl,
        'dd->d': npyfuncs.np_real_minimum_impl,
        'FF->F': npyfuncs.np_complex_minimum_impl,
        'DD->D': npyfuncs.np_complex_minimum_impl,
    }

    db[np.fmax] = {
        '??->?': npyfuncs.np_logical_or_impl,
        'bb->b': npyfuncs.np_int_smax_impl,
        'BB->B': npyfuncs.np_int_umax_impl,
        'hh->h': npyfuncs.np_int_smax_impl,
        'HH->H': npyfuncs.np_int_umax_impl,
        'ii->i': npyfuncs.np_int_smax_impl,
        'II->I': npyfuncs.np_int_umax_impl,
        'll->l': npyfuncs.np_int_smax_impl,
        'LL->L': npyfuncs.np_int_umax_impl,
        'qq->q': npyfuncs.np_int_smax_impl,
        'QQ->Q': npyfuncs.np_int_umax_impl,
        'ff->f': npyfuncs.np_real_fmax_impl,
        'dd->d': npyfuncs.np_real_fmax_impl,
        'FF->F': npyfuncs.np_complex_fmax_impl,
        'DD->D': npyfuncs.np_complex_fmax_impl,
    }

    db[np.fmin] = {
        '??->?': npyfuncs.np_logical_and_impl,
        'bb->b': npyfuncs.np_int_smin_impl,
        'BB->B': npyfuncs.np_int_umin_impl,
        'hh->h': npyfuncs.np_int_smin_impl,
        'HH->H': npyfuncs.np_int_umin_impl,
        'ii->i': npyfuncs.np_int_smin_impl,
        'II->I': npyfuncs.np_int_umin_impl,
        'll->l': npyfuncs.np_int_smin_impl,
        'LL->L': npyfuncs.np_int_umin_impl,
        'qq->q': npyfuncs.np_int_smin_impl,
        'QQ->Q': npyfuncs.np_int_umin_impl,
        'ff->f': npyfuncs.np_real_fmin_impl,
        'dd->d': npyfuncs.np_real_fmin_impl,
        'FF->F': npyfuncs.np_complex_fmin_impl,
        'DD->D': npyfuncs.np_complex_fmin_impl,
    }

    db[np.bitwise_and] = {
        '??->?': numbers.int_and_impl,
        'bb->b': numbers.int_and_impl,
        'BB->B': numbers.int_and_impl,
        'hh->h': numbers.int_and_impl,
        'HH->H': numbers.int_and_impl,
        'ii->i': numbers.int_and_impl,
        'II->I': numbers.int_and_impl,
        'll->l': numbers.int_and_impl,
        'LL->L': numbers.int_and_impl,
        'qq->q': numbers.int_and_impl,
        'QQ->Q': numbers.int_and_impl,
    }

    db[np.bitwise_or] = {
        '??->?': numbers.int_or_impl,
        'bb->b': numbers.int_or_impl,
        'BB->B': numbers.int_or_impl,
        'hh->h': numbers.int_or_impl,
        'HH->H': numbers.int_or_impl,
        'ii->i': numbers.int_or_impl,
        'II->I': numbers.int_or_impl,
        'll->l': numbers.int_or_impl,
        'LL->L': numbers.int_or_impl,
        'qq->q': numbers.int_or_impl,
        'QQ->Q': numbers.int_or_impl,
    }

    db[np.bitwise_xor] = {
        '??->?': numbers.int_xor_impl,
        'bb->b': numbers.int_xor_impl,
        'BB->B': numbers.int_xor_impl,
        'hh->h': numbers.int_xor_impl,
        'HH->H': numbers.int_xor_impl,
        'ii->i': numbers.int_xor_impl,
        'II->I': numbers.int_xor_impl,
        'll->l': numbers.int_xor_impl,
        'LL->L': numbers.int_xor_impl,
        'qq->q': numbers.int_xor_impl,
        'QQ->Q': numbers.int_xor_impl,
    }

    db[np.invert] = {
        '?->?': numbers.int_invert_impl,
        'b->b': numbers.int_invert_impl,
        'B->B': numbers.int_invert_impl,
        'h->h': numbers.int_invert_impl,
        'H->H': numbers.int_invert_impl,
        'i->i': numbers.int_invert_impl,
        'I->I': numbers.int_invert_impl,
        'l->l': numbers.int_invert_impl,
        'L->L': numbers.int_invert_impl,
        'q->q': numbers.int_invert_impl,
        'Q->Q': numbers.int_invert_impl,
    }

    db[np.left_shift] = {
        'bb->b': numbers.int_shl_impl,
        'BB->B': numbers.int_shl_impl,
        'hh->h': numbers.int_shl_impl,
        'HH->H': numbers.int_shl_impl,
        'ii->i': numbers.int_shl_impl,
        'II->I': numbers.int_shl_impl,
        'll->l': numbers.int_shl_impl,
        'LL->L': numbers.int_shl_impl,
        'qq->q': numbers.int_shl_impl,
        'QQ->Q': numbers.int_shl_impl,
    }

    db[np.right_shift] = {
        'bb->b': numbers.int_shr_impl,
        'BB->B': numbers.int_shr_impl,
        'hh->h': numbers.int_shr_impl,
        'HH->H': numbers.int_shr_impl,
        'ii->i': numbers.int_shr_impl,
        'II->I': numbers.int_shr_impl,
        'll->l': numbers.int_shr_impl,
        'LL->L': numbers.int_shr_impl,
        'qq->q': numbers.int_shr_impl,
        'QQ->Q': numbers.int_shr_impl,
    }

    db[np.log] = {
        'f->f': np_real_log_impl,
        'd->d': np_real_log_impl,
        'F->F': npyfuncs.np_complex_log_impl,
        'D->D': npyfuncs.np_complex_log_impl,
    }

    db[np.log2] = {
        'f->f': np_real_log2_impl,
        'd->d': np_real_log2_impl,
        'F->F': npyfuncs.np_complex_log2_impl,
        'D->D': npyfuncs.np_complex_log2_impl,
    }

    db[np.log10] = {
        'f->f': np_real_log10_impl,
        'd->d': np_real_log10_impl,
        'F->F': npyfuncs.np_complex_log10_impl,
        'D->D': npyfuncs.np_complex_log10_impl,
    }

    return db
