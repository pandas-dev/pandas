"""This file contains information on how to translate different ufuncs
into numba. It is a database of different ufuncs and how each of its
loops maps to a function that implements the inner kernel of that ufunc
(the inner kernel being the per-element function).

Use the function get_ufunc_info to get the information related to the
ufunc
"""


import numpy as np
import sys

# this is lazily initialized to avoid circular imports
IS_WIN32 = sys.platform.startswith('win32')
numpy_version = tuple(map(int, np.__version__.split('.')[:2]))
_ufunc_db = None


def _lazy_init_db():
    global _ufunc_db

    if _ufunc_db is None:
        _ufunc_db = {}
        _fill_ufunc_db(_ufunc_db)


def get_ufuncs():
    """obtain a list of supported ufuncs in the db"""
    _lazy_init_db()
    return _ufunc_db.keys()


def get_ufunc_info(ufunc_key):
    """get the lowering information for the ufunc with key ufunc_key.

    The lowering information is a dictionary that maps from a numpy
    loop string (as given by the ufunc types attribute) to a function
    that handles code generation for a scalar version of the ufunc
    (that is, generates the "per element" operation").

    raises a KeyError if the ufunc is not in the ufunc_db
    """
    _lazy_init_db()
    return _ufunc_db[ufunc_key]


def _fill_ufunc_db(ufunc_db):
    # some of these imports would cause a problem of circular
    # imports if done at global scope when importing the numba
    # module.
    from numba.np import npyfuncs
    from numba.np.math import cmathimpl, mathimpl, numbers
    from numba.np.numpy_support import numpy_version

    ufunc_db[np.isnat] = {
        # datetime & timedelta
        'M->?': npyfuncs.np_datetime_isnat_impl,
        'm->?': npyfuncs.np_datetime_isnat_impl,
    }

    ufunc_db[np.negative] = {
        '?->?': numbers.int_invert_impl,
        'b->b': numbers.int_negate_impl,
        'B->B': numbers.int_negate_impl,
        'h->h': numbers.int_negate_impl,
        'H->H': numbers.int_negate_impl,
        'i->i': numbers.int_negate_impl,
        'I->I': numbers.int_negate_impl,
        'l->l': numbers.int_negate_impl,
        'L->L': numbers.int_negate_impl,
        'q->q': numbers.int_negate_impl,
        'Q->Q': numbers.int_negate_impl,
        'f->f': numbers.real_negate_impl,
        'd->d': numbers.real_negate_impl,
        'F->F': numbers.complex_negate_impl,
        'D->D': numbers.complex_negate_impl,
    }

    ufunc_db[np.positive] = {
        '?->?': numbers.int_positive_impl,
        'b->b': numbers.int_positive_impl,
        'B->B': numbers.int_positive_impl,
        'h->h': numbers.int_positive_impl,
        'H->H': numbers.int_positive_impl,
        'i->i': numbers.int_positive_impl,
        'I->I': numbers.int_positive_impl,
        'l->l': numbers.int_positive_impl,
        'L->L': numbers.int_positive_impl,
        'q->q': numbers.int_positive_impl,
        'Q->Q': numbers.int_positive_impl,
        'f->f': numbers.real_positive_impl,
        'd->d': numbers.real_positive_impl,
        'F->F': numbers.complex_positive_impl,
        'D->D': numbers.complex_positive_impl,
    }

    ufunc_db[np.absolute] = {
        '?->?': numbers.int_abs_impl,
        'b->b': numbers.int_abs_impl,
        'B->B': numbers.uint_abs_impl,
        'h->h': numbers.int_abs_impl,
        'H->H': numbers.uint_abs_impl,
        'i->i': numbers.int_abs_impl,
        'I->I': numbers.uint_abs_impl,
        'l->l': numbers.int_abs_impl,
        'L->L': numbers.uint_abs_impl,
        'q->q': numbers.int_abs_impl,
        'Q->Q': numbers.uint_abs_impl,
        'f->f': numbers.real_abs_impl,
        'd->d': numbers.real_abs_impl,
        'F->f': numbers.complex_abs_impl,
        'D->d': numbers.complex_abs_impl,
    }

    ufunc_db[np.sign] = {
        'b->b': numbers.int_sign_impl,
        'B->B': numbers.int_sign_impl,
        'h->h': numbers.int_sign_impl,
        'H->H': numbers.int_sign_impl,
        'i->i': numbers.int_sign_impl,
        'I->I': numbers.int_sign_impl,
        'l->l': numbers.int_sign_impl,
        'L->L': numbers.int_sign_impl,
        'q->q': numbers.int_sign_impl,
        'Q->Q': numbers.int_sign_impl,
        'f->f': numbers.real_sign_impl,
        'd->d': numbers.real_sign_impl,
        'F->F': npyfuncs.np_complex_sign_impl,
        'D->D': npyfuncs.np_complex_sign_impl,
    }

    ufunc_db[np.add] = {
        '??->?': numbers.int_or_impl,
        'bb->b': numbers.int_add_impl,
        'BB->B': numbers.int_add_impl,
        'hh->h': numbers.int_add_impl,
        'HH->H': numbers.int_add_impl,
        'ii->i': numbers.int_add_impl,
        'II->I': numbers.int_add_impl,
        'll->l': numbers.int_add_impl,
        'LL->L': numbers.int_add_impl,
        'qq->q': numbers.int_add_impl,
        'QQ->Q': numbers.int_add_impl,
        'ff->f': numbers.real_add_impl,
        'dd->d': numbers.real_add_impl,
        'FF->F': numbers.complex_add_impl,
        'DD->D': numbers.complex_add_impl,
    }

    ufunc_db[np.subtract] = {
        '??->?': numbers.int_xor_impl,
        'bb->b': numbers.int_sub_impl,
        'BB->B': numbers.int_sub_impl,
        'hh->h': numbers.int_sub_impl,
        'HH->H': numbers.int_sub_impl,
        'ii->i': numbers.int_sub_impl,
        'II->I': numbers.int_sub_impl,
        'll->l': numbers.int_sub_impl,
        'LL->L': numbers.int_sub_impl,
        'qq->q': numbers.int_sub_impl,
        'QQ->Q': numbers.int_sub_impl,
        'ff->f': numbers.real_sub_impl,
        'dd->d': numbers.real_sub_impl,
        'FF->F': numbers.complex_sub_impl,
        'DD->D': numbers.complex_sub_impl,
    }

    ufunc_db[np.multiply] = {
        '??->?': numbers.int_and_impl,
        'bb->b': numbers.int_mul_impl,
        'BB->B': numbers.int_mul_impl,
        'hh->h': numbers.int_mul_impl,
        'HH->H': numbers.int_mul_impl,
        'ii->i': numbers.int_mul_impl,
        'II->I': numbers.int_mul_impl,
        'll->l': numbers.int_mul_impl,
        'LL->L': numbers.int_mul_impl,
        'qq->q': numbers.int_mul_impl,
        'QQ->Q': numbers.int_mul_impl,
        'ff->f': numbers.real_mul_impl,
        'dd->d': numbers.real_mul_impl,
        'FF->F': numbers.complex_mul_impl,
        'DD->D': numbers.complex_mul_impl,
    }

    if np.divide != np.true_divide:
        ufunc_db[np.divide] = {
            'bb->b': npyfuncs.np_int_sdiv_impl,
            'BB->B': npyfuncs.np_int_udiv_impl,
            'hh->h': npyfuncs.np_int_sdiv_impl,
            'HH->H': npyfuncs.np_int_udiv_impl,
            'ii->i': npyfuncs.np_int_sdiv_impl,
            'II->I': npyfuncs.np_int_udiv_impl,
            'll->l': npyfuncs.np_int_sdiv_impl,
            'LL->L': npyfuncs.np_int_udiv_impl,
            'qq->q': npyfuncs.np_int_sdiv_impl,
            'QQ->Q': npyfuncs.np_int_udiv_impl,
            'ff->f': npyfuncs.np_real_div_impl,
            'dd->d': npyfuncs.np_real_div_impl,
            'FF->F': npyfuncs.np_complex_div_impl,
            'DD->D': npyfuncs.np_complex_div_impl,
        }

    ufunc_db[np.true_divide] = {
        'bb->d': npyfuncs.np_int_truediv_impl,
        'BB->d': npyfuncs.np_int_truediv_impl,
        'hh->d': npyfuncs.np_int_truediv_impl,
        'HH->d': npyfuncs.np_int_truediv_impl,
        'ii->d': npyfuncs.np_int_truediv_impl,
        'II->d': npyfuncs.np_int_truediv_impl,
        'll->d': npyfuncs.np_int_truediv_impl,
        'LL->d': npyfuncs.np_int_truediv_impl,
        'qq->d': npyfuncs.np_int_truediv_impl,
        'QQ->d': npyfuncs.np_int_truediv_impl,
        'ff->f': npyfuncs.np_real_div_impl,
        'dd->d': npyfuncs.np_real_div_impl,
        'FF->F': npyfuncs.np_complex_div_impl,
        'DD->D': npyfuncs.np_complex_div_impl,
    }

    ufunc_db[np.floor_divide] = {
        'bb->b': npyfuncs.np_int_sdiv_impl,
        'BB->B': npyfuncs.np_int_udiv_impl,
        'hh->h': npyfuncs.np_int_sdiv_impl,
        'HH->H': npyfuncs.np_int_udiv_impl,
        'ii->i': npyfuncs.np_int_sdiv_impl,
        'II->I': npyfuncs.np_int_udiv_impl,
        'll->l': npyfuncs.np_int_sdiv_impl,
        'LL->L': npyfuncs.np_int_udiv_impl,
        'qq->q': npyfuncs.np_int_sdiv_impl,
        'QQ->Q': npyfuncs.np_int_udiv_impl,
        'ff->f': npyfuncs.np_real_floor_div_impl,
        'dd->d': npyfuncs.np_real_floor_div_impl,
    }

    ufunc_db[np.remainder] = {
        'bb->b': npyfuncs.np_int_srem_impl,
        'BB->B': npyfuncs.np_int_urem_impl,
        'hh->h': npyfuncs.np_int_srem_impl,
        'HH->H': npyfuncs.np_int_urem_impl,
        'ii->i': npyfuncs.np_int_srem_impl,
        'II->I': npyfuncs.np_int_urem_impl,
        'll->l': npyfuncs.np_int_srem_impl,
        'LL->L': npyfuncs.np_int_urem_impl,
        'qq->q': npyfuncs.np_int_srem_impl,
        'QQ->Q': npyfuncs.np_int_urem_impl,
        'ff->f': npyfuncs.np_real_mod_impl,
        'dd->d': npyfuncs.np_real_mod_impl,
    }

    ufunc_db[np.divmod] = {
        'bb->bb': npyfuncs.np_int_sdivrem_impl,
        'BB->BB': npyfuncs.np_int_udivrem_impl,
        'hh->hh': npyfuncs.np_int_sdivrem_impl,
        'HH->HH': npyfuncs.np_int_udivrem_impl,
        'ii->ii': npyfuncs.np_int_sdivrem_impl,
        'II->II': npyfuncs.np_int_udivrem_impl,
        'll->ll': npyfuncs.np_int_sdivrem_impl,
        'LL->LL': npyfuncs.np_int_udivrem_impl,
        'qq->qq': npyfuncs.np_int_sdivrem_impl,
        'QQ->QQ': npyfuncs.np_int_udivrem_impl,
        'ff->ff': npyfuncs.np_real_divmod_impl,
        'dd->dd': npyfuncs.np_real_divmod_impl,
    }

    ufunc_db[np.fmod] = {
        'bb->b': npyfuncs.np_int_fmod_impl,
        'BB->B': npyfuncs.np_int_fmod_impl,
        'hh->h': npyfuncs.np_int_fmod_impl,
        'HH->H': npyfuncs.np_int_fmod_impl,
        'ii->i': npyfuncs.np_int_fmod_impl,
        'II->I': npyfuncs.np_int_fmod_impl,
        'll->l': npyfuncs.np_int_fmod_impl,
        'LL->L': npyfuncs.np_int_fmod_impl,
        'qq->q': npyfuncs.np_int_fmod_impl,
        'QQ->Q': npyfuncs.np_int_fmod_impl,
        'ff->f': npyfuncs.np_real_fmod_impl,
        'dd->d': npyfuncs.np_real_fmod_impl,
    }

    ufunc_db[np.logaddexp] = {
        'ff->f': npyfuncs.np_real_logaddexp_impl,
        'dd->d': npyfuncs.np_real_logaddexp_impl,
    }

    ufunc_db[np.logaddexp2] = {
        'ff->f': npyfuncs.np_real_logaddexp2_impl,
        'dd->d': npyfuncs.np_real_logaddexp2_impl,
    }

    ufunc_db[np.power] = {
        'bb->b': numbers.int_power_impl,
        'BB->B': numbers.int_power_impl,
        'hh->h': numbers.int_power_impl,
        'HH->H': numbers.int_power_impl,
        'ii->i': numbers.int_power_impl,
        'II->I': numbers.int_power_impl,
        'll->l': numbers.int_power_impl,
        'LL->L': numbers.int_power_impl,
        'qq->q': numbers.int_power_impl,
        'QQ->Q': numbers.int_power_impl,
        # XXX we would like to use `int_power_impl` for real ** integer
        # as well (for better performance), but the current ufunc typing
        # rules forbid that
        'ff->f': numbers.real_power_impl,
        'dd->d': numbers.real_power_impl,
        'FF->F': npyfuncs.np_complex_power_impl,
        'DD->D': npyfuncs.np_complex_power_impl,
    }

    ufunc_db[np.float_power] = {
        'ff->f': npyfuncs.real_float_power_impl,
        'dd->d': npyfuncs.real_float_power_impl,
        'FF->F': npyfuncs.np_complex_float_power_impl,
        'DD->D': npyfuncs.np_complex_float_power_impl,
    }

    ufunc_db[np.gcd] = {
        'bb->b': npyfuncs.np_gcd_impl,
        'BB->B': npyfuncs.np_gcd_impl,
        'hh->h': npyfuncs.np_gcd_impl,
        'HH->H': npyfuncs.np_gcd_impl,
        'ii->i': npyfuncs.np_gcd_impl,
        'II->I': npyfuncs.np_gcd_impl,
        'll->l': npyfuncs.np_gcd_impl,
        'LL->L': npyfuncs.np_gcd_impl,
        'qq->q': npyfuncs.np_gcd_impl,
        'QQ->Q': npyfuncs.np_gcd_impl,
    }

    ufunc_db[np.lcm] = {
        'bb->b': npyfuncs.np_lcm_impl,
        'BB->B': npyfuncs.np_lcm_impl,
        'hh->h': npyfuncs.np_lcm_impl,
        'HH->H': npyfuncs.np_lcm_impl,
        'ii->i': npyfuncs.np_lcm_impl,
        'II->I': npyfuncs.np_lcm_impl,
        'll->l': npyfuncs.np_lcm_impl,
        'LL->L': npyfuncs.np_lcm_impl,
        'qq->q': npyfuncs.np_lcm_impl,
        'QQ->Q': npyfuncs.np_lcm_impl,
    }

    ufunc_db[np.rint] = {
        'f->f': npyfuncs.np_real_rint_impl,
        'd->d': npyfuncs.np_real_rint_impl,
        'F->F': npyfuncs.np_complex_rint_impl,
        'D->D': npyfuncs.np_complex_rint_impl,
    }

    ufunc_db[np.conjugate] = {
        'b->b': numbers.real_conjugate_impl,
        'B->B': numbers.real_conjugate_impl,
        'h->h': numbers.real_conjugate_impl,
        'H->H': numbers.real_conjugate_impl,
        'i->i': numbers.real_conjugate_impl,
        'I->I': numbers.real_conjugate_impl,
        'l->l': numbers.real_conjugate_impl,
        'L->L': numbers.real_conjugate_impl,
        'q->q': numbers.real_conjugate_impl,
        'Q->Q': numbers.real_conjugate_impl,
        'f->f': numbers.real_conjugate_impl,
        'd->d': numbers.real_conjugate_impl,
        'F->F': numbers.complex_conjugate_impl,
        'D->D': numbers.complex_conjugate_impl,
    }

    ufunc_db[np.exp] = {
        'f->f': npyfuncs.np_real_exp_impl,
        'd->d': npyfuncs.np_real_exp_impl,
        'F->F': npyfuncs.np_complex_exp_impl,
        'D->D': npyfuncs.np_complex_exp_impl,
    }

    ufunc_db[np.exp2] = {
        'f->f': npyfuncs.np_real_exp2_impl,
        'd->d': npyfuncs.np_real_exp2_impl,
        'F->F': npyfuncs.np_complex_exp2_impl,
        'D->D': npyfuncs.np_complex_exp2_impl,
    }

    ufunc_db[np.log] = {
        'f->f': npyfuncs.np_real_log_impl,
        'd->d': npyfuncs.np_real_log_impl,
        'F->F': npyfuncs.np_complex_log_impl,
        'D->D': npyfuncs.np_complex_log_impl,
    }

    ufunc_db[np.log2] = {
        'f->f': npyfuncs.np_real_log2_impl,
        'd->d': npyfuncs.np_real_log2_impl,
        'F->F': npyfuncs.np_complex_log2_impl,
        'D->D': npyfuncs.np_complex_log2_impl,
    }

    ufunc_db[np.log10] = {
        'f->f': npyfuncs.np_real_log10_impl,
        'd->d': npyfuncs.np_real_log10_impl,
        'F->F': npyfuncs.np_complex_log10_impl,
        'D->D': npyfuncs.np_complex_log10_impl,
    }

    ufunc_db[np.expm1] = {
        'f->f': npyfuncs.np_real_expm1_impl,
        'd->d': npyfuncs.np_real_expm1_impl,
        'F->F': npyfuncs.np_complex_expm1_impl,
        'D->D': npyfuncs.np_complex_expm1_impl,
    }

    ufunc_db[np.log1p] = {
        'f->f': npyfuncs.np_real_log1p_impl,
        'd->d': npyfuncs.np_real_log1p_impl,
        'F->F': npyfuncs.np_complex_log1p_impl,
        'D->D': npyfuncs.np_complex_log1p_impl,
    }

    ufunc_db[np.sqrt] = {
        'f->f': npyfuncs.np_real_sqrt_impl,
        'd->d': npyfuncs.np_real_sqrt_impl,
        'F->F': npyfuncs.np_complex_sqrt_impl,
        'D->D': npyfuncs.np_complex_sqrt_impl,
    }

    ufunc_db[np.square] = {
        'b->b': npyfuncs.np_int_square_impl,
        'B->B': npyfuncs.np_int_square_impl,
        'h->h': npyfuncs.np_int_square_impl,
        'H->H': npyfuncs.np_int_square_impl,
        'i->i': npyfuncs.np_int_square_impl,
        'I->I': npyfuncs.np_int_square_impl,
        'l->l': npyfuncs.np_int_square_impl,
        'L->L': npyfuncs.np_int_square_impl,
        'q->q': npyfuncs.np_int_square_impl,
        'Q->Q': npyfuncs.np_int_square_impl,
        'f->f': npyfuncs.np_real_square_impl,
        'd->d': npyfuncs.np_real_square_impl,
        'F->F': npyfuncs.np_complex_square_impl,
        'D->D': npyfuncs.np_complex_square_impl,
    }

    ufunc_db[np.cbrt] = {
        'f->f': npyfuncs.np_real_cbrt_impl,
        'd->d': npyfuncs.np_real_cbrt_impl,
    }

    ufunc_db[np.reciprocal] = {
        'b->b': npyfuncs.np_int_reciprocal_impl,
        'B->B': npyfuncs.np_int_reciprocal_impl,
        'h->h': npyfuncs.np_int_reciprocal_impl,
        'H->H': npyfuncs.np_int_reciprocal_impl,
        'i->i': npyfuncs.np_int_reciprocal_impl,
        'I->I': npyfuncs.np_int_reciprocal_impl,
        'l->l': npyfuncs.np_int_reciprocal_impl,
        'L->L': npyfuncs.np_int_reciprocal_impl,
        'q->q': npyfuncs.np_int_reciprocal_impl,
        'Q->Q': npyfuncs.np_int_reciprocal_impl,
        'f->f': npyfuncs.np_real_reciprocal_impl,
        'd->d': npyfuncs.np_real_reciprocal_impl,
        'F->F': npyfuncs.np_complex_reciprocal_impl,
        'D->D': npyfuncs.np_complex_reciprocal_impl,
    }

    ufunc_db[np.sin] = {
        'f->f': npyfuncs.np_real_sin_impl,
        'd->d': npyfuncs.np_real_sin_impl,
        'F->F': npyfuncs.np_complex_sin_impl,
        'D->D': npyfuncs.np_complex_sin_impl,
    }

    ufunc_db[np.cos] = {
        'f->f': npyfuncs.np_real_cos_impl,
        'd->d': npyfuncs.np_real_cos_impl,
        'F->F': npyfuncs.np_complex_cos_impl,
        'D->D': npyfuncs.np_complex_cos_impl,
    }

    tan_impl = cmathimpl.tan_impl

    ufunc_db[np.tan] = {
        'f->f': npyfuncs.np_real_tan_impl,
        'd->d': npyfuncs.np_real_tan_impl,
        'F->F': tan_impl,
        'D->D': tan_impl,
    }

    arcsin_impl = cmathimpl.asin_impl

    ufunc_db[np.arcsin] = {
        'f->f': npyfuncs.np_real_asin_impl,
        'd->d': npyfuncs.np_real_asin_impl,
        'F->F': arcsin_impl,
        'D->D': arcsin_impl,
    }

    ufunc_db[np.arccos] = {
        'f->f': npyfuncs.np_real_acos_impl,
        'd->d': npyfuncs.np_real_acos_impl,
        'F->F': cmathimpl.acos_impl,
        'D->D': cmathimpl.acos_impl,
    }

    arctan_impl = cmathimpl.atan_impl

    ufunc_db[np.arctan] = {
        'f->f': npyfuncs.np_real_atan_impl,
        'd->d': npyfuncs.np_real_atan_impl,
        'F->F': arctan_impl,
        'D->D': arctan_impl,
    }

    ufunc_db[np.arctan2] = {
        'ff->f': npyfuncs.np_real_atan2_impl,
        'dd->d': npyfuncs.np_real_atan2_impl,
    }

    ufunc_db[np.hypot] = {
        'ff->f': npyfuncs.np_real_hypot_impl,
        'dd->d': npyfuncs.np_real_hypot_impl,
    }

    ufunc_db[np.sinh] = {
        'f->f': npyfuncs.np_real_sinh_impl,
        'd->d': npyfuncs.np_real_sinh_impl,
        'F->F': npyfuncs.np_complex_sinh_impl,
        'D->D': npyfuncs.np_complex_sinh_impl,
    }

    ufunc_db[np.cosh] = {
        'f->f': npyfuncs.np_real_cosh_impl,
        'd->d': npyfuncs.np_real_cosh_impl,
        'F->F': npyfuncs.np_complex_cosh_impl,
        'D->D': npyfuncs.np_complex_cosh_impl,
    }

    ufunc_db[np.tanh] = {
        'f->f': npyfuncs.np_real_tanh_impl,
        'd->d': npyfuncs.np_real_tanh_impl,
        'F->F': npyfuncs.np_complex_tanh_impl,
        'D->D': npyfuncs.np_complex_tanh_impl,
    }

    arcsinh_impl = cmathimpl.asinh_impl

    ufunc_db[np.arcsinh] = {
        'f->f': npyfuncs.np_real_asinh_impl,
        'd->d': npyfuncs.np_real_asinh_impl,
        'F->F': arcsinh_impl,
        'D->D': arcsinh_impl,
    }

    ufunc_db[np.arccosh] = {
        'f->f': npyfuncs.np_real_acosh_impl,
        'd->d': npyfuncs.np_real_acosh_impl,
        'F->F': npyfuncs.np_complex_acosh_impl,
        'D->D': npyfuncs.np_complex_acosh_impl,
    }

    arctanh_impl = cmathimpl.atanh_impl

    ufunc_db[np.arctanh] = {
        'f->f': npyfuncs.np_real_atanh_impl,
        'd->d': npyfuncs.np_real_atanh_impl,
        'F->F': arctanh_impl,
        'D->D': arctanh_impl,
    }

    ufunc_db[np.deg2rad] = {
        'f->f': mathimpl.radians_float_impl,
        'd->d': mathimpl.radians_float_impl,
    }

    ufunc_db[np.radians] = ufunc_db[np.deg2rad]

    ufunc_db[np.rad2deg] = {
        'f->f': mathimpl.degrees_float_impl,
        'd->d': mathimpl.degrees_float_impl,
    }

    ufunc_db[np.degrees] = ufunc_db[np.rad2deg]

    ufunc_db[np.floor] = {
        'f->f': npyfuncs.np_real_floor_impl,
        'd->d': npyfuncs.np_real_floor_impl,
    }
    if numpy_version >= (2, 1):
        ufunc_db[np.floor].update({
            '?->?': numbers.identity_impl,
            'b->b': numbers.identity_impl,
            'B->B': numbers.identity_impl,
            'h->h': numbers.identity_impl,
            'H->H': numbers.identity_impl,
            'i->i': numbers.identity_impl,
            'I->I': numbers.identity_impl,
            'l->l': numbers.identity_impl,
            'L->L': numbers.identity_impl,
            'q->q': numbers.identity_impl,
            'Q->Q': numbers.identity_impl,
        })

    ufunc_db[np.ceil] = {
        'f->f': npyfuncs.np_real_ceil_impl,
        'd->d': npyfuncs.np_real_ceil_impl,
    }
    if numpy_version >= (2, 1):
        ufunc_db[np.ceil].update({
            '?->?': numbers.identity_impl,
            'b->b': numbers.identity_impl,
            'B->B': numbers.identity_impl,
            'h->h': numbers.identity_impl,
            'H->H': numbers.identity_impl,
            'i->i': numbers.identity_impl,
            'I->I': numbers.identity_impl,
            'l->l': numbers.identity_impl,
            'L->L': numbers.identity_impl,
            'q->q': numbers.identity_impl,
            'Q->Q': numbers.identity_impl,
        })

    ufunc_db[np.trunc] = {
        'f->f': npyfuncs.np_real_trunc_impl,
        'd->d': npyfuncs.np_real_trunc_impl,
    }
    if numpy_version >= (2, 1):
        ufunc_db[np.trunc].update({
            '?->?': numbers.identity_impl,
            'b->b': numbers.identity_impl,
            'B->B': numbers.identity_impl,
            'h->h': numbers.identity_impl,
            'H->H': numbers.identity_impl,
            'i->i': numbers.identity_impl,
            'I->I': numbers.identity_impl,
            'l->l': numbers.identity_impl,
            'L->L': numbers.identity_impl,
            'q->q': numbers.identity_impl,
            'Q->Q': numbers.identity_impl,
        })

    ufunc_db[np.fabs] = {
        'f->f': npyfuncs.np_real_fabs_impl,
        'd->d': npyfuncs.np_real_fabs_impl,
    }

    # logical ufuncs
    ufunc_db[np.greater] = {
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
        ufunc_db[np.greater].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('>'),
            'Qq->?': numbers.int_unsigned_signed_cmp('>')})

    ufunc_db[np.greater_equal] = {
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
        ufunc_db[np.greater_equal].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('>='),
            'Qq->?': numbers.int_unsigned_signed_cmp('>=')})

    ufunc_db[np.less] = {
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
        ufunc_db[np.less].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('<'),
            'Qq->?': numbers.int_unsigned_signed_cmp('<')})

    ufunc_db[np.less_equal] = {
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
        ufunc_db[np.less_equal].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('<='),
            'Qq->?': numbers.int_unsigned_signed_cmp('<=')})

    ufunc_db[np.not_equal] = {
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
        ufunc_db[np.not_equal].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('!='),
            'Qq->?': numbers.int_unsigned_signed_cmp('!=')})

    ufunc_db[np.equal] = {
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
        ufunc_db[np.equal].update({
            'qQ->?': numbers.int_signed_unsigned_cmp('=='),
            'Qq->?': numbers.int_unsigned_signed_cmp('==')})

    ufunc_db[np.logical_and] = {
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

    ufunc_db[np.logical_or] = {
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

    ufunc_db[np.logical_xor] = {
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

    ufunc_db[np.logical_not] = {
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

    ufunc_db[np.maximum] = {
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

    ufunc_db[np.minimum] = {
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

    ufunc_db[np.fmax] = {
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

    ufunc_db[np.fmin] = {
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

    # misc floating functions
    ufunc_db[np.isnan] = {
        'f->?': npyfuncs.np_real_isnan_impl,
        'd->?': npyfuncs.np_real_isnan_impl,
        'F->?': npyfuncs.np_complex_isnan_impl,
        'D->?': npyfuncs.np_complex_isnan_impl,
        # int8
        'b->?': npyfuncs.np_int_isnan_impl,
        'B->?': npyfuncs.np_int_isnan_impl,
        # int16
        'h->?': npyfuncs.np_int_isnan_impl,
        'H->?': npyfuncs.np_int_isnan_impl,
        # int32
        'i->?': npyfuncs.np_int_isnan_impl,
        'I->?': npyfuncs.np_int_isnan_impl,
        # int64
        'l->?': npyfuncs.np_int_isnan_impl,
        'L->?': npyfuncs.np_int_isnan_impl,
        # intp
        'q->?': npyfuncs.np_int_isnan_impl,
        'Q->?': npyfuncs.np_int_isnan_impl,
        # boolean
        '?->?': npyfuncs.np_int_isnan_impl,
        # datetime & timedelta
        'm->?': npyfuncs.np_datetime_isnat_impl,
        'M->?': npyfuncs.np_datetime_isnat_impl,
    }

    ufunc_db[np.isinf] = {
        'f->?': npyfuncs.np_real_isinf_impl,
        'd->?': npyfuncs.np_real_isinf_impl,
        'F->?': npyfuncs.np_complex_isinf_impl,
        'D->?': npyfuncs.np_complex_isinf_impl,
        # int8
        'b->?': npyfuncs.np_int_isinf_impl,
        'B->?': npyfuncs.np_int_isinf_impl,
        # int16
        'h->?': npyfuncs.np_int_isinf_impl,
        'H->?': npyfuncs.np_int_isinf_impl,
        # int32
        'i->?': npyfuncs.np_int_isinf_impl,
        'I->?': npyfuncs.np_int_isinf_impl,
        # int64
        'l->?': npyfuncs.np_int_isinf_impl,
        'L->?': npyfuncs.np_int_isinf_impl,
        # intp
        'q->?': npyfuncs.np_int_isinf_impl,
        'Q->?': npyfuncs.np_int_isinf_impl,
        # boolean
        '?->?': npyfuncs.np_int_isinf_impl,
        # datetime & timedelta
        'm->?': npyfuncs.np_int_isinf_impl,
        'M->?': npyfuncs.np_int_isinf_impl,
    }

    ufunc_db[np.isfinite] = {
        'f->?': npyfuncs.np_real_isfinite_impl,
        'd->?': npyfuncs.np_real_isfinite_impl,
        'F->?': npyfuncs.np_complex_isfinite_impl,
        'D->?': npyfuncs.np_complex_isfinite_impl,
        # int8
        'b->?': npyfuncs.np_int_isfinite_impl,
        'B->?': npyfuncs.np_int_isfinite_impl,
        # int16
        'h->?': npyfuncs.np_int_isfinite_impl,
        'H->?': npyfuncs.np_int_isfinite_impl,
        # int32
        'i->?': npyfuncs.np_int_isfinite_impl,
        'I->?': npyfuncs.np_int_isfinite_impl,
        # int64
        'l->?': npyfuncs.np_int_isfinite_impl,
        'L->?': npyfuncs.np_int_isfinite_impl,
        # intp
        'q->?': npyfuncs.np_int_isfinite_impl,
        'Q->?': npyfuncs.np_int_isfinite_impl,
        # boolean
        '?->?': npyfuncs.np_int_isfinite_impl,
        # datetime & timedelta
        'M->?': npyfuncs.np_datetime_isfinite_impl,
        'm->?': npyfuncs.np_datetime_isfinite_impl,
    }

    ufunc_db[np.signbit] = {
        'f->?': npyfuncs.np_real_signbit_impl,
        'd->?': npyfuncs.np_real_signbit_impl,
    }

    ufunc_db[np.copysign] = {
        'ff->f': npyfuncs.np_real_copysign_impl,
        'dd->d': npyfuncs.np_real_copysign_impl,
    }

    ufunc_db[np.nextafter] = {
        'ff->f': npyfuncs.np_real_nextafter_impl,
        'dd->d': npyfuncs.np_real_nextafter_impl,
    }

    ufunc_db[np.spacing] = {
        'f->f': npyfuncs.np_real_spacing_impl,
        'd->d': npyfuncs.np_real_spacing_impl,
    }

    ufunc_db[np.ldexp] = {
        'fi->f': npyfuncs.np_real_ldexp_impl,
        'fl->f': npyfuncs.np_real_ldexp_impl,
        'di->d': npyfuncs.np_real_ldexp_impl,
        'dl->d': npyfuncs.np_real_ldexp_impl,
    }
    if numpy_version >= (2, 0) and IS_WIN32:
        ufunc_db[np.ldexp]['fq->f'] = ufunc_db[np.ldexp].pop('fl->f')
        ufunc_db[np.ldexp]['dq->d'] = ufunc_db[np.ldexp].pop('dl->d')

    # bit twiddling functions
    ufunc_db[np.bitwise_and] = {
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

    ufunc_db[np.bitwise_or] = {
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

    ufunc_db[np.bitwise_xor] = {
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

    ufunc_db[np.invert] = {  # aka np.bitwise_not
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

    ufunc_db[np.left_shift] = {
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

    ufunc_db[np.right_shift] = {
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

    # Inject datetime64 support
    from numba.np import npdatetime
    ufunc_db[np.negative].update({
        'm->m': npdatetime.timedelta_neg_impl,
    })
    ufunc_db[np.positive].update({
        'm->m': npdatetime.timedelta_pos_impl,
    })
    ufunc_db[np.absolute].update({
        'm->m': npdatetime.timedelta_abs_impl,
    })
    ufunc_db[np.sign].update({
        'm->m': npdatetime.timedelta_sign_impl,
    })
    ufunc_db[np.add].update({
        'mm->m': npdatetime.timedelta_add_impl,
        'Mm->M': npdatetime.datetime_plus_timedelta,
        'mM->M': npdatetime.timedelta_plus_datetime,
    })
    ufunc_db[np.subtract].update({
        'mm->m': npdatetime.timedelta_sub_impl,
        'Mm->M': npdatetime.datetime_minus_timedelta,
        'MM->m': npdatetime.datetime_minus_datetime,
    })
    ufunc_db[np.multiply].update({
        'mq->m': npdatetime.timedelta_times_number,
        'md->m': npdatetime.timedelta_times_number,
        'qm->m': npdatetime.number_times_timedelta,
        'dm->m': npdatetime.number_times_timedelta,
    })
    if np.divide != np.true_divide:
        ufunc_db[np.divide].update({
            'mq->m': npdatetime.timedelta_over_number,
            'md->m': npdatetime.timedelta_over_number,
            'mm->d': npdatetime.timedelta_over_timedelta,
        })
    ufunc_db[np.true_divide].update({
        'mq->m': npdatetime.timedelta_over_number,
        'md->m': npdatetime.timedelta_over_number,
        'mm->d': npdatetime.timedelta_over_timedelta,
    })
    ufunc_db[np.floor_divide].update({
        'mq->m': npdatetime.timedelta_over_number,
        'md->m': npdatetime.timedelta_over_number,
    })

    ufunc_db[np.floor_divide].update({
        'mm->q': npdatetime.timedelta_floor_div_timedelta,
    })

    ufunc_db[np.equal].update({
        'MM->?': npdatetime.datetime_eq_datetime_impl,
        'mm->?': npdatetime.timedelta_eq_timedelta_impl,
    })
    ufunc_db[np.not_equal].update({
        'MM->?': npdatetime.datetime_ne_datetime_impl,
        'mm->?': npdatetime.timedelta_ne_timedelta_impl,
    })
    ufunc_db[np.less].update({
        'MM->?': npdatetime.datetime_lt_datetime_impl,
        'mm->?': npdatetime.timedelta_lt_timedelta_impl,
    })
    ufunc_db[np.less_equal].update({
        'MM->?': npdatetime.datetime_le_datetime_impl,
        'mm->?': npdatetime.timedelta_le_timedelta_impl,
    })
    ufunc_db[np.greater].update({
        'MM->?': npdatetime.datetime_gt_datetime_impl,
        'mm->?': npdatetime.timedelta_gt_timedelta_impl,
    })
    ufunc_db[np.greater_equal].update({
        'MM->?': npdatetime.datetime_ge_datetime_impl,
        'mm->?': npdatetime.timedelta_ge_timedelta_impl,
    })
    ufunc_db[np.maximum].update({
        'MM->M': npdatetime.datetime_maximum_impl,
        'mm->m': npdatetime.timedelta_maximum_impl,
    })
    ufunc_db[np.minimum].update({
        'MM->M': npdatetime.datetime_minimum_impl,
        'mm->m': npdatetime.timedelta_minimum_impl,
    })
    # there is no difference for datetime/timedelta in maximum/fmax
    # and minimum/fmin
    ufunc_db[np.fmax].update({
        'MM->M': npdatetime.datetime_fmax_impl,
        'mm->m': npdatetime.timedelta_fmax_impl,
    })
    ufunc_db[np.fmin].update({
        'MM->M': npdatetime.datetime_fmin_impl,
        'mm->m': npdatetime.timedelta_fmin_impl,
    })

    ufunc_db[np.remainder].update({
        'mm->m': npdatetime.timedelta_mod_timedelta,
    })
