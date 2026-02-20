// -*- c-mode -*-
/*********************************************************************
  Numexpr - Fast numerical array expression evaluator for NumPy.

      License: MIT
      Author:  See AUTHORS.txt

  See LICENSE.txt for details about copyright and rights to use.
**********************************************************************/

/* These #if blocks make it easier to query this file, without having
   to define every row function before #including it. */
#ifndef FUNC_FF
#define ELIDE_FUNC_FF
#define FUNC_FF(...)
#endif
FUNC_FF(FUNC_SQRT_FF,    "sqrt_ff",     sqrtf,  sqrtf2,  vsSqrt)
FUNC_FF(FUNC_SIN_FF,     "sin_ff",      sinf,   sinf2,   vsSin)
FUNC_FF(FUNC_COS_FF,     "cos_ff",      cosf,   cosf2,   vsCos)
FUNC_FF(FUNC_TAN_FF,     "tan_ff",      tanf,   tanf2,   vsTan)
FUNC_FF(FUNC_ARCSIN_FF,  "arcsin_ff",   asinf,  asinf2,  vsAsin)
FUNC_FF(FUNC_ARCCOS_FF,  "arccos_ff",   acosf,  acosf2,  vsAcos)
FUNC_FF(FUNC_ARCTAN_FF,  "arctan_ff",   atanf,  atanf2,  vsAtan)
FUNC_FF(FUNC_SINH_FF,    "sinh_ff",     sinhf,  sinhf2,  vsSinh)
FUNC_FF(FUNC_COSH_FF,    "cosh_ff",     coshf,  coshf2,  vsCosh)
FUNC_FF(FUNC_TANH_FF,    "tanh_ff",     tanhf,  tanhf2,  vsTanh)
FUNC_FF(FUNC_ARCSINH_FF, "arcsinh_ff",  asinhf, asinhf2, vsAsinh)
FUNC_FF(FUNC_ARCCOSH_FF, "arccosh_ff",  acoshf, acoshf2, vsAcosh)
FUNC_FF(FUNC_ARCTANH_FF, "arctanh_ff",  atanhf, atanhf2, vsAtanh)
FUNC_FF(FUNC_LOG_FF,     "log_ff",      logf,   logf2,   vsLn)
FUNC_FF(FUNC_LOG1P_FF,   "log1p_ff",    log1pf, log1pf2, vsLog1p)
FUNC_FF(FUNC_LOG10_FF,   "log10_ff",    log10f, log10f2, vsLog10)
FUNC_FF(FUNC_LOG2_FF,    "log2_ff",     log2f,  log2f2,  vsLog2)
FUNC_FF(FUNC_EXP_FF,     "exp_ff",      expf,   expf2,   vsExp)
FUNC_FF(FUNC_EXPM1_FF,   "expm1_ff",    expm1f, expm1f2, vsExpm1)
FUNC_FF(FUNC_ABS_FF,     "absolute_ff", fabsf,  fabsf2,  vsAbs)
FUNC_FF(FUNC_CONJ_FF,    "conjugate_ff",fconjf, fconjf2, vsConj)
FUNC_FF(FUNC_CEIL_FF,    "ceil_ff",     ceilf,  ceilf2,  vsCeil)
FUNC_FF(FUNC_FLOOR_FF,   "floor_ff",    floorf, floorf2, vsFloor)
FUNC_FF(FUNC_TRUNC_FF,   "trunc_ff",    truncf, truncf2, vsTrunc)
FUNC_FF(FUNC_SIGN_FF,   "sign_ff",    signf, signf2, vsSign)
//rint rounds to nearest even integer, matching NumPy (round doesn't)
FUNC_FF(FUNC_ROUND_FF,   "round_ff",    rintf,  rintf2,  vsRint)
FUNC_FF(FUNC_FF_LAST,    NULL,          NULL,   NULL,    NULL)
#ifdef ELIDE_FUNC_FF
#undef ELIDE_FUNC_FF
#undef FUNC_FF
#endif

#ifndef FUNC_FFF
#define ELIDE_FUNC_FFF
#define FUNC_FFF(...)
#endif
FUNC_FFF(FUNC_FMOD_FFF,    "fmod_fff",    fmodf,  fmodf2,  vsfmod)
FUNC_FFF(FUNC_ARCTAN2_FFF, "arctan2_fff", atan2f, atan2f2, vsAtan2)
FUNC_FFF(FUNC_HYPOT_FFF, "hypot_fff", hypotf, hypotf2, vsHypot)
FUNC_FFF(FUNC_NEXTAFTER_FFF, "nextafter_fff", nextafterf, nextafterf2, vsNextAfter)
FUNC_FFF(FUNC_COPYSIGN_FFF, "copysign_fff", copysignf, copysignf2, vsCopySign)
FUNC_FFF(FUNC_MAXIMUM_FFF,   "maximum_fff",    fmaxf_, fmaxf2, vsFmax_)
FUNC_FFF(FUNC_MINIMUM_FFF,   "minimum_fff",    fminf_, fminf2, vsFmin_)
FUNC_FFF(FUNC_FFF_LAST,    NULL,          NULL,   NULL,    NULL)
#ifdef ELIDE_FUNC_FFF
#undef ELIDE_FUNC_FFF
#undef FUNC_FFF
#endif

#ifndef FUNC_DD
#define ELIDE_FUNC_DD
#define FUNC_DD(...)
#endif
FUNC_DD(FUNC_SQRT_DD,    "sqrt_dd",     sqrt,  vdSqrt)
FUNC_DD(FUNC_SIN_DD,     "sin_dd",      sin,   vdSin)
FUNC_DD(FUNC_COS_DD,     "cos_dd",      cos,   vdCos)
FUNC_DD(FUNC_TAN_DD,     "tan_dd",      tan,   vdTan)
FUNC_DD(FUNC_ARCSIN_DD,  "arcsin_dd",   asin,  vdAsin)
FUNC_DD(FUNC_ARCCOS_DD,  "arccos_dd",   acos,  vdAcos)
FUNC_DD(FUNC_ARCTAN_DD,  "arctan_dd",   atan,  vdAtan)
FUNC_DD(FUNC_SINH_DD,    "sinh_dd",     sinh,  vdSinh)
FUNC_DD(FUNC_COSH_DD,    "cosh_dd",     cosh,  vdCosh)
FUNC_DD(FUNC_TANH_DD,    "tanh_dd",     tanh,  vdTanh)
FUNC_DD(FUNC_ARCSINH_DD, "arcsinh_dd",  asinh, vdAsinh)
FUNC_DD(FUNC_ARCCOSH_DD, "arccosh_dd",  acosh, vdAcosh)
FUNC_DD(FUNC_ARCTANH_DD, "arctanh_dd",  atanh, vdAtanh)
FUNC_DD(FUNC_LOG_DD,     "log_dd",      log,   vdLn)
FUNC_DD(FUNC_LOG1P_DD,   "log1p_dd",    log1p, vdLog1p)
FUNC_DD(FUNC_LOG10_DD,   "log10_dd",    log10, vdLog10)
FUNC_DD(FUNC_LOG2_DD,   "log2_dd",      log2,  vdLog2)
FUNC_DD(FUNC_EXP_DD,     "exp_dd",      exp,   vdExp)
FUNC_DD(FUNC_EXPM1_DD,   "expm1_dd",    expm1, vdExpm1)
FUNC_DD(FUNC_ABS_DD,     "absolute_dd", fabs,  vdAbs)
FUNC_DD(FUNC_CONJ_DD,    "conjugate_dd",fconj, vdConj)
FUNC_DD(FUNC_CEIL_DD,    "ceil_dd",     ceil,  vdCeil)
FUNC_DD(FUNC_FLOOR_DD,   "floor_dd",    floor, vdFloor)
FUNC_DD(FUNC_TRUNC_DD,   "trunc_dd",    trunc, vdTrunc)
FUNC_DD(FUNC_SIGN_DD,    "sign_dd",     sign,  vdSign)
//rint rounds to nearest even integer, matching NumPy (round doesn't)
FUNC_DD(FUNC_ROUND_DD,   "round_dd",    rint,  vdRint)
FUNC_DD(FUNC_DD_LAST,    NULL,          NULL,  NULL)
#ifdef ELIDE_FUNC_DD
#undef ELIDE_FUNC_DD
#undef FUNC_DD
#endif

// double -> boolean functions
#ifndef FUNC_BD
#define ELIDE_FUNC_BD
#define FUNC_BD(...)
#endif
FUNC_BD(FUNC_ISNAN_BD,   "isnan_bd",    isnand, vdIsnan)
FUNC_BD(FUNC_ISFINITE_BD, "isfinite_bd", isfinited, vdIsfinite)
FUNC_BD(FUNC_ISINF_BD, "isinf_bd", isinfd, vdIsinf)
FUNC_BD(FUNC_SIGNBIT_BD, "signbit_bd",  signbit, vdSignBit)
FUNC_BD(FUNC_BD_LAST,    NULL,          NULL,  NULL)
#ifdef ELIDE_FUNC_BD
#undef ELIDE_FUNC_BD
#undef FUNC_BD
#endif

// float -> boolean functions (C99 defines the same function for all types)
#ifndef FUNC_BF
#define ELIDE_FUNC_BF
#define FUNC_BF(...)
#endif // use wrappers as there is name collision with isnanf in std
FUNC_BF(FUNC_ISNAN_BF,   "isnan_bf",    isnanf_,  isnanf2, vsIsnan)
FUNC_BF(FUNC_ISFINITE_BF, "isfinite_bf", isfinitef_, isfinitef2, vsIsfinite)
FUNC_BF(FUNC_ISINF_BF, "isinf_bf", isinff_, isinff2, vsIsinf)
FUNC_BF(FUNC_SIGNBIT_BF, "signbit_bf", signbitf, signbitf2, vsSignBit)
FUNC_BF(FUNC_BF_LAST,    NULL,            NULL,        NULL,    NULL)
#ifdef ELIDE_FUNC_BF
#undef ELIDE_FUNC_BF
#undef FUNC_BF
#endif

#ifndef FUNC_DDD
#define ELIDE_FUNC_DDD
#define FUNC_DDD(...)
#endif
FUNC_DDD(FUNC_FMOD_DDD,    "fmod_ddd",    fmod,  vdfmod)
FUNC_DDD(FUNC_ARCTAN2_DDD, "arctan2_ddd", atan2, vdAtan2)
FUNC_DDD(FUNC_HYPOT_DDD, "hypot_ddd", hypot, vdHypot)
FUNC_DDD(FUNC_NEXTAFTER_DDD, "nextafter_ddd", nextafter, vdNextAfter)
FUNC_DDD(FUNC_COPYSIGN_DDD, "copysign_ddd", copysign, vdCopySign)
FUNC_DDD(FUNC_MAXIMUM_DDD, "maximum_ddd",  fmaxd,  vdFmax_)
FUNC_DDD(FUNC_MINIMUM_DDD, "minimum_ddd",  fmind,  vdFmin_)
FUNC_DDD(FUNC_DDD_LAST,    NULL,          NULL,  NULL)
#ifdef ELIDE_FUNC_DDD
#undef ELIDE_FUNC_DDD
#undef FUNC_DDD
#endif

#ifndef FUNC_CC
#define ELIDE_FUNC_CC
#define FUNC_CC(...)
#endif
FUNC_CC(FUNC_SQRT_CC,    "sqrt_cc",     nc_sqrt,   vzSqrt)
FUNC_CC(FUNC_SIN_CC,     "sin_cc",      nc_sin,    vzSin)
FUNC_CC(FUNC_COS_CC,     "cos_cc",      nc_cos,    vzCos)
FUNC_CC(FUNC_TAN_CC,     "tan_cc",      nc_tan,    vzTan)
FUNC_CC(FUNC_ARCSIN_CC,  "arcsin_cc",   nc_asin,   vzAsin)
FUNC_CC(FUNC_ARCCOS_CC,  "arccos_cc",   nc_acos,   vzAcos)
FUNC_CC(FUNC_ARCTAN_CC,  "arctan_cc",   nc_atan,   vzAtan)
FUNC_CC(FUNC_SINH_CC,    "sinh_cc",     nc_sinh,   vzSinh)
FUNC_CC(FUNC_COSH_CC,    "cosh_cc",     nc_cosh,   vzCosh)
FUNC_CC(FUNC_TANH_CC,    "tanh_cc",     nc_tanh,   vzTanh)
FUNC_CC(FUNC_ARCSINH_CC, "arcsinh_cc",  nc_asinh,  vzAsinh)
FUNC_CC(FUNC_ARCCOSH_CC, "arccosh_cc",  nc_acosh,  vzAcosh)
FUNC_CC(FUNC_ARCTANH_CC, "arctanh_cc",  nc_atanh,  vzAtanh)
FUNC_CC(FUNC_LOG_CC,     "log_cc",      nc_log,    vzLn)
FUNC_CC(FUNC_LOG1P_CC,   "log1p_cc",    nc_log1p,  vzLog1p)
FUNC_CC(FUNC_LOG10_CC,   "log10_cc",    nc_log10,  vzLog10)
FUNC_CC(FUNC_LOG2_CC,    "log2_cc",     nc_log2,   vzLog2)
FUNC_CC(FUNC_EXP_CC,     "exp_cc",      nc_exp,    vzExp)
FUNC_CC(FUNC_EXPM1_CC,   "expm1_cc",    nc_expm1,  vzExpm1)
FUNC_CC(FUNC_ABS_CC,     "absolute_cc", nc_abs,    vzAbs_)
FUNC_CC(FUNC_CONJ_CC,    "conjugate_cc",nc_conj,   vzConj)
FUNC_CC(FUNC_SIGN_CC,   "sign_cc",    nc_sign,   vzSign)
// rint rounds to nearest even integer, matches NumPy behaviour (round doesn't)
FUNC_CC(FUNC_ROUND_CC,   "round_cc",    nc_rint,   vzRint)
FUNC_CC(FUNC_CC_LAST,    NULL,          NULL,      NULL)
#ifdef ELIDE_FUNC_CC
#undef ELIDE_FUNC_CC
#undef FUNC_CC
#endif

#ifndef FUNC_CCC
#define ELIDE_FUNC_CCC
#define FUNC_CCC(...)
#endif
FUNC_CCC(FUNC_POW_CCC,   "pow_ccc", nc_pow)
FUNC_CCC(FUNC_CCC_LAST,  NULL,      NULL)
#ifdef ELIDE_FUNC_CCC
#undef ELIDE_FUNC_CCC
#undef FUNC_CCC
#endif

// complex -> boolean functions
#ifndef FUNC_BC
#define ELIDE_FUNC_BC
#define FUNC_BC(...)
#endif // use wrappers as there is name collision with isnanf in std
FUNC_BC(FUNC_ISNAN_BC,   "isnan_bc",    nc_isnan, vzIsnan)
FUNC_BC(FUNC_ISFINITE_BC, "isfinite_bc", nc_isfinite, vzIsfinite)
FUNC_BC(FUNC_ISINF_BC, "isinf_bc", nc_isinf, vzIsinf)
FUNC_BC(FUNC_BC_LAST,    NULL,            NULL,        NULL)
#ifdef ELIDE_FUNC_BC
#undef ELIDE_FUNC_BC
#undef FUNC_BC
#endif

// int -> int functions
#ifndef FUNC_II
#define ELIDE_FUNC_II
#define FUNC_II(...)
#endif
FUNC_II(FUNC_SIGN_II,   "sign_ii",  signi, viSign)
FUNC_II(FUNC_ROUND_II,  "round_ii", rinti, viRint)
FUNC_II(FUNC_ABS_II,  "absolute_ii", fabsi, viFabs)
FUNC_II(FUNC_II_LAST,  NULL,      NULL, NULL)
#ifdef ELIDE_FUNC_II
#undef ELIDE_FUNC_II
#undef FUNC_II
#endif

#ifndef FUNC_LL
#define ELIDE_FUNC_LL
#define FUNC_LL(...)
#endif
FUNC_LL(FUNC_SIGN_LL,  "sign_ll",   signl,   vlSign)
FUNC_LL(FUNC_ROUND_LL, "round_ll",  rintl,   vlRint)
FUNC_LL(FUNC_ABS_LL,  "absolute_ll", fabsl, vlFabs)
FUNC_LL(FUNC_LL_LAST,  NULL,      NULL, NULL)
#ifdef ELIDE_FUNC_LL
#undef ELIDE_FUNC_LL
#undef FUNC_LL
#endif
