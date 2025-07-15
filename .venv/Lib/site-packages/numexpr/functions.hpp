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
FUNC_FF(FUNC_EXP_FF,     "exp_ff",      expf,   expf2,   vsExp)
FUNC_FF(FUNC_EXPM1_FF,   "expm1_ff",    expm1f, expm1f2, vsExpm1)
FUNC_FF(FUNC_ABS_FF,     "absolute_ff", fabsf,  fabsf2,  vsAbs)
FUNC_FF(FUNC_CONJ_FF,    "conjugate_ff",fconjf, fconjf2, vsConj)
FUNC_FF(FUNC_CEIL_FF,    "ceil_ff",     ceilf,  ceilf2,  vsCeil)
FUNC_FF(FUNC_FLOOR_FF,   "floor_ff",    floorf, floorf2, vsFloor)
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
FUNC_DD(FUNC_EXP_DD,     "exp_dd",      exp,   vdExp)
FUNC_DD(FUNC_EXPM1_DD,   "expm1_dd",    expm1, vdExpm1)
FUNC_DD(FUNC_ABS_DD,     "absolute_dd", fabs,  vdAbs)
FUNC_DD(FUNC_CONJ_DD,    "conjugate_dd",fconj, vdConj)
FUNC_DD(FUNC_CEIL_DD,    "ceil_dd",     ceil,  vdCeil)
FUNC_DD(FUNC_FLOOR_DD,   "floor_dd",    floor, vdFloor)
FUNC_DD(FUNC_DD_LAST,    NULL,          NULL,  NULL)
#ifdef ELIDE_FUNC_DD
#undef ELIDE_FUNC_DD
#undef FUNC_DD
#endif

#ifndef FUNC_DDD
#define ELIDE_FUNC_DDD
#define FUNC_DDD(...)
#endif
FUNC_DDD(FUNC_FMOD_DDD,    "fmod_ddd",    fmod,  vdfmod)
FUNC_DDD(FUNC_ARCTAN2_DDD, "arctan2_ddd", atan2, vdAtan2)
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
FUNC_CC(FUNC_EXP_CC,     "exp_cc",      nc_exp,    vzExp)
FUNC_CC(FUNC_EXPM1_CC,   "expm1_cc",    nc_expm1,  vzExpm1)
FUNC_CC(FUNC_ABS_CC,     "absolute_cc", nc_abs,    vzAbs_)
FUNC_CC(FUNC_CONJ_CC,    "conjugate_cc",nc_conj,   vzConj)
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
