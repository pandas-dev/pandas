#ifndef NUMEXPR_MSVC_FUNCTION_STUBS_HPP
#define NUMEXPR_MSVC_FUNCTION_STUBS_HPP

/*********************************************************************
  Numexpr - Fast numerical array expression evaluator for NumPy.

      License: MIT
      Author:  See AUTHORS.txt

  See LICENSE.txt for details about copyright and rights to use.
**********************************************************************/

/* Declare stub functions for MSVC.  It turns out that single precision
   definitions in <math.h> are actually #define'd and are not usable
   as function pointers :-/ */

#if _MSC_VER < 1400  // 1310 == MSVC 7.1
/* Apparently, single precision functions are not included in MSVC 7.1 */

#define sqrtf(x)    ((float)sqrt((double)(x)))
#define sinf(x)    ((float)sin((double)(x)))
#define cosf(x)    ((float)cos((double)(x)))
#define tanf(x)    ((float)tan((double)(x)))
#define asinf(x)    ((float)asin((double)(x)))
#define acosf(x)    ((float)acos((double)(x)))
#define atanf(x)    ((float)atan((double)(x)))
#define sinhf(x)    ((float)sinh((double)(x)))
#define coshf(x)    ((float)cosh((double)(x)))
#define tanhf(x)    ((float)tanh((double)(x)))
#define asinhf(x)    ((float)asinh((double)(x)))
#define acoshf(x)    ((float)acosh((double)(x)))
#define atanhf(x)    ((float)atanh((double)(x)))
#define logf(x)    ((float)log((double)(x)))
#define log1pf(x)    ((float)log1p((double)(x)))
#define log10f(x)    ((float)log10((double)(x)))
#define expf(x)    ((float)exp((double)(x)))
#define expm1f(x)    ((float)expm1((double)(x)))
#define fabsf(x)    ((float)fabs((double)(x)))
#define fmodf(x, y)    ((float)fmod((double)(x), (double)(y)))
#define atan2f(x, y)    ((float)atan2((double)(x), (double)(y)))
#define ceilf(x)    ((float)ceil((double)(x)))

/* The next are directly called from interp_body.cpp */
#define powf(x, y)    ((float)pow((double)(x), (double)(y)))
#define floorf(x)    ((float)floor((double)(x)))

#endif  // _MSC_VER < 1400


/* Now the actual stubs */

inline float sqrtf2(float x) {
    return sqrtf(x);
}

inline float sinf2(float x) {
    return sinf(x);
}

inline float cosf2(float x) {
    return cosf(x);
}

inline float tanf2(float x) {
    return tanf(x);
}

inline float asinf2(float x) {
    return asinf(x);
}

inline float acosf2(float x) {
    return acosf(x);
}

inline float atanf2(float x) {
    return atanf(x);
}

inline float sinhf2(float x) {
    return sinhf(x);
}

inline float coshf2(float x) {
    return coshf(x);
}

inline float tanhf2(float x) {
    return tanhf(x);
}

inline float asinhf2(float x) {
    return asinhf(x);
}

inline float acoshf2(float x) {
    return acoshf(x);
}

inline float atanhf2(float x) {
    return atanhf(x);
}

inline float logf2(float x) {
    return logf(x);
}

inline float log1pf2(float x) {
    return log1pf(x);
}

inline float log10f2(float x) {
    return log10f(x);
}

inline float expf2(float x) {
    return expf(x);
}

inline float expm1f2(float x) {
    return expm1f(x);
}

inline float fabsf2(float x) {
    return fabsf(x);
}

inline float fmodf2(float x, float y) {
    return fmodf(x, y);
}

inline float atan2f2(float x, float y) {
    return atan2f(x, y);
}

// Needed for allowing the internal casting in numexpr machinery for
// conjugate operations
inline float fconjf2(float x) {
    return x;
}

inline float ceilf2(float x) {
    return ceilf(x);
}

inline float floorf2(float x) {
    return floorf(x);
}

#endif // NUMEXPR_MSVC_FUNCTION_STUBS_HPP
