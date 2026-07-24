#include <float.h>  // for _finite, _isnan on MSVC

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

/* Due to casting problems (normally return ints not bools, easiest to define
non-overloaded wrappers for these functions) */
// MSVC version: use global ::isfinite / ::isnan
inline bool isfinitef_(float x) { return !!::_finite(x); }   // MSVC has _finite
inline bool isnanf_(float x)    { return !!::_isnan(x); }    // MSVC has _isnan
inline bool isfinited(double x) { return !!::_finite(x); }
inline bool isnand(double x)    { return !!::_isnan(x); }
inline bool isinfd(double x) { return !!::isinf(x); }
inline bool isinff_(float x)    { return !!::isinf(x); }

// To handle overloading of fmax/fmin in cmath and match NumPy behaviour for NaNs
inline double fmaxd(double x, double y)    { return (isnand(x) | isnand(y))? NAN : fmax(x, y); }
inline double fmind(double x, double y)    { return (isnand(x) | isnand(y))? NAN : fmin(x, y); }


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
    #define log2f(x)    ((float)log2((double)(x)))
    #define expf(x)    ((float)exp((double)(x)))
    #define expm1f(x)    ((float)expm1((double)(x)))
    #define fabsf(x)    ((float)fabs((double)(x)))
    #define fmodf(x, y)    ((float)fmod((double)(x), (double)(y)))
    #define atan2f(x, y)    ((float)atan2((double)(x), (double)(y)))
    #define hypotf(x, y)    ((float)hypot((double)(x), (double)(y)))
    #define copysignf(x, y)    ((float)copysign((double)(x), (double)(y)))
    #define nextafterf(x, y)    ((float)nextafter((double)(x), (double)(y)))
    #define ceilf(x)    ((float)ceil((double)(x)))
    #define hypotf(x)    ((float)hypot((double)(x)))
    #define rintf(x)    ((float)rint((double)(x)))
    #define truncf(x)    ((float)trunc((double)(x)))


    /* The next are directly called from interp_body.cpp */
    #define powf(x, y)    ((float)pow((double)(x), (double)(y)))
    #define floorf(x)    ((float)floor((double)(x)))

    #define fmaxf_(x, y)    ((float)fmaxd((double)(x), (double)(y))) // define fmaxf_ since fmaxf doesn't exist for early MSVC
    #define fminf_(x, y)    ((float)fmind((double)(x), (double)(y)))
#else
    inline float fmaxf_(float x, float y)    { return (isnanf_(x) | isnanf_(y))? NAN : fmaxf(x, y); }
    inline float fminf_(float x, float y)    { return (isnanf_(x) | isnanf_(y))? NAN : fminf(x, y); }
#endif // _MSC_VER < 1400


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

inline float log2f2(float x) {
    return log2f(x);
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

inline float hypotf2(float x, float y) {
    return hypotf(x, y);
}

inline float nextafterf2(float x, float y) {
    return nextafterf(x, y);
}

inline float copysignf2(float x, float y) {
    return copysignf(x, y);
}

inline float fmaxf2(float x, float y) {
    return fmaxf_(x, y);
}

inline float fminf2(float x, float y) {
    return fminf_(x, y);
}


// Boolean output functions
inline bool isnanf2(float x) {
    return isnanf_(x);
}

inline bool isfinitef2(float x) {
    return isfinitef_(x);
}

inline bool isinff2(float x) {
    return isinff_(x);
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

inline float rintf2(float x) {
    return rintf(x);
}

inline float truncf2(float x) {
    return truncf(x);
}

inline bool signbitf2(float x) {
    return signbitf(x);
}

#endif // NUMEXPR_MSVC_FUNCTION_STUBS_HPP
