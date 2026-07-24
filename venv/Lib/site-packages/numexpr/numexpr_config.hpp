#ifndef NUMEXPR_CONFIG_HPP
#define NUMEXPR_CONFIG_HPP

// x86 platform works with unaligned reads and writes
// MW: I have seen exceptions to this when the compiler chooses to use aligned SSE
#if (defined(NPY_CPU_X86) || defined(NPY_CPU_AMD64))
#  define USE_UNALIGNED_ACCESS 1
#endif

// #ifdef SCIPY_MKL_H
// #define USE_VML
// #endif

#ifdef USE_VML
/* The values below have been tuned for a Skylake processor (E3-1245 v5 @ 3.50GHz) */
#define BLOCK_SIZE1 1024
#else
/* The values below have been tuned for a Skylake processor (E3-1245 v5 @ 3.50GHz) */
#define BLOCK_SIZE1 1024
#endif

// The default threadpool size. It's prefer that the user set this via an
// environment variable, "NUMEXPR_MAX_THREADS"
#define DEFAULT_MAX_THREADS 64

// Remove dependence on NPY_MAXARGS, which would be a runtime constant instead of compiletime
// constant. If numpy raises NPY_MAXARGS, we should notice and raise this as well
#define NE_MAXARGS 64

#if defined(_WIN32)
  #include "win32/pthread.h"
  #include <process.h>
  #define getpid _getpid
#else
  #include <pthread.h>
  #include "unistd.h"
#endif

#ifdef USE_VML
#include "mkl_vml.h"
#include "mkl_service.h"
#endif
#include <cmath>
//no single precision version of signbit in C++ standard
inline bool signbitf(float x) { return signbit((double)x); }

#ifdef _WIN32
  #ifndef __MINGW32__
    #include "missing_posix_functions.hpp"
  #endif
  #include "msvc_function_stubs.hpp"
#else
/* GCC/Clang version: use std:: (can't use it for windows)
   msvc_function_stubs contains windows alternatives */
/* Due to casting problems (normally return ints not bools, easiest to define
   non-overloaded wrappers for these functions) */
inline bool isfinitef_(float x) { return !!std::isfinite(x); }
inline bool isnanf_(float x)    { return !!std::isnan(x); }
inline bool isfinited(double x) { return !!std::isfinite(x); }
inline bool isnand(double x)    { return !!std::isnan(x); }
inline bool isinff_(float x) { return !!std::isinf(x); }
inline bool isinfd(double x)    { return !!std::isinf(x); }

// To handle overloading of fmax/fmin in cmath and match NumPy behaviour for NaNs
inline double fmaxd(double x, double y)    { return (isnand(x) | isnand(y))? NAN : fmax(x, y); }
inline double fmind(double x, double y)    { return (isnand(x) | isnand(y))? NAN : fmin(x, y); }
inline float fmaxf_(float x, float y)    { return (isnanf_(x) | isnanf_(y))? NAN : fmaxf(x, y); }
inline float fminf_(float x, float y)    { return (isnanf_(x) | isnanf_(y))? NAN : fminf(x, y); }
#endif

#endif // NUMEXPR_CONFIG_HPP
