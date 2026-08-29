/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#pragma once

#include <string.h>

#if defined(_MSC_VER)
#  define strncasecmp(s1, s2, n) _strnicmp(s1, s2, n)
#endif

// GH-23516 - works around locale perf issues
// from MUSL libc, licence at LICENSES/MUSL_LICENSE
#define isdigit_ascii(c) (((unsigned)(c) - '0') < 10u)
#define getdigit_ascii(c, default)                                             \
  (isdigit_ascii(c) ? ((int)((c) - '0')) : default)
#define isspace_ascii(c) (((c) == ' ') || (((unsigned)(c) - '\t') < 5))
#define toupper_ascii(c) ((((unsigned)(c) - 'a') < 26) ? ((c) & 0x5f) : (c))
#define tolower_ascii(c) ((((unsigned)(c) - 'A') < 26) ? ((c) | 0x20) : (c))

#ifdef __has_attribute
#  if __has_attribute(__fallthrough__)
#    define PD_FALLTHROUGH __attribute__((__fallthrough__))
#  endif
#endif

#ifndef PD_FALLTHROUGH
#  define PD_FALLTHROUGH                                                       \
    do {                                                                       \
    } while (0) /* fallthrough */
#endif

#if defined(_WIN32)
#  ifndef ENABLE_INTSAFE_SIGNED_FUNCTIONS
#    define ENABLE_INTSAFE_SIGNED_FUNCTIONS
#  endif
#  include <intsafe.h>
#  ifdef __cplusplus
// _Generic (below) is C-only. C++ TUs get inline overloads instead; only the
// overload whose result-pointer type matches exactly is viable, so overload
// resolution reproduces the _Generic dispatch, including the argument
// conversions implied by the selected intsafe prototype.
static inline HRESULT pd_checked_add(short a, short b, short *res) {
  return ShortAdd(a, b, res);
}
static inline HRESULT pd_checked_add(unsigned short a, unsigned short b,
                                     unsigned short *res) {
  return UShortAdd(a, b, res);
}
static inline HRESULT pd_checked_add(int a, int b, int *res) {
  return IntAdd(a, b, res);
}
static inline HRESULT pd_checked_add(unsigned int a, unsigned int b,
                                     unsigned int *res) {
  return UIntAdd(a, b, res);
}
static inline HRESULT pd_checked_add(long a, long b, long *res) {
  return LongAdd(a, b, res);
}
static inline HRESULT pd_checked_add(unsigned long a, unsigned long b,
                                     unsigned long *res) {
  return ULongAdd(a, b, res);
}
static inline HRESULT pd_checked_add(long long a, long long b, long long *res) {
  return LongLongAdd(a, b, res);
}
static inline HRESULT pd_checked_add(unsigned long long a, unsigned long long b,
                                     unsigned long long *res) {
  return ULongLongAdd(a, b, res);
}
static inline HRESULT pd_checked_sub(short a, short b, short *res) {
  return ShortSub(a, b, res);
}
static inline HRESULT pd_checked_sub(unsigned short a, unsigned short b,
                                     unsigned short *res) {
  return UShortSub(a, b, res);
}
static inline HRESULT pd_checked_sub(int a, int b, int *res) {
  return IntSub(a, b, res);
}
static inline HRESULT pd_checked_sub(unsigned int a, unsigned int b,
                                     unsigned int *res) {
  return UIntSub(a, b, res);
}
static inline HRESULT pd_checked_sub(long a, long b, long *res) {
  return LongSub(a, b, res);
}
static inline HRESULT pd_checked_sub(unsigned long a, unsigned long b,
                                     unsigned long *res) {
  return ULongSub(a, b, res);
}
static inline HRESULT pd_checked_sub(long long a, long long b, long long *res) {
  return LongLongSub(a, b, res);
}
static inline HRESULT pd_checked_sub(unsigned long long a, unsigned long long b,
                                     unsigned long long *res) {
  return ULongLongSub(a, b, res);
}
static inline HRESULT pd_checked_mul(short a, short b, short *res) {
  return ShortMult(a, b, res);
}
static inline HRESULT pd_checked_mul(unsigned short a, unsigned short b,
                                     unsigned short *res) {
  return UShortMult(a, b, res);
}
static inline HRESULT pd_checked_mul(int a, int b, int *res) {
  return IntMult(a, b, res);
}
static inline HRESULT pd_checked_mul(unsigned int a, unsigned int b,
                                     unsigned int *res) {
  return UIntMult(a, b, res);
}
static inline HRESULT pd_checked_mul(long a, long b, long *res) {
  return LongMult(a, b, res);
}
static inline HRESULT pd_checked_mul(unsigned long a, unsigned long b,
                                     unsigned long *res) {
  return ULongMult(a, b, res);
}
static inline HRESULT pd_checked_mul(long long a, long long b, long long *res) {
  return LongLongMult(a, b, res);
}
static inline HRESULT pd_checked_mul(unsigned long long a, unsigned long long b,
                                     unsigned long long *res) {
  return ULongLongMult(a, b, res);
}
#    define checked_add(a, b, res) pd_checked_add((a), (b), (res))
#    define checked_sub(a, b, res) pd_checked_sub((a), (b), (res))
#    define checked_mul(a, b, res) pd_checked_mul((a), (b), (res))
#  else
#    define checked_add(a, b, res)                                             \
      _Generic((res),                                                          \
          int *: IntAdd,                                                       \
          unsigned int *: UIntAdd,                                             \
          long *: LongAdd,                                                     \
          unsigned long *: ULongAdd,                                           \
          long long *: LongLongAdd,                                            \
          unsigned long long *: ULongLongAdd,                                  \
          short *: ShortAdd,                                                   \
          unsigned short *: UShortAdd)(a, b, res)

#    define checked_sub(a, b, res)                                             \
      _Generic((res),                                                          \
          int *: IntSub,                                                       \
          unsigned int *: UIntSub,                                             \
          long *: LongSub,                                                     \
          unsigned long *: ULongSub,                                           \
          long long *: LongLongSub,                                            \
          unsigned long long *: ULongLongSub,                                  \
          short *: ShortSub,                                                   \
          unsigned short *: UShortSub)(a, b, res)

#    define checked_mul(a, b, res)                                             \
      _Generic((res),                                                          \
          int *: IntMult,                                                      \
          unsigned int *: UIntMult,                                            \
          long *: LongMult,                                                    \
          unsigned long *: ULongMult,                                          \
          long long *: LongLongMult,                                           \
          unsigned long long *: ULongLongMult,                                 \
          short *: ShortMult,                                                  \
          unsigned short *: UShortMult)(a, b, res)
#  endif

#elif (defined(__has_builtin) && __has_builtin(__builtin_add_overflow)) ||     \
    __GNUC__ > 7
#  define checked_add(a, b, res) __builtin_add_overflow(a, b, res)
#  define checked_sub(a, b, res) __builtin_sub_overflow(a, b, res)
#  define checked_mul(a, b, res) __builtin_mul_overflow(a, b, res)
#else
_Static_assert(0,
               "Overflow checking not detected; please try a newer compiler");
#endif
