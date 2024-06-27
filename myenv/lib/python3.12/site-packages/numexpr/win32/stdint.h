/* ISO C9x  7.18  Integer types <stdint.h>
 * Based on ISO/IEC SC22/WG14 9899 Committee draft (SC22 N2794)
 *
 *  THIS SOFTWARE IS NOT COPYRIGHTED
 *
 *  Contributor: Danny Smith <danny_r_smith_2001@yahoo.co.nz>
 *
 *  This source code is offered for use in the public domain. You may
 *  use, modify or distribute it freely.
 *
 *  This code is distributed in the hope that it will be useful but
 *  WITHOUT ANY WARRANTY. ALL WARRANTIES, EXPRESS OR IMPLIED ARE HEREBY
 *  DISCLAIMED. This includes but is not limited to warranties of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *  Date: 2000-12-02
 *
 * mwb: This was modified in the following ways:
 *
 *      - make it compatible with Visual C++ 6 (which uses 
 *          non-standard keywords and suffixes for 64-bit types)
 *      - some environments need stddef.h included (for wchar stuff?)
 *      - handle the fact that Microsoft's limits.h header defines
 *          SIZE_MAX
 *      - make corrections for SIZE_MAX, INTPTR_MIN, INTPTR_MAX, UINTPTR_MAX,
 *          PTRDIFF_MIN, PTRDIFF_MAX, SIG_ATOMIC_MIN, and SIG_ATOMIC_MAX
 *          to be 64-bit aware.
 */


#ifndef _STDINT_H
#define _STDINT_H
#define __need_wint_t
#define __need_wchar_t
#include <wchar.h>
#include <stddef.h>

#if _MSC_VER && (_MSC_VER < 1300)
/* using MSVC 6 or earlier - no "long long" type, but might have _int64 type */
#define __STDINT_LONGLONG           __int64
#define __STDINT_LONGLONG_SUFFIX    i64
#else
#define __STDINT_LONGLONG           long long
#define __STDINT_LONGLONG_SUFFIX    LL
#endif

#if !defined( PASTE)
#define PASTE2( x, y) x##y
#define PASTE( x, y)  PASTE2( x, y)
#endif /* PASTE */


/* 7.18.1.1  Exact-width integer types */
typedef signed char int8_t;
typedef unsigned char   uint8_t;
typedef short  int16_t;
typedef unsigned short  uint16_t;
typedef int  int32_t;
typedef unsigned   uint32_t;
typedef __STDINT_LONGLONG  int64_t;
typedef unsigned __STDINT_LONGLONG   uint64_t;

/* 7.18.1.2  Minimum-width integer types */
typedef signed char int_least8_t;
typedef unsigned char   uint_least8_t;
typedef short  int_least16_t;
typedef unsigned short  uint_least16_t;
typedef int  int_least32_t;
typedef unsigned   uint_least32_t;
typedef __STDINT_LONGLONG  int_least64_t;
typedef unsigned __STDINT_LONGLONG   uint_least64_t;

/*  7.18.1.3  Fastest minimum-width integer types 
 *  Not actually guaranteed to be fastest for all purposes
 *  Here we use the exact-width types for 8 and 16-bit ints. 
 */
typedef char int_fast8_t;
typedef unsigned char uint_fast8_t;
typedef short  int_fast16_t;
typedef unsigned short  uint_fast16_t;
typedef int  int_fast32_t;
typedef unsigned  int  uint_fast32_t;
typedef __STDINT_LONGLONG  int_fast64_t;
typedef unsigned __STDINT_LONGLONG   uint_fast64_t;

/* 7.18.1.4  Integer types capable of holding object pointers */
#ifndef _INTPTR_T_DEFINED
#define _INTPTR_T_DEFINED
#ifdef _WIN64
typedef __STDINT_LONGLONG intptr_t
#else
typedef int intptr_t;
#endif /* _WIN64 */
#endif /* _INTPTR_T_DEFINED */

#ifndef _UINTPTR_T_DEFINED
#define _UINTPTR_T_DEFINED
#ifdef _WIN64
typedef unsigned __STDINT_LONGLONG uintptr_t
#else
typedef unsigned int uintptr_t;
#endif /* _WIN64 */
#endif /* _UINTPTR_T_DEFINED */

/* 7.18.1.5  Greatest-width integer types */
typedef __STDINT_LONGLONG  intmax_t;
typedef unsigned __STDINT_LONGLONG   uintmax_t;

/* 7.18.2  Limits of specified-width integer types */
#if !defined ( __cplusplus) || defined (__STDC_LIMIT_MACROS)

/* 7.18.2.1  Limits of exact-width integer types */
#define INT8_MIN (-128) 
#define INT16_MIN (-32768)
#define INT32_MIN (-2147483647 - 1)
#define INT64_MIN  (PASTE( -9223372036854775807, __STDINT_LONGLONG_SUFFIX) - 1)

#define INT8_MAX 127
#define INT16_MAX 32767
#define INT32_MAX 2147483647
#define INT64_MAX (PASTE( 9223372036854775807, __STDINT_LONGLONG_SUFFIX))

#define UINT8_MAX 0xff /* 255U */
#define UINT16_MAX 0xffff /* 65535U */
#define UINT32_MAX 0xffffffff  /* 4294967295U */
#define UINT64_MAX (PASTE( 0xffffffffffffffffU, __STDINT_LONGLONG_SUFFIX)) /* 18446744073709551615ULL */

/* 7.18.2.2  Limits of minimum-width integer types */
#define INT_LEAST8_MIN INT8_MIN
#define INT_LEAST16_MIN INT16_MIN
#define INT_LEAST32_MIN INT32_MIN
#define INT_LEAST64_MIN INT64_MIN

#define INT_LEAST8_MAX INT8_MAX
#define INT_LEAST16_MAX INT16_MAX
#define INT_LEAST32_MAX INT32_MAX
#define INT_LEAST64_MAX INT64_MAX

#define UINT_LEAST8_MAX UINT8_MAX
#define UINT_LEAST16_MAX UINT16_MAX
#define UINT_LEAST32_MAX UINT32_MAX
#define UINT_LEAST64_MAX UINT64_MAX

/* 7.18.2.3  Limits of fastest minimum-width integer types */
#define INT_FAST8_MIN INT8_MIN
#define INT_FAST16_MIN INT16_MIN
#define INT_FAST32_MIN INT32_MIN
#define INT_FAST64_MIN INT64_MIN

#define INT_FAST8_MAX INT8_MAX
#define INT_FAST16_MAX INT16_MAX
#define INT_FAST32_MAX INT32_MAX
#define INT_FAST64_MAX INT64_MAX

#define UINT_FAST8_MAX UINT8_MAX
#define UINT_FAST16_MAX UINT16_MAX
#define UINT_FAST32_MAX UINT32_MAX
#define UINT_FAST64_MAX UINT64_MAX

/* 7.18.2.4  Limits of integer types capable of holding
    object pointers */ 
#ifdef _WIN64
#define INTPTR_MIN INT64_MIN
#define INTPTR_MAX INT64_MAX
#define UINTPTR_MAX UINT64_MAX
#else
#define INTPTR_MIN INT32_MIN
#define INTPTR_MAX INT32_MAX
#define UINTPTR_MAX UINT32_MAX
#endif /* _WIN64 */

/* 7.18.2.5  Limits of greatest-width integer types */
#define INTMAX_MIN INT64_MIN
#define INTMAX_MAX INT64_MAX
#define UINTMAX_MAX UINT64_MAX

/* 7.18.3  Limits of other integer types */
#define PTRDIFF_MIN INTPTR_MIN
#define PTRDIFF_MAX INTPTR_MAX

#define SIG_ATOMIC_MIN INTPTR_MIN
#define SIG_ATOMIC_MAX INTPTR_MAX

/* we need to check for SIZE_MAX already defined because MS defines it in limits.h */
#ifndef SIZE_MAX
#define SIZE_MAX UINTPTR_MAX
#endif

#ifndef WCHAR_MIN  /* also in wchar.h */ 
#define WCHAR_MIN 0
#define WCHAR_MAX ((wchar_t)-1) /* UINT16_MAX */
#endif

/*
 * wint_t is unsigned short for compatibility with MS runtime
 */
#define WINT_MIN 0
#define WINT_MAX ((wint_t)-1) /* UINT16_MAX */

#endif /* !defined ( __cplusplus) || defined __STDC_LIMIT_MACROS */


/* 7.18.4  Macros for integer constants */
#if !defined ( __cplusplus) || defined (__STDC_CONSTANT_MACROS)

/* 7.18.4.1  Macros for minimum-width integer constants

    Accoding to Douglas Gwyn <gwyn@arl.mil>:
	"This spec was changed in ISO/IEC 9899:1999 TC1; in ISO/IEC
	9899:1999 as initially published, the expansion was required
	to be an integer constant of precisely matching type, which
	is impossible to accomplish for the shorter types on most
	platforms, because C99 provides no standard way to designate
	an integer constant with width less than that of type int.
	TC1 changed this to require just an integer constant
	*expression* with *promoted* type."
*/

#define INT8_C(val) ((int8_t) + (val))
#define UINT8_C(val) ((uint8_t) + (val##U))
#define INT16_C(val) ((int16_t) + (val))
#define UINT16_C(val) ((uint16_t) + (val##U))

#define INT32_C(val) val##L
#define UINT32_C(val) val##UL
#define INT64_C(val) (PASTE( val, __STDINT_LONGLONG_SUFFIX))
#define UINT64_C(val)(PASTE( PASTE( val, U), __STDINT_LONGLONG_SUFFIX))

/* 7.18.4.2  Macros for greatest-width integer constants */
#define INTMAX_C(val)  INT64_C(val)
#define UINTMAX_C(val) UINT64_C(val)

#endif  /* !defined ( __cplusplus) || defined __STDC_CONSTANT_MACROS */

#endif
