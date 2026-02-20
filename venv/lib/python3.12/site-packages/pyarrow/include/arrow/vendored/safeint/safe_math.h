// Licensed under the MIT License.
// Copyright David LeBlanc - dcl@dleblanc.net

#if !defined SAFE_MATH_H
#define SAFE_MATH_H

#if defined SAFEINT_HPP
#error use either the C++ SafeInt, or safe_math, not both
#endif

// C wants a prototype, if all warnings enabled
// #if !defined SAFE_MATH_FAIL_DEFINED
// static inline void safe_math_fail(const char* msg);
// #endif

#include "safe_math_impl.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*
	The following functions are defined in safe_math_impl.h:

	// All check_cast functions return 0 if safe, non-zero if unsafe
	// Casting test to int8
	int check_cast_int8_int32(int32_t in)
	int check_cast_int8_uint32(uint32_t in)
	int check_cast_int8_int64(int64_t in)
	int check_cast_int8_uint64(uint64_t in)

	// Casting to int64
	int check_cast_int16_int32(int32_t in)
	int check_cast_int16_uint32(uint32_t in)
	int check_cast_int16_int64(int64_t in)
	int check_cast_int16_uint64(uint64_t in)

	// Casting to int32
	int check_cast_int32_uint32(uint32_t in)
	int check_cast_int32_int64(int64_t in)
	int check_cast_int32_uint64(uint64_t in)
	int check_cast_int64_uint64(uint64_t in)

	// Casting to uint8
	int check_cast_uint8_int32(int32_t in)
	int check_cast_uint8_uint32(uint32_t in)
	int check_cast_uint8_int64(int64_t in)
	int check_cast_uint8_uint64(uint64_t in)

	// Casting to uint16
	int check_cast_uint16_int32(int32_t in)
	int check_cast_uint16_uint32(uint32_t in)
	int check_cast_uint16_int64(int64_t in)
	int check_cast_uint16_uint64(uint64_t in)

	// Casting to uint32
	int check_cast_uint32_int32(int32_t in)
	int check_cast_uint32_int64(int64_t in)
	int check_cast_uint32_uint64(uint64_t in)

	// Casting to uint64
	int check_cast_uint64_int64(int64_t in)

	// safe_cast functions all abort on failure

	// Casting to int8
	int8_t safe_cast_int8_int32(int32_t in)
	int8_t safe_cast_int8_uint32(uint32_t in)
	int8_t safe_cast_int8_int64(int64_t in)
	int8_t safe_cast_int8_uint64(uint64_t in)

	// Casting to int16
	int16_t safe_cast_int16_int32(int32_t in)
	int16_t safe_cast_int16_uint32(uint32_t in)
	int16_t safe_cast_int16_int64(int64_t in)
	int16_t safe_cast_int16_uint64(uint64_t in)

	// Casting to int32
	int32_t safe_cast_int32_uint32(uint32_t in)
	int32_t safe_cast_int32_int64(int64_t in)
	int32_t safe_cast_int32_uint64(uint64_t in)

	// Casting to int64
	int64_t safe_cast_int64_uint64(uint64_t in)

	// Casting to uint8
	uint8_t safe_cast_uint8_int32(int32_t in)
	uint8_t safe_cast_uint8_uint32(uint32_t in)
	uint8_t safe_cast_uint8_int64(int64_t in)
	uint8_t safe_cast_uint8_uint64(uint64_t in)

	// Casting to uint16
	uint16_t safe_cast_uint16_int32(int32_t in)
	uint16_t safe_cast_uint16_uint32(uint32_t in)
	uint16_t safe_cast_uint16_int64(int64_t in)
	uint16_t safe_cast_uint16_uint64(uint64_t in)

	// Casting to uint32
	uint32_t safe_cast_uint32_int32(int32_t in)
	uint32_t safe_cast_uint32_int64(int64_t in)
	uint32_t safe_cast_uint32_uint64(uint64_t in)

	// Casting to uint64
	uint64_t safe_cast_uint64_int64(int64_t in)

	// Only 32-bit or larger types are supported for 
	// addition, subtraction, multiplication and division

	// If smaller types are needed, either wrap the result
	// in a safe_cast, or pass the smaller type in as a
	// 32-bit type of the same signedness

	// Addition functions, all of these abort on failure
	// For all of the below, there are also non-aborting versions
	// that have the signature of:
	//
	// bool check_op_intXX_intYY(intXX a, intYY b, intXX* ret)

	int32_t safe_add_int32_int32(int32_t a, int32_t b)
	int32_t safe_add_int32_uint32(int32_t a, uint32_t b)
	int32_t safe_add_int32_int64(int32_t a, int64_t b)
	int32_t safe_add_int32_uint64(int32_t a, uint64_t b)

	uint32_t safe_add_uint32_int32(uint32_t a, int32_t b)
	uint32_t safe_add_uint32_uint32(uint32_t a, uint32_t b)
	uint32_t safe_add_uint32_int64(uint32_t a, int64_t b)
	uint32_t safe_add_uint32_uint64(uint32_t a, uint64_t b)

	int64_t safe_add_int64_int32(int64_t a, int32_t b)
	int64_t safe_add_int64_uint32(int64_t a, uint32_t b)
	int64_t safe_add_int64_int64(int64_t a, int64_t b)
	int64_t safe_add_int64_uint64(int64_t a, uint64_t b)

	uint64_t safe_add_uint64_int32(uint64_t a, int32_t b)
	uint64_t safe_add_uint64_uint32(uint64_t a, uint32_t b)
	uint64_t safe_add_uint64_int64(uint64_t a, int64_t b)
	uint64_t safe_add_uint64_uint64(uint64_t a, uint64_t b)

	// Multiplication
	int32_t safe_div_int32_int32(int32_t a, int32_t b)
	int32_t safe_div_int32_uint32(int32_t a, uint32_t b)
	int32_t safe_div_int32_int64(int32_t a, int64_t b)
	int32_t safe_div_int32_uint64(int32_t a, uint64_t b)

	uint32_t safe_div_uint32_int32(uint32_t a, int32_t b)
	uint32_t safe_div_uint32_uint32(uint32_t a, uint32_t b)
	uint32_t safe_div_uint32_int64(uint32_t a, int64_t b)
	uint32_t safe_div_uint32_uint64(uint32_t a, uint64_t b)

	int64_t safe_div_int64_int32(int64_t a, int32_t b)
	int64_t safe_div_int64_uint32(int64_t a, uint32_t b)
	int64_t safe_div_int64_int64(int64_t a, int64_t b)
	int64_t safe_div_int64_uint64(int64_t a, uint64_t b)

	uint64_t safe_div_uint64_int32(uint64_t a, int32_t b)
	uint64_t safe_div_uint64_uint32(uint64_t a, uint32_t b)
	uint64_t safe_div_uint64_int64(uint64_t a, int64_t b)
	uint64_t safe_div_uint64_uint64(uint64_t a, uint64_t b)

	// Division
	int32_t safe_div_int32_int32(int32_t a, int32_t b)
	int32_t safe_div_int32_uint32(int32_t a, uint32_t b)
	int32_t safe_div_int32_int64(int32_t a, int64_t b)
	int32_t safe_div_int32_uint64(int32_t a, uint64_t b)

	uint32_t safe_div_uint32_int32(uint32_t a, int32_t b)
	uint32_t safe_div_uint32_uint32(uint32_t a, uint32_t b)
	uint32_t safe_div_uint32_int64(uint32_t a, int64_t b)
	uint32_t safe_div_uint32_uint64(uint32_t a, uint64_t b)

	int64_t safe_div_int64_int32(int64_t a, int32_t b)
	int64_t safe_div_int64_uint32(int64_t a, uint32_t b)
	int64_t safe_div_int64_int64(int64_t a, int64_t b)
	int64_t safe_div_int64_uint64(int64_t a, uint64_t b)

	uint64_t safe_div_uint64_int32(uint64_t a, int32_t b)
	uint64_t safe_div_uint64_uint32(uint64_t a, uint32_t b)
	uint64_t safe_div_uint64_int64(uint64_t a, int64_t b)
	uint64_t safe_div_uint64_uint64(uint64_t a, uint64_t b)

	// Subtraction
	int32_t safe_sub_int32_int32(int32_t a, int32_t b)
	int32_t safe_sub_int32_uint32(int32_t a, uint32_t b)
	int32_t safe_sub_int32_int64(int32_t a, int64_t b)
	int32_t safe_sub_int32_uint64(int32_t a, uint64_t b)

	uint32_t safe_sub_uint32_int32(uint32_t a, int32_t b)
	uint32_t safe_sub_uint32_uint32(uint32_t a, uint32_t b)
	uint32_t safe_sub_uint32_int64(uint32_t a, int64_t b)
	uint32_t safe_sub_uint32_uint64(uint32_t a, uint64_t b)

	int64_t safe_sub_int64_int32(int64_t a, int32_t b)
	int64_t safe_sub_int64_uint32(int64_t a, uint32_t b)
	int64_t safe_sub_int64_int64(int64_t a, int64_t b)
	int64_t safe_sub_int64_uint64(int64_t a, uint64_t b)

	uint64_t safe_sub_uint64_int32(uint64_t a, int32_t b)
	uint64_t safe_sub_uint64_uint32(uint64_t a, uint32_t b)
	uint64_t safe_sub_uint64_int64(uint64_t a, int64_t b)
	uint64_t safe_sub_uint64_uint64(uint64_t a, uint64_t b)
*/

// Do some sorting out of standard types and sizes

#if CHAR_MIN != 0
#define SAFE_MATH_SIGNED_CHAR 1
#else
#define SAFE_MATH_SIGNED_CHAR 0
#endif

#if LONG_MAX == LLONG_MAX
#define SAFE_MATH_LONG 64
#else
#define SAFE_MATH_LONG 32
#endif

// Not going to support odd sizes of things
extern char SAFE_MATH_CHECK_SHORT_IS_16[1 / ((sizeof(short)-2) ? 0 : 1)];
extern char SAFE_MATH_CHECK_INT_IS_32[1 / ((sizeof(int) - 4) ? 0 : 1)];

// In order to help keep people from making mistakes by 
// incorrectly guessing which types match which of the intXX types,
// make some functions.

// Cast to char, char might be signed or unsigned
#if SAFE_MATH_SIGNED_CHAR
static inline char safe_cast_char_int(int in) { return safe_cast_int8_int32(in); }
static inline char safe_cast_char_uint(unsigned int in) { return safe_cast_int8_uint32(in); }

static inline int check_cast_char_int(int in) { return safe_cast_int8_int32(in); }
static inline int check_cast_char_uint(unsigned int in) { return safe_cast_int8_uint32(in); }

#if SAFE_MATH_LONG == 64
static inline char safe_cast_char_long(long in) { return safe_cast_int8_int64(in); }
static inline int check_cast_char_long(long in) { return check_cast_int8_int64(in); }
#else
static inline char safe_cast_char_long(long in) { return safe_cast_int8_int32(in); }
static inline int check_cast_char_long(long in) { return check_cast_int8_int32(in); }
#endif

static inline char safe_cast_char_longlong(long long in) { return safe_cast_int8_int64(in); }
static inline char safe_cast_char_ulonglong(unsigned long long in) { return safe_cast_int8_uint64(in); }

static inline int check_cast_char_longlong(long long in) { return check_cast_int8_int64(in); }
static inline int check_cast_char_ulonglong(unsigned long long in) { return check_cast_int8_uint64(in); }
#else
static inline char safe_cast_char_int(int in) { return safe_cast_uint8_int32(in); }
static inline char safe_cast_char_uint(unsigned int in) { return safe_cast_uint8_uint32(in); }

static inline int check_cast_char_int(int in) { return check_cast_uint8_int32(in); }
static inline int check_cast_char_uint(unsigned int in) { return check_cast_uint8_uint32(in); }

#if SAFE_MATH_LONG == 64
static inline char safe_cast_char_long(long in) { return safe_cast_uint8_int64(in); }
static inline int check_cast_char_long(long in) { return check_cast_uint8_int64(in); }
#else
static inline char safe_cast_char_long(long in) { return safe_cast_uint8_int32(in); }
static inline int check_cast_char_long(long in) { return check_cast_uint8_int32(in); }
#endif

static inline char safe_cast_char_longlong(long long in) { return safe_cast_uint8_int64(in); }
static inline char safe_cast_char_ulonglong(unsigned long long in) { return safe_cast_uint8_uint64(in); }

static inline int check_cast_char_longlong(long long in) { return check_cast_uint8_int64(in); }
static inline int check_cast_char_ulonglong(unsigned long long in) { return check_cast_uint8_uint64(in); }
#endif

// Signed char
static inline signed char safe_cast_schar_int(int in) { return safe_cast_int8_int32(in); }
static inline signed char safe_cast_schar_uint(unsigned int in) { return safe_cast_int8_uint32(in); }

static inline int check_cast_schar_int(int in) { return check_cast_int8_int32(in); }
static inline int check_cast_schar_uint(unsigned int in) { return check_cast_int8_uint32(in); }

#if SAFE_MATH_LONG == 64
static inline signed char safe_cast_schar_long(long in) { return safe_cast_int8_int64(in); }
static inline int check_cast_schar_long(long in) { return check_cast_int8_int64(in); }
#else
static inline signed char safe_cast_schar_long(long in) { return safe_cast_int8_int32(in); }
static inline int check_cast_schar_long(long in) { return check_cast_int8_int32(in); }
#endif

static inline signed char safe_cast_schar_longlong(long long in) { return safe_cast_int8_int64(in); }
static inline signed char safe_cast_schar_ulonglong(unsigned long long in) { return safe_cast_int8_uint64(in); }

static inline int check_cast_schar_longlong(long long in) { return check_cast_int8_int64(in); }
static inline int check_cast_schar_ulonglong(unsigned long long in) { return check_cast_int8_uint64(in); }

// Unsigned char
static inline unsigned char safe_cast_uchar_int(int in) { return safe_cast_uint8_int32(in); }
static inline unsigned char safe_cast_uchar_uint(unsigned int in) { return safe_cast_uint8_uint32(in); }

static inline int check_cast_uchar_int(int in) { return check_cast_uint8_int32(in); }
static inline int check_cast_uchar_uint(unsigned int in) { return check_cast_uint8_uint32(in); }

#if SAFE_MATH_LONG == 64
static inline unsigned char safe_cast_uchar_long(long in) { return safe_cast_uint8_int64(in); }
static inline int check_cast_uchar_long(long in) { return check_cast_uint8_int64(in); }
#else
static inline unsigned char safe_cast_uchar_long(long in) { return safe_cast_uint8_int32(in); }
static inline int check_cast_uchar_long(long in) { return check_cast_uint8_int32(in); }
#endif

static inline unsigned char safe_cast_uchar_longlong(long long in) { return safe_cast_uint8_int64(in); }
static inline unsigned char safe_cast_uchar_ulonglong(unsigned long long in) { return safe_cast_uint8_uint64(in); }

static inline int check_cast_uchar_longlong(long long in) { return check_cast_uint8_int64(in); }
static inline int check_cast_uchar_ulonglong(unsigned long long in) { return check_cast_uint8_uint64(in); }

// 16-bit signed casting
static inline short safe_cast_short_int(int in) { return safe_cast_int16_int32(in); }
static inline short safe_cast_short_uint(unsigned int in) { return safe_cast_int16_uint32(in); }

static inline int check_cast_short_int(int in) { return check_cast_int16_int32(in); }
static inline int check_cast_short_uint(unsigned int in) { return check_cast_int16_uint32(in); }
#if SAFE_MATH_LONG == 64
static inline short safe_cast_short_long(long in) { return safe_cast_int16_int64(in); }
static inline int check_cast_short_long(long in) { return check_cast_int16_int64(in); }
#else
static inline short safe_cast_short_long(long in) { return safe_cast_int16_int32(in); }
static inline int check_cast_short_long(long in) { return check_cast_int16_int32(in); }
#endif

static inline short safe_cast_short_longlong(long long in) { return safe_cast_int16_int64(in); }
static inline short safe_cast_short_ulonglong(unsigned long long in) { return safe_cast_int16_uint64(in); }

static inline int check_cast_short_longlong(long long in) { return check_cast_int16_int64(in); }
static inline int check_cast_short_ulonglong(unsigned long long in) { return check_cast_int16_uint64(in); }

// 16-bit unsigned casting
static inline unsigned short safe_cast_ushort_int(int in) { return safe_cast_uint16_int32(in); }
static inline unsigned short safe_cast_ushort_uint(unsigned int in) { return safe_cast_uint16_uint32(in); }

static inline int check_cast_ushort_int(int in) { return check_cast_uint16_int32(in); }
static inline int check_cast_ushort_uint(unsigned int in) { return check_cast_uint16_uint32(in); }
#if SAFE_MATH_LONG == 64
static inline unsigned short safe_cast_ushort_long(long in) { return safe_cast_uint16_int64(in); }
static inline int check_cast_ushort_long(long in) { return check_cast_uint16_int64(in); }
#else
static inline unsigned short safe_cast_ushort_long(long in) { return safe_cast_uint16_int32(in); }
static inline int check_cast_ushort_long(long in) { return check_cast_uint16_int32(in); }
#endif

static inline unsigned short safe_cast_ushort_longlong(long long in) { return safe_cast_uint16_int64(in); }
static inline unsigned short safe_cast_ushort_ulonglong(unsigned long long in) { return safe_cast_uint16_uint64(in); }

static inline int check_cast_ushort_longlong(long long in) { return check_cast_uint16_int64(in); }
static inline int check_cast_ushort_ulonglong(unsigned long long in) { return check_cast_uint16_uint64(in); }

// Cast to int
static inline int safe_cast_int_uint(unsigned int in) { return safe_cast_int32_uint32(in); }
static inline int check_cast_int_uint(unsigned int in) { return check_cast_int32_uint32(in); }

#if SAFE_MATH_LONG == 64
static inline int safe_cast_int_long(long in) { return safe_cast_int32_int64(in); }
static inline int safe_cast_int_ulong(unsigned long in) { return safe_cast_int32_uint64(in); }

static inline int check_cast_int_long(long in) { return check_cast_int32_int64(in); }
static inline int check_cast_int_ulong(unsigned long in) { return check_cast_int32_uint64(in); }
#else
static inline int safe_cast_int_long(long in) { return in; }
static inline int safe_cast_int_ulong(unsigned long in) { return safe_cast_int32_uint32(in); }

static inline int check_cast_int_long(long in) { (void)in;  return 0; }
static inline int check_cast_int_ulong(unsigned long in) { return check_cast_int32_uint32(in); }
#endif

static inline int safe_cast_int_longlong(long long in) { return safe_cast_int32_int64(in); }
static inline int safe_cast_int_ulonglong(unsigned long long in) { return safe_cast_int32_uint64(in); }

static inline int check_cast_int_longlong(long long in) { return check_cast_int32_int64(in); }
static inline int check_cast_int_ulonglong(unsigned long long in) { return check_cast_int32_uint64(in); }

// Cast to unsigned int
static inline unsigned int safe_cast_uint_int(int in) { return safe_cast_uint32_int32(in); }
static inline int check_cast_uint_int(int in) { return check_cast_uint32_int32(in); }
#if SAFE_MATH_LONG == 64
static inline unsigned int safe_cast_uint_long(long in) { return safe_cast_uint32_int64(in); }
static inline unsigned int safe_cast_uint_ulong(unsigned long in) { return safe_cast_uint32_uint64(in); }

static inline int check_cast_uint_long(long in) { return check_cast_uint32_int64(in); }
static inline int check_cast_uint_ulong(unsigned long in) { return check_cast_uint32_uint64(in); }
#else
static inline unsigned int safe_cast_uint_long(long in) { return safe_cast_uint32_int32(in); }
static inline unsigned int safe_cast_uint_ulong(unsigned long in) { return in; }

static inline int check_cast_uint_long(long in) { return check_cast_uint32_int32(in); }
static inline int check_cast_uint_ulong(unsigned long in) { (void)in; return 0; }
#endif

static inline unsigned int safe_cast_uint_longlong(long long in) { return safe_cast_uint32_int64(in); }
static inline unsigned int safe_cast_uint_ulonglong(unsigned long long in) { return safe_cast_uint32_uint64(in); }

static inline int check_cast_uint_longlong(long long in) { return check_cast_uint32_int64(in); }
static inline int check_cast_uint_ulonglong(unsigned long long in) { return check_cast_uint32_uint64(in); }

// Cast to long
// Also have to keep parity in the case of different compilations
// of the same code.
#if SAFE_MATH_LONG == 64
static inline long safe_cast_long_ulong(unsigned long in) { return safe_cast_int64_uint64(in); }
static inline long safe_cast_long_longlong(long long in) { return in; }
static inline long safe_cast_long_ulonglong(unsigned long long in) { return safe_cast_int64_uint64(in); }

static inline int check_cast_long_ulong(unsigned long in) { return check_cast_int64_uint64(in); }
static inline int check_cast_long_longlong(long long in) { (void)in; return 0; }
static inline int check_cast_long_ulonglong(unsigned long long in) { return check_cast_int64_uint64(in); }

static inline unsigned long safe_cast_ulong_long(long in) { return safe_cast_uint64_int64(in); }
static inline unsigned long safe_cast_ulong_ulonglong(unsigned long long in) { return in; }
static inline unsigned long safe_cast_ulong_longlong(long long in) { return safe_cast_uint64_int64(in); }

static inline int check_cast_ulong_long(long in) { return check_cast_uint64_int64(in); }
static inline int check_cast_ulong_ulonglong(unsigned long long in) { (void)in; return 0; }
static inline int check_cast_ulong_longlong(long long in) { return check_cast_uint64_int64(in); }
#else
static inline long safe_cast_long_ulong(unsigned long in) { return safe_cast_int32_uint32(in); }
static inline long safe_cast_long_longlong(long long in) { return safe_cast_int32_int64(in); }
static inline long safe_cast_long_ulonglong(unsigned long long in) { return safe_cast_int32_uint64(in); }

static inline int check_cast_long_ulong(unsigned long in) { return check_cast_int32_uint32(in); }
static inline int check_cast_long_longlong(long long in) { return check_cast_int32_int64(in); }
static inline int check_cast_long_ulonglong(unsigned long long in) { return check_cast_int32_uint64(in); }

static inline unsigned long safe_cast_ulong_long(long in) { return safe_cast_uint32_int32(in); }
static inline unsigned long safe_cast_ulong_ulonglong(unsigned long long in) { return safe_cast_uint32_uint64(in); }
static inline unsigned long safe_cast_ulong_longlong(long long in) { return safe_cast_uint32_int64(in); }

static inline int check_cast_ulong_long(long in) { return check_cast_uint32_int32(in); }
static inline int check_cast_ulong_ulonglong(unsigned long long in) { return check_cast_uint32_uint64(in); }
static inline int check_cast_ulong_longlong(long long in) { return check_cast_uint32_int64(in); }
#endif

// And long long
static inline long long safe_cast_longlong_ulonglong(unsigned long long in) { return safe_cast_int64_uint64(in); }
static inline unsigned long long safe_cast_ulonglong_longlong(long long in) { return safe_cast_uint64_int64(in); }

static inline int check_cast_longlong_ulonglong(unsigned long long in) { return check_cast_int64_uint64(in); }
static inline int check_cast_ulonglong_longlong(long long in) { return check_cast_uint64_int64(in); }

// Addition
static inline int safe_add_int_int(int a, int b) { return safe_add_int32_int32(a, b); }
static inline int safe_add_int_uint(int a, unsigned int b) { return safe_add_int32_uint32(a, b); }
static inline int safe_add_int_longlong(int a, long long b) { return safe_add_int32_int64(a, b); }
static inline int safe_add_int_ulonglong(int a, unsigned long long b) { return safe_add_int32_uint64(a, b); }

static inline bool check_add_int_int(int a, int b, int* ret) { return check_add_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_add_int_uint(int a, unsigned int b, int* ret) { return check_add_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_add_int_longlong(int a, long long b, int* ret) { return check_add_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_add_int_ulonglong(int a, unsigned long long b, int* ret) { return check_add_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned int safe_add_uint_int(unsigned int a, int b) { return safe_add_uint32_int32(a, b); }
static inline unsigned int safe_add_uint_uint(unsigned int a, unsigned int b) { return safe_add_uint32_uint32(a, b); }
static inline unsigned int safe_add_uint_longlong(unsigned int a, long long b) { return safe_add_uint32_int64(a, b); }
static inline unsigned int safe_add_uint_ulonglong(unsigned int a, unsigned long long b) { return safe_add_uint32_uint64(a, b); }

static inline bool check_add_uint_int(unsigned int a, int b, unsigned int* ret) { return check_add_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_add_uint_uint(unsigned int a, unsigned int b, unsigned int* ret) { return check_add_uint32_uint32(a, b, (uint32_t*)ret); }
static inline bool check_add_uint_longlong(unsigned int a, long long b, unsigned int* ret) { return check_add_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_add_uint_ulonglong(unsigned int a, unsigned long long b, unsigned int* ret) { return check_add_uint32_uint64(a, b, (uint32_t*)ret); }

#if SAFE_MATH_LONG == 64
static inline int safe_add_int_long(int a, long b) { return safe_add_int32_int64(a, b); }
static inline int safe_add_int_ulong(int a, unsigned long b) { return safe_add_int32_uint64(a, b); }

static inline bool check_add_int_long(int a, long b, int* ret) { return check_add_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_add_int_ulong(int a, unsigned long b, int* ret) { return check_add_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned int safe_add_uint_long(unsigned int a, long b) { return safe_add_uint32_int64(a, b); }
static inline unsigned int safe_add_uint_ulong(unsigned int a, unsigned long b) { return safe_add_uint32_uint64(a, b); }

static inline bool check_add_uint_long(unsigned int a, long b, unsigned int* ret) { return check_add_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_add_uint_ulong(unsigned int a, unsigned long b, unsigned int* ret) { return check_add_uint32_uint64(a, b, (uint32_t*)ret); }

static inline long safe_add_long_int(long a, int b) { return safe_add_int64_int32(a, b); }
static inline long safe_add_long_uint(long a, unsigned int b) { return safe_add_int64_uint32(a, b); }
static inline long safe_add_long_long(long a, long b) { return safe_add_int64_int64(a, b); }
static inline long safe_add_long_ulong(long a, unsigned long b) { return safe_add_int64_uint64(a, b); }
static inline long safe_add_long_longlong(long a, long long b) { return safe_add_int64_int64(a, b); }
static inline long safe_add_long_ulonglong(long a, unsigned long long b) { return safe_add_int64_uint64(a, b); }

static inline bool check_add_long_int(long a, int b, long* ret) { return check_add_int64_int32(a, b, (int64_t*)ret); }
static inline bool check_add_long_uint(long a, unsigned int b, long* ret) { return check_add_int64_uint32(a, b, (int64_t*)ret); }
static inline bool check_add_long_long(long a, long b, long* ret) { return check_add_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_add_long_ulong(long a, unsigned long b, long* ret) { return check_add_int64_uint64(a, b, (int64_t*)ret); }
static inline bool check_add_long_longlong(long a, long long b, long* ret) { return check_add_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_add_long_ulonglong(long a, unsigned long long b, long* ret) { return check_add_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long safe_add_ulong_int(unsigned long a, int b) { return safe_add_uint64_int32(a, b); }
static inline unsigned long safe_add_ulong_uint(unsigned long a, unsigned int b) { return safe_add_uint64_uint32(a, b); }
static inline unsigned long safe_add_ulong_long(unsigned long a, long b) { return safe_add_uint64_int64(a, b); }
static inline unsigned long safe_add_ulong_ulong(unsigned long a, unsigned long b) { return safe_add_uint64_uint64(a, b); }
static inline unsigned long safe_add_ulong_longlong(unsigned long a, long long b) { return safe_add_uint64_int64(a, b); }
static inline unsigned long safe_add_ulong_ulonglong(unsigned long a, unsigned long long b) { return safe_add_uint64_uint64(a, b); }

static inline bool check_add_ulong_int(unsigned long a, int b, unsigned long* ret) { return check_add_uint64_int32(a, b, (uint64_t*)ret); }
static inline bool check_add_ulong_uint(unsigned long a, unsigned int b, unsigned long* ret) { return check_add_uint64_uint32(a, b, (uint64_t*)ret); }
static inline bool check_add_ulong_long(unsigned long a, long b, unsigned long* ret) { return check_add_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_add_ulong_ulong(unsigned long a, unsigned long b, unsigned long* ret) { return check_add_uint64_uint64(a, b, (uint64_t*)ret); }
static inline bool check_add_ulong_longlong(unsigned long a, long long b, unsigned long* ret) { return check_add_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_add_ulong_ulonglong(unsigned long a, unsigned long long b, unsigned long* ret) { return check_add_uint64_uint64(a, b, (uint64_t*)ret); }

static inline long long safe_add_longlong_long(long long a, long b) { return safe_add_int64_int64(a, b); }
static inline long long safe_add_longlong_ulong(long long a, unsigned long b) { return safe_add_int64_uint64(a, b); }

static inline bool check_add_longlong_long(long long a, long b, long long* ret) { return check_add_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_add_longlong_ulong(long long a, unsigned long b, long long* ret) { return check_add_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long long safe_add_ulonglong_long(unsigned long long a, long b) { return safe_add_uint64_int64(a, b); }
static inline unsigned long long safe_add_ulonglong_ulong(unsigned long long a, unsigned long b) { return safe_add_uint64_uint64(a, b); }

static inline bool check_add_ulonglong_long(unsigned long long a, long b, unsigned long long* ret) { return check_add_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_add_ulonglong_ulong(unsigned long long a, unsigned long b, unsigned long long* ret) { return check_add_uint64_uint64(a, b, (uint64_t*)ret); }
#else
static inline int safe_add_int_long(int a, long b) { return safe_add_int32_int32(a, b); }
static inline int safe_add_int_ulong(int a, unsigned long b) { return safe_add_int32_uint32(a, b); }

static inline bool check_add_int_long(int a, long b, int* ret) { return check_add_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_add_int_ulong(int a, unsigned long b, int* ret) { return check_add_int32_uint32(a, b, (int32_t*)ret); }

static inline unsigned int safe_add_uint_long(unsigned int a, long b) { return safe_add_uint32_int32(a, b); }
static inline unsigned int safe_add_uint_ulong(unsigned int a, unsigned long b) { return safe_add_uint32_uint32(a, b); }

static inline bool check_add_uint_long(unsigned int a, long b, unsigned int* ret) { return check_add_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_add_uint_ulong(unsigned int a, unsigned long b, unsigned int* ret) { return check_add_uint32_uint32(a, b, (uint32_t*)ret); }

static inline long safe_add_long_int(long a, int b) { return safe_add_int32_int32(a, b); }
static inline long safe_add_long_uint(long a, unsigned int b) { return safe_add_int32_uint32(a, b); }
static inline long safe_add_long_long(long a, long b) { return safe_add_int32_int32(a, b); }
static inline long safe_add_long_ulong(long a, unsigned long b) { return safe_add_int32_uint32(a, b); }
static inline long safe_add_long_longlong(long a, long long b) { return safe_add_int32_int64(a, b); }
static inline long safe_add_long_ulonglong(long a, unsigned long long b) { return safe_add_int32_uint64(a, b); }

static inline bool check_add_long_int(long a, int b, long* ret) { return check_add_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_add_long_uint(long a, unsigned int b, long* ret) { return check_add_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_add_long_long(long a, long b, long* ret) { return check_add_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_add_long_ulong(long a, unsigned long b, long* ret) { return check_add_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_add_long_longlong(long a, long long b, long* ret) { return check_add_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_add_long_ulonglong(long a, unsigned long long b, long* ret) { return check_add_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned long safe_add_ulong_int(unsigned long a, int b) { return safe_add_uint32_int32(a, b); }
static inline unsigned long safe_add_ulong_uint(unsigned long a, unsigned int b) { return safe_add_uint32_uint32(a, b); }
static inline unsigned long safe_add_ulong_long(unsigned long a, long b) { return safe_add_uint32_int32(a, b); }
static inline unsigned long safe_add_ulong_ulong(unsigned long a, unsigned long b) { return safe_add_uint32_uint32(a, b); }
static inline unsigned long safe_add_ulong_longlong(unsigned long a, long long b) { return safe_add_uint32_int64(a, b); }
static inline unsigned long safe_add_ulong_ulonglong(unsigned long a, unsigned long long b) { return safe_add_uint32_uint64(a, b); }

static inline bool check_add_ulong_int(unsigned long a, int b, unsigned long* ret) { return check_add_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_add_ulong_uint(unsigned long a, unsigned int b, unsigned long* ret) { return check_add_uint32_uint32(a, b, (uint32_t*)ret); }
static inline bool check_add_ulong_long(unsigned long a, long b, unsigned long* ret) { return check_add_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_add_ulong_ulong(unsigned long a, unsigned long b, unsigned long* ret) { return check_add_uint32_uint32(a, b, (uint32_t*)ret); }
static inline bool check_add_ulong_longlong(unsigned long a, long long b, unsigned long* ret) { return check_add_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_add_ulong_ulonglong(unsigned long a, unsigned long long b, unsigned long* ret) { return check_add_uint32_uint64(a, b, (uint32_t*)ret); }

static inline long long safe_add_longlong_long(long long a, long b) { return safe_add_int64_int32(a, b); }
static inline long long safe_add_longlong_ulong(long long a, unsigned long b) { return safe_add_int64_uint32(a, b); }

static inline bool check_add_longlong_long(long long a, long b, long long* ret) { return check_add_int64_int32(a, b, ret); }
static inline bool check_add_longlong_ulong(long long a, unsigned long b, long long* ret) { return check_add_int64_uint32(a, b, ret); }

static inline unsigned long long safe_add_ulonglong_long(unsigned long long a, long b) { return safe_add_uint64_int32(a, b); }
static inline unsigned long long safe_add_ulonglong_ulong(unsigned long long a, unsigned long b) { return safe_add_uint64_uint32(a, b); }

static inline bool check_add_ulonglong_long(unsigned long long a, long b, unsigned long long* ret) { return check_add_uint64_int32(a, b, ret); }
static inline bool check_add_ulonglong_ulong(unsigned long long a, unsigned long b, unsigned long long* ret) { return check_add_uint64_uint32(a, b, ret); }
#endif

static inline long long safe_add_longlong_int(long long a, int b) { return safe_add_int64_int32(a, b); }
static inline long long safe_add_longlong_uint(long long a, unsigned int b) { return safe_add_int64_uint32(a, b); }
static inline long long safe_add_longlong_longlong(long long a, long long b) { return safe_add_int64_int64(a, b); }
static inline long long safe_add_longlong_ulonglong(long long a, unsigned long long b) { return safe_add_int64_uint64(a, b); }

static inline bool check_add_longlong_int(long long a, int b, long long* ret) { return check_add_int64_int32(a, b, (int64_t*)ret); }
static inline bool check_add_longlong_uint(long long a, unsigned int b, long long* ret) { return check_add_int64_uint32(a, b, (int64_t*)ret); }
static inline bool check_add_longlong_longlong(long long a, long long b, long long* ret) { return check_add_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_add_longlong_ulonglong(long long a, unsigned long long b, long long* ret) { return check_add_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long long safe_add_ulonglong_int(unsigned long long a, int b) { return safe_add_uint64_int32(a, b); }
static inline unsigned long long safe_add_ulonglong_uint(unsigned long long a, unsigned int b) { return safe_add_uint64_uint32(a, b); }
static inline unsigned long long safe_add_ulonglong_longlong(unsigned long long a, long long b) { return safe_add_uint64_int64(a, b); }
static inline unsigned long long safe_add_ulonglong_ulonglong(unsigned long long a, unsigned long long b) { return safe_add_uint64_uint64(a, b); }

static inline bool check_add_ulonglong_int(unsigned long long a, int b, unsigned long long* ret) { return check_add_uint64_int32(a, b, (uint64_t*)ret); }
static inline bool check_add_ulonglong_uint(unsigned long long a, unsigned int b, unsigned long long* ret) { return check_add_uint64_uint32(a, b, (uint64_t*)ret); }
static inline bool check_add_ulonglong_longlong(unsigned long long a, long long b, unsigned long long* ret) { return check_add_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_add_ulonglong_ulonglong(unsigned long long a, unsigned long long b, unsigned long long* ret) { return check_add_uint64_uint64(a, b, (uint64_t*)ret); }

// Multiplication
static inline int safe_mul_int_int(int a, int b) { return safe_mul_int32_int32(a, b); }
static inline int safe_mul_int_uint(int a, unsigned int b) { return safe_mul_int32_uint32(a, b); }
static inline int safe_mul_int_longlong(int a, long long b) { return safe_mul_int32_int64(a, b); }
static inline int safe_mul_int_ulonglong(int a, unsigned long long b) { return safe_mul_int32_uint64(a, b); }

static inline bool check_mul_int_int(int a, int b, int* ret) { return check_mul_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_mul_int_uint(int a, unsigned int b, int* ret) { return check_mul_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_mul_int_longlong(int a, long long b, int* ret) { return check_mul_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_mul_int_ulonglong(int a, unsigned long long b, int* ret) { return check_mul_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned int safe_mul_uint_int(unsigned int a, int b) { return safe_mul_uint32_int32(a, b); }
static inline unsigned int safe_mul_uint_uint(unsigned int a, unsigned int b) { return safe_mul_uint32_uint32(a, b); }
static inline unsigned int safe_mul_uint_longlong(unsigned int a, long long b) { return safe_mul_uint32_int64(a, b); }
static inline unsigned int safe_mul_uint_ulonglong(unsigned int a, unsigned long long b) { return safe_mul_uint32_uint64(a, b); }

static inline bool check_mul_uint_int(unsigned int a, int b, unsigned int* ret) { return check_mul_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_mul_uint_uint(unsigned int a, unsigned int b, unsigned int* ret) { return check_mul_uint32_uint32(a, b, (uint32_t*)ret); }
static inline bool check_mul_uint_longlong(unsigned int a, long long b, unsigned int* ret) { return check_mul_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_mul_uint_ulonglong(unsigned int a, unsigned long long b, unsigned int* ret) { return check_mul_uint32_uint64(a, b, (uint32_t*)ret); }

#if SAFE_MATH_LONG == 64
static inline int safe_mul_int_long(int a, long b) { return safe_mul_int32_int64(a, b); }
static inline int safe_mul_int_ulong(int a, unsigned long b) { return safe_mul_int32_uint64(a, b); }

static inline bool check_mul_int_long(int a, long b, int* ret) { return check_mul_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_mul_int_ulong(int a, unsigned long b, int* ret) { return check_mul_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned int safe_mul_uint_long(unsigned int a, long b) { return safe_mul_uint32_int64(a, b); }
static inline unsigned int safe_mul_uint_ulong(unsigned int a, unsigned long b) { return safe_mul_uint32_uint64(a, b); }

static inline bool check_mul_uint_long(unsigned int a, long b, unsigned int* ret) { return check_mul_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_mul_uint_ulong(unsigned int a, unsigned long b, unsigned int* ret) { return check_mul_uint32_uint64(a, b, (uint32_t*)ret); }

static inline long safe_mul_long_int(long a, int b) { return safe_mul_int64_int32(a, b); }
static inline long safe_mul_long_uint(long a, unsigned int b) { return safe_mul_int64_uint32(a, b); }
static inline long safe_mul_long_long(long a, long b) { return safe_mul_int64_int64(a, b); }
static inline long safe_mul_long_ulong(long a, unsigned long b) { return safe_mul_int64_uint64(a, b); }
static inline long safe_mul_long_longlong(long a, long long b) { return safe_mul_int64_int64(a, b); }
static inline long safe_mul_long_ulonglong(long a, unsigned long long b) { return safe_mul_int64_uint64(a, b); }

static inline bool check_mul_long_int(long a, int b, long* ret) { return check_mul_int64_int32(a, b, (int64_t*)ret); }
static inline bool check_mul_long_uint(long a, unsigned int b, long* ret) { return check_mul_int64_uint32(a, b, (int64_t*)ret); }
static inline bool check_mul_long_long(long a, long b, long* ret) { return check_mul_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_mul_long_ulong(long a, unsigned long b, long* ret) { return check_mul_int64_uint64(a, b, (int64_t*)ret); }
static inline bool check_mul_long_longlong(long a, long long b, long* ret) { return check_mul_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_mul_long_ulonglong(long a, unsigned long long b, long* ret) { return check_mul_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long safe_mul_ulong_int(unsigned long a, int b) { return safe_mul_uint64_int32(a, b); }
static inline unsigned long safe_mul_ulong_uint(unsigned long a, unsigned int b) { return safe_mul_uint64_uint32(a, b); }
static inline unsigned long safe_mul_ulong_long(unsigned long a, long b) { return safe_mul_uint64_int64(a, b); }
static inline unsigned long safe_mul_ulong_ulong(unsigned long a, unsigned long b) { return safe_mul_uint64_uint64(a, b); }
static inline unsigned long safe_mul_ulong_longlong(unsigned long a, long long b) { return safe_mul_uint64_int64(a, b); }
static inline unsigned long safe_mul_ulong_ulonglong(unsigned long a, unsigned long long b) { return safe_mul_uint64_uint64(a, b); }

static inline bool check_mul_ulong_int(unsigned long a, int b, unsigned long* ret) { return check_mul_uint64_int32(a, b, (uint64_t*)ret); }
static inline bool check_mul_ulong_uint(unsigned long a, unsigned int b, unsigned long* ret) { return check_mul_uint64_uint32(a, b, (uint64_t*)ret); }
static inline bool check_mul_ulong_long(unsigned long a, long b, unsigned long* ret) { return check_mul_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_mul_ulong_ulong(unsigned long a, unsigned long b, unsigned long* ret) { return check_mul_uint64_uint64(a, b, (uint64_t*)ret); }
static inline bool check_mul_ulong_longlong(unsigned long a, long long b, unsigned long* ret) { return check_mul_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_mul_ulong_ulonglong(unsigned long a, unsigned long long b, unsigned long* ret) { return check_mul_uint64_uint64(a, b,(uint64_t*)ret); }

static inline long long safe_mul_longlong_long(long long a, long b) { return safe_mul_int64_int64(a, b); }
static inline long long safe_mul_longlong_ulong(long long a, unsigned long b) { return safe_mul_int64_uint64(a, b); }

static inline bool check_mul_longlong_long(long long a, long b, long long* ret) { return check_mul_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_mul_longlong_ulong(long long a, unsigned long b, long long* ret) { return check_mul_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long long safe_mul_ulonglong_long(unsigned long long a, long b) { return safe_mul_uint64_int64(a, b); }
static inline unsigned long long safe_mul_ulonglong_ulong(unsigned long long a, unsigned long b) { return safe_mul_uint64_uint64(a, b); }

static inline bool check_mul_ulonglong_long(unsigned long long a, long b, unsigned long long* ret) { return check_mul_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_mul_ulonglong_ulong(unsigned long long a, unsigned long b, unsigned long long* ret) { return check_mul_uint64_uint64(a, b, (uint64_t*)ret); }
#else
static inline int safe_mul_int_long(int a, long b) { return safe_mul_int32_int32(a, b); }
static inline int safe_mul_int_ulong(int a, unsigned long b) { return safe_mul_int32_uint32(a, b); }

static inline bool check_mul_int_long(int a, long b, int* ret) { return check_mul_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_mul_int_ulong(int a, unsigned long b, int* ret) { return check_mul_int32_uint32(a, b, (int32_t*)ret); }

static inline unsigned int safe_mul_uint_long(unsigned int a, long b) { return safe_mul_uint32_int32(a, b); }
static inline unsigned int safe_mul_uint_ulong(unsigned int a, unsigned long b) { return safe_mul_uint32_uint32(a, b); }

static inline bool check_mul_uint_long(unsigned int a, long b, unsigned int* ret) { return check_mul_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_mul_uint_ulong(unsigned int a, unsigned long b, unsigned int* ret) { return check_mul_uint32_uint32(a, b, (uint32_t*)ret); }

static inline long safe_mul_long_int(long a, int b) { return safe_mul_int32_int32(a, b); }
static inline long safe_mul_long_uint(long a, unsigned int b) { return safe_mul_int32_uint32(a, b); }
static inline long safe_mul_long_long(long a, long b) { return safe_mul_int32_int32(a, b); }
static inline long safe_mul_long_ulong(long a, unsigned long b) { return safe_mul_int32_uint32(a, b); }
static inline long safe_mul_long_longlong(long a, long long b) { return safe_mul_int32_int64(a, b); }
static inline long safe_mul_long_ulonglong(long a, unsigned long long b) { return safe_mul_int32_uint64(a, b); }

static inline bool check_mul_long_int(long a, int b, long* ret) { return check_mul_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_mul_long_uint(long a, unsigned int b, long* ret) { return check_mul_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_mul_long_long(long a, long b, long* ret) { return check_mul_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_mul_long_ulong(long a, unsigned long b, long* ret) { return check_mul_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_mul_long_longlong(long a, long long b, long* ret) { return check_mul_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_mul_long_ulonglong(long a, unsigned long long b, long* ret) { return check_mul_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned long safe_mul_ulong_int(unsigned long a, int b) { return safe_mul_uint32_int32(a, b); }
static inline unsigned long safe_mul_ulong_uint(unsigned long a, unsigned int b) { return safe_mul_uint32_uint32(a, b); }
static inline unsigned long safe_mul_ulong_long(unsigned long a, long b) { return safe_mul_uint32_int32(a, b); }
static inline unsigned long safe_mul_ulong_ulong(unsigned long a, unsigned long b) { return safe_mul_uint32_uint32(a, b); }
static inline unsigned long safe_mul_ulong_longlong(unsigned long a, long long b) { return safe_mul_uint32_int64(a, b); }
static inline unsigned long safe_mul_ulong_ulonglong(unsigned long a, unsigned long long b) { return safe_mul_uint32_uint64(a, b); }

static inline bool check_mul_ulong_int(unsigned long a, int b, unsigned long* ret) { return check_mul_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_mul_ulong_uint(unsigned long a, unsigned int b, unsigned long* ret) { return check_mul_uint32_uint32(a, b, (uint32_t*)ret); }
static inline bool check_mul_ulong_long(unsigned long a, long b, unsigned long* ret) { return check_mul_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_mul_ulong_ulong(unsigned long a, unsigned long b, unsigned long* ret) { return check_mul_uint32_uint32(a, b, (uint32_t*)ret); }
static inline bool check_mul_ulong_longlong(unsigned long a, long long b, unsigned long* ret) { return check_mul_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_mul_ulong_ulonglong(unsigned long a, unsigned long long b, unsigned long* ret) { return check_mul_uint32_uint64(a, b, (uint32_t*)ret); }

static inline long long safe_mul_longlong_long(long long a, long b) { return safe_mul_int64_int32(a, b); }
static inline long long safe_mul_longlong_ulong(long long a, unsigned long b) { return safe_mul_int64_uint32(a, b); }

static inline bool check_mul_longlong_long(long long a, long b, long long* ret) { return check_mul_int64_int32(a, b, ret); }
static inline bool check_mul_longlong_ulong(long long a, unsigned long b, long long* ret) { return check_mul_int64_uint64(a, b, ret); }

static inline unsigned long long safe_mul_ulonglong_long(unsigned long long a, long b) { return safe_mul_uint64_int32(a, b); }
static inline unsigned long long safe_mul_ulonglong_ulong(unsigned long long a, unsigned long b) { return safe_mul_uint64_uint32(a, b); }

static inline bool check_mul_ulonglong_long(unsigned long long a, long b, unsigned long long* ret) { return check_mul_uint64_int32(a, b, ret); }
static inline bool check_mul_ulonglong_ulong(unsigned long long a, unsigned long b, unsigned long long* ret) { return check_mul_uint64_uint32(a, b, ret); }
#endif

static inline long long safe_mul_longlong_int(long long a, int b) { return safe_mul_int64_int32(a, b); }
static inline long long safe_mul_longlong_uint(long long a, unsigned int b) { return safe_mul_int64_uint32(a, b); }
static inline long long safe_mul_longlong_longlong(long long a, long long b) { return safe_mul_int64_int64(a, b); }
static inline long long safe_mul_longlong_ulonglong(long long a, unsigned long long b) { return safe_mul_int64_uint64(a, b); }

static inline bool check_mul_longlong_int(long long a, int b, long long* ret) { return check_mul_int64_int32(a, b, (int64_t*)ret); }
static inline bool check_mul_longlong_uint(long long a, unsigned int b, long long* ret) { return check_mul_int64_uint32(a, b, (int64_t*)ret); }
static inline bool check_mul_longlong_longlong(long long a, long long b, long long* ret) { return check_mul_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_mul_longlong_ulonglong(long long a, unsigned long long b, long long* ret) { return check_mul_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long long safe_mul_ulonglong_int(unsigned long long a, int b) { return safe_mul_uint64_int32(a, b); }
static inline unsigned long long safe_mul_ulonglong_uint(unsigned long long a, unsigned int b) { return safe_mul_uint64_uint32(a, b); }
static inline unsigned long long safe_mul_ulonglong_longlong(unsigned long long a, long long b) { return safe_mul_uint64_int64(a, b); }
static inline unsigned long long safe_mul_ulonglong_ulonglong(unsigned long long a, unsigned long long b) { return safe_mul_uint64_uint64(a, b); }

static inline bool check_mul_ulonglong_int(unsigned long long a, int b, unsigned long long* ret) { return check_mul_uint64_int32(a, b, (uint64_t*)ret); }
static inline bool check_mul_ulonglong_uint(unsigned long long a, unsigned int b, unsigned long long* ret) { return check_mul_uint64_uint32(a, b, (uint64_t*)ret); }
static inline bool check_mul_ulonglong_longlong(unsigned long long a, long long b, unsigned long long* ret) { return check_mul_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_mul_ulonglong_ulonglong(unsigned long long a, unsigned long long b, unsigned long long* ret) { return check_mul_uint64_uint64(a, b, (uint64_t*)ret); }

// Subtraction
static inline int safe_sub_int_int(int a, int b) { return safe_sub_int32_int32(a, b); }
static inline int safe_sub_int_uint(int a, unsigned int b) { return safe_sub_int32_uint32(a, b); }
static inline int safe_sub_int_longlong(int a, long long b) { return safe_sub_int32_int64(a, b); }
static inline int safe_sub_int_ulonglong(int a, unsigned long long b) { return safe_sub_int32_uint64(a, b); }

static inline bool check_sub_int_int(int a, int b, int* ret) { return check_sub_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_sub_int_uint(int a, unsigned int b, int* ret) { return check_sub_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_sub_int_longlong(int a, long long b, int* ret) { return check_sub_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_sub_int_ulonglong(int a, unsigned long long b, int* ret) { return check_sub_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned int safe_sub_uint_int(unsigned int a, int b) { return safe_sub_uint32_int32(a, b); }
static inline unsigned int safe_sub_uint_uint(unsigned int a, unsigned int b) { return safe_sub_uint32_uint32(a, b); }
static inline unsigned int safe_sub_uint_longlong(unsigned int a, long long b) { return safe_sub_uint32_int64(a, b); }
static inline unsigned int safe_sub_uint_ulonglong(unsigned int a, unsigned long long b) { return safe_sub_uint32_uint64(a, b); }

static inline bool check_sub_uint_int(unsigned int a, int b, unsigned int* ret) { return check_sub_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_sub_uint_uint(unsigned int a, unsigned int b, unsigned int* ret) { return check_sub_uint32_uint32(a, b, (uint32_t*)ret); }
static inline bool check_sub_uint_longlong(unsigned int a, long long b, unsigned int* ret) { return check_sub_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_sub_uint_ulonglong(unsigned int a, unsigned long long b, unsigned int* ret) { return check_sub_uint32_uint64(a, b, (uint32_t*)ret); }

#if SAFE_MATH_LONG == 64
static inline int safe_sub_int_long(int a, long b) { return safe_sub_int32_int64(a, b); }
static inline int safe_sub_int_ulong(int a, unsigned long b) { return safe_sub_int32_uint64(a, b); }

static inline bool check_sub_int_long(int a, long b, int* ret) { return check_sub_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_sub_int_ulong(int a, unsigned long b, int* ret) { return check_sub_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned int safe_sub_uint_long(unsigned int a, long b) { return safe_sub_uint32_int64(a, b); }
static inline unsigned int safe_sub_uint_ulong(unsigned int a, unsigned long b) { return safe_sub_uint32_uint64(a, b); }

static inline bool check_sub_uint_long(unsigned int a, long b, unsigned int* ret) { return check_sub_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_sub_uint_ulong(unsigned int a, unsigned long b, unsigned int* ret) { return check_sub_uint32_uint64(a, b, (uint32_t*)ret); }

static inline long safe_sub_long_int(long a, int b) { return safe_sub_int64_int32(a, b); }
static inline long safe_sub_long_uint(long a, unsigned int b) { return safe_sub_int64_uint32(a, b); }
static inline long safe_sub_long_long(long a, long b) { return safe_sub_int64_int64(a, b); }
static inline long safe_sub_long_ulong(long a, unsigned long b) { return safe_sub_int64_uint64(a, b); }
static inline long safe_sub_long_longlong(long a, long long b) { return safe_sub_int64_int64(a, b); }
static inline long safe_sub_long_ulonglong(long a, unsigned long long b) { return safe_sub_int64_uint64(a, b); }

static inline bool check_sub_long_int(long a, int b, long* ret) { return check_sub_int64_int32(a, b, (int64_t*)ret); }
static inline bool check_sub_long_uint(long a, unsigned int b, long* ret) { return check_sub_int64_uint32(a, b, (int64_t*)ret); }
static inline bool check_sub_long_long(long a, long b, long* ret) { return check_sub_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_sub_long_ulong(long a, unsigned long b, long* ret) { return check_sub_int64_uint64(a, b, (int64_t*)ret); }
static inline bool check_sub_long_longlong(long a, long long b, long* ret) { return check_sub_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_sub_long_ulonglong(long a, unsigned long long b, long* ret) { return check_sub_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long safe_sub_ulong_int(unsigned long a, int b) { return safe_sub_uint64_int32(a, b); }
static inline unsigned long safe_sub_ulong_uint(unsigned long a, unsigned int b) { return safe_sub_uint64_uint32(a, b); }
static inline unsigned long safe_sub_ulong_long(unsigned long a, long b) { return safe_sub_uint64_int64(a, b); }
static inline unsigned long safe_sub_ulong_ulong(unsigned long a, unsigned long b) { return safe_sub_uint64_uint64(a, b); }
static inline unsigned long safe_sub_ulong_longlong(unsigned long a, long long b) { return safe_sub_uint64_int64(a, b); }
static inline unsigned long safe_sub_ulong_ulonglong(unsigned long a, unsigned long long b) { return safe_sub_uint64_uint64(a, b); }

static inline bool check_sub_ulong_int(unsigned long a, int b, unsigned long* ret) { return check_sub_uint64_int32(a, b, (uint64_t*)ret); }
static inline bool check_sub_ulong_uint(unsigned long a, unsigned int b, unsigned long* ret) { return check_sub_uint64_uint32(a, b, (uint64_t*)ret); }
static inline bool check_sub_ulong_long(unsigned long a, long b, unsigned long* ret) { return check_sub_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_sub_ulong_ulong(unsigned long a, unsigned long b, unsigned long* ret) { return check_sub_uint64_uint64(a, b, (uint64_t*)ret); }
static inline bool check_sub_ulong_longlong(unsigned long a, long long b, unsigned long* ret) { return check_sub_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_sub_ulong_ulonglong(unsigned long a, unsigned long long b, unsigned long* ret) { return check_sub_uint64_uint64(a, b, (uint64_t*)ret); }

static inline long long safe_sub_longlong_long(long long a, long b) { return safe_sub_int64_int64(a, b); }
static inline long long safe_sub_longlong_ulong(long long a, unsigned long b) { return safe_sub_int64_uint64(a, b); }

static inline bool check_sub_longlong_long(long long a, long b, long long* ret) { return check_sub_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_sub_longlong_ulong(long long a, unsigned long b, long long* ret) { return check_sub_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long long safe_sub_ulonglong_long(unsigned long long a, long b) { return safe_sub_uint64_int64(a, b); }
static inline unsigned long long safe_sub_ulonglong_ulong(unsigned long long a, unsigned long b) { return safe_sub_uint64_uint64(a, b); }

static inline bool check_sub_ulonglong_long(unsigned long long a, long b, unsigned long long* ret) { return check_sub_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_sub_ulonglong_ulong(unsigned long long a, unsigned long b, unsigned long long* ret) { return check_sub_uint64_uint64(a, b, (uint64_t*)ret); }
#else
static inline int safe_sub_int_long(int a, long b) { return safe_sub_int32_int32(a, b); }
static inline int safe_sub_int_ulong(int a, unsigned long b) { return safe_sub_int32_uint32(a, b); }

static inline bool check_sub_int_long(int a, long b, int* ret) { return check_sub_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_sub_int_ulong(int a, unsigned long b, int* ret) { return check_sub_int32_uint32(a, b, (int32_t*)ret); }

static inline unsigned int safe_sub_uint_long(unsigned int a, long b) { return safe_sub_uint32_int32(a, b); }
static inline unsigned int safe_sub_uint_ulong(unsigned int a, unsigned long b) { return safe_sub_uint32_uint32(a, b); }

static inline bool check_sub_uint_long(unsigned int a, long b, unsigned int* ret) { return check_sub_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_sub_uint_ulong(unsigned int a, unsigned long b, unsigned int* ret) { return check_sub_uint32_uint32(a, b, (uint32_t*)ret); }

static inline long safe_sub_long_int(long a, int b) { return safe_sub_int32_int32(a, b); }
static inline long safe_sub_long_uint(long a, unsigned int b) { return safe_sub_int32_uint32(a, b); }
static inline long safe_sub_long_long(long a, long b) { return safe_sub_int32_int32(a, b); }
static inline long safe_sub_long_ulong(long a, unsigned long b) { return safe_sub_int32_uint32(a, b); }
static inline long safe_sub_long_longlong(long a, long long b) { return safe_sub_int32_int64(a, b); }
static inline long safe_sub_long_ulonglong(long a, unsigned long long b) { return safe_sub_int32_uint64(a, b); }

static inline bool check_sub_long_int(long a, int b, long* ret) { return check_sub_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_sub_long_uint(long a, unsigned int b, long* ret) { return check_sub_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_sub_long_long(long a, long b, long* ret) { return check_sub_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_sub_long_ulong(long a, unsigned long b, long* ret) { return check_sub_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_sub_long_longlong(long a, long long b, long* ret) { return check_sub_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_sub_long_ulonglong(long a, unsigned long long b, long* ret) { return check_sub_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned long safe_sub_ulong_int(unsigned long a, int b) { return safe_sub_uint32_int32(a, b); }
static inline unsigned long safe_sub_ulong_uint(unsigned long a, unsigned int b) { return safe_sub_uint32_uint32(a, b); }
static inline unsigned long safe_sub_ulong_long(unsigned long a, long b) { return safe_sub_uint32_int32(a, b); }
static inline unsigned long safe_sub_ulong_ulong(unsigned long a, unsigned long b) { return safe_sub_uint32_uint32(a, b); }
static inline unsigned long safe_sub_ulong_longlong(unsigned long a, long long b) { return safe_sub_uint32_int64(a, b); }
static inline unsigned long safe_sub_ulong_ulonglong(unsigned long a, unsigned long long b) { return safe_sub_uint32_uint64(a, b); }

static inline bool check_sub_ulong_int(unsigned long a, int b, unsigned long* ret) { return check_sub_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_sub_ulong_uint(unsigned long a, unsigned int b, unsigned long* ret) { return check_sub_uint32_uint32(a, b, (uint32_t*)ret); }
static inline bool check_sub_ulong_long(unsigned long a, long b, unsigned long* ret) { return check_sub_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_sub_ulong_ulong(unsigned long a, unsigned long b, unsigned long* ret) { return check_sub_uint32_uint32(a, b, (uint32_t*)ret); }
static inline bool check_sub_ulong_longlong(unsigned long a, long long b, unsigned long* ret) { return check_sub_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_sub_ulong_ulonglong(unsigned long a, unsigned long long b, unsigned long* ret) { return check_sub_uint32_uint64(a, b, (uint32_t*)ret); }

static inline long long safe_sub_longlong_long(long long a, long b) { return safe_sub_int64_int32(a, b); }
static inline long long safe_sub_longlong_ulong(long long a, unsigned long b) { return safe_sub_int64_uint32(a, b); }

static inline bool check_sub_longlong_long(long long a, long b, long long* ret) { return check_sub_int64_int32(a, b, ret); }
static inline bool check_sub_longlong_ulong(long long a, unsigned long b, long long* ret) { return check_sub_int64_uint32(a, b, ret); }

static inline unsigned long long safe_sub_ulonglong_long(unsigned long long a, long b) { return safe_sub_uint64_int32(a, b); }
static inline unsigned long long safe_sub_ulonglong_ulong(unsigned long long a, unsigned long b) { return safe_sub_uint64_uint32(a, b); }

static inline bool check_sub_ulonglong_long(unsigned long long a, long b, unsigned long long* ret) { return check_sub_uint64_int32(a, b, ret); }
static inline bool check_sub_ulonglong_ulong(unsigned long long a, unsigned long b, unsigned long long* ret) { return check_sub_uint64_uint64(a, b, ret); }
#endif

static inline long long safe_sub_longlong_int(long long a, int b) { return safe_sub_int64_int32(a, b); }
static inline long long safe_sub_longlong_uint(long long a, unsigned int b) { return safe_sub_int64_uint32(a, b); }
static inline long long safe_sub_longlong_longlong(long long a, long long b) { return safe_sub_int64_int64(a, b); }
static inline long long safe_sub_longlong_ulonglong(long long a, unsigned long long b) { return safe_sub_int64_uint64(a, b); }

static inline bool check_sub_longlong_int(long long a, int b, long long* ret) { return check_sub_int64_int32(a, b, (int64_t*)ret); }
static inline bool check_sub_longlong_uint(long long a, unsigned int b, long long* ret) { return check_sub_int64_uint32(a, b, (int64_t*)ret); }
static inline bool check_sub_longlong_longlong(long long a, long long b, long long* ret) { return check_sub_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_sub_longlong_ulonglong(long long a, unsigned long long b, long long* ret) { return check_sub_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long long safe_sub_ulonglong_int(unsigned long long a, int b) { return safe_sub_uint64_int32(a, b); }
static inline unsigned long long safe_sub_ulonglong_uint(unsigned long long a, unsigned int b) { return safe_sub_uint64_uint32(a, b); }
static inline unsigned long long safe_sub_ulonglong_longlong(unsigned long long a, long long b) { return safe_sub_uint64_int64(a, b); }
static inline unsigned long long safe_sub_ulonglong_ulonglong(unsigned long long a, unsigned long long b) { return safe_sub_uint64_uint64(a, b); }

static inline bool check_sub_ulonglong_int(unsigned long long a, int b, unsigned long long* ret) { return check_sub_uint64_int32(a, b, (uint64_t*)ret); }
static inline bool check_sub_ulonglong_uint(unsigned long long a, unsigned int b, unsigned long long* ret) { return check_sub_uint64_uint32(a, b, (uint64_t*)ret); }
static inline bool check_sub_ulonglong_longlong(unsigned long long a, long long b, unsigned long long* ret) { return check_sub_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_sub_ulonglong_ulonglong(unsigned long long a, unsigned long long b, unsigned long long* ret) { return check_sub_uint64_uint64(a, b, (uint64_t*)ret); }

// Division
static inline int safe_div_int_int(int a, int b) { return safe_div_int32_int32(a, b); }
static inline int safe_div_int_uint(int a, unsigned int b) { return safe_div_int32_uint32(a, b); }
static inline int safe_div_int_longlong(int a, long long b) { return safe_div_int32_int64(a, b); }
static inline int safe_div_int_ulonglong(int a, unsigned long long b) { return safe_div_int32_uint64(a, b); }

static inline bool check_div_int_int(int a, int b, int* ret) { return check_div_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_div_int_uint(int a, unsigned int b, int* ret) { return check_div_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_div_int_longlong(int a, long long b, int* ret) { return check_div_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_div_int_ulonglong(int a, unsigned long long b, int* ret) { return check_div_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned int safe_div_uint_int(unsigned int a, int b) { return safe_div_uint32_int32(a, b); }
static inline unsigned int safe_div_uint_uint(unsigned int a, unsigned int b) { return safe_div_uint32_uint32(a, b); }
static inline unsigned int safe_div_uint_longlong(unsigned int a, long long b) { return safe_div_uint32_int64(a, b); }
static inline unsigned int safe_div_uint_ulonglong(unsigned int a, unsigned long long b) { return safe_div_uint32_uint64(a, b); }

static inline bool check_div_uint_int(unsigned int a, int b, unsigned int* ret) { return check_div_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_div_uint_uint(unsigned int a, unsigned int b, unsigned int* ret) { return check_div_uint32_uint32(a, b, (uint32_t*)ret); }
static inline bool check_div_uint_longlong(unsigned int a, long long b, unsigned int* ret) { return check_div_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_div_uint_ulonglong(unsigned int a, unsigned long long b, unsigned int* ret) { return check_div_uint32_uint64(a, b, (uint32_t*)ret); }

#if SAFE_MATH_LONG == 64
static inline int safe_div_int_long(int a, long b) { return safe_div_int32_int64(a, b); }
static inline int safe_div_int_ulong(int a, unsigned long b) { return safe_div_int32_uint64(a, b); }

static inline bool check_div_int_long(int a, long b, int* ret) { return check_div_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_div_int_ulong(int a, unsigned long b, int* ret) { return check_div_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned int safe_div_uint_long(unsigned int a, long b) { return safe_div_uint32_int64(a, b); }
static inline unsigned int safe_div_uint_ulong(unsigned int a, unsigned long b) { return safe_div_uint32_uint64(a, b); }

static inline bool check_div_uint_long(unsigned int a, long b, unsigned int* ret) { return check_div_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_div_uint_ulong(unsigned int a, unsigned long b, unsigned int* ret) { return check_div_uint32_uint64(a, b, (uint32_t*)ret); }

static inline long safe_div_long_int(long a, int b) { return safe_div_int64_int32(a, b); }
static inline long safe_div_long_uint(long a, unsigned int b) { return safe_div_int64_uint32(a, b); }
static inline long safe_div_long_long(long a, long b) { return safe_div_int64_int64(a, b); }
static inline long safe_div_long_ulong(long a, unsigned long b) { return safe_div_int64_uint64(a, b); }
static inline long safe_div_long_longlong(long a, long long b) { return safe_div_int64_int64(a, b); }
static inline long safe_div_long_ulonglong(long a, unsigned long long b) { return safe_div_int64_uint64(a, b); }

static inline bool check_div_long_int(long a, int b, long* ret) { return check_div_int64_int32(a, b, (int64_t*)ret); }
static inline bool check_div_long_uint(long a, unsigned int b, long* ret) { return check_div_int64_uint32(a, b, (int64_t*)ret); }
static inline bool check_div_long_long(long a, long b, long* ret) { return check_div_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_div_long_ulong(long a, unsigned long b, long* ret) { return check_div_int64_uint64(a, b, (int64_t*)ret); }
static inline bool check_div_long_longlong(long a, long long b, long* ret) { return check_div_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_div_long_ulonglong(long a, unsigned long long b, long* ret) { return check_div_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long safe_div_ulong_int(unsigned long a, int b) { return safe_div_uint64_int32(a, b); }
static inline unsigned long safe_div_ulong_uint(unsigned long a, unsigned int b) { return safe_div_uint64_uint32(a, b); }
static inline unsigned long safe_div_ulong_long(unsigned long a, long b) { return safe_div_uint64_int64(a, b); }
static inline unsigned long safe_div_ulong_ulong(unsigned long a, unsigned long b) { return safe_div_uint64_uint64(a, b); }
static inline unsigned long safe_div_ulong_longlong(unsigned long a, long long b) { return safe_div_uint64_int64(a, b); }
static inline unsigned long safe_div_ulong_ulonglong(unsigned long a, unsigned long long b) { return safe_div_uint64_uint64(a, b); }

static inline bool check_div_ulong_int(unsigned long a, int b, unsigned long* ret) { return check_div_uint64_int32(a, b, (uint64_t*)ret); }
static inline bool check_div_ulong_uint(unsigned long a, unsigned int b, unsigned long* ret) { return check_div_uint64_uint32(a, b, (uint64_t*)ret); }
static inline bool check_div_ulong_long(unsigned long a, long b, unsigned long* ret) { return check_div_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_div_ulong_ulong(unsigned long a, unsigned long b, unsigned long* ret) { return check_div_uint64_uint64(a, b, (uint64_t*)ret); }
static inline bool check_div_ulong_longlong(unsigned long a, long long b, unsigned long* ret) { return check_div_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_div_ulong_ulonglong(unsigned long a, unsigned long long b, unsigned long* ret) { return check_div_uint64_uint64(a, b, (uint64_t*)ret); }

static inline long long safe_div_longlong_long(long long a, long b) { return safe_div_int64_int64(a, b); }
static inline long long safe_div_longlong_ulong(long long a, unsigned long b) { return safe_div_int64_uint64(a, b); }

static inline bool check_div_longlong_long(long long a, long b, long long* ret) { return check_div_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_div_longlong_ulong(long long a, unsigned long b, long long* ret) { return check_div_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long long safe_div_ulonglong_long(unsigned long long a, long b) { return safe_div_uint64_int64(a, b); }
static inline unsigned long long safe_div_ulonglong_ulong(unsigned long long a, unsigned long b) { return safe_div_uint64_uint64(a, b); }

static inline bool check_div_ulonglong_long(unsigned long long a, long b, unsigned long long* ret) { return check_div_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_div_ulonglong_ulong(unsigned long long a, unsigned long b, unsigned long long* ret) { return check_div_uint64_uint64(a, b, (uint64_t*)ret); }
#else
static inline int safe_div_int_long(int a, long b) { return safe_div_int32_int32(a, b); }
static inline int safe_div_int_ulong(int a, unsigned long b) { return safe_div_int32_uint32(a, b); }

static inline bool check_div_int_long(int a, long b, int* ret) { return check_div_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_div_int_ulong(int a, unsigned long b, int* ret) { return check_div_int32_uint32(a, b, (int32_t*)ret); }

static inline unsigned int safe_div_uint_long(unsigned int a, long b) { return safe_div_uint32_int32(a, b); }
static inline unsigned int safe_div_uint_ulong(unsigned int a, unsigned long b) { return safe_div_uint32_uint32(a, b); }

static inline bool check_div_uint_long(unsigned int a, long b, unsigned int* ret) { return check_div_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_div_uint_ulong(unsigned int a, unsigned long b, unsigned int* ret) { return check_div_uint32_uint32(a, b, (uint32_t*)ret); }

static inline long safe_div_long_int(long a, int b) { return safe_div_int32_int32(a, b); }
static inline long safe_div_long_uint(long a, unsigned int b) { return safe_div_int32_uint32(a, b); }
static inline long safe_div_long_long(long a, long b) { return safe_div_int32_int32(a, b); }
static inline long safe_div_long_ulong(long a, unsigned long b) { return safe_div_int32_uint32(a, b); }
static inline long safe_div_long_longlong(long a, long long b) { return safe_div_int32_int64(a, b); }
static inline long safe_div_long_ulonglong(long a, unsigned long long b) { return safe_div_int32_uint64(a, b); }

static inline bool check_div_long_int(long a, int b, long* ret) { return check_div_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_div_long_uint(long a, unsigned int b, long* ret) { return check_div_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_div_long_long(long a, long b, long* ret) { return check_div_int32_int32(a, b, (int32_t*)ret); }
static inline bool check_div_long_ulong(long a, unsigned long b, long* ret) { return check_div_int32_uint32(a, b, (int32_t*)ret); }
static inline bool check_div_long_longlong(long a, long long b, long* ret) { return check_div_int32_int64(a, b, (int32_t*)ret); }
static inline bool check_div_long_ulonglong(long a, unsigned long long b, long* ret) { return check_div_int32_uint64(a, b, (int32_t*)ret); }

static inline unsigned long safe_div_ulong_int(unsigned long a, int b) { return safe_div_uint32_int32(a, b); }
static inline unsigned long safe_div_ulong_uint(unsigned long a, unsigned int b) { return safe_div_uint32_uint32(a, b); }
static inline unsigned long safe_div_ulong_long(unsigned long a, long b) { return safe_div_uint32_int32(a, b); }
static inline unsigned long safe_div_ulong_ulong(unsigned long a, unsigned long b) { return safe_div_uint32_uint32(a, b); }
static inline unsigned long safe_div_ulong_longlong(unsigned long a, long long b) { return safe_div_uint32_int64(a, b); }
static inline unsigned long safe_div_ulong_ulonglong(unsigned long a, unsigned long long b) { return safe_div_uint32_uint64(a, b); }

static inline bool check_div_ulong_int(unsigned long a, int b, unsigned long* ret) { return check_div_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_div_ulong_uint(unsigned long a, unsigned int b, unsigned long* ret) { return check_div_uint32_uint32(a, b, (uint32_t*)ret); }
static inline bool check_div_ulong_long(unsigned long a, long b, unsigned long* ret) { return check_div_uint32_int32(a, b, (uint32_t*)ret); }
static inline bool check_div_ulong_ulong(unsigned long a, unsigned long b, unsigned long* ret) { return check_div_uint32_uint64(a, b, (uint32_t*)ret); }
static inline bool check_div_ulong_longlong(unsigned long a, long long b, unsigned long* ret) { return check_div_uint32_int64(a, b, (uint32_t*)ret); }
static inline bool check_div_ulong_ulonglong(unsigned long a, unsigned long long b, unsigned long* ret) { return check_div_uint32_uint64(a, b, (uint32_t*)ret); }

static inline long long safe_div_longlong_long(long long a, long b) { return safe_div_int64_int32(a, b); }
static inline long long safe_div_longlong_ulong(long long a, unsigned long b) { return safe_div_int64_uint32(a, b); }

static inline bool check_div_longlong_long(long long a, long b, long long* ret) { return check_div_int64_int32(a, b, ret); }
static inline bool check_div_longlong_ulong(long long a, unsigned long b, long long* ret) { return check_div_int64_uint32(a, b, ret); }

static inline unsigned long long safe_div_ulonglong_long(unsigned long long a, long b) { return safe_div_uint64_int32(a, b); }
static inline unsigned long long safe_div_ulonglong_ulong(unsigned long long a, unsigned long b) { return safe_div_uint64_uint32(a, b); }

static inline bool check_div_ulonglong_long(unsigned long long a, long b, unsigned long long* ret) { return check_div_uint64_int32(a, b, ret); }
static inline bool check_div_ulonglong_ulong(unsigned long long a, unsigned long b, unsigned long long* ret) { return check_div_uint64_uint32(a, b, ret); }
#endif

static inline long long safe_div_longlong_int(long long a, int b) { return safe_div_int64_int32(a, b); }
static inline long long safe_div_longlong_uint(long long a, unsigned int b) { return safe_div_int64_uint32(a, b); }
static inline long long safe_div_longlong_longlong(long long a, long long b) { return safe_div_int64_int64(a, b); }
static inline long long safe_div_longlong_ulonglong(long long a, unsigned long long b) { return safe_div_int64_uint64(a, b); }

static inline bool check_div_longlong_int(long long a, int b, long long* ret) { return check_div_int64_int32(a, b, (int64_t*)ret); }
static inline bool check_div_longlong_uint(long long a, unsigned int b, long long* ret) { return check_div_int64_uint32(a, b, (int64_t*)ret); }
static inline bool check_div_longlong_longlong(long long a, long long b, long long* ret) { return check_div_int64_int64(a, b, (int64_t*)ret); }
static inline bool check_div_longlong_ulonglong(long long a, unsigned long long b, long long* ret) { return check_div_int64_uint64(a, b, (int64_t*)ret); }

static inline unsigned long long safe_div_ulonglong_int(unsigned long long a, int b) { return safe_div_uint64_int32(a, b); }
static inline unsigned long long safe_div_ulonglong_uint(unsigned long long a, unsigned int b) { return safe_div_uint64_uint32(a, b); }
static inline unsigned long long safe_div_ulonglong_longlong(unsigned long long a, long long b) { return safe_div_uint64_int64(a, b); }
static inline unsigned long long safe_div_ulonglong_ulonglong(unsigned long long a, unsigned long long b) { return safe_div_uint64_uint64(a, b); }

static inline bool check_div_ulonglong_int(unsigned long long a, int b, unsigned long long* ret) { return check_div_uint64_int32(a, b, (uint64_t*)ret); }
static inline bool check_div_ulonglong_uint(unsigned long long a, unsigned int b, unsigned long long* ret) { return check_div_uint64_uint32(a, b, (uint64_t*)ret); }
static inline bool check_div_ulonglong_longlong(unsigned long long a, long long b, unsigned long long* ret) { return check_div_uint64_int64(a, b, (uint64_t*)ret); }
static inline bool check_div_ulonglong_ulonglong(unsigned long long a, unsigned long long b, unsigned long long* ret) { return check_div_uint64_uint64(a, b, (uint64_t*)ret); }

#ifdef __cplusplus
}
#endif

#endif
