/*
 * Swiss Table Type Specializations for pandas
 *
 * Defines hash functions and equality for all pandas numeric types.
 *
 * Hash Function Summary:
 * ----------------------
 * 1. swiss_mix64(x) - Fast single-round Fibonacci hash for integers
 *    Used by: int8, int16, int32, int64, uint8, uint16, uint32, uint64
 *    Cost: 1 multiply + 1 xor-shift
 *    Properties: Fast, adequate distribution for typical pandas workloads
 *
 * 2. swiss_splitmix64(x) - High-quality finalizer for floats
 *    Used by: float32, float64, complex64, complex128
 *    Cost: 2 multiplies + 3 xor-shifts
 *    Properties: Excellent avalanche, handles pathological float bit patterns
 *
 * Special Value Handling:
 * -----------------------
 * - NaN: All NaN values hash to the same value and compare equal
 * - Zero: +0.0 and -0.0 hash to the same value (they compare equal via ==)
 * - Complex NaN: If either component is NaN, the complex value is treated as NaN
 */

#ifndef PANDAS_SWISSTABLE_TYPES_H
#define PANDAS_SWISSTABLE_TYPES_H

#include "swisstable_generic.h"
#include <stdint.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Hash Function 1: swiss_mix64 - Fast Fibonacci Hash
 *
 * Single multiply + xor-shift mixing for integer types.
 * The golden ratio constant (2^64 / phi) spreads bits upward,
 * then xor-shift mixes high bits back into low bits.
 *
 * This provides adequate H1/H2 independence for Swiss Tables:
 * - H1 (index) uses low bits of hash
 * - H2 (control byte) uses top 7 bits of hash
 *
 * Used by: All integer types (int8 through int64, uint8 through uint64)
 * ============================================================================ */
static inline uint64_t swiss_mix64(uint64_t x) {
    x *= 0x9e3779b97f4a7c15ULL;  /* 2^64 / phi (golden ratio constant) */
    x ^= x >> 33;                 /* Mix high bits back into low bits */
    return x;
}

/* Integer hash functions - all use swiss_mix64 */
static inline uint64_t swiss_hash_int64(int64_t key) {
    return swiss_mix64((uint64_t)key);
}

static inline uint64_t swiss_hash_uint64(uint64_t key) {
    return swiss_mix64(key);
}

static inline uint64_t swiss_hash_int32(int32_t key) {
    return swiss_mix64((uint64_t)(uint32_t)key);
}

static inline uint64_t swiss_hash_uint32(uint32_t key) {
    return swiss_mix64((uint64_t)key);
}

static inline uint64_t swiss_hash_int16(int16_t key) {
    return swiss_mix64((uint64_t)(uint16_t)key);
}

static inline uint64_t swiss_hash_uint16(uint16_t key) {
    return swiss_mix64((uint64_t)key);
}

static inline uint64_t swiss_hash_int8(int8_t key) {
    return swiss_mix64((uint64_t)(uint8_t)key);
}

static inline uint64_t swiss_hash_uint8(uint8_t key) {
    return swiss_mix64((uint64_t)key);
}

/* ============================================================================
 * Hash Function 2: swiss_splitmix64 - High-Quality Finalizer
 *
 * Full splitmix64 finalizer with excellent avalanche properties.
 * More expensive than swiss_mix64 but necessary for floating point types
 * where bit patterns can be pathological (e.g., sequential floats have
 * very similar mantissa bits).
 *
 * Used by: float32, float64, complex64, complex128
 * ============================================================================ */
static inline uint64_t swiss_splitmix64(uint64_t x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

/* ============================================================================
 * Float Helpers (bit reinterpretation)
 * ============================================================================ */

static inline uint64_t swiss_asuint64(double key) {
    uint64_t val;
    memcpy(&val, &key, sizeof(double));
    return val;
}

static inline uint32_t swiss_asuint32(float key) {
    uint32_t val;
    memcpy(&val, &key, sizeof(float));
    return val;
}

/* ============================================================================
 * Float Hash Functions (with NaN and zero handling)
 * Using splitmix64 for better avalanche properties.
 * ============================================================================ */

/* Use splitmix64 output for special values to get good distribution */
#define SWISS_ZERO_HASH swiss_splitmix64(0)
#define SWISS_NAN_HASH swiss_splitmix64(0x7ff8000000000000ULL)

static inline uint64_t swiss_hash_float64(double val) {
    /* +0.0 and -0.0 should have the same hash */
    if (val == 0.0) {
        return SWISS_ZERO_HASH;
    }
    /* All NaN values should have the same hash */
    if (val != val) {  /* NaN check: NaN != NaN */
        return SWISS_NAN_HASH;
    }
    uint64_t as_int = swiss_asuint64(val);
    return swiss_splitmix64(as_int);
}

static inline uint64_t swiss_hash_float32(float val) {
    if (val == 0.0f) {
        return SWISS_ZERO_HASH;
    }
    if (val != val) {
        return SWISS_NAN_HASH;
    }
    uint32_t as_int = swiss_asuint32(val);
    /* Expand to 64-bit then use splitmix64 */
    uint64_t x = (uint64_t)as_int * 0x9e3779b97f4a7c15ULL;
    return swiss_splitmix64(x);
}

/* ============================================================================
 * Complex Number Types and Hash Functions
 * ============================================================================ */

typedef struct {
    float real;
    float imag;
} swiss_complex64_t;

typedef struct {
    double real;
    double imag;
} swiss_complex128_t;

/* For complex, combine the hashes with rotation to avoid symmetry issues */
static inline uint64_t swiss_hash_complex64(swiss_complex64_t val) {
    uint64_t h1 = swiss_hash_float32(val.real);
    uint64_t h2 = swiss_hash_float32(val.imag);
    /* Rotate h2 to avoid (a,b) and (b,a) colliding */
    return h1 ^ ((h2 << 32) | (h2 >> 32));
}

static inline uint64_t swiss_hash_complex128(swiss_complex128_t val) {
    uint64_t h1 = swiss_hash_float64(val.real);
    uint64_t h2 = swiss_hash_float64(val.imag);
    return h1 ^ ((h2 << 32) | (h2 >> 32));
}

/* ============================================================================
 * Integer Equality Functions (simple comparison)
 * ============================================================================ */

static inline bool swiss_equal_int64(int64_t a, int64_t b) { return a == b; }
static inline bool swiss_equal_uint64(uint64_t a, uint64_t b) { return a == b; }
static inline bool swiss_equal_int32(int32_t a, int32_t b) { return a == b; }
static inline bool swiss_equal_uint32(uint32_t a, uint32_t b) { return a == b; }
static inline bool swiss_equal_int16(int16_t a, int16_t b) { return a == b; }
static inline bool swiss_equal_uint16(uint16_t a, uint16_t b) { return a == b; }
static inline bool swiss_equal_int8(int8_t a, int8_t b) { return a == b; }
static inline bool swiss_equal_uint8(uint8_t a, uint8_t b) { return a == b; }

/* ============================================================================
 * Float Equality Functions (NaN == NaN semantics)
 * ============================================================================ */

/* NaN values are considered equal (for groupby, unique, etc.) */
static inline bool swiss_equal_float64(double a, double b) {
    return (a == b) || (a != a && b != b);  /* NaN == NaN */
}

static inline bool swiss_equal_float32(float a, float b) {
    return (a == b) || (a != a && b != b);
}

static inline bool swiss_equal_complex64(swiss_complex64_t a, swiss_complex64_t b) {
    return swiss_equal_float32(a.real, b.real) && swiss_equal_float32(a.imag, b.imag);
}

static inline bool swiss_equal_complex128(swiss_complex128_t a, swiss_complex128_t b) {
    return swiss_equal_float64(a.real, b.real) && swiss_equal_float64(a.imag, b.imag);
}

/* ============================================================================
 * Swiss Table Type Instantiations
 *
 * All tables map keys to size_t values (for factorize/indexing operations)
 * ============================================================================ */

/* Integer types */
SWISSTABLE_INIT(int64, int64_t, size_t, swiss_hash_int64, swiss_equal_int64)
SWISSTABLE_INIT(uint64, uint64_t, size_t, swiss_hash_uint64, swiss_equal_uint64)
SWISSTABLE_INIT(int32, int32_t, size_t, swiss_hash_int32, swiss_equal_int32)
SWISSTABLE_INIT(uint32, uint32_t, size_t, swiss_hash_uint32, swiss_equal_uint32)
SWISSTABLE_INIT(int16, int16_t, size_t, swiss_hash_int16, swiss_equal_int16)
SWISSTABLE_INIT(uint16, uint16_t, size_t, swiss_hash_uint16, swiss_equal_uint16)
SWISSTABLE_INIT(int8, int8_t, size_t, swiss_hash_int8, swiss_equal_int8)
SWISSTABLE_INIT(uint8, uint8_t, size_t, swiss_hash_uint8, swiss_equal_uint8)

/* Floating point types */
SWISSTABLE_INIT(float64, double, size_t, swiss_hash_float64, swiss_equal_float64)
SWISSTABLE_INIT(float32, float, size_t, swiss_hash_float32, swiss_equal_float32)

/* Complex types */
SWISSTABLE_INIT(complex64, swiss_complex64_t, size_t, swiss_hash_complex64, swiss_equal_complex64)
SWISSTABLE_INIT(complex128, swiss_complex128_t, size_t, swiss_hash_complex128, swiss_equal_complex128)

#endif  /* PANDAS_SWISSTABLE_TYPES_H */
