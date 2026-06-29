/*
 * Swiss Table Hash Map Implementation for pandas
 *
 * Based on Google's Abseil Swiss Tables design:
 * https://abseil.io/about/design/swisstables
 *
 * Algorithm: Open addressing with SIMD-accelerated quadratic probing
 *
 * Key features:
 * - Control bytes store 7-bit hash fragments (h2)
 * - SIMD compares 16 control bytes in parallel (SSE2 on x86, NEON on ARM)
 * - Load factor: 81.25% (vs khash's 77%)
 * - Cache-efficient metadata scanning
 *
 * License: BSD-3-Clause (pandas)
 */

#ifndef PANDAS_SWISSTABLE_CORE_H
#define PANDAS_SWISSTABLE_CORE_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>  // for malloc, free (used by SWISSTABLE_MALLOC/FREE defaults)
#include <string.h>  // for memcpy, memset
#include <assert.h>  // for assert in countr_zero

// MSVC intrinsics for _BitScanForward
#if defined(_MSC_VER)
#include <intrin.h>
#endif

// SIMD support - detect architecture
#if defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
    #define SWISSTABLE_SIMD_SSE2
    #include <emmintrin.h>  // SSE2
#elif defined(__aarch64__)
    // AArch64 only - vqtbl1q_u8 not available on 32-bit ARM NEON
    #define SWISSTABLE_SIMD_NEON
    #include <arm_neon.h>
#else
    // Fallback to portable implementation (includes 32-bit ARM)
    #define SWISSTABLE_SIMD_PORTABLE
#endif

// Control byte values
#define CTRL_EMPTY   ((int8_t)0x80)  // -128: Empty slot
#define CTRL_DELETED ((int8_t)0xFE)  // -2: Deleted slot (tombstone)
// 0x00-0x7F: Occupied slot, stores swiss_h2(hash)

// Extract H2 (top 7 bits) from hash
// Returns int8_t in range 0x00-0x7F (always non-negative)
// Using top bits for H2 while H1 (index) uses low bits reduces correlation
// between the two, improving probe sequence independence.
// All hash functions must produce well-mixed 64-bit output for this to work.
static inline int8_t swiss_h2(uint64_t hash) {
    return (int8_t)(hash >> 57);  // Top 7 bits
}

// Memory allocation macros (default to malloc/free; can be overridden by pandas).
#ifndef SWISSTABLE_MALLOC
#define SWISSTABLE_MALLOC(size) malloc(size)
#endif

#ifndef SWISSTABLE_FREE
#define SWISSTABLE_FREE(ptr) free(ptr)
#endif

// Group: 16 control bytes for SIMD operations
typedef struct {
#ifdef SWISSTABLE_SIMD_SSE2
    __m128i ctrl;
#elif defined(SWISSTABLE_SIMD_NEON)
    int8x16_t ctrl;
#else
    int8_t ctrl[16];
#endif
} Group;

// Load 16 control bytes into SIMD register
static inline Group group_load(const int8_t *ctrl) {
    Group g;
#ifdef SWISSTABLE_SIMD_SSE2
    g.ctrl = _mm_loadu_si128((const __m128i*)ctrl);
#elif defined(SWISSTABLE_SIMD_NEON)
    g.ctrl = vld1q_s8(ctrl);
#else
    memcpy(g.ctrl, ctrl, 16);
#endif
    return g;
}

// Helper function to convert NEON comparison result to bitmask
// Optimized for AArch64 using 64-bit multiply pack approach
#ifdef SWISSTABLE_SIMD_NEON
static inline uint16_t neon_movemask(uint8x16_t v) {
    // v is 0xFF/0x00 per byte. Extract MSB -> 0x01/0x00 per byte.
    uint8x16_t msb = vshrq_n_u8(v, 7);

    // Two 64-bit lanes, each lane holds 8 bytes (little-endian on AArch64).
    uint64x2_t lanes = vreinterpretq_u64_u8(msb);
    uint64_t lo = vgetq_lane_u64(lanes, 0);
    uint64_t hi = vgetq_lane_u64(lanes, 1);

    // Pack 8 bytes into an 8-bit mask using magic multiply.
    // Magic constant 0x0102040810204080 packs bits correctly for little-endian.
    // Each byte position i in the 64-bit word gets multiplied and shifted
    // such that bit i ends up in the high byte.
    const uint64_t magic = 0x0102040810204080ULL;
    uint8_t mask_lo = (uint8_t)((lo * magic) >> 56);
    uint8_t mask_hi = (uint8_t)((hi * magic) >> 56);

    return (uint16_t)mask_lo | ((uint16_t)mask_hi << 8);
}
#endif

// Match control bytes against h2, returns bitmask
static inline uint16_t group_match(Group g, int8_t h2) {
#ifdef SWISSTABLE_SIMD_SSE2
    __m128i target = _mm_set1_epi8(h2);
    __m128i matches = _mm_cmpeq_epi8(g.ctrl, target);
    return (uint16_t)_mm_movemask_epi8(matches);
#elif defined(SWISSTABLE_SIMD_NEON)
    int8x16_t target = vdupq_n_s8(h2);
    uint8x16_t matches = vceqq_s8(g.ctrl, target);
    // Fast-reject: if no 0xFF bytes, skip expensive movemask
    if (vmaxvq_u8(matches) == 0) {
        return 0;
    }
    return neon_movemask(matches);
#else
    uint16_t mask = 0;
    for (int i = 0; i < 16; i++) {
        if (g.ctrl[i] == h2) {
            mask |= (1 << i);
        }
    }
    return mask;
#endif
}

// Match EMPTY control bytes
static inline uint16_t group_match_empty(Group g) {
    return group_match(g, CTRL_EMPTY);
}

// Match EMPTY or DELETED control bytes
static inline uint16_t group_match_empty_or_deleted(Group g) {
#ifdef SWISSTABLE_SIMD_SSE2
    __m128i ctrl_bytes = g.ctrl;
    __m128i msb_mask = _mm_set1_epi8((char)0x80);
    __m128i result = _mm_and_si128(ctrl_bytes, msb_mask);
    __m128i cmp = _mm_cmpeq_epi8(result, msb_mask);
    return (uint16_t)_mm_movemask_epi8(cmp);
#elif defined(SWISSTABLE_SIMD_NEON)
    // Empty: 0x80 (-128), Deleted: 0xFE (-2)
    // Both have MSB set (< 0 when signed)
    uint8x16_t msb_set = vcltq_s8(g.ctrl, vdupq_n_s8(0));
    // Fast-reject: if no matches, skip expensive movemask
    if (vmaxvq_u8(msb_set) == 0) {
        return 0;
    }
    return neon_movemask(msb_set);
#else
    uint16_t mask = 0;
    for (int i = 0; i < 16; i++) {
        if (g.ctrl[i] == CTRL_EMPTY || g.ctrl[i] == CTRL_DELETED) {
            mask |= (1 << i);
        }
    }
    return mask;
#endif
}

// Check if control byte indicates occupied slot
static inline bool is_full(int8_t ctrl) {
    return ctrl >= 0;  // 0x00-0x7F
}

// Count trailing zeros (position of first set bit)
// Precondition: x must be non-zero (calling with x == 0 is undefined behavior)
static inline int countr_zero(uint16_t x) {
    assert(x != 0);  // Verify precondition
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctz(x);
#elif defined(_MSC_VER)
    unsigned long index;
    _BitScanForward(&index, x);
    return (int)index;
#else
    int count = 0;
    while ((x & 1) == 0 && count < 16) {
        x >>= 1;
        count++;
    }
    return count;
#endif
}

// Probe sequence for quadratic probing
typedef struct {
    size_t mask;
    size_t offset;
    size_t index;
} ProbeSeq;

static inline void probe_seq_init(ProbeSeq *seq, uint64_t hash, size_t mask) {
    seq->mask = mask;
    seq->offset = hash & mask;
    seq->index = 0;
}

static inline size_t probe_seq_offset(const ProbeSeq *seq) {
    return seq->offset;
}

static inline void probe_seq_next(ProbeSeq *seq) {
    seq->index++;
    seq->offset = (seq->offset + seq->index) & seq->mask;
}

// Normalize capacity to power of 2, minimum 16
static inline size_t normalize_capacity(size_t n) {
    if (n < 16) return 16;

    // Round up to next power of 2
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
#if SIZE_MAX > UINT32_MAX
    n |= n >> 32;
#endif
    n++;

    return n;
}

// Calculate growth_left from capacity (81.25% load factor)
static inline size_t capacity_to_growth(size_t capacity) {
    if (capacity == 0) return 0;
    // 13/16 = 81.25% load factor
    // Multiply first to preserve precision (though capacity is always power-of-2)
    return capacity * 13 / 16;
}

#endif  // PANDAS_SWISSTABLE_CORE_H
