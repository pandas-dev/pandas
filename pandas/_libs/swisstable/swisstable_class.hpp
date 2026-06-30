/*
 * Swiss Table C++17 Template Implementation for pandas
 *
 * Based on Google's Abseil Swiss Tables design:
 * https://abseil.io/about/design/swisstables
 *
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <new>
#include <limits>
#include <type_traits>

// MSVC intrinsics for _BitScanForward
#if defined(_MSC_VER)
#include <intrin.h>
#endif

// SIMD support - detect architecture
#if defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
    #define SWISSTABLE_SIMD_SSE2
    #include <emmintrin.h>  // SSE2
    #define GROUP_WIDTH 16
#elif defined(__aarch64__)
    // AArch64 only - vqtbl1q_u8 not available on 32-bit ARM NEON
    #define SWISSTABLE_SIMD_NEON
    #include <arm_neon.h>
    #include <arm_acle.h> // for crc
    #define GROUP_WIDTH 16
#else
    // Fallback to portable implementation (includes 32-bit ARM)
    #define SWISSTABLE_SIMD_PORTABLE
    #define GROUP_WIDTH 16
#endif

// =============================================================================
// ssize_t definition for Windows compatibility
// =============================================================================
#if defined(_MSC_VER)
typedef intptr_t ssize_t;
#else
#include <sys/types.h>
#endif

#if defined(_MSC_VER)
#define SWISS_FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define SWISS_FORCE_INLINE inline __attribute__((always_inline))
#else
#define SWISS_FORCE_INLINE inline
#endif

namespace pandas::swisstable {

using ctrl_t = int8_t;

// Control byte values
#define CTRL_EMPTY   ((ctrl_t)0x80)  // -128: Empty slot
#define CTRL_DELETED ((ctrl_t)0xFE)  // -2: Deleted slot (tombstone)
// 0x00-0x7F: Occupied slot, stores swiss_h2(hash)

// Extract H2 (top 7 bits) from hash
// Returns ctrl_t in range 0x00-0x7F (always non-negative)
// Using top bits for H2 while H1 (index) uses low bits reduces correlation
// between the two, improving probe sequence independence.
// All hash functions must produce well-mixed 64-bit output for this to work.
inline ctrl_t swiss_h2(uint64_t hash) noexcept {
    return (ctrl_t)(hash >> 57);  // Top 7 bits, masked to ensure non-negative
}

// Memory allocation macros (default to malloc/free; can be overridden by pandas).
#ifndef SWISSTABLE_MALLOC
#define SWISSTABLE_MALLOC(size) malloc(size)
#endif

#ifndef SWISSTABLE_FREE
#define SWISSTABLE_FREE(ptr) free(ptr)
#endif

// =============================================================================
// Probe Sequence (POD struct with forced inline methods)
// =============================================================================
struct ProbeSeq {
    size_t mask;
    size_t offset;
    size_t index;

    ProbeSeq(uint64_t hash, size_t m) noexcept
    {
        mask = m;
        offset = hash & m;
        index = 0;
    }

    inline void next() noexcept
    {
        index += GROUP_WIDTH;
        offset = (offset + index) & mask;
    }
};

// =============================================================================
// SIMD Group Operations (fully inlined)
// =============================================================================
class Group {
public:
    inline static Group load(const ctrl_t *ctrl) noexcept
    {
        Group g;
#if defined(SWISSTABLE_SIMD_NEON)
        g.ctrl_ = vld1q_s8(ctrl);
#elif defined(SWISSTABLE_SIMD_SSE2)
        g.ctrl_ = _mm_loadu_si128(reinterpret_cast<const __m128i *>(ctrl));
#else
        std::memcpy(g.ctrl_, ctrl, GROUP_WIDTH);
#endif
        return g;
    }

    inline uint16_t match(ctrl_t h2) const noexcept
    {
#if defined(SWISSTABLE_SIMD_NEON)
        int8x16_t target = vdupq_n_s8(h2);
        uint8x16_t matches = vceqq_s8(ctrl_, target);
        return neon_movemask(matches);
#elif defined(SWISSTABLE_SIMD_SSE2)
        __m128i target = _mm_set1_epi8(h2);
        __m128i matches = _mm_cmpeq_epi8(ctrl_, target);
        return (uint16_t)_mm_movemask_epi8(matches);
#else
        uint16_t mask = 0;
        for (int i = 0; i < GROUP_WIDTH; i++) {
            if (ctrl_[i] == h2)
                mask |= (1 << i);
        }
        return mask;
#endif
    }

    inline uint16_t match_empty() const noexcept
    {
        return match(CTRL_EMPTY);
    }

    inline uint8_t match_any_empty() const noexcept
    {
#if defined(SWISSTABLE_SIMD_NEON)
        return vminvq_s8(ctrl_) == CTRL_EMPTY;
#else
        return match(CTRL_EMPTY) != 0;
#endif
    }

private:
#if defined(SWISSTABLE_SIMD_NEON)
    static inline uint16_t neon_movemask(uint8x16_t v) noexcept
    {
        const uint8x8_t weights = {
            1,2,4,8,16,32,64,128
        };

        uint8x8_t bits0 = vand_u8(vget_low_u8(v), weights);
        uint8x8_t bits1 = vand_u8(vget_high_u8(v), weights);
        return (uint16_t)vaddv_u8(bits0) | (uint16_t)vaddv_u8(bits1) << 8;
    }
#endif

private:
#if defined(SWISSTABLE_SIMD_NEON)
    int8x16_t ctrl_;
#elif defined(SWISSTABLE_SIMD_SSE2)
    __m128i ctrl_;
#else
    ctrl_t ctrl_[GROUP_WIDTH];
#endif
};

// Check if control byte indicates occupied slot
inline bool is_full(ctrl_t ctrl) noexcept {
    return ctrl >= 0;  // 0x00-0x7F
}

// Count trailing zeros (position of first set bit)
// Precondition: x must be non-zero (calling with x == 0 is undefined behavior)
inline int countr_zero(uint16_t x) noexcept
{
    assert(x != 0);  // Verify precondition
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctz(x);
#elif defined(_MSC_VER)
    unsigned long index;
    _BitScanForward(&index, x);
    return (int)index;
#else
    int count = 0;
    while ((x & 1) == 0 && count < GROUP_WIDTH) {
        x >>= 1;
        count++;
    }
    return count;
#endif
}

// Normalize capacity to power of 2, minimum GROUP_WIDTH
inline size_t normalize_capacity(size_t n) noexcept
{
    if (n < GROUP_WIDTH) return GROUP_WIDTH;

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
inline size_t capacity_to_growth(size_t capacity) noexcept
{
    if (capacity == 0) return 0;
    // 13/16 = 81.25% load factor
    // Multiply first to preserve precision (though capacity is always power-of-2)
    return capacity * 13 / GROUP_WIDTH;
}

// =============================================================================
// Complex Types (custom layout for performance, not std::complex)
// =============================================================================
// Note: CTRL_EMPTY and CTRL_DELETED macros are defined in swisstable_core.h
struct Complex64 {
    float real;
    float imag;
};

struct Complex128 {
    double real;
    double imag;
};

using swiss_complex64_t = Complex64;
using swiss_complex128_t = Complex128;

// =============================================================================
// Hash Strategy - Default (Fibonacci hash for integers)
// =============================================================================
template <typename Key>
struct Hash {
    static inline uint64_t hash(Key key) noexcept
    {
        static_assert(std::is_integral_v<Key>, "Hash only defined for integral types");
        uint64_t x = static_cast<uint64_t>(key);
#ifdef __ARM_FEATURE_CRC32
        return __crc32d(0x9E3779B1, x) * 0xbf58476d1ce4e5b9ULL;
#else
        // Fibonacci hash (golden ratio constant)
        x *= 0x9e3779b97f4a7c15ULL;
        x ^= x >> 33;
        return x;
#endif
    }
};

inline uint64_t swiss_splitmix64(uint64_t x) noexcept
{
#ifdef __ARM_FEATURE_CRC32
    return __crc32d(0x9E3779B1, x) * 0xbf58476d1ce4e5b9ULL;
#else
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
#endif
}

// =============================================================================
// Hash Strategy - Float32 (with NaN handling)
// =============================================================================
template <>
struct Hash<float> {
    static inline uint64_t hash(float val) noexcept
    {
        if (val == 0.0f) {
            // +0.0 and -0.0 should have the same hash
            return swiss_splitmix64(0);
        }
        if (val != val) {                            // NaN check
            return swiss_splitmix64(0x7fc00000ULL);  // Quiet NaN
        }
        uint32_t as_int;
        std::memcpy(&as_int, &val, sizeof(float));
        return swiss_splitmix64(as_int);
    }
};

// =============================================================================
// Hash Strategy - Float64 (with NaN handling)
// =============================================================================
template <>
struct Hash<double> {
    static inline uint64_t hash(double val) noexcept
    {
        if (val == 0.0) {
            return swiss_splitmix64(0);
        }
        if (val != val) {                                    // NaN check
            return swiss_splitmix64(0x7ff8000000000000ULL);  // Quiet NaN
        }
        uint64_t as_int;
        std::memcpy(&as_int, &val, sizeof(double));
        return swiss_splitmix64(as_int);
    }
};

// =============================================================================
// Hash Strategy - Complex64
// =============================================================================
template <>
struct Hash<Complex64> {
    static inline uint64_t hash(Complex64 val) noexcept
    {
        uint64_t h1 = Hash<float>::hash(val.real);
        uint64_t h2 = Hash<float>::hash(val.imag);
        // Rotate h2 to avoid (a,b) and (b,a) colliding
        return h1 ^ ((h2 << 32) | (h2 >> 32));
    }
};

// =============================================================================
// Hash Strategy - Complex128
// =============================================================================
template <>
struct Hash<Complex128> {
    static inline uint64_t hash(Complex128 val) noexcept
    {
        uint64_t h1 = Hash<double>::hash(val.real);
        uint64_t h2 = Hash<double>::hash(val.imag);
        return h1 ^ ((h2 << 32) | (h2 >> 32));
    }
};

// =============================================================================
// Equal Strategy - Default (simple equality)
// =============================================================================
template <typename Key>
struct Equal {
    static inline bool equal(Key a, Key b) noexcept
    {
        return a == b;
    }
};

// =============================================================================
// Equal Strategy - Float32 (NaN == NaN semantics)
// =============================================================================
template <>
struct Equal<float> {
    static inline bool equal(float a, float b) noexcept
    {
        return (a == b) || (a != a && b != b);  // NaN == NaN
    }
};

// =============================================================================
// Equal Strategy - Float64 (NaN == NaN semantics)
// =============================================================================
template <>
struct Equal<double> {
    static inline bool equal(double a, double b) noexcept
    {
        return (a == b) || (a != a && b != b);  // NaN == NaN
    }
};

// =============================================================================
// Equal Strategy - Complex64
// =============================================================================
template <>
struct Equal<Complex64> {
    static inline bool equal(Complex64 a, Complex64 b) noexcept
    {
        return Equal<float>::equal(a.real, b.real) && Equal<float>::equal(a.imag, b.imag);
    }
};

// =============================================================================
// Equal Strategy - Complex128
// =============================================================================
template <>
struct Equal<Complex128> {
    static inline bool equal(Complex128 a, Complex128 b) noexcept
    {
        return Equal<double>::equal(a.real, b.real) && Equal<double>::equal(a.imag, b.imag);
    }
};

// NaN Traits
template <typename Key>
struct NaNTraits {
    static inline Key NaN() noexcept
    {
        if constexpr (std::is_same_v<Key, swiss_complex64_t>) {
            return {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
        } else if constexpr (std::is_same_v<Key, swiss_complex128_t>) {
            return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
        } else {
            return std::numeric_limits<Key>::quiet_NaN();
        }
    }

    static inline bool IsNaN(Key val) noexcept
    {
        if constexpr (std::is_integral_v<Key>) {
            return false;
        } else if constexpr (std::is_same_v<Key, swiss_complex64_t> || std::is_same_v<Key, swiss_complex128_t>) {
            return std::isnan(val.real) || std::isnan(val.imag);
        } else {
            return std::isnan(val);
        }
    }

    static inline bool AreEqual(Key val1, Key val2) noexcept
    {
        return Equal<Key>::equal(val1, val2);
    }
};

// =============================================================================
// SwissTable Class Template
// =============================================================================
template <typename Key, typename Value, typename HashFn = Hash<Key>, typename EqualFn = Equal<Key>>
class SwissTable {
public:
    // =========================================================================
    // Constructor/Destructor
    // =========================================================================

    SwissTable() noexcept
        : capacity_(0),
          mask_(0),
          size_(0),
          growth_left_(0),
          ctrl_(nullptr),
          keys_(nullptr),
          vals_(nullptr),
          alloc_(nullptr)
    {
    }

    ~SwissTable() noexcept
    {
        destroy();
    }

    SwissTable(const SwissTable &) = delete;

    SwissTable(SwissTable &&other) noexcept
    {
        std::swap(capacity_, other.capacity_);
        std::swap(mask_, other.mask_);
        std::swap(size_, other.size_);
        std::swap(growth_left_, other.growth_left_);
        std::swap(ctrl_, other.ctrl_);
        std::swap(keys_, other.keys_);
        std::swap(vals_, other.vals_);
        std::swap(alloc_, other.alloc_);
    }

    // =========================================================================
    // Initialization
    // =========================================================================

    int init_with_capacity(size_t capacity) noexcept
    {
        destroy();
        if (capacity > 0) {
            capacity = normalize_capacity((capacity * GROUP_WIDTH + 12) / 13);  // load factory
            size_t keys_offset, vals_offset;
            size_t alloc_size = calc_alloc_size(capacity, &keys_offset, &vals_offset);

            alloc_ = SWISSTABLE_MALLOC(alloc_size);
            if (!alloc_) {
                return 0;
            }

            // Set up pointers into the single allocation
            ctrl_ = static_cast<ctrl_t *>(alloc_);
            keys_ = static_cast<Key *>(static_cast<void *>(static_cast<char *>(alloc_) + keys_offset));
            vals_ = static_cast<Value *>(static_cast<void *>(static_cast<char *>(alloc_) + vals_offset));

            std::memset(ctrl_, CTRL_EMPTY, capacity + GROUP_WIDTH);
            capacity_ = capacity;
            mask_ = capacity - 1;
            growth_left_ = capacity_to_growth(capacity);
            size_ = 0;
        }
        return 1;
    }

    // =========================================================================
    // Destruction
    // =========================================================================

    void destroy() noexcept
    {
        if (alloc_) {
            SWISSTABLE_FREE(alloc_);
            alloc_ = nullptr;
            ctrl_ = nullptr;
            keys_ = nullptr;
            vals_ = nullptr;
            capacity_ = 0;
            size_ = 0;
            growth_left_ = 0;
            mask_ = 0;
        }
    }

    // =========================================================================
    // Clear all elements but keep capacity
    // =========================================================================

    void clear() noexcept
    {
        if (capacity_ == 0) {
            return;
        }
        std::memset(ctrl_, CTRL_EMPTY, capacity_ + GROUP_WIDTH);
        size_ = 0;
        growth_left_ = capacity_to_growth(capacity_);
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    size_t size() const noexcept
    {
        return size_;
    }
    size_t capacity() const noexcept
    {
        return capacity_;
    }
    size_t growth_left() const noexcept
    {
        return growth_left_;
    }
    size_t mask() const noexcept
    {
        return mask_;
    }

    ctrl_t *ctrl() noexcept
    {
        return ctrl_;
    }
    const ctrl_t *ctrl() const noexcept
    {
        return ctrl_;
    }
    Key *keys() noexcept
    {
        return keys_;
    }
    const Key *keys() const noexcept
    {
        return keys_;
    }
    Value *vals() noexcept
    {
        return vals_;
    }
    const Value *vals() const noexcept
    {
        return vals_;
    }

    // =========================================================================
    // Find position of key (returns capacity_ if not found)
    // =========================================================================
    inline size_t find(const Key &key) const noexcept
    {
        if (capacity_ == 0) {
            return capacity_;
        }

        uint64_t hash = HashFn::hash(key);
        ctrl_t h2 = swiss_h2(hash);
        size_t index = hash & mask_;
        ctrl_t ctrl = ctrl_[index];

        // Fast path: first slot is empty (key not in table) - COMMON in sparse tables
        if (ctrl == CTRL_EMPTY) {
            return capacity_;
        }

        // Fast path: first slot matches - COMMON in cache-friendly access
        if (ctrl == h2 && EqualFn::equal(keys_[index], key)) {
            return index;
        }

        // SIMD path: scan groups
        ProbeSeq seq(hash, mask_);
        bool first_group = true;

        while (true) {
            size_t offset = seq.offset;
            Group g = Group::load(&ctrl_[offset]);
            uint16_t match_mask = g.match(h2);

            // Skip already-checked first slot (index) on first group
            if (first_group) {
                match_mask &= ~1;
                first_group = false;
            }

            while (match_mask != 0) {
                int bit = countr_zero(match_mask);
                size_t idx = (offset + bit) & mask_;
                if (EqualFn::equal(keys_[idx], key)) {
                    return idx;
                }
                match_mask &= match_mask - 1;
            }

            if (g.match_any_empty() != 0) {
                return capacity_;
            }

            seq.next();
        }
    }

    // =========================================================================
    // Insert or update key-value pair
    // Returns: 0 if key already existed (updated), 1 if newly inserted, -1 on error
    // =========================================================================
    inline int insert(const Key &key, const Value &val) noexcept
    {
        if (growth_left_ == 0) {
            size_t new_size = size_ == 0 ? GROUP_WIDTH : size_ * 2;
            if (!resize(new_size)) {
                return -1;
            }
        }

        return _insert(key, val);
    }

    // =========================================================================
    // Insert key only (no value) - for unique without return_inverse
    // Returns: 1 if newly inserted, 0 if key exists, -1 on memory error
    // =========================================================================
    inline int insert_key_only(const Key &key) noexcept
    {
        if (growth_left_ == 0) {
            size_t new_size = size_ == 0 ? GROUP_WIDTH : size_ * 2;
            if (!resize(new_size)) {
                return -1;
            }
        }

        return _insert_key_only(key);
    }

    // =========================================================================
    // Get value by key (returns true if found, false otherwise)
    // =========================================================================
    inline bool get(const Key &key, Value *val_out) const noexcept
    {
        size_t index = find(key);
        if (index == capacity_) {
            return false;
        }
        *val_out = vals_[index];
        return true;
    }

    // =========================================================================
    // Check if key exists
    // =========================================================================
    inline bool contains(const Key &key) const noexcept
    {
        return find(key) != capacity_;
    }

    // =========================================================================
    // Iteration support - get next valid index
    // Start with index = 0, call repeatedly until returns capacity
    // =========================================================================
    inline size_t iter_next(size_t index) const noexcept
    {
        while (index < capacity_) {
            if (is_full(ctrl_[index])) {
                return index;
            }
            index++;
        }
        return capacity_;
    }

    // =========================================================================
    // Get key at index (for iteration)
    // =========================================================================
    inline Key key_at(size_t index) const noexcept
    {
        return keys_[index];
    }

    // =========================================================================
    // Get value at index (for iteration)
    // =========================================================================
    inline Value val_at(size_t index) const noexcept
    {
        return vals_[index];
    }

    // =========================================================================
    // Increment value if key exists, insert with value 1 if not
    // Returns: 1 if newly inserted, 0 if incremented, -1 on error
    // =========================================================================
    inline int increment(const Key &key, Value *idx_out) noexcept
    {
        if (growth_left_ == 0) {
            size_t new_size = size_ == 0 ? GROUP_WIDTH : size_ * 2;
            if (!resize(new_size)) {
                return -1;
            }
        }

        return _increment(key, idx_out);
    }

    // =========================================================================
    // Insert if key doesn't exist, return existing value if it does
    // Returns: 1 if newly inserted, 0 if key exists (existing value in *val_out)
    // Returns: -1 on memory allocation error
    // =========================================================================
    inline int insert_if_absent(const Key &key, Value new_val, Value *val_out) noexcept
    {
        if (growth_left_ == 0) {
            size_t new_size = size_ == 0 ? GROUP_WIDTH : size_ * 2;
            if (!resize(new_size)) {
                return -1;
            }
        }

        return _insert_if_absent(key, new_val, val_out);
    }

    // =========================================================================
    // Reserve capacity for at least 'want' more insertions without resize
    // Returns: true on success, false on memory allocation failure
    // =========================================================================
    bool reserve(size_t want) noexcept
    {
        if (capacity_ == 0) {
            return resize(want);
        }
        if (growth_left_ >= want) {
            return true;  // Already have enough capacity
        }
        // Need to grow: ensure we can fit current size + want
        size_t need_total = size_ + want;
        return resize(need_total);
    }

    // =========================================================================
    // Batch factorize - optimized for 256-element batches
    // Returns: number of unique values on success, -1 on memory error
    // =========================================================================
    int64_t factorize_batch(const Key *keys, size_t n, Key *uniques_out, int64_t *labels_out,
        int64_t na_sentinel, Key nan_key, int64_t count_prior = 0) noexcept
    {
        // Reserve capacity upfront
        if (!reserve(n)) {
            return -1;
        }

        size_t count = static_cast<size_t>(count_prior);
        size_t i = 0;
        Value found_idx0, found_idx1, found_idx2, found_idx3;
        int ret0, ret1, ret2, ret3;

        for (; i + 4 <= n; i += 4) {
            const Key &k0 = keys[i + 0];
            const Key &k1 = keys[i + 1];
            const Key &k2 = keys[i + 2];
            const Key &k3 = keys[i + 3];

            if (EqualFn::equal(nan_key, k0)) {
                labels_out[i + 0] = na_sentinel;
            } else {
                ret0 = _insert_if_absent(k0, static_cast<Value>(count), &found_idx0);
                if (ret0 == 0) {
                    labels_out[i + 0] = static_cast<int64_t>(found_idx0);
                } else {
                    uniques_out[count] = k0;
                    labels_out[i + 0] = static_cast<int64_t>(count);
                    count++;
                }
            }

            if (EqualFn::equal(nan_key, k1)) {
                labels_out[i + 1] = na_sentinel;
            } else {
                ret1 = _insert_if_absent(k1, static_cast<Value>(count), &found_idx1);
                if (ret1 == 0) {
                    labels_out[i + 1] = static_cast<int64_t>(found_idx1);
                } else {
                    uniques_out[count] = k1;
                    labels_out[i + 1] = static_cast<int64_t>(count);
                    count++;
                }
            }

            if (EqualFn::equal(nan_key, k2)) {
                labels_out[i + 2] = na_sentinel;
            } else {
                ret2 = _insert_if_absent(k2, static_cast<Value>(count), &found_idx2);
                if (ret2 == 0) {
                    labels_out[i + 2] = static_cast<int64_t>(found_idx2);
                } else {
                    uniques_out[count] = k2;
                    labels_out[i + 2] = static_cast<int64_t>(count);
                    count++;
                }
            }

            if (EqualFn::equal(nan_key, k3)) {
                labels_out[i + 3] = na_sentinel;
            } else {
                ret3 = _insert_if_absent(k3, static_cast<Value>(count), &found_idx3);
                if (ret3 == 0) {
                    labels_out[i + 3] = static_cast<int64_t>(found_idx3);
                } else {
                    uniques_out[count] = k3;
                    labels_out[i + 3] = static_cast<int64_t>(count);
                    count++;
                }
            }
        }

        for (; i < n; i++) {
            const Key &key = keys[i];
            if (EqualFn::equal(nan_key, key)) {
                labels_out[i] = na_sentinel;
            } else {
                ret0 = _insert_if_absent(key, static_cast<Value>(count), &found_idx0);
                if (ret0 == 0) {
                    labels_out[i] = static_cast<int64_t>(found_idx0);
                } else {
                    uniques_out[count] = key;
                    labels_out[i] = static_cast<int64_t>(count);
                    count++;
                }
            }
        }

        return static_cast<int64_t>(count);
    }

    // =========================================================================
    // Batch value_count - optimized for 256-element batches
    // Returns: number of unique values on success, -1 on memory error
    // IMPORTANT: Assumes table is freshly initialized (no tombstones)
    // =========================================================================
    int64_t value_count_batch(const Key *keys, size_t n, Key *keys_out, size_t *indices_out, bool dropna) noexcept
    {
        // Reserve capacity upfront
        if (!reserve(n)) {
            return -1;
        }

        size_t count = 0;
        size_t i = 0;
        Value found_idx0, found_idx1, found_idx2, found_idx3;
        int ret0, ret1, ret2, ret3;

        for (; i + 4 <= n; i += 4) {
            const Key &k0 = keys[i + 0];
            const Key &k1 = keys[i + 1];
            const Key &k2 = keys[i + 2];
            const Key &k3 = keys[i + 3];

            if (!(dropna && NaNTraits<Key>::IsNaN(k0))) {
                ret0 = _increment(k0, &found_idx0);
                if (ret0 == 1) {
                    keys_out[count] = k0;
                    indices_out[count] = found_idx0;
                    count++;
                }
            }

            if (!(dropna && NaNTraits<Key>::IsNaN(k1))) {
                ret1 = _increment(k1, &found_idx1);
                if (ret1 == 1) {
                    keys_out[count] = k1;
                    indices_out[count] = found_idx1;
                    count++;
                }
            }

            if (!(dropna && NaNTraits<Key>::IsNaN(k2))) {
                ret2 = _increment(k2, &found_idx2);
                if (ret2 == 1) {
                    keys_out[count] = k2;
                    indices_out[count] = found_idx2;
                    count++;
                }
            }

            if (!(dropna && NaNTraits<Key>::IsNaN(k3))) {
                ret3 = _increment(k3, &found_idx3);
                if (ret3 == 1) {
                    keys_out[count] = k3;
                    indices_out[count] = found_idx3;
                    count++;
                }
            }
        }

        for (; i < n; i++) {
            const Key &key = keys[i];
            if (!(dropna && NaNTraits<Key>::IsNaN(key))) {
                ret0 = _increment(key, &found_idx0);
                if (ret0 == 1) {
                    keys_out[count] = key;
                    indices_out[count] = found_idx0;
                    count++;
                }
            }
        }

        return static_cast<int64_t>(count);
    }

    // =========================================================================
    // Integration Methods - for pandas internal operations
    // =========================================================================

    // ------------------------------------------------------------------------
    // map_locations - Build table mapping values to their array positions
    //
    // For each key in keys[0..n-1], stores (key -> index) in the table.
    // If a key appears multiple times, the last occurrence wins.
    //
    // Returns: 0 on success, -1 on memory allocation failure
    // Used in: libindex, safe_sort
    // ------------------------------------------------------------------------
    int map_locations(const Key *keys, size_t n) noexcept
    {
        // Reserve capacity upfront
        if (!reserve(n)) {
            return -1;
        }

        // Insert each key with its index as value
        for (size_t i = 0; i < n; i++) {
            int ret = insert(keys[i], i);
            if (ret == -1) {
                return -1;
            }
        }

        return 0;
    }

    // ------------------------------------------------------------------------
    // lookup_batch - Look up array keys in table, returning positions
    //
    // For each key in keys[0..n-1], stores the corresponding value (position)
    // in locs[i]. If key not found, locs[i] = -1.
    //
    // Used in: safe_sort, IndexEngine.get_indexer
    // ------------------------------------------------------------------------
    void lookup_batch(const Key *keys, size_t n, int64_t *locs) const noexcept
    {
        if (capacity_ == 0) {
            std::fill(locs, locs + n, -1);
            return;
        }

        size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            size_t idx0 = find(keys[i + 0]);
            locs[i + 0] = (idx0 != capacity_) ? static_cast<int64_t>(vals_[idx0]) : -1;
            size_t idx1 = find(keys[i + 1]);
            locs[i + 1] = (idx1 != capacity_) ? static_cast<int64_t>(vals_[idx1]) : -1;
            size_t idx2 = find(keys[i + 2]);
            locs[i + 2] = (idx2 != capacity_) ? static_cast<int64_t>(vals_[idx2]) : -1;
            size_t idx3 = find(keys[i + 3]);
            locs[i + 3] = (idx3 != capacity_) ? static_cast<int64_t>(vals_[idx3]) : -1;
        }

        for (; i < n; i++) {
            size_t idx = find(keys[i]);
            locs[i] = (idx != capacity_) ? static_cast<int64_t>(vals_[idx]) : -1;
        }
    }

    // ------------------------------------------------------------------------
    // build_set - Build a set from keys (no values, just keys)
    //
    // Inserts each key into the table using insert_key_only (value ignored).
    //
    // Returns: 0 on success, -1 on memory allocation failure
    // Used in: ismember, unique operations
    // ------------------------------------------------------------------------
    int build_set(const Key *keys, size_t n) noexcept
    {
        // Reserve capacity upfront
        if (!reserve(n)) {
            return -1;
        }
        size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            (void)_insert_key_only(keys[i + 0]);
            (void)_insert_key_only(keys[i + 1]);
            (void)_insert_key_only(keys[i + 2]);
            (void)_insert_key_only(keys[i + 3]);
        }
        for (; i < n; i++) {
            (void)_insert_key_only(keys[i]);
        }
        return 0;
    }

    // ------------------------------------------------------------------------
    // contains_batch - Batch membership test
    //
    // For each key in keys[0..n-1], sets result[i] = 1 if key exists, else 0.
    //
    // Used in: ismember operations
    // ------------------------------------------------------------------------
    void contains_batch(const Key *keys, size_t n, uint8_t *result) const noexcept
    {
        if (capacity_ == 0) {
            std::fill(result, result + n, 0);
            return;
        }

        size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            result[i + 0] = (find(keys[i + 0]) != capacity_) ? 1 : 0;
            result[i + 1] = (find(keys[i + 1]) != capacity_) ? 1 : 0;
            result[i + 2] = (find(keys[i + 2]) != capacity_) ? 1 : 0;
            result[i + 3] = (find(keys[i + 3]) != capacity_) ? 1 : 0;
        }

        for (; i < n; i++) {
            result[i] = (find(keys[i]) != capacity_) ? 1 : 0;
        }
    }

    int64_t unique_batch(const Key *keys, size_t n, Key *uniques_out) noexcept
    {
        // Reserve capacity upfront
        if (!reserve(n)) {
            return -1;
        }
        size_t count = 0;
        size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            if (_insert_key_only(keys[i + 0]) > 0) {
                uniques_out[count++] = keys[i + 0];
            }
            if (_insert_key_only(keys[i + 1]) > 0) {
                uniques_out[count++] = keys[i + 1];
            }
            if (_insert_key_only(keys[i + 2]) > 0) {
                uniques_out[count++] = keys[i + 2];
            }
            if (_insert_key_only(keys[i + 3]) > 0) {
                uniques_out[count++] = keys[i + 3];
            }
        }
        for (; i < n; i++) {
            if (_insert_key_only(keys[i]) > 0) {
                uniques_out[count++] = keys[i];
            }
        }
        return static_cast<int64_t>(count);
    }

    int64_t unique_with_inverse(const Key *keys, size_t n, Key *uniques_out, int64_t *labels_out,
        int64_t count_prior = 0) noexcept
    {
        // Reserve capacity upfront
        if (!reserve(n)) {
            return -1;
        }

        size_t count = static_cast<size_t>(count_prior);
        size_t i = 0;
        Value found_idx0, found_idx1, found_idx2, found_idx3;
        int ret0, ret1, ret2, ret3;

        for (; i + 4 <= n; i += 4) {
            const Key &k0 = keys[i + 0];
            const Key &k1 = keys[i + 1];
            const Key &k2 = keys[i + 2];
            const Key &k3 = keys[i + 3];

            ret0 = _insert_if_absent(k0, static_cast<Value>(count), &found_idx0);
            if (ret0 == 0) {
                labels_out[i + 0] = static_cast<int64_t>(found_idx0);
            } else {
                labels_out[i + 0] = static_cast<int64_t>(count);
                uniques_out[count] = k0;
                count++;
            }

            ret1 = _insert_if_absent(k1, static_cast<Value>(count), &found_idx1);
            if (ret1 == 0) {
                labels_out[i + 1] = static_cast<int64_t>(found_idx1);
            } else {
                labels_out[i + 1] = static_cast<int64_t>(count);
                uniques_out[count] = k1;
                count++;
            }

            ret2 = _insert_if_absent(k2, static_cast<Value>(count), &found_idx2);
            if (ret2 == 0) {
                labels_out[i + 2] = static_cast<int64_t>(found_idx2);
            } else {
                labels_out[i + 2] = static_cast<int64_t>(count);
                uniques_out[count] = k2;
                count++;
            }

            ret3 = _insert_if_absent(k3, static_cast<Value>(count), &found_idx3);
            if (ret3 == 0) {
                labels_out[i + 3] = static_cast<int64_t>(found_idx3);
            } else {
                labels_out[i + 3] = static_cast<int64_t>(count);
                uniques_out[count] = k3;
                count++;
            }
        }

        for (; i < n; i++) {
            const Key &key = keys[i];
            /* First, try to find the key (fast path for duplicates) */
            ret0 = _insert_if_absent(key, static_cast<Value>(count), &found_idx0);
            if (ret0 == 0) {
                /* Key exists - this is the common case for duplicate-heavy data */
                labels_out[i] = static_cast<int64_t>(found_idx0);
            } else {
                labels_out[i] = static_cast<int64_t>(count);
                uniques_out[count] = key;
                count++;
            }
        }

        return static_cast<int64_t>(count);
    }

    int duplicated_keep_batch(const Key *keys, size_t n, uint8_t keep_first, uint8_t *result) noexcept
    {
        // Reserve capacity upfront
        if (!reserve(n)) {
            return -1;
        }
        if (keep_first) {
            size_t i = 0;
            for (; i + 4 <= n; i += 4) {
                result[i + 0] = (_insert_key_only(keys[i + 0]) == 0) ? 1 : 0;
                result[i + 1] = (_insert_key_only(keys[i + 1]) == 0) ? 1 : 0;
                result[i + 2] = (_insert_key_only(keys[i + 2]) == 0) ? 1 : 0;
                result[i + 3] = (_insert_key_only(keys[i + 3]) == 0) ? 1 : 0;
            }
            for (; i < n; i++) {
                result[i] = (_insert_key_only(keys[i]) == 0) ? 1 : 0;
            }
        } else {
            size_t i = n;
            for (; i >= 4; i -= 4) {
                result[i - 1] = (_insert_key_only(keys[i - 1]) == 0) ? 1 : 0;
                result[i - 2] = (_insert_key_only(keys[i - 2]) == 0) ? 1 : 0;
                result[i - 3] = (_insert_key_only(keys[i - 3]) == 0) ? 1 : 0;
                result[i - 4] = (_insert_key_only(keys[i - 4]) == 0) ? 1 : 0;
            }
            for (; i > 0; --i) {
                result[i - 1] = (_insert_key_only(keys[i - 1]) == 0) ? 1 : 0;
            }
        }
        return 0;
    }

    int duplicated_false_batch(const Key *keys, size_t n, uint8_t *result) noexcept
    {
        // Reserve capacity upfront
        if (!reserve(n)) {
            return -1;
        }

        size_t i = 0;
        Value found_idx0, found_idx1, found_idx2, found_idx3;
        int ret0, ret1, ret2, ret3;

        for (; i + 4 <= n; i += 4) {
            const Key &k0 = keys[i + 0];
            const Key &k1 = keys[i + 1];
            const Key &k2 = keys[i + 2];
            const Key &k3 = keys[i + 3];

            ret0 = _insert_if_absent(k0, static_cast<Value>(i + 0), &found_idx0);
            if (ret0 == 0) {
                result[i + 0] = 1;
                result[found_idx0] = 1;
            } else {
                result[i + 0] = 0;
            }

            ret1 = _insert_if_absent(k1, static_cast<Value>(i + 1), &found_idx1);
            if (ret1 == 0) {
                result[i + 1] = 1;
                result[found_idx1] = 1;
            } else {
                result[i + 1] = 0;
            }

            ret2 = _insert_if_absent(k2, static_cast<Value>(i + 2), &found_idx2);
            if (ret2 == 0) {
                result[i + 2] = 1;
                result[found_idx2] = 1;
            } else {
                result[i + 2] = 0;
            }

            ret3 = _insert_if_absent(k3, static_cast<Value>(i + 3), &found_idx3);
            if (ret3 == 0) {
                result[i + 3] = 1;
                result[found_idx3] = 1;
            } else {
                result[i + 3] = 0;
            }
        }

        for (; i < n; ++i) {
            /* First, try to find the key (fast path for duplicates) */
            ret0 = _insert_if_absent(keys[i], static_cast<Value>(i), &found_idx0);
            if (ret0 == 0) {
                /* Key exists - this is the common case for duplicate-heavy data */
                result[i] = 1;
                result[found_idx0] = 1;
            } else {
                result[i] = 0;
            }
        }

        return 0;
    }

private:
    SWISS_FORCE_INLINE int _insert_key_only(const Key &key) noexcept
    {
        uint64_t hash = HashFn::hash(key);
        ctrl_t h2 = swiss_h2(hash);
        size_t index = hash & mask_;
        ctrl_t c0 = ctrl_[index];

        // Fast path: check if first slot is empty
        if (c0 == CTRL_EMPTY) {
            set_ctrl(index, h2);
            keys_[index] = key;
            size_++;
            growth_left_--;
            return 1;
        }

        // Fast path: check if first slot has matching key
        if (c0 == h2 && EqualFn::equal(keys_[index], key)) {
            return 0;
        }

        // Slow path: SIMD group operations
        ProbeSeq seq(hash, mask_);
        bool first_group = true;

        while (true) {
            size_t offset = seq.offset;
            Group g = Group::load(&ctrl_[offset]);
            uint16_t match_mask = g.match(h2);

            // Skip already-checked first slot (index) on first group
            if (first_group) {
                match_mask &= ~1;
                first_group = false;
            }

            while (match_mask != 0) {
                int bit = countr_zero(match_mask);
                size_t idx = (offset + bit) & mask_;
                if (EqualFn::equal(keys_[idx], key)) {
                    return 0;
                }
                match_mask &= match_mask - 1;
            }

            uint16_t empty_mask = g.match_empty();
            if (empty_mask != 0) {
                int bit = countr_zero(empty_mask);
                size_t idx = (offset + bit) & mask_;
                set_ctrl(idx, h2);
                keys_[idx] = key;
                size_++;
                growth_left_--;
                return 1;
            }

            seq.next();
        }
    }

    SWISS_FORCE_INLINE int _insert(const Key &key, const Value &val) noexcept
    {
        uint64_t hash = HashFn::hash(key);
        ctrl_t h2 = swiss_h2(hash);
        size_t index = hash & mask_;
        ctrl_t c0 = ctrl_[index];

        // Fast path: first slot is empty (common case for sparse tables)
        if (c0 == CTRL_EMPTY) {
            set_ctrl(index, h2);
            keys_[index] = key;
            vals_[index] = val;
            size_++;
            growth_left_--;
            return 1;
        }

        // Fast path: first slot has matching key
        if (c0 == h2 && EqualFn::equal(keys_[index], key)) {
            vals_[index] = val;
            return 0;
        }

        // Slow path: need to check for existing key or find another empty slot
        ProbeSeq seq(hash, mask_);
        bool first_group = true;

        while (true) {
            size_t offset = seq.offset;

            Group g = Group::load(&ctrl_[offset]);
            uint16_t match_mask = g.match(h2);

            // Skip already-checked first slot (index) on first group
            if (first_group) {
                match_mask &= ~1;
                first_group = false;
            }

            while (match_mask != 0) {
                int bit = countr_zero(match_mask);
                size_t idx = (offset + bit) & mask_;
                if (EqualFn::equal(keys_[idx], key)) {
                    vals_[idx] = val;
                    return 0;
                }
                match_mask &= match_mask - 1;
            }

            uint16_t empty_mask = g.match_empty();
            if (empty_mask != 0) {
                int bit = countr_zero(empty_mask);
                size_t idx = (offset + bit) & mask_;
                set_ctrl(idx, h2);
                keys_[idx] = key;
                vals_[idx] = val;
                size_++;
                growth_left_--;
                return 1;
            }

            seq.next();
        }
    }

    SWISS_FORCE_INLINE int _insert_if_absent(const Key &key, Value new_val, Value *val_out) noexcept
    {
        uint64_t hash = HashFn::hash(key);
        ctrl_t h2 = swiss_h2(hash);
        size_t index = hash & mask_;
        ctrl_t c0 = ctrl_[index];

        // Fast path: first slot is empty
        if (c0 == CTRL_EMPTY) {
            set_ctrl(index, h2);
            keys_[index] = key;
            vals_[index] = new_val;
            size_++;
            growth_left_--;
            return 1;
        }

        // Fast path: first slot matches
        if (c0 == h2 && EqualFn::equal(keys_[index], key)) {
            *val_out = vals_[index];
            return 0;
        }

        // SIMD search
        ProbeSeq seq(hash, mask_);
        bool first_group = true;

        while (true) {
            size_t offset = seq.offset;
            Group g = Group::load(&ctrl_[offset]);
            uint16_t match_mask = g.match(h2);

            if (first_group) {
                match_mask &= ~1;
                first_group = false;
            }

            while (match_mask != 0) {
                int bit = countr_zero(match_mask);
                size_t idx = (offset + bit) & mask_;
                if (EqualFn::equal(keys_[idx], key)) {
                    *val_out = vals_[idx];
                    return 0;
                }
                match_mask &= match_mask - 1;
            }

            uint16_t empty_mask = g.match_empty();
            if (empty_mask != 0) {
                int bit = countr_zero(empty_mask);
                size_t idx = (offset + bit) & mask_;
                set_ctrl(idx, h2);
                keys_[idx] = key;
                vals_[idx] = new_val;
                size_++;
                growth_left_--;
                return 1;
            }

            seq.next();
        }
    }

    SWISS_FORCE_INLINE int _increment(const Key &key, Value *idx_out) noexcept
    {
        uint64_t hash = HashFn::hash(key);
        ctrl_t h2 = swiss_h2(hash);
        size_t index = hash & mask_;
        ctrl_t c0 = ctrl_[index];

        // Fast path: first slot is empty
        if (c0 == CTRL_EMPTY) {
            set_ctrl(index, h2);
            keys_[index] = key;
            vals_[index] = 1;
            *idx_out = index;
            size_++;
            growth_left_--;
            return 1;
        }

        // Fast path: first slot matches - increment
        if (c0 == h2 && EqualFn::equal(keys_[index], key)) {
            vals_[index]++;
            *idx_out = index;
            return 0;
        }

        // Slow path
        ProbeSeq seq(hash, mask_);
        bool first_group = true;

        while (true) {
            size_t offset = seq.offset;
            Group g = Group::load(&ctrl_[offset]);
            uint16_t match_mask = g.match(h2);

            if (first_group) {
                match_mask &= ~1;
                first_group = false;
            }

            while (match_mask != 0) {
                int bit = countr_zero(match_mask);
                size_t idx = (offset + bit) & mask_;
                if (EqualFn::equal(keys_[idx], key)) {
                    vals_[idx]++;
                    *idx_out = idx;
                    return 0;
                }
                match_mask &= match_mask - 1;
            }

            uint16_t empty_mask = g.match_empty();
            if (empty_mask != 0) {
                int bit = countr_zero(empty_mask);
                size_t idx = (offset + bit) & mask_;
                set_ctrl(idx, h2);
                keys_[idx] = key;
                vals_[idx] = 1;
                *idx_out = idx;
                size_++;
                growth_left_--;
                return 1;
            }

            seq.next();
        }
    }

    // =========================================================================
    // Helper Functions
    // =========================================================================

    bool resize(size_t new_capacity) noexcept
    {
        SwissTable new_table;
        if (!new_table.init_with_capacity(new_capacity)) {
            return false;
        }

        // Rehash all existing elements
        if (size_ > 0) {
            for (size_t i = 0; i < capacity_; i++) {
                if (is_full(ctrl_[i])) {
                    auto &key = keys_[i];
                    Value val = vals_[i];
                    uint64_t hash = HashFn::hash(key);
                    ctrl_t h2 = swiss_h2(hash);
                    ProbeSeq seq(hash, new_table.mask_);

                    while (true) {
                        size_t offset = seq.offset;
                        Group g = Group::load(&new_table.ctrl_[offset]);
                        uint16_t empty_mask = g.match_empty();

                        if (empty_mask != 0) {
                            int bit = countr_zero(empty_mask);
                            size_t index = (offset + bit) & new_table.mask_;
                            // Use set_ctrl for correct sentinel byte handling
                            new_table.set_ctrl(index, h2);
                            new_table.keys_[index] = key;
                            new_table.vals_[index] = val;
                            new_table.size_++;
                            new_table.growth_left_--;
                            break;
                        }

                        seq.next();
                    }
                }
            }
        }

        // Swap with new table (single free for old allocation)
        std::swap(capacity_, new_table.capacity_);
        std::swap(mask_, new_table.mask_);
        std::swap(size_, new_table.size_);
        std::swap(growth_left_, new_table.growth_left_);
        std::swap(ctrl_, new_table.ctrl_);
        std::swap(keys_, new_table.keys_);
        std::swap(vals_, new_table.vals_);
        std::swap(alloc_, new_table.alloc_);
        return true;
    }

    // ------------------------------------------------------------------------
    // set_ctrl - Control byte setting with sentinel cloning
    //
    // Memory layout (matches C implementation):
    //   ctrl[0..capacity-1]: Actual control bytes
    //   ctrl[capacity..capacity+15]: Cloned bytes (copy of ctrl[0..15])
    //
    // Note: The sentinel byte is at ctrl[capacity] which is the first cloned byte.
    // The C implementation sets ctrl[capacity + index] for index < GROUP_WIDTH.
    // ------------------------------------------------------------------------
    inline void set_ctrl(size_t index, ctrl_t value) noexcept
    {
        ctrl_[index] = value;
        // Clone to sentinel area for index < GROUP_WIDTH
        if (index < GROUP_WIDTH) {
            ctrl_[capacity_ + index] = value;
        }
    }

    static inline size_t calc_alloc_size(size_t capacity, size_t *keys_offset_out, size_t *vals_offset_out) noexcept
    {
        size_t ctrl_size = capacity + GROUP_WIDTH;
        // Align keys to max of Key size and GROUP_WIDTH bytes for SIMD
        size_t key_align = sizeof(Key) > GROUP_WIDTH ? sizeof(Key) : GROUP_WIDTH;
        size_t keys_offset = (ctrl_size + key_align - 1) & ~(key_align - 1);
        size_t keys_size = capacity * sizeof(Key);
        // Align vals to Value size
        size_t vals_offset = (keys_offset + keys_size + sizeof(Value) - 1) & ~(sizeof(Value) - 1);
        size_t vals_size = capacity * sizeof(Value);

        if (keys_offset_out)
            *keys_offset_out = keys_offset;
        if (vals_offset_out)
            *vals_offset_out = vals_offset;

        return vals_offset + vals_size;
    }

private:
    // =========================================================================
    // Member Variables
    // =========================================================================
    size_t capacity_;     // Must be power of 2, >= GROUP_WIDTH
    size_t mask_;         // capacity - 1 (for fast modulo)
    size_t size_;         // Number of elements
    size_t growth_left_;  // Elements until resize
    ctrl_t *ctrl_;        // -> ctrl bytes within alloc
    Key *keys_;           // -> keys within alloc
    Value *vals_;         // -> vals within alloc
    void *alloc_;         // Single allocation block
};

}  // namespace pandas::swisstable