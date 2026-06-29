/*
 * Type-Generic Swiss Table Implementation for pandas
 *
 * Usage:
 *   // Define a Swiss Table for int64_t keys -> size_t values:
 *   SWISSTABLE_INIT(int64, int64_t, size_t, hash_int64, key_equal_int64)
 *
 *   // This creates:
 *   //   - SwissTable_int64 (the struct type)
 *   //   - swiss_int64_init(), swiss_int64_destroy(), etc.
 *
 * Similar to khash's KHASH_INIT macro system.
 *
 * Memory Layout (single allocation for cache efficiency):
 *   [ctrl bytes: capacity + 16] [padding] [keys: capacity] [vals: capacity]
 *
 * This reduces cache misses from 3 per lookup to 1-2 since ctrl/keys/vals
 * are in contiguous memory.
 */

#ifndef PANDAS_SWISSTABLE_GENERIC_H
#define PANDAS_SWISSTABLE_GENERIC_H

#include "swisstable_core.h"
#include <stdlib.h>

/* ssize_t for factorize_batch - needed for pandas intp_t compatibility */
#ifdef _MSC_VER
  /* Windows: ssize_t is not standard, use intptr_t */
  typedef intptr_t ssize_t;
#else
  #include <sys/types.h>  /* POSIX ssize_t */
#endif

/* Calculate aligned offset for keys array after ctrl bytes */
/* Ensures proper alignment for key_t type */
static inline size_t swiss_align_offset(size_t offset, size_t alignment) {
    return (offset + alignment - 1) & ~(alignment - 1);
}

/*
 * SWISSTABLE_INIT - Generate a complete Swiss Table implementation
 *
 * Parameters:
 *   name      - Suffix for function/type names (e.g., "int64" -> SwissTable_int64)
 *   key_t     - Key type (e.g., int64_t)
 *   val_t     - Value type (e.g., size_t)
 *   hash_fn   - Hash function: uint64_t hash_fn(key_t key)
 *   equal_fn  - Equality function: bool equal_fn(key_t a, key_t b)
 */
#define SWISSTABLE_INIT(name, key_t, val_t, hash_fn, equal_fn) \
    SWISSTABLE_INIT_IMPL(name, key_t, val_t, hash_fn, equal_fn)

#define SWISSTABLE_INIT_IMPL(name, key_t, val_t, hash_fn, equal_fn) \
\
/* Table structure - single allocation for better cache locality */ \
typedef struct { \
    size_t capacity;      /* Must be power of 2, >= 16 */ \
    size_t size;          /* Number of elements */ \
    size_t deleted;       /* Number of tombstones (deleted slots) */ \
    size_t growth_left;   /* Elements until resize */ \
    size_t mask;          /* capacity - 1 (for fast modulo) */ \
    void *alloc;          /* Single allocation block */ \
    int8_t *ctrl;         /* -> ctrl bytes within alloc */ \
    key_t *keys;          /* -> keys within alloc */ \
    val_t *vals;          /* -> vals within alloc */ \
} SwissTable_##name; \
\
/* Calculate total allocation size and offsets for single-block layout */ \
/* Layout: [ctrl: capacity+16] [padding] [keys: capacity] [vals: capacity] */ \
static inline size_t swiss_##name##_calc_alloc_size(size_t capacity, \
        size_t *keys_offset_out, size_t *vals_offset_out) { \
    size_t ctrl_size = capacity + 16; \
    /* Align keys to max of key_t and 16 bytes for SIMD */ \
    size_t key_align = sizeof(key_t) > 16 ? sizeof(key_t) : 16; \
    size_t keys_offset = swiss_align_offset(ctrl_size, key_align); \
    size_t keys_size = capacity * sizeof(key_t); \
    /* Align vals to val_t size */ \
    size_t vals_offset = swiss_align_offset(keys_offset + keys_size, sizeof(val_t)); \
    size_t vals_size = capacity * sizeof(val_t); \
    if (keys_offset_out) *keys_offset_out = keys_offset; \
    if (vals_offset_out) *vals_offset_out = vals_offset; \
    return vals_offset + vals_size; \
} \
\
/* Initialize empty table */ \
static inline SwissTable_##name* swiss_##name##_init(void) { \
    SwissTable_##name *table = (SwissTable_##name*)SWISSTABLE_MALLOC(sizeof(SwissTable_##name)); \
    if (!table) return NULL; \
    table->capacity = 0; \
    table->size = 0; \
    table->deleted = 0; \
    table->growth_left = 0; \
    table->mask = 0; \
    table->alloc = NULL; \
    table->ctrl = NULL; \
    table->keys = NULL; \
    table->vals = NULL; \
    return table; \
} \
\
/* Initialize with specific capacity - single allocation */ \
static inline SwissTable_##name* swiss_##name##_init_with_capacity(size_t capacity) { \
    SwissTable_##name *table = swiss_##name##_init(); \
    if (!table) return NULL; \
    if (capacity > 0) { \
        capacity = normalize_capacity(capacity); \
        size_t keys_offset, vals_offset; \
        size_t alloc_size = swiss_##name##_calc_alloc_size(capacity, &keys_offset, &vals_offset); \
        table->alloc = SWISSTABLE_MALLOC(alloc_size); \
        if (!table->alloc) { \
            SWISSTABLE_FREE(table); \
            return NULL; \
        } \
        /* Set up pointers into the single allocation */ \
        table->ctrl = (int8_t*)table->alloc; \
        table->keys = (key_t*)((char*)table->alloc + keys_offset); \
        table->vals = (val_t*)((char*)table->alloc + vals_offset); \
        memset(table->ctrl, CTRL_EMPTY, capacity + 16); \
        table->capacity = capacity; \
        table->mask = capacity - 1; \
        table->growth_left = capacity_to_growth(capacity); \
    } \
    return table; \
} \
\
/* Destroy table and free memory - single free */ \
static inline void swiss_##name##_destroy(SwissTable_##name *table) { \
    if (!table) return; \
    SWISSTABLE_FREE(table->alloc); \
    SWISSTABLE_FREE(table); \
} \
\
/* Clear all elements but keep capacity */ \
static inline void swiss_##name##_clear(SwissTable_##name *table) { \
    if (!table || table->capacity == 0) return; \
    memset(table->ctrl, CTRL_EMPTY, table->capacity + 16); \
    table->size = 0; \
    table->deleted = 0; \
    table->growth_left = capacity_to_growth(table->capacity); \
} \
\
/* Find position of key (returns capacity if not found) */ \
static inline size_t swiss_##name##_find(const SwissTable_##name *table, key_t key) { \
    if (table->capacity == 0) { \
        return table->capacity; \
    } \
    uint64_t hash = hash_fn(key); \
    int8_t h2 = swiss_h2(hash); \
    size_t index = hash & table->mask; \
    int8_t ctrl = table->ctrl[index]; \
    /* Fast path: first slot is empty (key not in table) */ \
    if (ctrl == CTRL_EMPTY) { \
        return table->capacity; \
    } \
    /* Fast path: first slot matches */ \
    if (ctrl == h2 && equal_fn(table->keys[index], key)) { \
        return index; \
    } \
    /* SIMD path: scan groups */ \
    ProbeSeq seq; \
    probe_seq_init(&seq, hash, table->mask); \
    bool first_group = true; \
    while (true) { \
        size_t offset = probe_seq_offset(&seq); \
        Group g = group_load(&table->ctrl[offset]); \
        uint16_t match_mask = group_match(g, h2); \
        /* Skip already-checked first slot (index) on first group */ \
        /* Compute relative position within the group */ \
        if (first_group) { \
            size_t rel = (index - offset) & table->mask; \
            if (rel < 16) { \
                match_mask &= (uint16_t)~((uint16_t)1u << (uint16_t)rel); \
            } \
            first_group = false; \
        } \
        while (match_mask != 0) { \
            int bit = countr_zero(match_mask); \
            size_t idx = (offset + bit) & table->mask; \
            if (equal_fn(table->keys[idx], key)) { \
                return idx; \
            } \
            match_mask &= match_mask - 1; \
        } \
        if (group_match_empty(g) != 0) { \
            return table->capacity; \
        } \
        probe_seq_next(&seq); \
    } \
} \
\
/* Forward declaration for resize */ \
static inline bool swiss_##name##_resize(SwissTable_##name *table, size_t new_capacity); \
\
/* Insert or update key-value pair */ \
/* Returns: 0 if key already existed (updated), 1 if newly inserted, -1 on error */ \
static inline int swiss_##name##_insert(SwissTable_##name *table, key_t key, val_t val) { \
    if (table->growth_left == 0) { \
        size_t new_capacity = table->capacity == 0 ? 16 : table->capacity * 2; \
        if (!swiss_##name##_resize(table, new_capacity)) { \
            return -1; \
        } \
    } \
    uint64_t hash = hash_fn(key); \
    int8_t h2 = swiss_h2(hash); \
    size_t index = hash & table->mask; \
    int8_t c0 = table->ctrl[index]; \
    /* Fast path: check if first slot is empty (common case for sparse tables) */ \
    if (c0 == CTRL_EMPTY) { \
        table->ctrl[index] = h2; \
        table->keys[index] = key; \
        table->vals[index] = val; \
        table->size++; \
        table->growth_left--; \
        if (index < 16) { \
            table->ctrl[table->capacity + index] = h2; \
        } \
        return 1; \
    } \
    /* Fast path: check if first slot has matching key */ \
    if (c0 == h2 && equal_fn(table->keys[index], key)) { \
        table->vals[index] = val; \
        return 0; \
    } \
    /* Slow path: need to check for existing key or find another empty slot */ \
    ProbeSeq seq; \
    probe_seq_init(&seq, hash, table->mask); \
    bool first_group = true; \
    while (true) { \
        size_t offset = probe_seq_offset(&seq); \
        Group g = group_load(&table->ctrl[offset]); \
        uint16_t match_mask = group_match(g, h2); \
        /* Skip already-checked first slot (index) on first group */ \
        /* Compute relative position within the group */ \
        if (first_group) { \
            size_t rel = (index - offset) & table->mask; \
            if (rel < 16) { \
                match_mask &= (uint16_t)~((uint16_t)1u << (uint16_t)rel); \
            } \
            first_group = false; \
        } \
        while (match_mask != 0) { \
            int bit = countr_zero(match_mask); \
            size_t idx = (offset + bit) & table->mask; \
            if (equal_fn(table->keys[idx], key)) { \
                table->vals[idx] = val; \
                return 0; \
            } \
            match_mask &= match_mask - 1; \
        } \
        uint16_t empty_mask = group_match_empty_or_deleted(g); \
        if (empty_mask != 0) { \
            int bit = countr_zero(empty_mask); \
            size_t idx = (offset + bit) & table->mask; \
            int8_t old_ctrl = table->ctrl[idx]; \
            table->ctrl[idx] = h2; \
            table->keys[idx] = key; \
            table->vals[idx] = val; \
            table->size++; \
            if (old_ctrl == CTRL_EMPTY) { \
                table->growth_left--; \
            } else { \
                /* Was CTRL_DELETED - reusing tombstone */ \
                table->deleted--; \
                table->growth_left--; /* Still consuming a slot */ \
            } \
            if (idx < 16) { \
                table->ctrl[table->capacity + idx] = h2; \
            } \
            return 1; \
        } \
        probe_seq_next(&seq); \
    } \
} \
\
/* Resize table to new capacity - single allocation */ \
static inline bool swiss_##name##_resize(SwissTable_##name *table, size_t new_capacity) { \
    new_capacity = normalize_capacity(new_capacity); \
    SwissTable_##name *new_table = swiss_##name##_init_with_capacity(new_capacity); \
    if (!new_table) return false; \
    for (size_t i = 0; i < table->capacity; i++) { \
        if (is_full(table->ctrl[i])) { \
            key_t key = table->keys[i]; \
            val_t val = table->vals[i]; \
            uint64_t hash = hash_fn(key); \
            int8_t h2 = swiss_h2(hash); \
            ProbeSeq seq; \
            probe_seq_init(&seq, hash, new_capacity - 1); \
            while (true) { \
                size_t offset = probe_seq_offset(&seq); \
                Group g = group_load(&new_table->ctrl[offset]); \
                uint16_t empty_mask = group_match_empty(g); \
                if (empty_mask != 0) { \
                    int bit = countr_zero(empty_mask); \
                    size_t index = (offset + bit) & (new_capacity - 1); \
                    new_table->ctrl[index] = h2; \
                    new_table->keys[index] = key; \
                    new_table->vals[index] = val; \
                    new_table->size++; \
                    break; \
                } \
                probe_seq_next(&seq); \
            } \
        } \
    } \
    memcpy(&new_table->ctrl[new_capacity], new_table->ctrl, 16); \
    /* Single free for old allocation */ \
    void *old_alloc = table->alloc; \
    table->capacity = new_table->capacity; \
    table->size = new_table->size; \
    table->growth_left = capacity_to_growth(new_capacity) - new_table->size; \
    table->mask = table->capacity - 1; \
    table->alloc = new_table->alloc; \
    table->ctrl = new_table->ctrl; \
    table->keys = new_table->keys; \
    table->vals = new_table->vals; \
    SWISSTABLE_FREE(old_alloc); \
    SWISSTABLE_FREE(new_table); \
    return true; \
} \
\
/* Reserve capacity for at least 'want' more insertions without resize */ \
/* Call before tight loops to eliminate resize checks from hot path */ \
/* Returns: true on success, false on memory allocation failure */ \
static inline bool swiss_##name##_reserve(SwissTable_##name *table, size_t want) { \
    if (table->capacity == 0) { \
        /* Empty table: allocate enough for 'want' at 13/16 load factor */ \
        size_t cap = normalize_capacity((want * 16 + 12) / 13); \
        return swiss_##name##_resize(table, cap); \
    } \
    if (table->growth_left >= want) { \
        return true; /* Already have enough capacity */ \
    } \
    /* Need to grow: ensure we can fit current size + want */ \
    size_t need_total = table->size + want; \
    size_t cap = normalize_capacity((need_total * 16 + 12) / 13); \
    return swiss_##name##_resize(table, cap); \
} \
\
/* Delete key from table */ \
/* Returns: true if key was found and deleted, false if not found */ \
static inline bool swiss_##name##_delete(SwissTable_##name *table, key_t key) { \
    size_t index = swiss_##name##_find(table, key); \
    if (index == table->capacity) { \
        return false; \
    } \
    table->ctrl[index] = CTRL_DELETED; \
    if (index < 16) { \
        table->ctrl[table->capacity + index] = CTRL_DELETED; \
    } \
    table->size--; \
    table->deleted++; \
    table->growth_left++; /* Deleted slot can be reused for insert */ \
    return true; \
} \
\
/* Lookup value by key */ \
/* Returns: true if found (value written to *val_out), false if not found */ \
static inline bool swiss_##name##_get(const SwissTable_##name *table, key_t key, val_t *val_out) { \
    size_t index = swiss_##name##_find(table, key); \
    if (index == table->capacity) { \
        return false; \
    } \
    if (val_out) { \
        *val_out = table->vals[index]; \
    } \
    return true; \
} \
\
/* Check if key exists */ \
static inline bool swiss_##name##_contains(const SwissTable_##name *table, key_t key) { \
    return swiss_##name##_find(table, key) != table->capacity; \
} \
\
/* Get size (number of elements) */ \
static inline size_t swiss_##name##_size(const SwissTable_##name *table) { \
    return table->size; \
} \
\
/* Get capacity */ \
static inline size_t swiss_##name##_capacity(const SwissTable_##name *table) { \
    return table->capacity; \
} \
\
/* Iteration support - get next valid index */ \
/* Start with index = 0, call repeatedly until returns capacity */ \
static inline size_t swiss_##name##_iter_next(const SwissTable_##name *table, size_t index) { \
    while (index < table->capacity) { \
        if (is_full(table->ctrl[index])) { \
            return index; \
        } \
        index++; \
    } \
    return table->capacity; \
} \
\
/* Get key at index (for iteration) */ \
static inline key_t swiss_##name##_key_at(const SwissTable_##name *table, size_t index) { \
    return table->keys[index]; \
} \
\
/* Get value at index (for iteration) */ \
static inline val_t swiss_##name##_val_at(const SwissTable_##name *table, size_t index) { \
    return table->vals[index]; \
} \
\
/* Set value at index (for iteration-based updates) */ \
static inline void swiss_##name##_set_val_at(SwissTable_##name *table, size_t index, val_t val) { \
    table->vals[index] = val; \
} \
\
/* Insert key only (no value) - for unique without return_inverse */ \
/* Returns: 1 if newly inserted, 0 if key exists, -1 on memory error */ \
static inline int swiss_##name##_insert_key_only(SwissTable_##name *table, key_t key) { \
    if (table->growth_left == 0) { \
        size_t new_capacity = table->capacity == 0 ? 16 : table->capacity * 2; \
        if (!swiss_##name##_resize(table, new_capacity)) { \
            return -1; \
        } \
    } \
    uint64_t hash = hash_fn(key); \
    int8_t h2 = swiss_h2(hash); \
    size_t index = hash & table->mask; \
    int8_t c0 = table->ctrl[index]; \
    /* Fast path: check if first slot is empty */ \
    if (c0 == CTRL_EMPTY) { \
        table->ctrl[index] = h2; \
        table->keys[index] = key; \
        table->size++; \
        table->growth_left--; \
        if (index < 16) { \
            table->ctrl[table->capacity + index] = h2; \
        } \
        return 1; \
    } \
    /* Fast path: check if first slot has matching key */ \
    if (c0 == h2 && equal_fn(table->keys[index], key)) { \
        return 0; \
    } \
    /* Slow path: SIMD group operations */ \
    ProbeSeq seq; \
    probe_seq_init(&seq, hash, table->mask); \
    bool first_group = true; \
    while (true) { \
        size_t offset = probe_seq_offset(&seq); \
        Group g = group_load(&table->ctrl[offset]); \
        uint16_t match_mask = group_match(g, h2); \
        /* Skip already-checked first slot (index) on first group */ \
        /* Compute relative position within the group */ \
        if (first_group) { \
            size_t rel = (index - offset) & table->mask; \
            if (rel < 16) { \
                match_mask &= (uint16_t)~((uint16_t)1u << (uint16_t)rel); \
            } \
            first_group = false; \
        } \
        while (match_mask != 0) { \
            int bit = countr_zero(match_mask); \
            size_t idx = (offset + bit) & table->mask; \
            if (equal_fn(table->keys[idx], key)) { \
                return 0; \
            } \
            match_mask &= match_mask - 1; \
        } \
        uint16_t empty_mask = group_match_empty_or_deleted(g); \
        if (empty_mask != 0) { \
            int bit = countr_zero(empty_mask); \
            size_t idx = (offset + bit) & table->mask; \
            int8_t old_ctrl = table->ctrl[idx]; \
            table->ctrl[idx] = h2; \
            table->keys[idx] = key; \
            table->size++; \
            if (old_ctrl == CTRL_EMPTY) { \
                table->growth_left--; \
            } else { \
                /* Was CTRL_DELETED - reusing tombstone */ \
                table->deleted--; \
                table->growth_left--; /* Still consuming a slot */ \
            } \
            if (idx < 16) { \
                table->ctrl[table->capacity + idx] = h2; \
            } \
            return 1; \
        } \
        probe_seq_next(&seq); \
    } \
} \
\
/* Batch insert keys only (no values) - for building set from array */ \
/* Returns: 0 on success, -1 on memory error */ \
static inline int swiss_##name##_build_set(SwissTable_##name *table, const key_t *keys, size_t n) { \
    for (size_t i = 0; i < n; i++) { \
        if (swiss_##name##_insert_key_only(table, keys[i]) == -1) { \
            return -1; \
        } \
    } \
    return 0; \
} \
\
/* Insert if key doesn't exist, return existing value if it does */ \
/* Returns: 1 if newly inserted, 0 if key exists (existing value in *val_out) */ \
/* Returns: -1 on memory allocation error */ \
/* On insert, sets value to new_val; on existing, writes existing to *val_out */ \
static inline int swiss_##name##_insert_if_absent(SwissTable_##name *table, key_t key, val_t new_val, val_t *val_out) { \
    if (table->growth_left == 0) { \
        size_t new_capacity = table->capacity == 0 ? 16 : table->capacity * 2; \
        if (!swiss_##name##_resize(table, new_capacity)) { \
            return -1; \
        } \
    } \
    uint64_t hash = hash_fn(key); \
    int8_t h2 = swiss_h2(hash); \
    size_t index = hash & table->mask; \
    int8_t c0 = table->ctrl[index]; \
    /* Fast path: check if first slot is empty (common case for unique values) */ \
    if (c0 == CTRL_EMPTY) { \
        table->ctrl[index] = h2; \
        table->keys[index] = key; \
        table->vals[index] = new_val; \
        table->size++; \
        table->growth_left--; \
        if (index < 16) { \
            table->ctrl[table->capacity + index] = h2; \
        } \
        return 1; \
    } \
    /* Fast path: check if first slot has matching key */ \
    if (c0 == h2 && equal_fn(table->keys[index], key)) { \
        if (val_out) *val_out = table->vals[index]; \
        return 0; \
    } \
    /* SIMD search */ \
    ProbeSeq seq; \
    probe_seq_init(&seq, hash, table->mask); \
    bool first_group = true; \
    while (true) { \
        size_t offset = probe_seq_offset(&seq); \
        Group g = group_load(&table->ctrl[offset]); \
        uint16_t match_mask = group_match(g, h2); \
        /* Skip already-checked first slot (index) on first group */ \
        /* Compute relative position within the group */ \
        if (first_group) { \
            size_t rel = (index - offset) & table->mask; \
            if (rel < 16) { \
                match_mask &= (uint16_t)~((uint16_t)1u << (uint16_t)rel); \
            } \
            first_group = false; \
        } \
        while (match_mask != 0) { \
            int bit = countr_zero(match_mask); \
            size_t idx = (offset + bit) & table->mask; \
            if (equal_fn(table->keys[idx], key)) { \
                if (val_out) *val_out = table->vals[idx]; \
                return 0; \
            } \
            match_mask &= match_mask - 1; \
        } \
        uint16_t empty_mask = group_match_empty_or_deleted(g); \
        if (empty_mask != 0) { \
            int bit = countr_zero(empty_mask); \
            size_t idx = (offset + bit) & table->mask; \
            int8_t old_ctrl = table->ctrl[idx]; \
            table->ctrl[idx] = h2; \
            table->keys[idx] = key; \
            table->vals[idx] = new_val; \
            table->size++; \
            if (old_ctrl == CTRL_EMPTY) { \
                table->growth_left--; \
            } else { \
                /* Was CTRL_DELETED - reusing tombstone */ \
                table->deleted--; \
                table->growth_left--; /* Still consuming a slot */ \
            } \
            if (idx < 16) { \
                table->ctrl[table->capacity + idx] = h2; \
            } \
            return 1; \
        } \
        probe_seq_next(&seq); \
    } \
} \
\
/* Increment value if key exists, insert with value 1 if not */ \
/* Returns: 1 if newly inserted, 0 if incremented, -1 on error */ \
/* Specialized for value_count operations */ \
static inline int swiss_##name##_increment(SwissTable_##name *table, key_t key) { \
    if (table->growth_left == 0) { \
        size_t new_capacity = table->capacity == 0 ? 16 : table->capacity * 2; \
        if (!swiss_##name##_resize(table, new_capacity)) { \
            return -1; \
        } \
    } \
    uint64_t hash = hash_fn(key); \
    int8_t h2 = swiss_h2(hash); \
    size_t index = hash & table->mask; \
    int8_t c0 = table->ctrl[index]; \
    /* Fast path: check if first slot is empty (common case for new keys) */ \
    if (c0 == CTRL_EMPTY) { \
        table->ctrl[index] = h2; \
        table->keys[index] = key; \
        table->vals[index] = 1; \
        table->size++; \
        table->growth_left--; \
        if (index < 16) { \
            table->ctrl[table->capacity + index] = h2; \
        } \
        return 1; \
    } \
    /* Fast path: check if first slot has matching key */ \
    if (c0 == h2 && equal_fn(table->keys[index], key)) { \
        table->vals[index]++; \
        return 0; \
    } \
    /* Slow path: need to check for existing key or find another empty slot */ \
    ProbeSeq seq; \
    probe_seq_init(&seq, hash, table->mask); \
    bool first_group = true; \
    while (true) { \
        size_t offset = probe_seq_offset(&seq); \
        Group g = group_load(&table->ctrl[offset]); \
        uint16_t match_mask = group_match(g, h2); \
        /* Skip already-checked first slot (index) on first group */ \
        /* Compute relative position within the group */ \
        if (first_group) { \
            size_t rel = (index - offset) & table->mask; \
            if (rel < 16) { \
                match_mask &= (uint16_t)~((uint16_t)1u << (uint16_t)rel); \
            } \
            first_group = false; \
        } \
        while (match_mask != 0) { \
            int bit = countr_zero(match_mask); \
            size_t idx = (offset + bit) & table->mask; \
            if (equal_fn(table->keys[idx], key)) { \
                table->vals[idx]++; \
                return 0; \
            } \
            match_mask &= match_mask - 1; \
        } \
        uint16_t empty_mask = group_match_empty_or_deleted(g); \
        if (empty_mask != 0) { \
            int bit = countr_zero(empty_mask); \
            size_t idx = (offset + bit) & table->mask; \
            int8_t old_ctrl = table->ctrl[idx]; \
            table->ctrl[idx] = h2; \
            table->keys[idx] = key; \
            table->vals[idx] = 1; \
            table->size++; \
            if (old_ctrl == CTRL_EMPTY) { \
                table->growth_left--; \
            } else { \
                /* Was CTRL_DELETED - reusing tombstone */ \
                table->deleted--; \
                table->growth_left--; /* Still consuming a slot */ \
            } \
            if (idx < 16) { \
                table->ctrl[table->capacity + idx] = h2; \
            } \
            return 1; \
        } \
        probe_seq_next(&seq); \
    } \
} \
\
/* Batch insert for map_locations - inserts array of (key, index) pairs */ \
/* Returns: 0 on success, -1 on memory error */ \
/* */ \
/* Optimization: Batch hashing with reserve-upfront pattern. */ \
/* Each key maps to its array index as the value. */ \
static inline int swiss_##name##_map_locations(SwissTable_##name *table, const key_t *keys, size_t n) { \
    /* Reserve capacity upfront to avoid resize during batch processing. */ \
    if (!swiss_##name##_reserve(table, n)) { \
        return -1; \
    } \
    /* Batch size 256 */ \
    uint64_t hashes[256]; \
    int8_t h2s[256]; \
    size_t indices_arr[256]; \
    \
    size_t i = 0; \
    while (i + 256 <= n) { \
        /* Phase 1: Batch hash computation */ \
        for (size_t b = 0; b < 256; b++) { \
            hashes[b] = hash_fn(keys[i + b]); \
            h2s[b] = swiss_h2(hashes[b]); \
            indices_arr[b] = hashes[b] & table->mask; \
        } \
        /* Phase 2: Process each key */ \
        for (size_t b = 0; b < 256; b++) { \
            key_t key = keys[i + b]; \
            val_t val = (val_t)(i + b); \
            uint64_t hash = hashes[b]; \
            int8_t h2 = h2s[b]; \
            size_t index = indices_arr[b]; \
            int8_t c0 = table->ctrl[index]; \
            /* Fast path: first slot is empty */ \
            if (c0 == CTRL_EMPTY) { \
                table->ctrl[index] = h2; \
                table->keys[index] = key; \
                table->vals[index] = val; \
                table->size++; \
                table->growth_left--; \
                if (index < 16) { \
                    table->ctrl[table->capacity + index] = h2; \
                } \
                continue; \
            } \
            /* Fast path: first slot matches - update value */ \
            if (c0 == h2 && equal_fn(table->keys[index], key)) { \
                table->vals[index] = val; \
                continue; \
            } \
            /* Scalar probes for slots 1-2 */ \
            bool found_in_scalar = false; \
            for (size_t step = 1; step <= 2; step++) { \
                size_t j = (index + step) & table->mask; \
                int8_t c = table->ctrl[j]; \
                if (c == CTRL_EMPTY) { \
                    table->ctrl[j] = h2; \
                    table->keys[j] = key; \
                    table->vals[j] = val; \
                    table->size++; \
                    table->growth_left--; \
                    if (j < 16) { \
                        table->ctrl[table->capacity + j] = h2; \
                    } \
                    found_in_scalar = true; \
                    break; \
                } \
                if (c == h2 && equal_fn(table->keys[j], key)) { \
                    table->vals[j] = val; \
                    found_in_scalar = true; \
                    break; \
                } \
            } \
            if (found_in_scalar) continue; \
            /* SIMD search */ \
            ProbeSeq seq; \
            probe_seq_init(&seq, hash, table->mask); \
            bool first_group = true; \
            bool done = false; \
            while (!done) { \
                size_t offset = probe_seq_offset(&seq); \
                Group g = group_load(&table->ctrl[offset]); \
                uint16_t match_mask = group_match(g, h2); \
                if (first_group) { \
                    for (size_t step = 0; step <= 2; step++) { \
                        size_t slot = (index + step) & table->mask; \
                        size_t rel = (slot - offset) & table->mask; \
                        if (rel < 16) { \
                            match_mask &= (uint16_t)~((uint16_t)1u << (uint16_t)rel); \
                        } \
                    } \
                    first_group = false; \
                } \
                while (match_mask != 0) { \
                    int bit = countr_zero(match_mask); \
                    size_t idx = (offset + bit) & table->mask; \
                    if (equal_fn(table->keys[idx], key)) { \
                        table->vals[idx] = val; \
                        done = true; \
                        break; \
                    } \
                    match_mask &= match_mask - 1; \
                } \
                if (done) break; \
                uint16_t empty_mask = group_match_empty(g); \
                if (empty_mask != 0) { \
                    int bit = countr_zero(empty_mask); \
                    size_t idx = (offset + bit) & table->mask; \
                    table->ctrl[idx] = h2; \
                    table->keys[idx] = key; \
                    table->vals[idx] = val; \
                    table->size++; \
                    table->growth_left--; \
                    if (idx < 16) { \
                        table->ctrl[table->capacity + idx] = h2; \
                    } \
                    break; \
                } \
                probe_seq_next(&seq); \
            } \
        } \
        i += 256; \
    } \
    /* Handle remaining keys */ \
    for (; i < n; i++) { \
        if (swiss_##name##_insert(table, keys[i], (val_t)i) == -1) { \
            return -1; \
        } \
    } \
    return 0; \
} \
\
/* Batch lookup - looks up array of keys, writes results to locs array */ \
/* locs[i] = index if found, -1 if not found */ \
static inline void swiss_##name##_lookup_batch(const SwissTable_##name *table, const key_t *keys, size_t n, int64_t *locs) { \
    for (size_t i = 0; i < n; i++) { \
        size_t idx = swiss_##name##_find(table, keys[i]); \
        if (idx == table->capacity) { \
            locs[i] = -1; \
        } else { \
            locs[i] = (int64_t)table->vals[idx]; \
        } \
    } \
} \
\
/* Batch contains - checks if keys exist, writes results to result array */ \
static inline void swiss_##name##_contains_batch(const SwissTable_##name *table, const key_t *keys, size_t n, uint8_t *result) { \
    for (size_t i = 0; i < n; i++) { \
        result[i] = (swiss_##name##_find(table, keys[i]) != table->capacity) ? 1 : 0; \
    } \
} \
\
/* Batch unique - extracts unique values from input array */ \
/* uniques_out must have space for n elements */ \
/* Returns: number of unique values on success, -1 on memory error */ \
static inline int64_t swiss_##name##_unique(SwissTable_##name *table, const key_t *keys, size_t n, key_t *uniques_out) { \
    size_t count = 0; \
    for (size_t i = 0; i < n; i++) { \
        key_t key = keys[i]; \
        val_t idx; \
        int ret = swiss_##name##_insert_if_absent(table, key, (val_t)count, &idx); \
        if (ret == 1) { \
            uniques_out[count] = key; \
            count++; \
        } else if (ret == -1) { \
            return -1; \
        } \
    } \
    return (int64_t)count; \
} \
\
/* Batch unique with inverse - extracts unique values and labels */ \
/* uniques_out must have space for n elements */ \
/* labels_out must have space for n elements */ \
/* Returns: number of unique values on success, -1 on memory error */ \
/* Optimized: uses find-first approach like khash for duplicate-heavy workloads */ \
static inline int64_t swiss_##name##_unique_with_inverse(SwissTable_##name *table, const key_t *keys, size_t n, key_t *uniques_out, int64_t *labels_out) { \
    size_t count = 0; \
    for (size_t i = 0; i < n; i++) { \
        key_t key = keys[i]; \
        /* First, try to find the key (fast path for duplicates) */ \
        size_t found_idx = swiss_##name##_find(table, key); \
        if (found_idx != table->capacity) { \
            /* Key exists - this is the common case for duplicate-heavy data */ \
            labels_out[i] = (int64_t)table->vals[found_idx]; \
        } else { \
            /* Key doesn't exist - insert it */ \
            int ret = swiss_##name##_insert(table, key, (val_t)count); \
            if (ret == -1) { \
                return -1; \
            } \
            uniques_out[count] = key; \
            labels_out[i] = (int64_t)count; \
            count++; \
        } \
    } \
    return (int64_t)count; \
} \
\
/* Batch factorize - same as unique_with_inverse but with ssize_t labels */ \
/* This matches pandas' intp_t (ssize_t) label type directly */ \
/* Returns: number of unique values on success, -1 on memory error */ \
/* */ \
/* Optimization: Batch hashing with reserve-upfront pattern (same as value_count_batch). */ \
/* Separates hash computation from probing for better CPU pipelining. */ \
static inline ssize_t swiss_##name##_factorize_batch(SwissTable_##name *table, const key_t *keys, size_t n, key_t *uniques_out, ssize_t *labels_out) { \
    size_t count = 0; \
    /* Reserve capacity upfront to avoid resize during batch processing. */ \
    if (!swiss_##name##_reserve(table, n)) { \
        return -1; \
    } \
    /* Batch size 256: good balance between stack usage and pipelining benefits. */ \
    uint64_t hashes[256]; \
    int8_t h2s[256]; \
    size_t indices_arr[256]; \
    \
    size_t i = 0; \
    /* Process in batches of 256 */ \
    while (i + 256 <= n) { \
        /* Phase 1: Batch hash computation */ \
        for (size_t b = 0; b < 256; b++) { \
            hashes[b] = hash_fn(keys[i + b]); \
            h2s[b] = swiss_h2(hashes[b]); \
            indices_arr[b] = hashes[b] & table->mask; \
        } \
        /* Phase 2: Process each key with precomputed hash */ \
        for (size_t b = 0; b < 256; b++) { \
            key_t key = keys[i + b]; \
            uint64_t hash = hashes[b]; \
            int8_t h2 = h2s[b]; \
            size_t index = indices_arr[b]; \
            int8_t c0 = table->ctrl[index]; \
            /* Fast path: first slot is empty - insert new key */ \
            if (c0 == CTRL_EMPTY) { \
                table->ctrl[index] = h2; \
                table->keys[index] = key; \
                table->vals[index] = (val_t)count; \
                table->size++; \
                table->growth_left--; \
                if (index < 16) { \
                    table->ctrl[table->capacity + index] = h2; \
                } \
                uniques_out[count] = key; \
                labels_out[i + b] = (ssize_t)count; \
                count++; \
                continue; \
            } \
            /* Fast path: first slot matches - return existing label */ \
            if (c0 == h2 && equal_fn(table->keys[index], key)) { \
                labels_out[i + b] = (ssize_t)table->vals[index]; \
                continue; \
            } \
            /* Scalar probes for slots 1-2 */ \
            bool found_in_scalar = false; \
            for (size_t step = 1; step <= 2; step++) { \
                size_t j = (index + step) & table->mask; \
                int8_t c = table->ctrl[j]; \
                if (c == CTRL_EMPTY) { \
                    table->ctrl[j] = h2; \
                    table->keys[j] = key; \
                    table->vals[j] = (val_t)count; \
                    table->size++; \
                    table->growth_left--; \
                    if (j < 16) { \
                        table->ctrl[table->capacity + j] = h2; \
                    } \
                    uniques_out[count] = key; \
                    labels_out[i + b] = (ssize_t)count; \
                    count++; \
                    found_in_scalar = true; \
                    break; \
                } \
                if (c == h2 && equal_fn(table->keys[j], key)) { \
                    labels_out[i + b] = (ssize_t)table->vals[j]; \
                    found_in_scalar = true; \
                    break; \
                } \
            } \
            if (found_in_scalar) continue; \
            /* SIMD search */ \
            ProbeSeq seq; \
            probe_seq_init(&seq, hash, table->mask); \
            bool first_group = true; \
            bool done = false; \
            while (!done) { \
                size_t offset = probe_seq_offset(&seq); \
                Group g = group_load(&table->ctrl[offset]); \
                uint16_t match_mask = group_match(g, h2); \
                /* Skip already-checked slots on first group */ \
                if (first_group) { \
                    for (size_t step = 0; step <= 2; step++) { \
                        size_t slot = (index + step) & table->mask; \
                        size_t rel = (slot - offset) & table->mask; \
                        if (rel < 16) { \
                            match_mask &= (uint16_t)~((uint16_t)1u << (uint16_t)rel); \
                        } \
                    } \
                    first_group = false; \
                } \
                while (match_mask != 0) { \
                    int bit = countr_zero(match_mask); \
                    size_t idx = (offset + bit) & table->mask; \
                    if (equal_fn(table->keys[idx], key)) { \
                        labels_out[i + b] = (ssize_t)table->vals[idx]; \
                        done = true; \
                        break; \
                    } \
                    match_mask &= match_mask - 1; \
                } \
                if (done) break; \
                uint16_t empty_mask = group_match_empty(g); \
                if (empty_mask != 0) { \
                    int bit = countr_zero(empty_mask); \
                    size_t idx = (offset + bit) & table->mask; \
                    table->ctrl[idx] = h2; \
                    table->keys[idx] = key; \
                    table->vals[idx] = (val_t)count; \
                    table->size++; \
                    table->growth_left--; \
                    if (idx < 16) { \
                        table->ctrl[table->capacity + idx] = h2; \
                    } \
                    uniques_out[count] = key; \
                    labels_out[i + b] = (ssize_t)count; \
                    count++; \
                    break; \
                } \
                probe_seq_next(&seq); \
            } \
        } \
        i += 256; \
    } \
    /* Handle remaining keys */ \
    for (; i < n; i++) { \
        key_t key = keys[i]; \
        uint64_t hash = hash_fn(key); \
        int8_t h2 = swiss_h2(hash); \
        size_t index = hash & table->mask; \
        int8_t c0 = table->ctrl[index]; \
        /* Fast path: first slot is empty */ \
        if (c0 == CTRL_EMPTY) { \
            table->ctrl[index] = h2; \
            table->keys[index] = key; \
            table->vals[index] = (val_t)count; \
            table->size++; \
            table->growth_left--; \
            if (index < 16) { \
                table->ctrl[table->capacity + index] = h2; \
            } \
            uniques_out[count] = key; \
            labels_out[i] = (ssize_t)count; \
            count++; \
            continue; \
        } \
        /* Fast path: first slot matches */ \
        if (c0 == h2 && equal_fn(table->keys[index], key)) { \
            labels_out[i] = (ssize_t)table->vals[index]; \
            continue; \
        } \
        /* Scalar probes for slots 1-2 */ \
        bool found_in_scalar = false; \
        for (size_t step = 1; step <= 2; step++) { \
            size_t j = (index + step) & table->mask; \
            int8_t c = table->ctrl[j]; \
            if (c == CTRL_EMPTY) { \
                table->ctrl[j] = h2; \
                table->keys[j] = key; \
                table->vals[j] = (val_t)count; \
                table->size++; \
                table->growth_left--; \
                if (j < 16) { \
                    table->ctrl[table->capacity + j] = h2; \
                } \
                uniques_out[count] = key; \
                labels_out[i] = (ssize_t)count; \
                count++; \
                found_in_scalar = true; \
                break; \
            } \
            if (c == h2 && equal_fn(table->keys[j], key)) { \
                labels_out[i] = (ssize_t)table->vals[j]; \
                found_in_scalar = true; \
                break; \
            } \
        } \
        if (found_in_scalar) continue; \
        /* SIMD search */ \
        ProbeSeq seq; \
        probe_seq_init(&seq, hash, table->mask); \
        bool first_group = true; \
        bool done = false; \
        while (!done) { \
            size_t offset = probe_seq_offset(&seq); \
            Group g = group_load(&table->ctrl[offset]); \
            uint16_t match_mask = group_match(g, h2); \
            if (first_group) { \
                for (size_t step = 0; step <= 2; step++) { \
                    size_t slot = (index + step) & table->mask; \
                    size_t rel = (slot - offset) & table->mask; \
                    if (rel < 16) { \
                        match_mask &= (uint16_t)~((uint16_t)1u << (uint16_t)rel); \
                    } \
                } \
                first_group = false; \
            } \
            while (match_mask != 0) { \
                int bit = countr_zero(match_mask); \
                size_t idx = (offset + bit) & table->mask; \
                if (equal_fn(table->keys[idx], key)) { \
                    labels_out[i] = (ssize_t)table->vals[idx]; \
                    done = true; \
                    break; \
                } \
                match_mask &= match_mask - 1; \
            } \
            if (done) break; \
            uint16_t empty_mask = group_match_empty(g); \
            if (empty_mask != 0) { \
                int bit = countr_zero(empty_mask); \
                size_t idx = (offset + bit) & table->mask; \
                table->ctrl[idx] = h2; \
                table->keys[idx] = key; \
                table->vals[idx] = (val_t)count; \
                table->size++; \
                table->growth_left--; \
                if (idx < 16) { \
                    table->ctrl[table->capacity + idx] = h2; \
                } \
                uniques_out[count] = key; \
                labels_out[i] = (ssize_t)count; \
                count++; \
                break; \
            } \
            probe_seq_next(&seq); \
        } \
    } \
    return (ssize_t)count; \
} \
\
/* Batch value_count - counts occurrences, returns keys and counts */ \
/* keys_out must have space for n elements */ \
/* counts_out must have space for n elements */ \
/* indices_out must have space for n elements (stores table indices for direct access) */ \
/* Returns: number of unique values on success, -1 on memory error */ \
/* */ \
/* IMPORTANT: This function assumes the table is freshly initialized (no tombstones). */ \
/* The table must be created with init_with_capacity() and used ONLY for this call. */ \
/* This invariant allows optimizations: */ \
/*   - group_match_empty instead of group_match_empty_or_deleted (saves NEON movemask) */ \
/*   - Scalar probes can check only CTRL_EMPTY (not CTRL_DELETED) for inserts */ \
/* Callers: see swisstable.pyx value_count_int64() which creates a fresh table. */ \
/* */ \
/* Optimization: Batch hashing - compute hashes for BATCH_SIZE keys at once. */ \
/* This improves CPU pipelining by separating hash computation from probing. */ \
/* Inspired by Apache Arrow's vectorized Swiss Table implementation. */ \
static inline ssize_t swiss_##name##_value_count_batch(SwissTable_##name *table, const key_t *keys, size_t n, key_t *keys_out, size_t *indices_out) { \
    size_t count = 0; \
    /* Reserve capacity upfront to avoid resize during batch processing. */ \
    /* Worst case: all n keys are unique. Reserve ensures no resize needed. */ \
    /* This is critical because batch hashing precomputes indices using table->mask, */ \
    /* and a resize would invalidate all remaining precomputed indices in the batch. */ \
    if (!swiss_##name##_reserve(table, n)) { \
        return -1; \
    } \
    /* Batch size 256: good balance between stack usage (~6KB) and pipelining benefits. */ \
    /* Arrow uses 1024, but 256 is sufficient for most workloads. */ \
    uint64_t hashes[256]; \
    int8_t h2s[256]; \
    size_t indices_arr[256]; \
    \
    size_t i = 0; \
    /* Process in batches of 256 */ \
    while (i + 256 <= n) { \
        /* Phase 1: Batch hash computation - better CPU pipelining */ \
        for (size_t b = 0; b < 256; b++) { \
            hashes[b] = hash_fn(keys[i + b]); \
            h2s[b] = swiss_h2(hashes[b]); \
            indices_arr[b] = hashes[b] & table->mask; \
        } \
        /* Phase 2: Process each key with precomputed hash */ \
        for (size_t b = 0; b < 256; b++) { \
            key_t key = keys[i + b]; \
            uint64_t hash = hashes[b]; \
            int8_t h2 = h2s[b]; \
            size_t index = indices_arr[b]; \
            int8_t c0 = table->ctrl[index]; \
            /* Fast path: first slot is empty */ \
            if (c0 == CTRL_EMPTY) { \
                table->ctrl[index] = h2; \
                table->keys[index] = key; \
                table->vals[index] = 1; \
                table->size++; \
                table->growth_left--; \
                if (index < 16) { \
                    table->ctrl[table->capacity + index] = h2; \
                } \
                keys_out[count] = key; \
                indices_out[count] = index; \
                count++; \
                continue; \
            } \
            /* Fast path: first slot matches */ \
            if (c0 == h2 && equal_fn(table->keys[index], key)) { \
                table->vals[index]++; \
                continue; \
            } \
            /* Scalar probes for slots 1-2: avoids SIMD for clustered duplicates */ \
            bool found_in_scalar = false; \
            for (size_t step = 1; step <= 2; step++) { \
                size_t j = (index + step) & table->mask; \
                int8_t c = table->ctrl[j]; \
                if (c == CTRL_EMPTY) { \
                    table->ctrl[j] = h2; \
                    table->keys[j] = key; \
                    table->vals[j] = 1; \
                    table->size++; \
                    table->growth_left--; \
                    if (j < 16) { \
                        table->ctrl[table->capacity + j] = h2; \
                    } \
                    keys_out[count] = key; \
                    indices_out[count] = j; \
                    count++; \
                    found_in_scalar = true; \
                    break; \
                } \
                if (c == h2 && equal_fn(table->keys[j], key)) { \
                    table->vals[j]++; \
                    found_in_scalar = true; \
                    break; \
                } \
            } \
            if (found_in_scalar) continue; \
            /* SIMD search */ \
            ProbeSeq seq; \
            probe_seq_init(&seq, hash, table->mask); \
            bool first_group = true; \
            bool done = false; \
            while (!done) { \
                size_t offset = probe_seq_offset(&seq); \
                Group g = group_load(&table->ctrl[offset]); \
                uint16_t match_mask = group_match(g, h2); \
                /* Skip already-checked slots (index, index+1, index+2) on first group */ \
                if (first_group) { \
                    for (size_t step = 0; step <= 2; step++) { \
                        size_t slot = (index + step) & table->mask; \
                        size_t rel = (slot - offset) & table->mask; \
                        if (rel < 16) { \
                            match_mask &= (uint16_t)~((uint16_t)1u << (uint16_t)rel); \
                        } \
                    } \
                    first_group = false; \
                } \
                while (match_mask != 0) { \
                    int bit = countr_zero(match_mask); \
                    size_t idx = (offset + bit) & table->mask; \
                    if (equal_fn(table->keys[idx], key)) { \
                        table->vals[idx]++; \
                        done = true; \
                        break; \
                    } \
                    match_mask &= match_mask - 1; \
                } \
                if (done) break; \
                /* Use group_match_empty (not empty_or_deleted) since we never delete */ \
                uint16_t empty_mask = group_match_empty(g); \
                if (empty_mask != 0) { \
                    int bit = countr_zero(empty_mask); \
                    size_t idx = (offset + bit) & table->mask; \
                    table->ctrl[idx] = h2; \
                    table->keys[idx] = key; \
                    table->vals[idx] = 1; \
                    table->size++; \
                    table->growth_left--; \
                    if (idx < 16) { \
                        table->ctrl[table->capacity + idx] = h2; \
                    } \
                    keys_out[count] = key; \
                    indices_out[count] = idx; \
                    count++; \
                    break; \
                } \
                probe_seq_next(&seq); \
            } \
        } \
        i += 256; \
    } \
    /* Handle remaining keys (less than 256) - uses same logic as batch loop */ \
    for (; i < n; i++) { \
        key_t key = keys[i]; \
        uint64_t hash = hash_fn(key); \
        int8_t h2 = swiss_h2(hash); \
        size_t index = hash & table->mask; \
        int8_t c0 = table->ctrl[index]; \
        /* Fast path: first slot is empty */ \
        if (c0 == CTRL_EMPTY) { \
            table->ctrl[index] = h2; \
            table->keys[index] = key; \
            table->vals[index] = 1; \
            table->size++; \
            table->growth_left--; \
            if (index < 16) { \
                table->ctrl[table->capacity + index] = h2; \
            } \
            keys_out[count] = key; \
            indices_out[count] = index; \
            count++; \
            continue; \
        } \
        /* Fast path: first slot matches */ \
        if (c0 == h2 && equal_fn(table->keys[index], key)) { \
            table->vals[index]++; \
            continue; \
        } \
        /* Scalar probes for slots 1-2 */ \
        bool found_in_scalar = false; \
        for (size_t step = 1; step <= 2; step++) { \
            size_t j = (index + step) & table->mask; \
            int8_t c = table->ctrl[j]; \
            if (c == CTRL_EMPTY) { \
                table->ctrl[j] = h2; \
                table->keys[j] = key; \
                table->vals[j] = 1; \
                table->size++; \
                table->growth_left--; \
                if (j < 16) { \
                    table->ctrl[table->capacity + j] = h2; \
                } \
                keys_out[count] = key; \
                indices_out[count] = j; \
                count++; \
                found_in_scalar = true; \
                break; \
            } \
            if (c == h2 && equal_fn(table->keys[j], key)) { \
                table->vals[j]++; \
                found_in_scalar = true; \
                break; \
            } \
        } \
        if (found_in_scalar) continue; \
        /* SIMD search */ \
        ProbeSeq seq; \
        probe_seq_init(&seq, hash, table->mask); \
        bool first_group = true; \
        bool done = false; \
        while (!done) { \
            size_t offset = probe_seq_offset(&seq); \
            Group g = group_load(&table->ctrl[offset]); \
            uint16_t match_mask = group_match(g, h2); \
            /* Skip already-checked slots (index, index+1, index+2) on first group */ \
            if (first_group) { \
                for (size_t step = 0; step <= 2; step++) { \
                    size_t slot = (index + step) & table->mask; \
                    size_t rel = (slot - offset) & table->mask; \
                    if (rel < 16) { \
                        match_mask &= (uint16_t)~((uint16_t)1u << (uint16_t)rel); \
                    } \
                } \
                first_group = false; \
            } \
            while (match_mask != 0) { \
                int bit = countr_zero(match_mask); \
                size_t idx = (offset + bit) & table->mask; \
                if (equal_fn(table->keys[idx], key)) { \
                    table->vals[idx]++; \
                    done = true; \
                    break; \
                } \
                match_mask &= match_mask - 1; \
            } \
            if (done) break; \
            uint16_t empty_mask = group_match_empty(g); \
            if (empty_mask != 0) { \
                int bit = countr_zero(empty_mask); \
                size_t idx = (offset + bit) & table->mask; \
                table->ctrl[idx] = h2; \
                table->keys[idx] = key; \
                table->vals[idx] = 1; \
                table->size++; \
                table->growth_left--; \
                if (idx < 16) { \
                    table->ctrl[table->capacity + idx] = h2; \
                } \
                keys_out[count] = key; \
                indices_out[count] = idx; \
                count++; \
                break; \
            } \
            probe_seq_next(&seq); \
        } \
    } \
    return (ssize_t)count; \
}

#endif  /* PANDAS_SWISSTABLE_GENERIC_H */
