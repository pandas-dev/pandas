// CKHashTable - A simple hash table
// Copyright 2010-2020 Evan Miller (see LICENSE)

#ifndef PANDAS__LIBS_SRC_LIBRDATA_CKHASHTABLE_H_
#define PANDAS__LIBS_SRC_LIBRDATA_CKHASHTABLE_H_

#include <sys/types.h>
#include <stdint.h>

typedef struct ck_hash_entry_s {
    off_t       key_offset;
    size_t      key_length;
    const void *value;
} ck_hash_entry_t;

typedef struct ck_hash_table_s {
    size_t capacity;
    size_t count;
    ck_hash_entry_t *entries;
    char   *keys;
    size_t  keys_used;
    size_t  keys_capacity;
} ck_hash_table_t;

int ck_str_hash_insert(
    const char *key, const void *value, ck_hash_table_t *table
);
const void *ck_str_hash_lookup(const char *key, ck_hash_table_t *table);

int ck_str_n_hash_insert(
    const char *key, size_t keylen, const void *value, ck_hash_table_t *table
);
const void *ck_str_n_hash_lookup(
    const char *key, size_t keylen, ck_hash_table_t *table
);

int ck_float_hash_insert(
    float key, const void *value, ck_hash_table_t *table
);
const void *ck_float_hash_lookup(float key, ck_hash_table_t *table);

int ck_double_hash_insert(
    double key, const void *value, ck_hash_table_t *table
);
const void *ck_double_hash_lookup(double key, ck_hash_table_t *table);

ck_hash_table_t *ck_hash_table_init(
    size_t num_entries, size_t mean_key_length
);
void ck_hash_table_wipe(ck_hash_table_t *table);
int ck_hash_table_grow(ck_hash_table_t *table);
void ck_hash_table_free(ck_hash_table_t *table);
uint64_t ck_hash_str(const char *str, size_t keylen);

#endif  // PANDAS__LIBS_SRC_LIBRDATA_CKHASHTABLE_H_
