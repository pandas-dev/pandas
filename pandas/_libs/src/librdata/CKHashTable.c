// CKHashTable - A simple hash table
// Copyright 2010-2020 Evan Miller (see LICENSE)

#include "CKHashTable.h"

/*
 SipHash reference C implementation

 Copyright (c) 2012 Jean-Philippe Aumasson <jeanphilippe.aumasson@gmail.com>
 Copyright (c) 2012 Daniel J. Bernstein <djb@cr.yp.to>

 To the extent possible under law, the author(s) have dedicated all copyright
 and related and neighboring rights to this software to the public domain
 worldwide. This software is distributed without any warranty.

 You should have received a copy of the CC0 Public Domain Dedication along with
 this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
typedef uint64_t u64;
typedef uint32_t u32;
typedef uint8_t u8;


#define ROTL(x, b) (u64)( ((x) << (b)) | ( (x) >> (64 - (b))) )

#define U32TO8_LE(p, v)         \
(p)[0] = (u8)((v)      ); (p)[1] = (u8)((v) >>  8); \
(p)[2] = (u8)((v) >> 16); (p)[3] = (u8)((v) >> 24);

#define U64TO8_LE(p, v)         \
U32TO8_LE((p),     (u32)((v)      ));   \
U32TO8_LE((p) + 4, (u32)((v) >> 32));

#define U8TO64_LE(p) \
(((u64)((p)[0])) | \
((u64)((p)[1]) <<  8) | \
((u64)((p)[2]) << 16) | \
((u64)((p)[3]) << 24) | \
((u64)((p)[4]) << 32) | \
((u64)((p)[5]) << 40) | \
((u64)((p)[6]) << 48) | \
((u64)((p)[7]) << 56))

#define SIPROUND            \
do {              \
v0 += v1; v1=ROTL(v1, 13); v1 ^= v0; v0=ROTL(v0, 32); \
v2 += v3; v3=ROTL(v3, 16); v3 ^= v2;     \
v0 += v3; v3=ROTL(v3, 21); v3 ^= v0;     \
v2 += v1; v1=ROTL(v1, 17); v1 ^= v2; v2=ROTL(v2, 32); \
} while (0)

/* SipHash-1-2 */
static int siphash(
    unsigned char *out,
    const unsigned char *in,
    unsigned long long inlen,
    const unsigned char *k) {
    /* "somepseudorandomlygeneratedbytes" */
    u64 v0 = 0x736f6d6570736575ULL;
    u64 v1 = 0x646f72616e646f6dULL;
    u64 v2 = 0x6c7967656e657261ULL;
    u64 v3 = 0x7465646279746573ULL;
    u64 b;
    u64 k0 = U8TO64_LE(k);
    u64 k1 = U8TO64_LE(k + 8);
    u64 m;
    const u8 *end = in + inlen - ( inlen % sizeof( u64 ) );
    const int left = inlen & 7;
    b = ((u64)inlen) << 56;
    v3 ^= k1;
    v2 ^= k0;
    v1 ^= k1;
    v0 ^= k0;

    for ( ; in != end; in += 8 ) {
        m = U8TO64_LE(in);

        v3 ^= m;

        SIPROUND;

        v0 ^= m;
    }

    switch ( left ) {
        case 7: b |= ((u64)in[ 6])  << 48;

        case 6: b |= ((u64)in[ 5])  << 40;

        case 5: b |= ((u64)in[ 4])  << 32;

        case 4: b |= ((u64)in[ 3])  << 24;

        case 3: b |= ((u64)in[ 2])  << 16;

        case 2: b |= ((u64)in[ 1])  <<  8;

        case 1: b |= ((u64)in[ 0]); break;

        case 0: break;
    }

    v3 ^= b;

    SIPROUND;

    v0 ^= b;
    v2 ^= 0xff;

    SIPROUND;
    SIPROUND;

    b = v0 ^ v1 ^ v2  ^ v3;
    U64TO8_LE(out, b);
    return 0;
}

inline uint64_t ck_hash_str(const char *str, size_t keylen) {
    uint64_t hash;
    unsigned char k[16] = { 0 };
    siphash((unsigned char *)&hash, (const unsigned char *)str, keylen, k);
    return hash;
}

const void *ck_float_hash_lookup(float key, ck_hash_table_t *table) {
    return ck_str_n_hash_lookup((const char *)&key, sizeof(float), table);
}

int ck_float_hash_insert(
    float key,
    const void *value,
    ck_hash_table_t *table
) {
    return ck_str_n_hash_insert(
         (const char *)&key,
         sizeof(float),
         value,
         table);
}

const void *ck_double_hash_lookup(double key, ck_hash_table_t *table) {
    return ck_str_n_hash_lookup((const char *)&key, sizeof(double), table);
}

int ck_double_hash_insert(
    double key,
    const void *value,
    ck_hash_table_t *table
) {
    return ck_str_n_hash_insert(
        (const char *)&key,
        sizeof(double),
        value,
        table);
}

const void *ck_str_hash_lookup(const char *key, ck_hash_table_t *table) {
    size_t keylen = strlen(key);
    return ck_str_n_hash_lookup(key, keylen, table);
}

const void *ck_str_n_hash_lookup(
    const char *key,
    size_t keylen,
    ck_hash_table_t *table
) {
    if (table->count == 0)
        return NULL;

    if (keylen == 0)
        return NULL;

    uint64_t hash_key = ck_hash_str(key, keylen);
    hash_key %= table->capacity;
    uint64_t end = hash_key;
    do {
        char *this_key = &table->keys[table->entries[hash_key].key_offset];
        size_t this_keylen = table->entries[hash_key].key_length;
        if (this_keylen == 0)
            return NULL;
        if (this_keylen == keylen && memcmp(this_key, key, keylen) == 0) {
            return table->entries[hash_key].value;
        }
        hash_key++;
        hash_key %= table->capacity;
    } while (hash_key != end);
    return NULL;
}

int ck_str_hash_insert(
    const char *key,
    const void *value,
    ck_hash_table_t *table
) {
    size_t keylen = strlen(key);
    return ck_str_n_hash_insert(key, keylen, value, table);
}

static int ck_hash_insert_nocopy(
    off_t key_offset,
    size_t keylen,
    uint64_t hash_key,
    const void *value,
    ck_hash_table_t *table
) {
    if (table->capacity == 0)
        return 0;

    hash_key %= table->capacity;
    uint64_t end = (hash_key + table->capacity - 1) % table->capacity;
    while (hash_key != end) {
        ck_hash_entry_t *entry = &table->entries[hash_key];
        if (table->entries[hash_key].key_length == 0) {
            table->count++;
            entry->key_offset = key_offset;
            entry->key_length = keylen;
            entry->value = value;
            return 1;
        } else if (entry->key_length == keylen &&
                   entry->key_offset == key_offset) {
            entry->value = value;
            return 1;
        }
        hash_key++;
        hash_key %= table->capacity;
    }
    return 0;
}

int ck_str_n_hash_insert(
    const char *key,
    size_t keylen,
    const void *value,
    ck_hash_table_t *table
) {
    if (table->capacity == 0)
        return 0;

    if (keylen == 0)
        return 0;

    if (table->count >= 0.75 * table->capacity) {
        if (ck_hash_table_grow(table) == -1) {
            return 0;
        }
    }

    uint64_t hash_key = ck_hash_str(key, keylen);
    hash_key %= table->capacity;
    uint64_t end = hash_key;
    do {
        ck_hash_entry_t *entry = &table->entries[hash_key];
        char *this_key = &table->keys[entry->key_offset];
        if (entry->key_length == 0) {
            table->count++;
            while (table->keys_used + keylen > table->keys_capacity) {
                table->keys_capacity *= 2;
                table->keys = realloc(table->keys, table->keys_capacity);
            }
            memcpy(table->keys + table->keys_used, key, keylen);
            entry->key_offset = table->keys_used;
            entry->key_length = keylen;
            table->keys_used += keylen;
            entry->value = value;
            return 1;
        } else if (entry->key_length == keylen &&
                   memcmp(this_key, key, keylen) == 0) {
            table->entries[hash_key].value = value;
            return 1;
        }
        hash_key++;
        hash_key %= table->capacity;
    } while (hash_key != end);
    return 0;
}

ck_hash_table_t *ck_hash_table_init(
    size_t num_entries,
    size_t mean_key_length
) {
    ck_hash_table_t *table;
    if ((table = malloc(sizeof(ck_hash_table_t))) == NULL)
        return NULL;

    if ((table->keys = malloc(num_entries * mean_key_length)) == NULL) {
        free(table);
        return NULL;
    }
    table->keys_capacity = num_entries * mean_key_length;

    num_entries *= 2;

    if ((table->entries = malloc(
        num_entries * sizeof(ck_hash_entry_t))) == NULL
    ) {
        free(table->keys);
        free(table);
        return NULL;
    }
    table->capacity = num_entries;
    ck_hash_table_wipe(table);
    return table;
}

void ck_hash_table_free(ck_hash_table_t *table) {
    free(table->entries);
    if (table->keys)
        free(table->keys);
    free(table);
}

void ck_hash_table_wipe(ck_hash_table_t *table) {
    table->keys_used = 0;
    table->count = 0;
    memset(table->entries, 0, table->capacity * sizeof(ck_hash_entry_t));
}

int ck_hash_table_grow(ck_hash_table_t *table) {
    ck_hash_entry_t *old_entries = table->entries;
    uint64_t old_capacity = table->capacity;
    uint64_t new_capacity = 2 * table->capacity;
    if ((table->entries = calloc(
        new_capacity,
        sizeof(ck_hash_entry_t))) == NULL
    ) {
        return -1;
    }
    table->capacity = new_capacity;
    table->count = 0;
    for (unsigned int i = 0; i < old_capacity; i++) {
        if (old_entries[i].key_length != 0) {
            char *this_key = &table->keys[old_entries[i].key_offset];
            uint64_t hash_key = ck_hash_str(
                this_key,
                old_entries[i].key_length);
            if (!ck_hash_insert_nocopy(
                old_entries[i].key_offset,
                old_entries[i].key_length,
                hash_key,
                old_entries[i].value, table)
            )
                return -1;
        }
    }
    free(old_entries);
    return 0;
}
