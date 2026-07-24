/* Adapted from CPython3.7 Objects/set-common.h */
#include "cext.h"

#ifndef NUMBA_SET_COMMON_H
#define NUMBA_SET_COMMON_H

typedef struct {
    /* Uses Py_ssize_t instead of Py_hash_t to guarantee word size alignment */
    Py_ssize_t  hash;
    char        key[];
} NB_SetEntry;


typedef int (*set_key_comparator_t)(const char *lhs, const char *rhs);
typedef void (*set_refcount_op_t)(const void*);

typedef struct {
    set_key_comparator_t    key_equal;
    set_refcount_op_t       key_incref;
    set_refcount_op_t       key_decref;
} set_type_based_methods_table;

typedef struct {
   /* hash table size */
    Py_ssize_t      size;
    /* Usable size of the hash table.
       Also, size of the entries */
    Py_ssize_t      usable;
    /* hash table used entries */
    Py_ssize_t      nentries;
    /* Entry info
        - key_size is the sizeof key type
        - entry_size is key_size + alignment
    */
    Py_ssize_t      key_size, entry_size;
    /* Byte offset from indices to the first entry. */
    Py_ssize_t      entry_offset;

    /* Method table for type-dependent operations. */
    set_type_based_methods_table methods;

    /* hash table */
    char            indices[];
} NB_SetKeys;


typedef struct {
    /* num of elements in the hashtable */
    Py_ssize_t        used;
    NB_SetKeys      *keys;
} NB_Set;

typedef struct {
    /* parent set */
    NB_Set         *parent;
    /* parent keys object */
    NB_SetKeys     *parent_keys;
    /* set size */
    Py_ssize_t       size;
    /* iterator position; indicates the next position to read */
    Py_ssize_t       pos;
} NB_SetIter;

/* A test function for the set
Returns 0 for OK; 1 for failure.
*/
NUMBA_EXPORT_FUNC(int)
numba_test_set(void);

/* Allocate a new set
Parameters
- NB_Set **out
    Output for the new set.
- Py_ssize_t size
    Hashtable size. Must be power of two.
- Py_ssize_t key_size
    Size of a key entry.
*/
NUMBA_EXPORT_FUNC(int)
numba_set_new(NB_Set **out, Py_ssize_t key_size, Py_ssize_t size);

/* Allocate a new set with enough space to hold n_keys without resizing.
Parameters
- NB_Set **out
    Output for the new set.
- Py_ssize_t size
    The number of keys to fit without needing resize.
- Py_ssize_t key_size
    Size of a key entry.
*/
NUMBA_EXPORT_FUNC(int)
numba_set_new_sized(NB_Set** out, Py_ssize_t key_size, Py_ssize_t size);

/* Free a set */
NUMBA_EXPORT_FUNC(void)
numba_set_free(NB_Set *setp);

/* Returns length of a set */
NUMBA_EXPORT_FUNC(Py_ssize_t)
numba_set_length(NB_Set *setp);

/* Set the method table for type specific operations
*/
NUMBA_EXPORT_FUNC(void)
numba_set_set_method_table(NB_Set *setp, set_type_based_methods_table *methods);

/* Lookup a key

Parameters
- NB_Set *setp
    The set object.
- const char *key_bytes
    The key as a byte buffer.
- Py_hash_t hash
    The precomputed hash of the key.
*/
NUMBA_EXPORT_FUNC(Py_ssize_t)
numba_set_contains(NB_Set *setp, char *key_bytes, Py_hash_t hash);

/* Insert to the set

Parameters
- NB_Set *setp
    The set object.
- const char *key_bytes
    The key as a byte buffer.
- Py_hash_t hash
    The precomputed hash of key.

Returns
- < 0 for error
- 0 for ok
*/
NUMBA_EXPORT_FUNC(int)
numba_set_add(NB_Set *setp, const char *key_bytes, Py_hash_t hash);

/* Delete an entry from the set
Parameters
- NB_Set *setp
    The set
- char *key_bytes
    Output. The key as a byte buffer
- Py_hash_t hash
    Precomputed hash of the key to be deleted
*/
NUMBA_EXPORT_FUNC(int)
numba_set_discard(NB_Set *setp, char *key_bytes, Py_hash_t hash);

/* Returns the sizeof a set iterator
*/
NUMBA_EXPORT_FUNC(size_t)
numba_set_iter_sizeof(void);

/* Fill a NB_SetIter for a set to begin iteration
Parameters
- NB_SetIter *it
    Output.  Must points to memory of size at least `numba_set_iter_sizeof()`.
- NB_Set *setp
    The set to be iterated.
*/
NUMBA_EXPORT_FUNC(void)
numba_set_iter(NB_SetIter *it, NB_Set *setp);

/* Advance the iterator
Parameters
- NB_SetIter *it
    The iterator
- const char **key_ptr
    Output pointer for the key.  Points to data in the set.

Returns
- 0 for success; valid key_ptr
- ERR_ITER_EXHAUSTED for end of iterator.
- ERR_SET_MUTATED for detected set mutation.
*/
NUMBA_EXPORT_FUNC(int)
numba_set_iter_next(NB_SetIter *it, const char **key_ptr);


NUMBA_EXPORT_FUNC(void)
numba_set_dump(NB_Set *setp);

#endif
