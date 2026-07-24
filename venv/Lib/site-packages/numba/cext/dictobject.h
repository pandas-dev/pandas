/* Adapted from CPython3.7 Objects/dict-common.h */
#include "cext.h"

#ifndef NUMBA_DICT_COMMON_H
#define NUMBA_DICT_COMMON_H

typedef struct {
    /* Uses Py_ssize_t instead of Py_hash_t to guarantee word size alignment */
    Py_ssize_t  hash;
    char        keyvalue[];
} NB_DictEntry;


typedef int (*dict_key_comparator_t)(const char *lhs, const char *rhs);
typedef void (*dict_refcount_op_t)(const void*);


typedef struct {
    dict_key_comparator_t    key_equal;
    dict_refcount_op_t       key_incref;
    dict_refcount_op_t       key_decref;
    dict_refcount_op_t       value_incref;
    dict_refcount_op_t       value_decref;
} type_based_methods_table;


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
        - val_size is the sizeof value type
        - entry_size is key_size + val_size + alignment
    */
    Py_ssize_t      key_size, val_size, entry_size;
    /* Byte offset from indices to the first entry. */
    Py_ssize_t      entry_offset;

    /* Method table for type-dependent operations. */
    type_based_methods_table methods;

    /* hash table */
    char            indices[];
} NB_DictKeys;


typedef struct {
    /* num of elements in the hashtable */
    Py_ssize_t        used;
    NB_DictKeys      *keys;
} NB_Dict;


typedef struct {
    /* parent dictionary */
    NB_Dict         *parent;
    /* parent keys object */
    NB_DictKeys     *parent_keys;
    /* dict size */
    Py_ssize_t       size;
    /* iterator position; indicates the next position to read */
    Py_ssize_t       pos;
} NB_DictIter;



/* A test function for the dict
Returns 0 for OK; 1 for failure.
*/
NUMBA_EXPORT_FUNC(int)
numba_test_dict(void);

/* Allocate a new dict
Parameters
- NB_Dict **out
    Output for the new dictionary.
- Py_ssize_t size
    Hashtable size. Must be power of two.
- Py_ssize_t key_size
    Size of a key entry.
- Py_ssize_t val_size
    Size of a value entry.
*/
NUMBA_EXPORT_FUNC(int)
numba_dict_new(NB_Dict **out, Py_ssize_t size, Py_ssize_t key_size, Py_ssize_t val_size);

/* Allocate a new dict with enough space to hold n_keys without resizing.
Parameters
- NB_Dict **out
    Output for the new dictionary.
- Py_ssize_t n_keys
    The number of keys to fit without needing resize.
- Py_ssize_t key_size
    Size of a key entry.
- Py_ssize_t val_size
    Size of a value entry.
*/
NUMBA_EXPORT_FUNC(int)
numba_dict_new_sized(NB_Dict** out, Py_ssize_t n_keys, Py_ssize_t key_size, Py_ssize_t val_size);

/* Free a dict */
NUMBA_EXPORT_FUNC(void)
numba_dict_free(NB_Dict *d);

/* Returns length of a dict */
NUMBA_EXPORT_FUNC(Py_ssize_t)
numba_dict_length(NB_Dict *d);

/* Set the method table for type specific operations
*/
NUMBA_EXPORT_FUNC(void)
numba_dict_set_method_table(NB_Dict *d, type_based_methods_table *methods);

/* Lookup a key

Parameters
- NB_Dict *d
    The dictionary object.
- const char *key_bytes
    The key as a byte buffer.
- Py_hash_t hash
    The precomputed hash of the key.
- char *oldval_bytes
    An output parameter to store the associated value if the key is found.
    Must point to memory of sufficient size to store the value.
*/
NUMBA_EXPORT_FUNC(Py_ssize_t)
numba_dict_lookup(NB_Dict *d, const char *key_bytes, Py_hash_t hash, char *oldval_bytes);

/* Resize the dict to at least *minsize*.
*/
NUMBA_EXPORT_FUNC(int)
numba_dict_resize(NB_Dict *d, Py_ssize_t minsize);

/* Insert to the dict

Parameters
- NB_Dict *d
    The dictionary object.
- const char *key_bytes
    The key as a byte buffer.
- Py_hash_t hash
    The precomputed hash of key.
- const char *val_bytes
    The value as a byte buffer.
- char *oldval_bytes
    An output buffer to store the replaced value.
    Must point to memory of sufficient size to store the value.

Returns
- < 0 for error
- 0 for ok
- 1 for ok and oldval_bytes has a copy of the replaced value.
*/
NUMBA_EXPORT_FUNC(int)
numba_dict_insert(NB_Dict *d, const char *key_bytes, Py_hash_t hash, const char *val_bytes, char *oldval_bytes);

/* Same as numba_dict_insert() but oldval_bytes is not needed */
NUMBA_EXPORT_FUNC(int)
numba_dict_insert_ez(NB_Dict *d, const char *key_bytes, Py_hash_t hash, const char *val_bytes);

/* Delete an entry from the dict
Parameters
- NB_Dict *d
    The dictionary
- Py_hash_t hash
    Precomputed hash of the key to be deleted
- Py_ssize_t ix
    Precomputed entry index of the key to be deleted.
    Usually results of numba_dict_lookup().
*/
NUMBA_EXPORT_FUNC(int)
numba_dict_delitem(NB_Dict *d, Py_hash_t hash, Py_ssize_t ix);

/* Remove an item from the dict
Parameters
- NB_Dict *d
    The dictionary
- char *key_bytes
    Output. The key as a byte buffer
- char *val_bytes
    Output. The value as a byte buffer
*/
NUMBA_EXPORT_FUNC(int)
numba_dict_popitem(NB_Dict *d, char *key_bytes, char *val_bytes);

/* Returns the sizeof a dictionary iterator
*/
NUMBA_EXPORT_FUNC(size_t)
numba_dict_iter_sizeof(void);

/* Fill a NB_DictIter for a dictionary to begin iteration
Parameters
- NB_DictIter *it
    Output.  Must points to memory of size at least `numba_dict_iter_sizeof()`.
- NB_Dict *d
    The dictionary to be iterated.
*/
NUMBA_EXPORT_FUNC(void)
numba_dict_iter(NB_DictIter *it, NB_Dict *d);

/* Advance the iterator
Parameters
- NB_DictIter *it
    The iterator
- const char **key_ptr
    Output pointer for the key.  Points to data in the dictionary.
- const char **val_ptr
    Output pointer for the key.  Points to data in the dictionary.

Returns
- 0 for success; valid key_ptr and val_ptr
- ERR_ITER_EXHAUSTED for end of iterator.
- ERR_DICT_MUTATED for detected dictionary mutation.
*/
NUMBA_EXPORT_FUNC(int)
numba_dict_iter_next(NB_DictIter *it, const char **key_ptr, const char **val_ptr);


NUMBA_EXPORT_FUNC(void)
numba_dict_dump(NB_Dict *);

#endif
