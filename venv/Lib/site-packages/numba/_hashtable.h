/*
 * See _hashtable.c for more information about this file.
 */

#ifndef Py_HASHTABLE_H
#define Py_HASHTABLE_H

/* The whole API is private */
#ifndef Py_LIMITED_API

typedef struct _Py_slist_item_s {
    struct _Py_slist_item_s *next;
} _Py_slist_item_t;

typedef struct {
    _Py_slist_item_t *head;
} _Py_slist_t;

#define _Py_SLIST_ITEM_NEXT(ITEM) (((_Py_slist_item_t *)ITEM)->next)

#define _Py_SLIST_HEAD(SLIST) (((_Py_slist_t *)SLIST)->head)

typedef struct {
    /* used by _Numba_hashtable_t.buckets to link entries */
    _Py_slist_item_t _Py_slist_item;

    const void *key;
    Py_uhash_t key_hash;

    /* data follows */
} _Numba_hashtable_entry_t;

#define _Numba_HASHTABLE_ENTRY_DATA(ENTRY) \
        ((char *)(ENTRY) + sizeof(_Numba_hashtable_entry_t))

#define _Numba_HASHTABLE_ENTRY_DATA_AS_VOID_P(ENTRY) \
        (*(void **)_Numba_HASHTABLE_ENTRY_DATA(ENTRY))

#define _Numba_HASHTABLE_ENTRY_READ_DATA(TABLE, DATA, DATA_SIZE, ENTRY) \
    do { \
        assert((DATA_SIZE) == (TABLE)->data_size); \
        memcpy(DATA, _Numba_HASHTABLE_ENTRY_DATA(ENTRY), DATA_SIZE); \
    } while (0)

typedef Py_uhash_t (*_Numba_hashtable_hash_func) (const void *key);
typedef int (*_Numba_hashtable_compare_func) (const void *key, const _Numba_hashtable_entry_t *he);
typedef void* (*_Numba_hashtable_copy_data_func)(void *data);
typedef void (*_Numba_hashtable_free_data_func)(void *data);
typedef size_t (*_Numba_hashtable_get_data_size_func)(void *data);

typedef struct {
    /* allocate a memory block */
    void* (*malloc) (size_t size);

    /* release a memory block */
    void (*free) (void *ptr);
} _Numba_hashtable_allocator_t;

typedef struct {
    size_t num_buckets;
    size_t entries; /* Total number of entries in the table. */
    _Py_slist_t *buckets;
    size_t data_size;

    _Numba_hashtable_hash_func hash_func;
    _Numba_hashtable_compare_func compare_func;
    _Numba_hashtable_copy_data_func copy_data_func;
    _Numba_hashtable_free_data_func free_data_func;
    _Numba_hashtable_get_data_size_func get_data_size_func;
    _Numba_hashtable_allocator_t alloc;
} _Numba_hashtable_t;

/* hash and compare functions for integers and pointers */
extern "C" PyAPI_FUNC(Py_uhash_t) _Numba_hashtable_hash_ptr(const void *key);
extern "C" PyAPI_FUNC(Py_uhash_t) _Numba_hashtable_hash_int(const void *key);
extern "C" PyAPI_FUNC(int) _Numba_hashtable_compare_direct(const void *key, const _Numba_hashtable_entry_t *entry);

extern "C" PyAPI_FUNC(_Numba_hashtable_t *) _Numba_hashtable_new(
    size_t data_size,
    _Numba_hashtable_hash_func hash_func,
    _Numba_hashtable_compare_func compare_func);
extern "C" PyAPI_FUNC(_Numba_hashtable_t *) _Numba_hashtable_new_full(
    size_t data_size,
    size_t init_size,
    _Numba_hashtable_hash_func hash_func,
    _Numba_hashtable_compare_func compare_func,
    _Numba_hashtable_copy_data_func copy_data_func,
    _Numba_hashtable_free_data_func free_data_func,
    _Numba_hashtable_get_data_size_func get_data_size_func,
    _Numba_hashtable_allocator_t *allocator);
extern "C" PyAPI_FUNC(_Numba_hashtable_t *) _Numba_hashtable_copy(_Numba_hashtable_t *src);
extern "C" PyAPI_FUNC(void) _Numba_hashtable_clear(_Numba_hashtable_t *ht);
extern "C" PyAPI_FUNC(void) _Numba_hashtable_destroy(_Numba_hashtable_t *ht);

typedef int (*_Numba_hashtable_foreach_func) (_Numba_hashtable_entry_t *entry, void *arg);

extern "C" PyAPI_FUNC(int) _Numba_hashtable_foreach(
    _Numba_hashtable_t *ht,
    _Numba_hashtable_foreach_func func, void *arg);
extern "C" PyAPI_FUNC(size_t) _Numba_hashtable_size(_Numba_hashtable_t *ht);

extern "C" PyAPI_FUNC(_Numba_hashtable_entry_t*) _Numba_hashtable_get_entry(
    _Numba_hashtable_t *ht,
    const void *key);
extern "C" PyAPI_FUNC(int) _Numba_hashtable_set(
    _Numba_hashtable_t *ht,
    const void *key,
    void *data,
    size_t data_size);
extern "C" PyAPI_FUNC(int) _Numba_hashtable_get(
    _Numba_hashtable_t *ht,
    const void *key,
    void *data,
    size_t data_size);
extern "C" PyAPI_FUNC(int) _Numba_hashtable_pop(
    _Numba_hashtable_t *ht,
    const void *key,
    void *data,
    size_t data_size);
extern "C" PyAPI_FUNC(void) _Numba_hashtable_delete(
    _Numba_hashtable_t *ht,
    const void *key);

#define _Numba_HASHTABLE_SET(TABLE, KEY, DATA) \
    _Numba_hashtable_set(TABLE, KEY, &(DATA), sizeof(DATA))

#define _Numba_HASHTABLE_GET(TABLE, KEY, DATA) \
    _Numba_hashtable_get(TABLE, KEY, &(DATA), sizeof(DATA))

#endif   /* Py_LIMITED_API */

#endif
