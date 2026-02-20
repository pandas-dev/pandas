/* Adapted from CPython3.7 Include/listobject.h
 *
 * The exact commit-id of the relevant file is:
 *
 * https://github.com/python/cpython/blob/51ddab8dae056867f3595ab3400bffc93f67c8d4/Include/listobject.h
 *
 * WARNING:
 * Most interfaces listed here are exported (global), but they are not
 * supported, stable, or part of Numba's public API. These interfaces and their
 * underlying implementations may be changed or removed in future without
 * notice.
 * */

#ifndef NUMBA_LIST_H
#define NUMBA_LIST_H

#include "cext.h"

typedef void (*list_refcount_op_t)(const void*);

typedef struct {
    list_refcount_op_t       item_incref;
    list_refcount_op_t       item_decref;
} list_type_based_methods_table;

/* This is the struct for the Numba typed list. It is largely inspired by the
 * CPython list struct in listobject.h. In essence the list is a homogeneously
 * typed container that can grow and shrink upon insertion and deletion. This
 * means that appending an item to, or removing an item from, the end of the
 * list, this will have a O(1) amortized runtime. This matches the
 * behaviour of the CPython list type and it will grow with the same
 * increments.
 *
 * 'items' contains space for 'allocated' elements.  The number
 * currently in use is 'size'. The size in bytes of the items stored in the
 * list is given by 'item_size'.
 *
 * Invariants:
 *     0 <= size <= allocated
 *     len(list) == size
 *     item == NULL implies size == allocated == 0
 *
 * FIXME: list.sort() temporarily sets allocated to -1 to detect mutations.
 *
 * Items must normally not be NULL, except during construction when
 * the list is not yet visible outside the function that builds it.
 *
 * Additionally, this list has boolean member 'is_mutable' that can be used to
 * set a list as immutable. Two functions to query and set this member are
 * provided. Any attempt to mutate an immutable list will result in a status
 * of LIST_ERR_IMMUTABLE.
 *
 */
typedef struct {
    /* size of the list in items  */
    Py_ssize_t      size;
    /* size of the list items in bytes */
    Py_ssize_t      item_size;
    /* total allocated slots in items */
    Py_ssize_t allocated;
    /* is the list mutable */
    int is_mutable;
    /* method table for type-dependent operations */
    list_type_based_methods_table methods;
    /* array/pointer for items. Interpretation is governed by item_size */
    char  * items;
} NB_List;


typedef struct {
    /* parent list */
    NB_List         *parent;
    /* list size */
    Py_ssize_t       size;
    /* iterator position; indicates the next position to read */
    Py_ssize_t       pos;
} NB_ListIter;

NUMBA_GLOBAL_FUNC(void)
numba_list_set_method_table(NB_List *lp, list_type_based_methods_table *methods);

NUMBA_GLOBAL_FUNC(int)
numba_list_new(NB_List **out, Py_ssize_t item_size, Py_ssize_t allocated);

NUMBA_GLOBAL_FUNC(void)
numba_list_free(NB_List *lp);

NUMBA_GLOBAL_FUNC(char *)
numba_list_base_ptr(NB_List *lp);

NUMBA_GLOBAL_FUNC(Py_ssize_t)
numba_list_size_address(NB_List *lp);

NUMBA_GLOBAL_FUNC(Py_ssize_t)
numba_list_length(NB_List *lp);

NUMBA_GLOBAL_FUNC(Py_ssize_t)
numba_list_allocated(NB_List *lp);

NUMBA_GLOBAL_FUNC(int)
numba_list_is_mutable(NB_List *lp);

NUMBA_GLOBAL_FUNC(void)
numba_list_set_is_mutable(NB_List *lp, int is_mutable);

NUMBA_GLOBAL_FUNC(int)
numba_list_setitem(NB_List *lp, Py_ssize_t index, const char *item);

NUMBA_GLOBAL_FUNC(int)
numba_list_getitem(NB_List *lp, Py_ssize_t index, char *out);

NUMBA_GLOBAL_FUNC(int)
numba_list_append(NB_List *lp, const char *item);

NUMBA_GLOBAL_FUNC(int)
numba_list_resize(NB_List *lp, Py_ssize_t newsize);

NUMBA_GLOBAL_FUNC(int)
numba_list_delitem(NB_List *lp, Py_ssize_t index);

NUMBA_GLOBAL_FUNC(int)
numba_list_delete_slice(NB_List *lp,
                        Py_ssize_t start, Py_ssize_t stop, Py_ssize_t step);

NUMBA_GLOBAL_FUNC(size_t)
numba_list_iter_sizeof(void);

NUMBA_GLOBAL_FUNC(void)
numba_list_iter(NB_ListIter *it, NB_List *l);

NUMBA_GLOBAL_FUNC(int)
numba_list_iter_next(NB_ListIter *it, const char **item_ptr);

NUMBA_EXPORT_FUNC(int)
numba_test_list(void);

#endif
