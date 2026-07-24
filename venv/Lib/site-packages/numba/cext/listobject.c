#include "listobject.h"

/* This implements the C component of the Numba typed list. It is loosely
 * inspired by the list implementation of the cpython list with some parts
 * taken from the cpython slice implementation. The exact commit-id of the
 * relevant files are:
 *
 * https://github.com/python/cpython/blob/51ddab8dae056867f3595ab3400bffc93f67c8d4/Objects/listobject.c
 * https://github.com/python/cpython/blob/51ddab8dae056867f3595ab3400bffc93f67c8d4/Objects/sliceobject.c
 *
 * Algorithmically, this list is very similar to the cpython implementation so
 * it should have the same performance (Big-O) characteristics for accessing,
 * adding and removing elements/items. Specifically, it implements the same
 * algorithms for list overallocation and growth. However, it never deals with
 * PyObject types and instead must be typed with a type-size.  As a result, the
 * typed-list is type homogeneous and in contrast to the cpython version can
 * not store a mixture of arbitrarily typed objects. Reference counting via the
 * Numba Runtime (NRT) is supported and incrementing and decrementing functions
 * are store as part of the struct and can be setup from the compiler level.
 *
 * Importantly, only a very limited subset of the cpython c functions have been
 * ported over and the rest have been implemented (in Python) at the compiler
 * level using the c functions provided. Additionally, initialization of, and
 * iteration over, a ListIter is provided
 *
 * The following functions are implemented for the list:
 *
 * - Check valid index        valid_index
 * - Creation                 numba_list_new
 * - Deletion                 numba_list_free
 * - Accessing the length     numba_list_length
 * - Appending to the list    numba_list_append
 * - Getting an item          numba_list_setitem
 * - Setting an item          numba_list_getitem
 * - Resizing the list        numba_list_resize
 * - Deleting an item         numba_list_delitem
 * - Deleting a slice         numba_list_delete_slice
 *
 * As you can see, only a single function for slices is implemented. The rest
 * is all done entirely at the compiler level which then calls the c functions
 * to mutate the list accordingly. Since slicing allows for replace, insert and
 * delete operations over multiple items, we can simply implement those using
 * the basic functions above.
 *
 * The following additional functions are implemented for the list, these are
 * needed to make the list work within Numba.
 *
 * - Accessing the allocation numba_list_allocated
 * - Copying an item          copy_item
 * - Calling incref on item   list_incref_item
 * - Calling decref on item   list_decref_item
 * - Set method table         numba_list_set_method_table
 *
 * The following functions are implemented for the iterator:
 *
 * - Size of the iterator     numba_list_iter_size
 * - Initialization of iter   numba_list_iter
 * - Get next item from iter  numba_list_iter_next
 *
 * Two methods are provided to query and set the 'is_mutable':
 *
 * - Query                    numba_list_is_mutable
 * - Set                      numba_list_set_is_mutable
 *
 * Lastly a set of pure C level tests are provided which come in handy when
 * needing to use valgrind and friends.
 *
 */


/* Return status for the list functions.
 */
typedef enum {
    LIST_OK = 0,
    LIST_ERR_INDEX = -1,
    LIST_ERR_NO_MEMORY = -2,
    LIST_ERR_MUTATED = -3,
    LIST_ERR_ITER_EXHAUSTED = -4,
    LIST_ERR_IMMUTABLE = -5,
} ListStatus;

/* Copy an item from a list.
 *
 * lp: a list
 * dst: destination pointer
 * src: source pointer
 */
static void
copy_item(NB_List *lp, char *dst, const char *src){
    memcpy(dst, src, lp->item_size);
}

/* Increment a reference to an item in a list.
 *
 * lp: a list
 * item: the item to increment the reference for
 */
static void
list_incref_item(NB_List *lp, const char *item){
    if (lp->methods.item_incref) {
        lp->methods.item_incref(item);
    }
}

/* Decrement a reference to an item in a list.
 *
 * lp: a list
 * item: the item to decrement the reference for
 */
static void
list_decref_item(NB_List *lp, const char *item){
    if (lp->methods.item_decref) {
        lp->methods.item_decref(item);
    }
}

/* Setup the method table for a list.
 *
 * This function is used from the compiler level to initialize the internal
 * method table.
 *
 * lp: a list
 * methods: the methods table to set up
 */
void
numba_list_set_method_table(NB_List *lp, list_type_based_methods_table *methods)
{
    memcpy(&lp->methods, methods, sizeof(list_type_based_methods_table));
}

/* Check if a list index is valid.
 *
 * i: the index to check
 * limit: the size of a list
 *
 * Adapted from CPython's valid_index().
 *
 * FIXME: need to find a way to inline this, even for Python 2.7 on Windows
 */
static int
valid_index(Py_ssize_t i, Py_ssize_t limit){
    /* The cast to size_t lets us use just a single comparison
       to check whether i is in the range: 0 <= i < limit.

       See:  Section 14.2 "Bounds Checking" in the Agner Fog
       optimization manual found at:
       https://www.agner.org/optimize/optimizing_cpp.pdf
    */
    return (size_t) i < (size_t) limit;
}

/* Initialize a new list.
 *
 * out: pointer to hold an initialized list
 * item_size: the size in bytes of the items in the list
 * allocated: preallocation of the list in items
 *
 * This will allocate sufficient memory to hold the list structure and any
 * items if requested (allocated != 0). See _listobject.h for more information
 * on the NB_List struct.
 */
int
numba_list_new(NB_List **out, Py_ssize_t item_size, Py_ssize_t allocated){
    NB_List *lp;
    char *items;
    // allocate memory to hold the struct
    lp = malloc(aligned_size(sizeof(NB_List)));
    if (lp == NULL) {
        return LIST_ERR_NO_MEMORY;
    }
    // set up members
    lp->size = 0;
    lp->item_size = item_size;
    lp->allocated = allocated;
    lp->is_mutable = 1;
    // set method table to zero */
    memset(&lp->methods, 0x00, sizeof(list_type_based_methods_table));
    // allocate memory to hold items, if requested
    if (allocated != 0) {
        items = malloc(aligned_size(lp->item_size * allocated));
        // allocated was definitely not zero, if malloc returns NULL
        // this is definitely an error
        if (items == NULL) {
            // free previously allocated struct to avoid leaking memory
            free(lp);
            return LIST_ERR_NO_MEMORY;
        }
        lp->items = items;
    }
    else {
        // be explicit
        lp->items = NULL;
    }
    *out = lp;
    return LIST_OK;
}

/* Free the memory associated with a list.
 *
 * lp: a list
 */
void
numba_list_free(NB_List *lp) {
    // decref all items, if needed
    Py_ssize_t i;
    if (lp->methods.item_decref) {
        for (i = 0; i < lp->size; i++) {
            char *item = lp->items + lp->item_size * i;
            list_decref_item(lp, item);
        }
    }
    // free items and list
    if (lp->items != NULL) {
        free(lp->items);
    }
    free(lp);
}

/* Return the base pointer of the list items.
 */
char *
numba_list_base_ptr(NB_List *lp)
{
    return lp->items;
}

/* Return the address of the list size.
 */
Py_ssize_t
numba_list_size_address(NB_List *lp)
{
    return (Py_ssize_t)&lp->size;
}


/* Return the length of a list.
 *
 * lp: a list
 */
Py_ssize_t
numba_list_length(NB_List *lp) {
    return lp->size;
}

/* Return the current allocation of a list.
 *
 * lp: a list
 */
Py_ssize_t
numba_list_allocated(NB_List *lp) {
    return lp->allocated;
}

/* Return the mutability status of the list
 *
 * lp: a list
 *
 */
int
numba_list_is_mutable(NB_List *lp){
    return lp->is_mutable;
}

/* Set the is_mutable attribute
 *
 * lp: a list
 * is_mutable: an int, 0(False) or 1(True)
 *
 */
void
numba_list_set_is_mutable(NB_List *lp, int is_mutable){
    lp->is_mutable = is_mutable;
}

/* Set an item in a list.
 *
 * lp: a list
 * index: the index of the item to set (must be in range 0 <= index < len(list))
 * item: the item to set
 *
 * This assume there is already an element at the given index that will be
 * overwritten and thereby have its reference decremented.  DO NOT use this to
 * write to an unassigned location.
 */
int
numba_list_setitem(NB_List *lp, Py_ssize_t index, const char *item) {
    char *loc;
    // check for mutability
    if (!lp->is_mutable) {
        return LIST_ERR_IMMUTABLE;
    }
    // check index is valid
    // FIXME: this can be (and probably is) checked at the compiler level
    if (!valid_index(index, lp->size)) {
        return LIST_ERR_INDEX;
    }
    // set item at desired location
    loc = lp->items + lp-> item_size * index;
    list_decref_item(lp, loc);
    copy_item(lp, loc, item);
    list_incref_item(lp, loc);
    return LIST_OK;
}

/* Get an item from a list.
 *
 * lp: a list
 * index: the index of the item to get (must be in range 0 <= index < len(list))
 * out: a pointer to hold the item
 */
int
numba_list_getitem(NB_List *lp, Py_ssize_t index, char *out) {
    char *loc;
    // check index is valid
    // FIXME: this can be (and probably is) checked at the compiler level
    if (!valid_index(index, lp->size)) {
        return LIST_ERR_INDEX;
    }
    // get item at desired location
    loc = lp->items + lp->item_size * index;
    copy_item(lp, out, loc);
    return LIST_OK;
}

/* Append an item to the end of a list.
 *
 * lp: a list
 * item: the item to append.
 */
int
numba_list_append(NB_List *lp, const char *item) {
    char *loc;
    // check for mutability
    if (!lp->is_mutable) {
        return LIST_ERR_IMMUTABLE;
    }
    // resize by one, will change list size
    int result = numba_list_resize(lp, lp->size + 1);
    if(result < LIST_OK) {
        return result;
    }
    // insert item at index: original size before resize
    loc = lp->items + lp->item_size * (lp->size - 1);
    copy_item(lp, loc, item);
    list_incref_item(lp, loc);
    return LIST_OK;
}

/* Resize a list.
 *
 * lp: a list
 * newsize: the desired new size of the list.
 *
 * This will increase or decrease the size of the list, including reallocating
 * the required memory and increasing the total allocation (additional free
 * space to hold new items).
 *
 *
 * Adapted from CPython's list_resize().
 *
 * Ensure lp->items has room for at least newsize elements, and set
 * lp->size to newsize.  If newsize > lp->size on entry, the content
 * of the new slots at exit is undefined heap trash; it's the caller's
 * responsibility to overwrite them with sane values.
 * The number of allocated elements may grow, shrink, or stay the same.
 * Failure is impossible if newsize <= lp->allocated on entry, although
 * that partly relies on an assumption that the system realloc() never
 * fails when passed a number of bytes <= the number of bytes last
 * allocated (the C standard doesn't guarantee this, but it's hard to
 * imagine a realloc implementation where it wouldn't be true).
 * Note that lp->items may change, and even if newsize is less
 * than lp->size on entry.
 */
int
numba_list_resize(NB_List *lp, Py_ssize_t newsize) {
    char * items;
    // check for mutability
    if (!lp->is_mutable) {
        return LIST_ERR_IMMUTABLE;
    }
    size_t new_allocated, num_allocated_bytes;
    /* Bypass realloc() when a previous overallocation is large enough
       to accommodate the newsize.  If the newsize falls lower than half
       the allocated size, then proceed with the realloc() to shrink the list.
    */
    if (lp->allocated >= newsize && newsize >= (lp->allocated >> 1)) {
        assert(lp->items != NULL || newsize == 0);
        lp->size = newsize;
        return LIST_OK;
    }
    /* This over-allocates proportional to the list size, making room
     * for additional growth.  The over-allocation is mild, but is
     * enough to give linear-time amortized behavior over a long
     * sequence of appends() in the presence of a poorly-performing
     * system realloc().
     * The growth pattern is:  0, 4, 8, 16, 25, 35, 46, 58, 72, 88, ...
     * Note: new_allocated won't overflow because the largest possible value
     *       is PY_SSIZE_T_MAX * (9 / 8) + 6 which always fits in a size_t.
     */
    new_allocated = (size_t)newsize + (newsize >> 3) + (newsize < 9 ? 3 : 6);
    if (new_allocated > (size_t)PY_SSIZE_T_MAX / lp->item_size) {
        return LIST_ERR_NO_MEMORY;
    }

    if (newsize == 0)
        new_allocated = 0;
    num_allocated_bytes = new_allocated * lp->item_size;
    items = realloc(lp->items, aligned_size(num_allocated_bytes));
    // realloc may return NULL if requested size is 0
    if (num_allocated_bytes != 0 && items == NULL) {
        return LIST_ERR_NO_MEMORY;
    }
    lp->items = items;
    lp->size = newsize;
    lp->allocated = (Py_ssize_t)new_allocated;
    return LIST_OK;
}

/* Delete a single item.
 *
 * lp: a list
 * index: the index of the item to delete
 *        (must be in range 0 <= index < len(list))
 *
 * */
int
numba_list_delitem(NB_List *lp, Py_ssize_t index) {
    int result;
    char *loc, *new_loc;
    Py_ssize_t leftover_bytes;
    // check for mutability
    if (!lp->is_mutable) {
        return LIST_ERR_IMMUTABLE;
    }
    // check index is valid
    // FIXME: this can be (and probably is) checked at the compiler level
    if (!valid_index(index, lp->size)) {
        return LIST_ERR_INDEX;
    }
    // obtain item and decref if needed
    loc = lp->items + lp->item_size * index;
    list_decref_item(lp, loc);
    if (index != lp->size - 1) {
        // delitem from somewhere other than the end, incur the memory copy
        leftover_bytes = (lp->size - 1 - index) * lp->item_size;
        new_loc = lp->items + (lp->item_size * (index + 1));
        // use memmove instead of memcpy since we may be dealing with
        // overlapping regions of memory and the behaviour of memcpy is
        // undefined in such situation (C99).
        memmove(loc, new_loc, leftover_bytes);
    }
    // finally, shrink list by one
    result = numba_list_resize(lp, lp->size - 1);
    if(result < LIST_OK) {
         // Since we are decreasing the size, this should never happen
        return result;
    }
    return LIST_OK;

}

/* Delete a slice
 *
 * start: the start index of ths slice
 * stop: the stop index of the slice (not included)
 * step: the step to take
 *
 * This function assumes that the start and stop were clipped appropriately.
 * I.e. if step > 0 start >= 0 and stop <= len(l) and
 *      if step < 0 start <= length and stop >= -1
 *      step != 0 and no Python negative indexing allowed.
 *
 * This code was copied and edited from the relevant section in
 * list_ass_subscript from the cpython implementation, see the top of this file
 * for the exact source
 */
int
numba_list_delete_slice(NB_List *lp,
                        Py_ssize_t start, Py_ssize_t stop, Py_ssize_t step) {
    int result, i, slicelength, new_length;
    char *loc, *new_loc;
    Py_ssize_t leftover_bytes, cur, lim;
    // check for mutability
    if (!lp->is_mutable) {
        return LIST_ERR_IMMUTABLE;
    }
    // calculate the slicelength, taken from PySlice_AdjustIndices, see the top
    // of this file for the exact source
    if (step > 0) {
        slicelength = start < stop ? (stop - start - 1) / step + 1 : 0;
    } else {
        slicelength = stop < start ? (start - stop - 1) / -step + 1 : 0;
    }
    if (slicelength <= 0){
        return LIST_OK;
    }
    new_length = lp->size - slicelength;
    // reverse step and indices
    if (step < 0) {
        stop = start + 1;
        start = stop + step * (slicelength - 1) - 1;
        step = -step;
    }
    if (step == 1) {
        // decref if needed
        if (lp->methods.item_decref) {
            for (i = start ; i < stop ; i++){
                loc = lp->items + lp->item_size * i;
                lp->methods.item_decref(loc);
            }
        }
        // memmove items into place
        leftover_bytes = (lp->size - stop) * lp->item_size;
        loc = lp->items + lp->item_size * start;
        new_loc = lp->items + lp->item_size * stop;
        memmove(loc, new_loc, leftover_bytes);
    }
    else { // step != 1
        /* drawing pictures might help understand these for
         * loops. Basically, we memmove the parts of the
         * list that are *not* part of the slice: step-1
         * items for each item that is part of the slice,
         * and then tail end of the list that was not
         * covered by the slice
         *
         * */
        for (cur = start,   // index of item to be deleted
             i = 0;         // counter of total items deleted so far
             cur < stop;
             cur += step,
             i++) {
            lim = step - 1; // number of leftover items after deletion of item
            // clip limit, in case we are at the end of the slice, and there
            // are now less than step-1 items to be moved
            if (cur + step >= lp->size) {
                lim = lp->size - cur - 1;
            }
            // decref item being removed
            loc = lp->items + lp->item_size * cur;
            list_decref_item(lp, loc);
            /* memmove the aforementioned step-1 (or less) items
             * dst : index of deleted item minus total deleted sofar
             * src : index of deleted item plus one (next item)
             */
            memmove(lp->items + lp->item_size * (cur - i),
                    lp->items + lp->item_size * (cur + 1),
                    lim * lp->item_size);
        }
        // memmove tail of the list
        cur = start + slicelength * step;
        if (cur < lp->size) {
            memmove(lp->items + lp->item_size * (cur - slicelength),
                    lp->items + lp->item_size * cur,
                    (lp->size - cur) * lp->item_size);
        }
    }
    // resize to correct size
    result = numba_list_resize(lp, new_length);
    if(result < LIST_OK) {
        // Since we are decreasing the size, this should never happen
        return result;
    }
    return LIST_OK;
}


/* Return the size of the list iterator (NB_ListIter) struct.
 */
size_t
numba_list_iter_sizeof() {
    return sizeof(NB_ListIter);
}

/* Initialize a list iterator (NB_ListIter).
 *
 * it: an iterator
 * lp: a list to iterate over
 */
void
numba_list_iter(NB_ListIter *it, NB_List *lp) {
    // set members of iterator
    it->parent = lp;
    it->size = lp->size;
    it->pos = 0;
}

/* Obtain the next item from a list iterator.
 *
 * it: an iterator
 * item_ptr: pointer to hold the next item
 */
int
numba_list_iter_next(NB_ListIter *it, const char **item_ptr) {
    NB_List *lp;
    lp = it->parent;
    /* FIXME: Detect list mutation during iteration */
    if (lp->size != it->size) {
        return LIST_ERR_MUTATED;
    }
    // get next element
    if (it->pos < lp->size) {
        *item_ptr = lp->items + lp->item_size * it->pos++;
        return LIST_OK;
    }else{
        return LIST_ERR_ITER_EXHAUSTED;
    }
}


#define CHECK(CASE) {                                                   \
    if ( !(CASE) ) {                                                    \
        printf("'%s' failed file %s:%d\n", #CASE, __FILE__, __LINE__);   \
        return -1;                                                       \
    }                                                                   \
}

/* Basic C based tests for the list.
 */
int
numba_test_list(void) {
    NB_List *lp = NULL;
    int status, i;
    Py_ssize_t it_count;
    const char *it_item = NULL;
    NB_ListIter iter;
    char got_item[4] = "\x00\x00\x00\x00";
    const char *test_items_1 = NULL, *test_items_2 = NULL;
    char *test_items_3 = NULL;
    puts("test_list");


    status = numba_list_new(&lp, 4, 0);
    CHECK(status == LIST_OK);
    CHECK(lp->item_size == 4);
    CHECK(lp->size == 0);
    CHECK(lp->allocated == 0);
    CHECK(lp->is_mutable == 1);

    // flip and check the is_mutable bit
    CHECK(numba_list_is_mutable(lp) == 1);
    numba_list_set_is_mutable(lp, 0);
    CHECK(numba_list_is_mutable(lp) == 0);
    numba_list_set_is_mutable(lp, 1);
    CHECK(numba_list_is_mutable(lp) == 1);

    // append 1st item, this will cause a realloc
    status = numba_list_append(lp, "abc");
    CHECK(status == LIST_OK);
    CHECK(lp->size == 1);
    CHECK(lp->allocated == 4);
    status = numba_list_getitem(lp, 0, got_item);
    CHECK(status == LIST_OK);
    CHECK(memcmp(got_item, "abc", 4) == 0);

    // append 2nd item
    status = numba_list_append(lp, "def");
    CHECK(status == LIST_OK);
    CHECK(lp->size == 2);
    CHECK(lp->allocated == 4);
    status = numba_list_getitem(lp, 1, got_item);
    CHECK(status == LIST_OK);
    CHECK(memcmp(got_item, "def", 4) == 0);

    // append 3rd item
    status = numba_list_append(lp, "ghi");
    CHECK(status == LIST_OK);
    CHECK(lp->size == 3);
    CHECK(lp->allocated == 4);
    status = numba_list_getitem(lp, 2, got_item);
    CHECK(status == LIST_OK);
    CHECK(memcmp(got_item, "ghi", 4) == 0);

    // append 4th item
    status = numba_list_append(lp, "jkl");
    CHECK(status == LIST_OK);
    CHECK(lp->size == 4);
    CHECK(lp->allocated == 4);
    status = numba_list_getitem(lp, 3, got_item);
    CHECK(status == LIST_OK);
    CHECK(memcmp(got_item, "jkl", 4) == 0);

    // append 5th item, this will cause another realloc
    status = numba_list_append(lp, "mno");
    CHECK(status == LIST_OK);
    CHECK(lp->size == 5);
    CHECK(lp->allocated == 8);
    status = numba_list_getitem(lp, 4, got_item);
    CHECK(status == LIST_OK);
    CHECK(memcmp(got_item, "mno", 4) == 0);

    // overwrite 1st item
    status = numba_list_setitem(lp, 0, "pqr");
    CHECK(status == LIST_OK);
    CHECK(lp->size == 5);
    CHECK(lp->allocated == 8);
    status = numba_list_getitem(lp, 0, got_item);
    CHECK(status == LIST_OK);
    CHECK(memcmp(got_item, "pqr", 4) == 0);

    // get and del 1st item, check item shift
    status = numba_list_getitem(lp, 0, got_item);
    status = numba_list_delitem(lp, 0);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 4);
    CHECK(lp->allocated == 8);
    CHECK(memcmp(got_item, "pqr", 4) == 0);
    CHECK(memcmp(lp->items, "def\x00ghi\x00jkl\x00mno\x00", 16) == 0);

    // get and del last (4th) item, no shift since only last item affected
    status = numba_list_getitem(lp, 3, got_item);
    status = numba_list_delitem(lp, 3);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 3);
    CHECK(lp->allocated == 6);  // this also shrinks the allocation
    CHECK(memcmp(got_item, "mno", 4) == 0);
    CHECK(memcmp(lp->items, "def\x00ghi\x00jkl\x00", 12) == 0);

    // flip and check the is_mutable member
    CHECK(numba_list_is_mutable(lp) == 1);
    numba_list_set_is_mutable(lp, 0);
    CHECK(numba_list_is_mutable(lp) == 0);

    // ensure that any attempts to mutate an immutable list fail
    CHECK(numba_list_setitem(lp, 0, "zzz") == LIST_ERR_IMMUTABLE);
    CHECK(numba_list_append(lp, "zzz") == LIST_ERR_IMMUTABLE);
    CHECK(numba_list_delitem(lp, 0) == LIST_ERR_IMMUTABLE);
    CHECK(numba_list_resize(lp, 23) == LIST_ERR_IMMUTABLE);
    CHECK(numba_list_delete_slice(lp, 0, 3, 1) == LIST_ERR_IMMUTABLE);

    // ensure that all attempts to query/read from and immutable list succeed
    CHECK(numba_list_length(lp) == 3);
    status = numba_list_getitem(lp, 0, got_item);
    CHECK(status == LIST_OK);
    CHECK(memcmp(got_item, "def", 4) == 0);

    // flip the is_mutable member back  and check
    numba_list_set_is_mutable(lp, 1);
    CHECK(numba_list_is_mutable(lp) == 1);

    // test iterator
    CHECK(lp->size > 0);
    numba_list_iter(&iter, lp);
    it_count = 0;
    CHECK(iter.parent == lp);
    CHECK(iter.pos == it_count);

    // current contents of list
    test_items_1 = "def\x00ghi\x00jkl\x00";
    while ( (status = numba_list_iter_next(&iter, &it_item)) == LIST_OK) {
        it_count += 1;
        CHECK(iter.pos == it_count); // check iterator position
        CHECK(it_item != NULL); // quick check item is non-null
        // go fishing in test_items_1
        CHECK(memcmp((const char *)test_items_1 + ((it_count - 1) * 4), it_item, 4) == 0);
    }

    CHECK(status == LIST_ERR_ITER_EXHAUSTED);
    CHECK(lp->size == it_count);

    // free existing list
    numba_list_free(lp);

    // test growth upon append and shrink during delitem
    status = numba_list_new(&lp, 1, 0);
    CHECK(status == LIST_OK);
    CHECK(lp->item_size == 1);
    CHECK(lp->size == 0);
    CHECK(lp->allocated == 0);

    // first, grow the list
    // Use exactly 17 elements, should go through the allocation pattern:
    // 0, 4, 8, 16, 25
    for (i = 0; i < 17 ; i++) {
        switch(i) {
            // Check the allocation before
            case 0:  CHECK(lp->allocated == 0); break;
            case 4:  CHECK(lp->allocated == 4); break;
            case 8:  CHECK(lp->allocated == 8); break;
            case 16: CHECK(lp->allocated == 16); break;
        }
        status = numba_list_append(lp, (const char*)&i);
        CHECK(status == LIST_OK);
        switch(i) {
            // Check that the growth happened accordingly
            case 0:  CHECK(lp->allocated == 4); break;
            case 4:  CHECK(lp->allocated == 8); break;
            case 8:  CHECK(lp->allocated == 16); break;
            case 16: CHECK(lp->allocated == 25); break;
        }
    }
    CHECK(lp->size == 17);

    // Check current contents of list
    test_items_2  = "\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f\x10";
    CHECK(memcmp(lp->items, test_items_2, 17) == 0);

    // Now, delete them again and check that list shrinks
    for (i = 17; i > 0 ; i--) {
        switch(i) {
             // Check the allocation before delitem
             case 17:  CHECK(lp->allocated == 25); break;
             case 12:  CHECK(lp->allocated == 25); break;
             case 9:   CHECK(lp->allocated == 18); break;
             case 6:   CHECK(lp->allocated == 12); break;
             case 4:   CHECK(lp->allocated == 8); break;
             case 3:   CHECK(lp->allocated == 6); break;
             case 2:   CHECK(lp->allocated == 5); break;
             case 1:   CHECK(lp->allocated == 4); break;
        }
        status = numba_list_getitem(lp, i-1, got_item);
        status = numba_list_delitem(lp, i-1);
        CHECK(status == LIST_OK);
        switch(i) {
             // Check that the shrink happened accordingly
             case 17:  CHECK(lp->allocated == 25); break;
             case 12:  CHECK(lp->allocated == 18); break;
             case 9:   CHECK(lp->allocated == 12); break;
             case 6:   CHECK(lp->allocated == 8); break;
             case 4:   CHECK(lp->allocated == 6); break;
             case 3:   CHECK(lp->allocated == 5); break;
             case 2:   CHECK(lp->allocated == 4); break;
             case 1:   CHECK(lp->allocated == 0); break;
        }
    }
    // free existing list
    numba_list_free(lp);


    // Setup list for testing delete_slice
    status = numba_list_new(&lp, 1, 0);
    CHECK(status == LIST_OK);
    CHECK(lp->item_size == 1);
    CHECK(lp->size == 0);
    CHECK(lp->allocated == 0);
    for (i = 0; i < 17 ; i++) {
        status = numba_list_append(lp, (const char*)&i);
        CHECK(status == LIST_OK);
    }
    CHECK(lp->size == 17);
    test_items_3 = "\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f\x10";
    CHECK(memcmp(lp->items, test_items_3, 17) == 0);

    // delete multiple elements from the middle
    status = numba_list_delete_slice(lp, 2, 5, 1);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 14);
    test_items_3  = "\x00\x01\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f\x10";
    CHECK(memcmp(lp->items, test_items_3, 14) == 0);

    // delete single element from start
    status = numba_list_delete_slice(lp, 0, 1, 1);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 13);
    test_items_3  = "\x01\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f\x10";
    CHECK(memcmp(lp->items, test_items_3, 13) == 0);

    // delete single element from end
    status = numba_list_delete_slice(lp, 12, 13, 1);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 12);
    test_items_3  = "\x01\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f";
    CHECK(memcmp(lp->items, test_items_3, 12) == 0);

    // delete single element from middle
    status = numba_list_delete_slice(lp, 4, 5, 1);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 11);
    test_items_3  = "\x01\x05\x06\x07\x09\x0a\x0b\x0c\x0d\x0e\x0f";
    CHECK(memcmp(lp->items, test_items_3, 11) == 0);

    // delete all elements except first and last
    status = numba_list_delete_slice(lp, 1, 10, 1);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 2);
    test_items_3  = "\x01\x0f";
    CHECK(memcmp(lp->items, test_items_3, 2) == 0);

    // delete all remaining elements
    status = numba_list_delete_slice(lp, 0, lp->size, 1);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 0);
    test_items_3  = "";
    CHECK(memcmp(lp->items, test_items_3, 0) == 0);

    // free existing list
    numba_list_free(lp);

    // Setup list for testing delete_slice with non unary step
    status = numba_list_new(&lp, 1, 0);
    CHECK(status == LIST_OK);
    CHECK(lp->item_size == 1);
    CHECK(lp->size == 0);
    CHECK(lp->allocated == 0);
    for (i = 0; i < 17 ; i++) {
        status = numba_list_append(lp, (const char*)&i);
        CHECK(status == LIST_OK);
    }
    CHECK(lp->size == 17);

    // delete all items with odd index
    status = numba_list_delete_slice(lp, 0, 17, 2);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 8);
    test_items_3 = "\x01\x03\x05\x07\x09\x0b\x0d\x0f";
    CHECK(memcmp(lp->items, test_items_3, 8) == 0);

    // delete with a step of 4, starting at index 1
    status = numba_list_delete_slice(lp, 1, 8, 4);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 6);
    test_items_3 = "\x01\x05\x07\x09\x0d\x0f";
    CHECK(memcmp(lp->items, test_items_3, 6) == 0);

    // delete with a step of 2, but finish before end of list
    status = numba_list_delete_slice(lp, 0, 4, 2);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 4);
    test_items_3 = "\x05\x09\x0d\x0f";
    CHECK(memcmp(lp->items, test_items_3, 4) == 0);

    // no-op on empty slice
    status = numba_list_delete_slice(lp, 0, 0, 1);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 4);
    test_items_3 = "\x05\x09\x0d\x0f";
    CHECK(memcmp(lp->items, test_items_3, 4) == 0);

    // no-op on empty slice, non-zero index
    status = numba_list_delete_slice(lp, 2, 2, 1);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 4);
    test_items_3 = "\x05\x09\x0d\x0f";
    CHECK(memcmp(lp->items, test_items_3, 4) == 0);

    // free list and return 0
    numba_list_free(lp);

    // Setup list for testing delete_slice with negative step
    status = numba_list_new(&lp, 1, 0);
    CHECK(status == LIST_OK);
    CHECK(lp->item_size == 1);
    CHECK(lp->size == 0);
    CHECK(lp->allocated == 0);
    for (i = 0; i < 17 ; i++) {
        status = numba_list_append(lp, (const char*)&i);
        CHECK(status == LIST_OK);
    }
    CHECK(lp->size == 17);

    // delete all items using unary negative slice
    status = numba_list_delete_slice(lp, 16, -1, -1);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 0);

    // refill list
    for (i = 0; i < 17 ; i++) {
        status = numba_list_append(lp, (const char*)&i);
        CHECK(status == LIST_OK);
    }

    // delete all items using unary negative slice
    // need to start at index of last item (16) and
    // go beyond first item, i.e. -1 in Cd
    status = numba_list_delete_slice(lp, 16, -1, -2);
    CHECK(status == LIST_OK);
    CHECK(lp->size == 8);
    test_items_3 = "\x01\x03\x05\x07\x09\x0b\x0d\x0f";
    CHECK(memcmp(lp->items, test_items_3, 8) == 0);

    // free list and return 0
    numba_list_free(lp);
    return 0;


}

#undef CHECK
