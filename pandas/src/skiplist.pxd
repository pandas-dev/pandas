cdef extern from "skiplist.h":
    ctypedef struct node_t:
        node_t **next
        int *width
        double value
        int is_nil
        int levels
        int ref_count

    ctypedef struct skiplist_t:
        node_t *head
        node_t **tmp_chain
        int *tmp_steps
        int size
        int maxlevels

    inline skiplist_t* skiplist_init(int) nogil
    inline void skiplist_destroy(skiplist_t*) nogil
    inline double skiplist_get(skiplist_t*, int, int*) nogil
    inline int skiplist_insert(skiplist_t*, double) nogil
    inline int skiplist_remove(skiplist_t*, double) nogil

