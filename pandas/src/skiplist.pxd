cdef extern from "skiplist.h":
    ctypedef struct node_t:
        double value
        int is_nil
        int levels
        node_t **next
        int *width
        int ref_count

    ctypedef struct skiplist_t:
        node_t *head
        int size, maxlevels
        node_t **tmp_chain
        int *tmp_steps

    inline skiplist_t* skiplist_init(int)
    inline void skiplist_destroy(skiplist_t*)
    inline double skiplist_get(skiplist_t*, int, int*)
    inline int skiplist_insert(skiplist_t*, double)
    inline int skiplist_remove(skiplist_t*, double)

