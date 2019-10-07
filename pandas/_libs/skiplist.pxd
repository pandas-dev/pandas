# -*- coding: utf-8 -*-
# See GH#27465 for reference on related-but-unused cython code

cdef extern from "src/skiplist.h":
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

    skiplist_t* skiplist_init(int) nogil
    void skiplist_destroy(skiplist_t*) nogil
    double skiplist_get(skiplist_t*, int, int*) nogil
    int skiplist_insert(skiplist_t*, double) nogil
    int skiplist_remove(skiplist_t*, double) nogil
