# -*- coding: utf-8 -*-
# cython: profile=False

from cython cimport Py_ssize_t

from numpy cimport double_t


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


# Note: Node is declared here so that IndexableSkiplist can be exposed;
# Node itself not intended to be exposed.
cdef class Node:
    cdef public:
        double_t value
        list next
        list width


cdef class IndexableSkiplist:
    cdef:
        Py_ssize_t size, maxlevels
        Node head

    cpdef get(self, Py_ssize_t i)
    cpdef insert(self, double value)
    cpdef remove(self, double value)
