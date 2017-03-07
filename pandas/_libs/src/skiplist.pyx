# Cython version of IndexableSkiplist, for implementing moving median
# with O(log n) updates
# Original author: Raymond Hettinger
# Original license: MIT
# Link: http://code.activestate.com/recipes/576930/

# Cython version: Wes McKinney

cdef extern from "numpy/arrayobject.h":

    void import_array()

cdef extern from "math.h":
    double log(double x)

# MSVC does not have log2!

cdef double Log2(double x):
    return log(x) / log(2.)

cimport numpy as np
from numpy cimport *
import numpy as np

from random import random

# initialize numpy
import_array()

# TODO: optimize this, make less messy

cdef class Node:
    cdef public:
        double_t value
        list next
        list width

    def __init__(self, double_t value, list next, list width):
        self.value = value
        self.next = next
        self.width = width

# Singleton terminator node
NIL = Node(np.inf, [], [])

cdef class IndexableSkiplist:
    """
    Sorted collection supporting O(lg n) insertion, removal, and
    lookup by rank.
    """
    cdef:
        Py_ssize_t size, maxlevels
        Node head

    def __init__(self, expected_size=100):
        self.size = 0
        self.maxlevels = int(1 + Log2(expected_size))
        self.head = Node(np.NaN, [NIL] * self.maxlevels, [1] * self.maxlevels)

    def __len__(self):
        return self.size

    def __getitem__(self, i):
        return self.get(i)

    cpdef get(self, Py_ssize_t i):
        cdef Py_ssize_t level
        cdef Node node

        node = self.head
        i += 1

        for level in range(self.maxlevels - 1, -1, -1):
            while node.width[level] <= i:
                i -= node.width[level]
                node = node.next[level]

        return node.value

    cpdef insert(self, double value):
        cdef Py_ssize_t level, steps, d
        cdef Node node, prevnode, newnode, next_at_level, tmp
        cdef list chain, steps_at_level

        # find first node on each level where node.next[levels].value > value
        chain = [None] * self.maxlevels
        steps_at_level = [0] * self.maxlevels
        node = self.head

        for level in range(self.maxlevels - 1, -1, -1):
            next_at_level = node.next[level]

            while next_at_level.value <= value:
                steps_at_level[level] = (steps_at_level[level] +
                                         node.width[level])
                node = next_at_level
                next_at_level = node.next[level]

            chain[level] = node

        # insert a link to the newnode at each level
        d = min(self.maxlevels, 1 - int(Log2(random())))
        newnode = Node(value, [None] * d, [None] * d)
        steps = 0

        for level in range(d):
            prevnode = chain[level]
            newnode.next[level] = prevnode.next[level]
            prevnode.next[level] = newnode
            newnode.width[level] = (prevnode.width[level] - steps)
            prevnode.width[level] = steps + 1
            steps += steps_at_level[level]

        for level in range(d, self.maxlevels):
            (<Node> chain[level]).width[level] += 1

        self.size += 1

    cpdef remove(self, double value):
        cdef Py_ssize_t level, d
        cdef Node node, prevnode, tmpnode, next_at_level
        cdef list chain

        # find first node on each level where node.next[levels].value >= value
        chain = [None] * self.maxlevels
        node = self.head

        for level in range(self.maxlevels - 1, -1, -1):
            next_at_level = node.next[level]
            while next_at_level.value < value:
                node = next_at_level
                next_at_level = node.next[level]

            chain[level] = node

        if value != (<Node> (<Node> (<Node> chain[0]).next)[0]).value:
            raise KeyError('Not Found')

        # remove one link at each level
        d = len((<Node> (<Node> (<Node> chain[0]).next)[0]).next)

        for level in range(d):
            prevnode = chain[level]
            tmpnode = prevnode.next[level]
            prevnode.width[level] += tmpnode.width[level] - 1
            prevnode.next[level] = tmpnode.next[level]

        for level in range(d, self.maxlevels):
            tmpnode = chain[level]
            tmpnode.width[level] -= 1

        self.size -= 1
