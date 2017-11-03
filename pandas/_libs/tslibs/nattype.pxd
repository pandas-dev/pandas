# -*- coding: utf-8 -*-
# cython: profile=False

from numpy cimport int64_t
cdef int64_t NPY_NAT

cdef bint _nat_scalar_rules[6]
