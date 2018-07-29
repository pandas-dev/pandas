# -*- coding: utf-8 -*-
# cython: profile=False

from numpy cimport int64_t
cdef int64_t NPY_NAT

cdef bint _nat_scalar_rules[6]

cdef bint checknull_with_nat(object val)
cdef bint is_null_datetimelike(object val)
