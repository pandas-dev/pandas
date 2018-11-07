# -*- coding: utf-8 -*-

from numpy cimport int64_t
cdef int64_t NPY_NAT

cdef bint _nat_scalar_rules[6]

cdef bint checknull_with_nat(object val)
cdef bint is_null_datetimelike(object val)

from cpython.datetime cimport datetime


cdef class _NaT(datetime):
    cdef readonly:
        int64_t value
        object freq

# By declaring NAT in the .pxd file and cimporting it into other cython
#  modules, we make `if thing is NaT` into a C pointer lookup and avoid a
#  runtime lookup of NaT in the module namespace.
cdef _NaT NAT
