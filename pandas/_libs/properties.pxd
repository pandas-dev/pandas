#!/usr/bin/env python
# -*- coding: utf-8 -*-
# cython: profile=False
# cython: linetrace=False
# distutils: define_macros=CYTHON_TRACE=0
# distutils: define_macros=CYTHON_TRACE_NOGIL=0

from cython cimport Py_ssize_t


cdef class cache_readonly(object):
    cdef readonly:
        object func, name, allow_setting


cdef class AxisProperty(object):
    cdef:
        Py_ssize_t axis
