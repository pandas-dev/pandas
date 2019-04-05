# -*- coding: utf-8 -*-

from cpython.datetime cimport datetime

from numpy cimport int64_t

cdef class _Timestamp(datetime):
    cdef readonly:
        int64_t value, nanosecond
        object freq
        list _date_attributes
    cpdef bint _get_start_end_field(self, str field)
    cpdef _get_date_name_field(self, object field, object locale)
    cdef int64_t _maybe_convert_value_to_local(self)
    cpdef to_datetime64(self)
    cdef _assert_tzawareness_compat(_Timestamp self, datetime other)
    cpdef datetime to_pydatetime(_Timestamp self, bint warn=*)
    cdef bint _compare_outside_nanorange(_Timestamp self, datetime other,
                                         int op) except -1
