from cpython.datetime cimport tzinfo
from numpy cimport int64_t


cdef int64_t tz_convert_utc_to_tzlocal(int64_t utc_val, tzinfo tz)
cdef int64_t _tz_convert_tzlocal_fromutc(int64_t val, tzinfo tz, bint *fold)
cpdef int64_t tz_convert_single(int64_t val, object tz1, object tz2)


cdef class Localizer:
    cdef readonly:
        tzinfo tz

    cdef int64_t localize(self, int64_t value, Py_ssize_t i)
    cdef initialize(self, const int64_t[:] values)

cdef Localizer get_localizer(tz)
