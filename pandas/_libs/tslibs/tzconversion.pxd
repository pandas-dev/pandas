from cpython.datetime cimport tzinfo
from numpy cimport int64_t


cdef int64_t tz_convert_utc_to_tzlocal(int64_t utc_val, tzinfo tz, bint* fold=*)
cpdef int64_t tz_convert_single(int64_t val, tzinfo tz1, tzinfo tz2)
