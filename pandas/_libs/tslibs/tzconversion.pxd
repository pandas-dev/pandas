from cpython.datetime cimport tzinfo
from numpy cimport int64_t


cdef int64_t _tz_convert_tzlocal_utc(int64_t val, tzinfo tz, bint to_utc=*)
