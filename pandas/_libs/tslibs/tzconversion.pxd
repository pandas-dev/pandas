from cpython.datetime cimport tzinfo
from numpy cimport int64_t
from pandas._libs.tslibs.conversion cimport _TSObject


cdef int64_t tz_convert_utc_to_tzlocal(int64_t utc_val, tzinfo tz)
cdef int64_t _tz_convert_tzlocal_utc(int64_t val, tzinfo tz, bint to_utc=*)
cdef void _tz_convert_utctsobject_to_local(_TSObject obj, tzinfo tz)
cpdef int64_t tz_convert_single(int64_t val, object tz1, object tz2)
