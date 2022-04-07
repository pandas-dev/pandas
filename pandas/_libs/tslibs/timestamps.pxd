from cpython.datetime cimport (
    datetime,
    tzinfo,
)
from numpy cimport int64_t

from pandas._libs.tslibs.base cimport ABCTimestamp
from pandas._libs.tslibs.np_datetime cimport npy_datetimestruct
from pandas._libs.tslibs.offsets cimport BaseOffset


cdef _Timestamp create_timestamp_from_ts(int64_t value,
                                         npy_datetimestruct dts,
                                         tzinfo tz, BaseOffset freq, bint fold)


cdef class _Timestamp(ABCTimestamp):
    cdef readonly:
        int64_t value, nanosecond
        BaseOffset _freq

    cdef bint _get_start_end_field(self, str field, freq)
    cdef _get_date_name_field(self, str field, object locale)
    cdef int64_t _maybe_convert_value_to_local(self)
    cdef bint _can_compare(self, datetime other)
    cpdef to_datetime64(self)
    cpdef datetime to_pydatetime(_Timestamp self, bint warn=*)
    cdef bint _compare_outside_nanorange(_Timestamp self, datetime other,
                                         int op) except -1
    cpdef void _set_freq(self, freq)
    cdef _warn_on_field_deprecation(_Timestamp self, freq, str field)

cdef int64_t normalize_i8_stamp(int64_t local_val) nogil
