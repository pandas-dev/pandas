from cpython.datetime cimport (
    datetime,
    tzinfo,
)
from numpy cimport int64_t

from pandas._libs.tslibs.np_datetime cimport NPY_DATETIMEUNIT


cdef bint parse_today_now(str val, int64_t* iresult, bint utc)


cdef class DatetimeParseState:
    cdef:
        bint found_tz
        bint found_naive
        bint creso_changed
        NPY_DATETIMEUNIT creso

    cdef bint update_creso(self, NPY_DATETIMEUNIT creso)
    cdef tzinfo process_datetime(self, datetime dt, tzinfo tz, bint utc_convert)
