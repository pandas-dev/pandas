from cython cimport Py_ssize_t
from numpy cimport (
    int32_t,
    int64_t,
)

ctypedef (int64_t, int32_t, int32_t) iso_calendar_t

cdef int dayofweek(int64_t y, int m, int d) noexcept nogil
cdef bint is_leapyear(int64_t year) noexcept nogil
cpdef int32_t get_days_in_month(int64_t year, Py_ssize_t month) noexcept nogil
cpdef int32_t get_week_of_year(int64_t year, int month, int day) noexcept nogil
cpdef iso_calendar_t get_iso_calendar(int64_t year, int month, int day) noexcept nogil
cpdef int32_t get_day_of_year(int64_t year, int month, int day) noexcept nogil
cpdef int get_lastbday(int64_t year, int month) noexcept nogil
cpdef int get_firstbday(int64_t year, int month) noexcept nogil

cdef dict c_MONTH_NUMBERS, MONTH_TO_CAL_NUM

cdef int32_t* month_offset
