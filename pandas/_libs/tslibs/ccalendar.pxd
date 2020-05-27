from cython cimport Py_ssize_t

from numpy cimport int64_t, int32_t

ctypedef (int32_t, int32_t, int32_t) iso_calendar_t

cdef int dayofweek(int y, int m, int d) nogil
cdef bint is_leapyear(int64_t year) nogil
cpdef int32_t get_days_in_month(int year, Py_ssize_t month) nogil
cpdef int32_t get_week_of_year(int year, int month, int day) nogil
cpdef iso_calendar_t get_iso_calendar(int year, int month, int day) nogil
cpdef int32_t get_day_of_year(int year, int month, int day) nogil

cdef int64_t DAY_NANOS
cdef int64_t HOUR_NANOS
cdef dict c_MONTH_NUMBERS
