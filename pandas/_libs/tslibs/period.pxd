from numpy cimport int64_t

from .np_datetime cimport (
    NPY_DATETIMEUNIT,
    npy_datetimestruct,
)


cdef bint is_period_object(object obj)
cdef int64_t get_period_ordinal_unchecked(
    npy_datetimestruct *dts, int freq
) noexcept nogil
cdef bint get_period_bounds(
    int freq, NPY_DATETIMEUNIT* unit, int* min_year, int* max_year
) noexcept nogil
