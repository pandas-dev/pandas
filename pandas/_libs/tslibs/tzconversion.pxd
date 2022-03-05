from cpython.datetime cimport tzinfo
from numpy cimport (
    int64_t,
    intp_t,
)

from pandas._libs.tslibs.timezones cimport Localizer


cdef int64_t tz_convert_utc_to_tzlocal(
    int64_t utc_val, tzinfo tz, bint* fold=*
) except? -1
cpdef int64_t tz_convert_from_utc_single(int64_t val, tzinfo tz)
cdef int64_t tz_localize_to_utc_single(
    int64_t val, tzinfo tz, object ambiguous=*, object nonexistent=*
) except? -1


cdef int64_t utc_val_to_local_val(Localizer info, int64_t utc_val, intp_t[:] pos, Py_ssize_t i)
