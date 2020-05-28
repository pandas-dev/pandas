from numpy cimport int64_t

# Exposed for tslib, not intended for outside use.
cpdef int64_t delta_to_nanoseconds(delta) except? -1
cdef convert_to_timedelta64(object ts, object unit)
cdef bint is_any_td_scalar(object obj)
