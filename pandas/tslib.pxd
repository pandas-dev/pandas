from numpy cimport ndarray, int64_t

cdef convert_to_tsobject(object, object, object, bint, bint)
cpdef convert_to_timedelta64(object, object)
cpdef object maybe_get_tz(object)
cdef bint _is_utc(object)
cdef bint _is_tzlocal(object)
cdef object _get_dst_info(object)
cdef bint _nat_scalar_rules[6]
cdef bint _check_all_nulls(obj)
