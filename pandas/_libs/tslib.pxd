from numpy cimport ndarray, int64_t

cdef convert_to_tsobject(object, object, object, bint, bint)
cpdef convert_to_timedelta64(object, object)
cdef bint _nat_scalar_rules[6]
cdef bint _check_all_nulls(obj)
