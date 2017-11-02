from numpy cimport ndarray, int64_t

from tslibs.conversion cimport convert_to_tsobject

cpdef convert_to_timedelta64(object, object)

cdef _to_i8(object val)
