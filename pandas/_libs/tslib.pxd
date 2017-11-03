from numpy cimport ndarray, int64_t

from tslibs.conversion cimport convert_to_tsobject
from tslibs.timedeltas cimport convert_to_timedelta64

cdef bint _check_all_nulls(obj)

cdef _to_i8(object val)
