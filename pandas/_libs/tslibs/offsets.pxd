from numpy cimport int64_t


cpdef to_offset(object obj, bint is_period=*)
cdef bint is_offset_object(object obj)
cdef bint is_tick_object(object obj)

cdef class DateOffset:
    cdef readonly:
        int64_t _n
        bint _normalize
        dict _cache

    # Needed for the relativedelta keywords passed to DateOffset itself,
    #  which are set with object.__setattr__ in __init__.
    cdef dict __dict__
