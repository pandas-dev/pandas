from numpy cimport ndarray
cimport numpy as cnp

cdef extern from "numpy_helper.h":
    inline int is_integer_object(object)
    inline int is_float_object(object)
    inline int is_bool_object(object)
    inline int is_string_object(object)
    inline int assign_value_1d (ndarray, Py_ssize_t, object) except -1

cpdef inline object get_value_at(ndarray arr, object loc):
    cdef:
        Py_ssize_t i
        void* data_ptr
    if is_float_object(loc):
        casted = int(loc)
        if casted == loc:
            loc = casted
    i = <Py_ssize_t> loc
    if i < 0:
        i += cnp.PyArray_SIZE(arr)
    data_ptr = cnp.PyArray_GETPTR1(arr, i)
    return cnp.PyArray_GETITEM(arr, data_ptr)

cpdef inline set_value_at(ndarray arr, object loc, object value):
    cdef:
        Py_ssize_t i
    if is_float_object(loc):
        casted = int(loc)
        if casted == loc:
            loc = casted
    i = <Py_ssize_t> loc
    if i < 0:
        i += cnp.PyArray_SIZE(arr)

    assign_value_1d(arr, i, value)
