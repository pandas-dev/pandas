from numpy cimport ndarray

cdef extern from "numpy_helper.h":
    inline int is_integer_object(object)
    inline int is_float_object(object)
    inline int assign_value_1d (ndarray, Py_ssize_t, object) except -1
