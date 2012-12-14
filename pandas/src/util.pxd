from numpy cimport ndarray
cimport numpy as cnp
cimport cpython

cdef extern from "numpy_helper.h":
    inline void set_array_owndata(ndarray ao)
    inline void set_array_not_contiguous(ndarray ao)

    inline int is_integer_object(object)
    inline int is_float_object(object)
    inline int is_complex_object(object)
    inline int is_bool_object(object)
    inline int is_string_object(object)
    inline int is_datetime64_object(object)
    inline int is_timedelta64_object(object)
    inline int assign_value_1d(ndarray, Py_ssize_t, object) except -1
    inline cnp.int64_t get_nat()
    inline object get_value_1d(ndarray, Py_ssize_t)
    inline int floatify(object, double*) except -1
    inline char *get_c_string(object)
    inline object char_to_string(char*)
    inline void transfer_object_column(char *dst, char *src, size_t stride,
                                       size_t length)
    object sarr_from_data(cnp.dtype, int length, void* data)

cdef inline object get_value_at(ndarray arr, object loc):
    cdef:
        Py_ssize_t i, sz
        void* data_ptr
    if is_float_object(loc):
        casted = int(loc)
        if casted == loc:
            loc = casted
    i = <Py_ssize_t> loc
    sz = cnp.PyArray_SIZE(arr)

    if i < 0 and sz > 0:
        i += sz
    elif i >= sz or sz == 0:
        raise IndexError('index out of bounds')

    return get_value_1d(arr, i)

cdef inline set_value_at(ndarray arr, object loc, object value):
    cdef:
        Py_ssize_t i, sz
    if is_float_object(loc):
        casted = int(loc)
        if casted == loc:
            loc = casted
    i = <Py_ssize_t> loc
    sz = cnp.PyArray_SIZE(arr)

    if i < 0:
        i += sz
    elif i >= sz:
        raise IndexError('index out of bounds')

    assign_value_1d(arr, i, value)

cdef inline int is_contiguous(ndarray arr):
    return cnp.PyArray_CHKFLAGS(arr, cnp.NPY_C_CONTIGUOUS)

cdef inline is_array(object o):
    return cnp.PyArray_Check(o)


cdef inline bint _checknull(object val):
    try:
        return val is None or (cpython.PyFloat_Check(val) and val != val)
    except ValueError:
        return False

cdef inline bint _checknull_old(object val):
    import numpy as np
    cdef double INF = <double> np.inf
    cdef double NEGINF = -INF
    try:
        return bool(val is None or val != val and val != INF
                    and val != NEGINF)
    except ValueError:
        return False

cdef inline bint _checknan(object val):
    return not cnp.PyArray_Check(val) and val != val
