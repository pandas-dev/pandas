from tslibs.util cimport *

from cython cimport Py_ssize_t

cimport numpy as cnp
from numpy cimport ndarray

cdef extern from "numpy/ndarraytypes.h":
    void PyArray_CLEARFLAGS(ndarray arr, int flags) nogil


cdef extern from "numpy/arrayobject.h":
    enum:
        NPY_ARRAY_C_CONTIGUOUS
        NPY_ARRAY_F_CONTIGUOUS


cdef extern from *:
    """
    // returns ASCII or UTF8 (py3) view on python str
    // python object owns memory, should not be freed
    static const char* get_c_string(PyObject* obj) {
    #if PY_VERSION_HEX >= 0x03000000
        return PyUnicode_AsUTF8(obj);
    #else
        return PyString_AsString(obj);
    #endif
    }
    """
    const char *get_c_string(object) except NULL


cdef extern from "src/numpy_helper.h":
    int assign_value_1d(ndarray, Py_ssize_t, object) except -1
    object get_value_1d(ndarray, Py_ssize_t)


cdef extern from "src/headers/stdint.h":
    enum: UINT8_MAX
    enum: UINT16_MAX
    enum: UINT32_MAX
    enum: UINT64_MAX
    enum: INT8_MIN
    enum: INT8_MAX
    enum: INT16_MIN
    enum: INT16_MAX
    enum: INT32_MAX
    enum: INT32_MIN
    enum: INT64_MAX
    enum: INT64_MIN


ctypedef fused numeric:
    cnp.int8_t
    cnp.int16_t
    cnp.int32_t
    cnp.int64_t

    cnp.uint8_t
    cnp.uint16_t
    cnp.uint32_t
    cnp.uint64_t

    cnp.float32_t
    cnp.float64_t


cdef inline void set_array_not_contiguous(ndarray ao) nogil:
    # Numpy>=1.8-compliant equivalent to:
    # ao->flags &= ~(NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_F_CONTIGUOUS);
    PyArray_CLEARFLAGS(ao,
                       (NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_F_CONTIGUOUS))


cdef inline Py_ssize_t validate_indexer(ndarray arr, object loc) except -1:
    """
    Cast the given indexer `loc` to an integer.  If it is negative, i.e. a
    python-style indexing-from-the-end indexer, translate it to a
    from-the-front indexer.  Raise if this is not possible.

    Parameters
    ----------
    arr : ndarray
    loc : object

    Returns
    -------
    idx : Py_ssize_t

    Raises
    ------
    IndexError
    """
    cdef:
        Py_ssize_t idx, size
        int casted

    if is_float_object(loc):
        casted = int(loc)
        if casted == loc:
            loc = casted

    idx = <Py_ssize_t>loc
    size = cnp.PyArray_SIZE(arr)

    if idx < 0 and size > 0:
        idx += size
    if idx >= size or size == 0 or idx < 0:
        raise IndexError('index out of bounds')

    return idx


cdef inline object get_value_at(ndarray arr, object loc):
    cdef:
        Py_ssize_t i

    i = validate_indexer(arr, loc)
    return get_value_1d(arr, i)


cdef inline set_value_at_unsafe(ndarray arr, object loc, object value):
    """Sets a value into the array without checking the writeable flag.

    This should be used when setting values in a loop, check the writeable
    flag above the loop and then eschew the check on each iteration.
    """
    cdef:
        Py_ssize_t i

    i = validate_indexer(arr, loc)
    assign_value_1d(arr, i, value)


cdef inline set_value_at(ndarray arr, object loc, object value):
    """Sets a value into the array after checking that the array is mutable.
    """
    if not cnp.PyArray_ISWRITEABLE(arr):
        raise ValueError('assignment destination is read-only')

    set_value_at_unsafe(arr, loc, value)
