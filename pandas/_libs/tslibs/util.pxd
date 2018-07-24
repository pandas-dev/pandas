from numpy cimport ndarray
cimport numpy as cnp
cnp.import_array()

cimport cpython
from cpython cimport PyTypeObject


cdef extern from "Python.h":
    # Note: importing extern-style allows us to declare these as nogil
    # functions, whereas `from cpython cimport` does not.
    bint PyUnicode_Check(object obj) nogil
    bint PyString_Check(object obj) nogil
    bint PyBool_Check(object obj) nogil
    bint PyFloat_Check(object obj) nogil
    bint PyComplex_Check(object obj) nogil
    bint PyObject_TypeCheck(object obj, PyTypeObject* type) nogil


cdef extern from "numpy/arrayobject.h":
    PyTypeObject PyFloatingArrType_Type

cdef extern from "numpy/ndarrayobject.h":
    PyTypeObject PyTimedeltaArrType_Type
    PyTypeObject PyDatetimeArrType_Type
    PyTypeObject PyComplexFloatingArrType_Type
    PyTypeObject PyBoolArrType_Type

    bint PyArray_IsIntegerScalar(obj) nogil
    bint PyArray_Check(obj) nogil

# --------------------------------------------------------------------
# Type Checking

cdef inline bint is_string_object(object obj) nogil:
    return PyString_Check(obj) or PyUnicode_Check(obj)


cdef inline bint is_integer_object(object obj) nogil:
    return not PyBool_Check(obj) and PyArray_IsIntegerScalar(obj)


cdef inline bint is_float_object(object obj) nogil:
    return (PyFloat_Check(obj) or
            (PyObject_TypeCheck(obj, &PyFloatingArrType_Type)))


cdef inline bint is_complex_object(object obj) nogil:
    return (PyComplex_Check(obj) or
            PyObject_TypeCheck(obj, &PyComplexFloatingArrType_Type))


cdef inline bint is_bool_object(object obj) nogil:
    return (PyBool_Check(obj) or
            PyObject_TypeCheck(obj, &PyBoolArrType_Type))


cdef inline bint is_timedelta64_object(object obj) nogil:
    return PyObject_TypeCheck(obj, &PyTimedeltaArrType_Type)


cdef inline bint is_datetime64_object(object obj) nogil:
    return PyObject_TypeCheck(obj, &PyDatetimeArrType_Type)

# --------------------------------------------------------------------

cdef extern from "../src/numpy_helper.h":
    void set_array_not_contiguous(ndarray ao)

    int assign_value_1d(ndarray, Py_ssize_t, object) except -1
    cnp.int64_t get_nat()
    object get_value_1d(ndarray, Py_ssize_t)
    const char *get_c_string(object) except NULL
    object char_to_string(char*)

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

cdef extern from "../src/headers/stdint.h":
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


cdef inline object get_value_at(ndarray arr, object loc):
    cdef:
        Py_ssize_t i, sz
        int casted

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


cdef inline set_value_at_unsafe(ndarray arr, object loc, object value):
    """Sets a value into the array without checking the writeable flag.

    This should be used when setting values in a loop, check the writeable
    flag above the loop and then eschew the check on each iteration.
    """
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

cdef inline set_value_at(ndarray arr, object loc, object value):
    """Sets a value into the array after checking that the array is mutable.
    """
    if not cnp.PyArray_ISWRITEABLE(arr):
        raise ValueError('assignment destination is read-only')

    set_value_at_unsafe(arr, loc, value)


cdef inline is_array(object o):
    return cnp.PyArray_Check(o)


cdef inline bint _checknull(object val):
    try:
        return val is None or (cpython.PyFloat_Check(val) and val != val)
    except ValueError:
        return False


cdef inline bint is_period_object(object val):
    return getattr(val, '_typ', '_typ') == 'period'


cdef inline bint is_offset_object(object val):
    """
    Check if an object is a DateOffset object.

    Parameters
    ----------
    val : object

    Returns
    -------
    is_date_offset : bool
    """
    return getattr(val, '_typ', None) == "dateoffset"
