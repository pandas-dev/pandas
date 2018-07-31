
from cpython cimport PyTypeObject

cdef extern from *:
    """
    PyObject* char_to_string(const char* data) {
    #if PY_VERSION_HEX >= 0x03000000
        return PyUnicode_FromString(data);
    #else
        return PyString_FromString(data);
    #endif
    }
    """
    object char_to_string(const char* data)


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
    ctypedef signed long long int64_t
    int _import_array() except -1

cdef extern from "numpy/ndarrayobject.h":
    PyTypeObject PyTimedeltaArrType_Type
    PyTypeObject PyDatetimeArrType_Type
    PyTypeObject PyComplexFloatingArrType_Type
    PyTypeObject PyBoolArrType_Type

    bint PyArray_IsIntegerScalar(obj) nogil
    bint PyArray_Check(obj) nogil

cdef extern from  "numpy/npy_common.h":
    int64_t NPY_MIN_INT64


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


cdef inline int64_t get_nat():
    return NPY_MIN_INT64


cdef inline int import_array() except -1:
    _import_array()


# --------------------------------------------------------------------
# Type Checking

cdef inline bint is_string_object(object obj) nogil:
    """
    Cython equivalent of `isinstance(val, compat.string_types)`

    Parameters
    ----------
    val : object

    Returns
    -------
    is_string : bool
    """
    return PyString_Check(obj) or PyUnicode_Check(obj)


cdef inline bint is_integer_object(object obj) nogil:
    """
    Cython equivalent of

    `isinstance(val, (int, long, np.integer)) and not isinstance(val, bool)`

    Parameters
    ----------
    val : object

    Returns
    -------
    is_integer : bool

    Notes
    -----
    This counts np.timedelta64 objects as integers.
    """
    return not PyBool_Check(obj) and PyArray_IsIntegerScalar(obj)


cdef inline bint is_float_object(object obj) nogil:
    """
    Cython equivalent of `isinstance(val, (float, np.complex_))`

    Parameters
    ----------
    val : object

    Returns
    -------
    is_float : bool
    """
    return (PyFloat_Check(obj) or
            (PyObject_TypeCheck(obj, &PyFloatingArrType_Type)))


cdef inline bint is_complex_object(object obj) nogil:
    """
    Cython equivalent of `isinstance(val, (complex, np.complex_))`

    Parameters
    ----------
    val : object

    Returns
    -------
    is_complex : bool
    """
    return (PyComplex_Check(obj) or
            PyObject_TypeCheck(obj, &PyComplexFloatingArrType_Type))


cdef inline bint is_bool_object(object obj) nogil:
    """
    Cython equivalent of `isinstance(val, (bool, np.bool_))`

    Parameters
    ----------
    val : object

    Returns
    -------
    is_bool : bool
    """
    return (PyBool_Check(obj) or
            PyObject_TypeCheck(obj, &PyBoolArrType_Type))


cdef inline bint is_timedelta64_object(object obj) nogil:
    """
    Cython equivalent of `isinstance(val, np.timedelta64)`

    Parameters
    ----------
    val : object

    Returns
    -------
    is_timedelta64 : bool
    """
    return PyObject_TypeCheck(obj, &PyTimedeltaArrType_Type)


cdef inline bint is_datetime64_object(object obj) nogil:
    """
    Cython equivalent of `isinstance(val, np.datetime64)`

    Parameters
    ----------
    val : object

    Returns
    -------
    is_datetime64 : bool
    """
    return PyObject_TypeCheck(obj, &PyDatetimeArrType_Type)


cdef inline bint is_array(object val):
    """
    Cython equivalent of `isinstance(val, np.ndarray)`

    Parameters
    ----------
    val : object

    Returns
    -------
    is_ndarray : bool
    """
    return PyArray_Check(val)


cdef inline bint is_period_object(object val):
    """
    Cython equivalent of `isinstance(val, pd.Period)`

    Parameters
    ----------
    val : object

    Returns
    -------
    is_period : bool
    """
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


cdef inline bint is_nan(object val):
    """
    Check if val is a Not-A-Number float, including float('NaN') and np.nan.

    Parameters
    ----------
    val : object

    Returns
    -------
    is_nan : bool
    """
    return isinstance(val, float) and val != val
