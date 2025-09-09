
from cpython.object cimport PyTypeObject


cdef extern from "Python.h":
    # Note: importing extern-style allows us to declare these as nogil
    # functions, whereas `from cpython cimport` does not.
    bint PyBool_Check(object obj) nogil
    bint PyFloat_Check(object obj) nogil
    bint PyComplex_Check(object obj) nogil
    bint PyObject_TypeCheck(object obj, PyTypeObject* type) nogil

    # Note that following functions can potentially raise an exception,
    # thus they cannot be declared 'nogil'.
    object PyUnicode_EncodeLocale(object obj, const char *errors) nogil
    object PyUnicode_DecodeLocale(const char *str, const char *errors) nogil


cimport numpy as cnp
from numpy cimport (
    PyArray_Check,
    float64_t,
    int64_t,
    is_timedelta64_object,
)


cdef extern from "numpy/arrayobject.h":
    PyTypeObject PyFloatingArrType_Type

cdef extern from "numpy/ndarrayobject.h":
    PyTypeObject PyComplexFloatingArrType_Type
    PyTypeObject PyBoolArrType_Type

    bint PyArray_IsIntegerScalar(obj) nogil

cdef extern from "numpy/npy_common.h":
    int64_t NPY_MIN_INT64


cdef inline int64_t get_nat() noexcept:
    return NPY_MIN_INT64


# --------------------------------------------------------------------
# Type Checking

cdef inline bint is_integer_object(object obj) noexcept:
    """
    Cython equivalent of

    `isinstance(val, (int, np.integer)) and not isinstance(val, (bool, np.timedelta64))`

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
    return (not PyBool_Check(obj) and isinstance(obj, (int, cnp.integer))
            and not is_timedelta64_object(obj))


cdef inline bint is_float_object(object obj) noexcept nogil:
    """
    Cython equivalent of `isinstance(val, (float, np.floating))`

    Parameters
    ----------
    val : object

    Returns
    -------
    is_float : bool
    """
    return (PyFloat_Check(obj) or
            (PyObject_TypeCheck(obj, &PyFloatingArrType_Type)))


cdef inline bint is_complex_object(object obj) noexcept nogil:
    """
    Cython equivalent of `isinstance(val, (complex, np.complexfloating))`

    Parameters
    ----------
    val : object

    Returns
    -------
    is_complex : bool
    """
    return (PyComplex_Check(obj) or
            PyObject_TypeCheck(obj, &PyComplexFloatingArrType_Type))


cdef inline bint is_bool_object(object obj) noexcept nogil:
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


cdef inline bint is_real_number_object(object obj) noexcept:
    return is_bool_object(obj) or is_integer_object(obj) or is_float_object(obj)


cdef inline bint is_array(object val) noexcept:
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


cdef inline bint is_nan(object val):
    """
    Check if val is a Not-A-Number float or complex, including
    float('NaN') and np.nan.

    Parameters
    ----------
    val : object

    Returns
    -------
    is_nan : bool
    """
    cdef float64_t fval
    if is_float_object(val):
        fval = val
        return fval != fval
    return is_complex_object(val) and val != val


cdef inline bytes string_encode_locale(str py_string):
    """As opposed to PyUnicode_Encode, use current system locale to encode."""
    return PyUnicode_EncodeLocale(py_string, NULL)


cdef inline object char_to_string_locale(const char* data):
    """As opposed to PyUnicode_FromString, use current system locale to decode."""
    return PyUnicode_DecodeLocale(data, NULL)
