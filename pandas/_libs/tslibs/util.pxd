cdef extern from "Python.h":
    # Note that following functions can potentially raise an exception,
    # thus they cannot be declared 'nogil'. Also PyUnicode_AsUTF8AndSize() can
    # potentially allocate memory inside in unlikely case of when underlying
    # unicode object was stored as non-utf8 and utf8 wasn't requested before.
    const char* PyUnicode_AsUTF8AndSize(object obj,
                                        Py_ssize_t* length) except NULL

    object PyUnicode_EncodeLocale(object obj, const char *errors) nogil
    object PyUnicode_DecodeLocale(const char *str, const char *errors) nogil


from numpy cimport int64_t


cdef extern from "pandas/type.h":
    bint is_timedelta64_object(object obj)
    bint is_integer_object(object obj)
    bint is_float_object(object obj)
    bint is_complex_object(object obj)
    bint is_bool_object(object obj)
    bint is_real_number_object(object obj)
    bint is_datetime64_object(object obj)
    bint is_array(object obj)
    bint is_nan(object obj)

cdef extern from "numpy/npy_common.h":
    int64_t NPY_MIN_INT64


cdef inline int64_t get_nat():
    return NPY_MIN_INT64


# --------------------------------------------------------------------
# Type Checking

cdef inline const char* get_c_string_buf_and_size(str py_string,
                                                  Py_ssize_t *length) except NULL:
    """
    Extract internal char* buffer of unicode or bytes object `py_string` with
    getting length of this internal buffer saved in `length`.

    Notes
    -----
    Python object owns memory, thus returned char* must not be freed.
    `length` can be NULL if getting buffer length is not needed.

    Parameters
    ----------
    py_string : str
    length : Py_ssize_t*

    Returns
    -------
    buf : const char*
    """
    return PyUnicode_AsUTF8AndSize(py_string, length)


cdef inline const char* get_c_string(str py_string) except NULL:
    return get_c_string_buf_and_size(py_string, NULL)


cdef inline bytes string_encode_locale(str py_string):
    """As opposed to PyUnicode_Encode, use current system locale to encode."""
    return PyUnicode_EncodeLocale(py_string, NULL)


cdef inline object char_to_string_locale(const char* data):
    """As opposed to PyUnicode_FromString, use current system locale to decode."""
    return PyUnicode_DecodeLocale(data, NULL)
