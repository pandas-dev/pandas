"""
The following are faster versions of struct.unpack that avoid the overhead of Python
function calls.

In the SAS7BDAT parser, they may be called up to (n_rows * n_cols) times.
"""
from cython cimport Py_ssize_t
from libc.stdint cimport (
    uint16_t,
    uint32_t,
    uint64_t,
)


def read_float_with_byteswap(bytes data, Py_ssize_t offset, bint byteswap):
    assert offset + 4 < len(data)
    cdef:
        const char *data_ptr = data
        float res = (<float*>(data_ptr + offset))[0]
    if byteswap:
        res = _byteswap_float(res)
    return res


def read_double_with_byteswap(bytes data, Py_ssize_t offset, bint byteswap):
    assert offset + 8 < len(data)
    cdef:
        const char *data_ptr = data
        double res = (<double*>(data_ptr + offset))[0]
    if byteswap:
        res = _byteswap_double(res)
    return res


def read_uint16_with_byteswap(bytes data, Py_ssize_t offset, bint byteswap):
    assert offset + 2 < len(data)
    cdef:
        const char *data_ptr = data
        uint16_t res = (<uint16_t *>(data_ptr + offset))[0]
    if byteswap:
        res = _byteswap2(res)
    return res


def read_uint32_with_byteswap(bytes data, Py_ssize_t offset, bint byteswap):
    assert offset + 4 < len(data)
    cdef:
        const char *data_ptr = data
        uint32_t res = (<uint32_t *>(data_ptr + offset))[0]
    if byteswap:
        res = _byteswap4(res)
    return res


def read_uint64_with_byteswap(bytes data, Py_ssize_t offset, bint byteswap):
    assert offset + 8 < len(data)
    cdef:
        const char *data_ptr = data
        uint64_t res = (<uint64_t *>(data_ptr + offset))[0]
    if byteswap:
        res = _byteswap8(res)
    return res


# Byteswapping

cdef extern from *:
    """
    #ifdef _MSC_VER
        #define _byteswap2 _byteswap_ushort
        #define _byteswap4 _byteswap_ulong
        #define _byteswap8 _byteswap_uint64
    #else
        #define _byteswap2 __builtin_bswap16
        #define _byteswap4 __builtin_bswap32
        #define _byteswap8 __builtin_bswap64
    #endif
    """
    uint16_t _byteswap2(uint16_t)
    uint32_t _byteswap4(uint32_t)
    uint64_t _byteswap8(uint64_t)


cdef float _byteswap_float(float num):
    cdef uint32_t *intptr = <uint32_t *>&num
    intptr[0] = _byteswap4(intptr[0])
    return num


cdef double _byteswap_double(double num):
    cdef uint64_t *intptr = <uint64_t *>&num
    intptr[0] = _byteswap8(intptr[0])
    return num
