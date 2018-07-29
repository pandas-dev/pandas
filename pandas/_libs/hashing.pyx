# cython: profile=False
# Translated from the reference implementation
# at https://github.com/veorq/SipHash

import cython
cimport numpy as cnp
import numpy as np
from numpy cimport ndarray, uint8_t, uint32_t, uint64_t

from util cimport _checknull
from cpython cimport (PyBytes_Check,
                      PyUnicode_Check)
from libc.stdlib cimport malloc, free

DEF cROUNDS = 2
DEF dROUNDS = 4


@cython.boundscheck(False)
def hash_object_array(ndarray[object] arr, object key, object encoding='utf8'):
    """
    Parameters
    ----------
    arr : 1-d object ndarray of objects
    key : hash key, must be 16 byte len encoded
    encoding : encoding for key & arr, default to 'utf8'

    Returns
    -------
    1-d uint64 ndarray of hashes

    Notes
    -----
    allowed values must be strings, or nulls
    mixed array types will raise TypeError

    """
    cdef:
        Py_ssize_t i, l, n
        ndarray[uint64_t] result
        bytes data, k
        uint8_t *kb
        uint64_t *lens
        char **vecs
        char *cdata
        object val

    k = <bytes>key.encode(encoding)
    kb = <uint8_t *>k
    if len(k) != 16:
        raise ValueError(
            'key should be a 16-byte string encoded, got {!r} (len {})'.format(
                k, len(k)))

    n = len(arr)

    # create an array of bytes
    vecs = <char **> malloc(n * sizeof(char *))
    lens = <uint64_t*> malloc(n * sizeof(uint64_t))

    cdef list datas = []
    for i in range(n):
        val = arr[i]
        if PyBytes_Check(val):
            data = <bytes>val
        elif PyUnicode_Check(val):
            data = <bytes>val.encode(encoding)
        elif _checknull(val):
            # null, stringify and encode
            data = <bytes>str(val).encode(encoding)

        else:
            raise TypeError("{} of type {} is not a valid type for hashing, "
                            "must be string or null".format(val, type(val)))

        l = len(data)
        lens[i] = l
        cdata = data

        # keep the references alive thru the end of the
        # function
        datas.append(data)
        vecs[i] = cdata

    result = np.empty(n, dtype=np.uint64)
    with nogil:
        for i in range(n):
            result[i] = low_level_siphash(<uint8_t *>vecs[i], lens[i], kb)

    free(vecs)
    free(lens)
    return result


cdef inline uint64_t _rotl(uint64_t x, uint64_t b) nogil:
    return (x << b) | (x >> (64 - b))


cdef inline void u32to8_le(uint8_t* p, uint32_t v) nogil:
    p[0] = <uint8_t>(v)
    p[1] = <uint8_t>(v >> 8)
    p[2] = <uint8_t>(v >> 16)
    p[3] = <uint8_t>(v >> 24)


cdef inline uint64_t u8to64_le(uint8_t* p) nogil:
    return (<uint64_t>p[0] |
            <uint64_t>p[1] << 8 |
            <uint64_t>p[2] << 16 |
            <uint64_t>p[3] << 24 |
            <uint64_t>p[4] << 32 |
            <uint64_t>p[5] << 40 |
            <uint64_t>p[6] << 48 |
            <uint64_t>p[7] << 56)


cdef inline void _sipround(uint64_t* v0, uint64_t* v1,
                           uint64_t* v2, uint64_t* v3) nogil:
    v0[0] += v1[0]
    v1[0] = _rotl(v1[0], 13)
    v1[0] ^= v0[0]
    v0[0] = _rotl(v0[0], 32)
    v2[0] += v3[0]
    v3[0] = _rotl(v3[0], 16)
    v3[0] ^= v2[0]
    v0[0] += v3[0]
    v3[0] = _rotl(v3[0], 21)
    v3[0] ^= v0[0]
    v2[0] += v1[0]
    v1[0] = _rotl(v1[0], 17)
    v1[0] ^= v2[0]
    v2[0] = _rotl(v2[0], 32)


cpdef uint64_t siphash(bytes data, bytes key) except? 0:
    if len(key) != 16:
        raise ValueError(
            'key should be a 16-byte bytestring, got {!r} (len {})'.format(
                key, len(key)))
    return low_level_siphash(data, len(data), key)


@cython.cdivision(True)
cdef uint64_t low_level_siphash(uint8_t* data, size_t datalen,
                                uint8_t* key) nogil:
    cdef uint64_t v0 = 0x736f6d6570736575ULL
    cdef uint64_t v1 = 0x646f72616e646f6dULL
    cdef uint64_t v2 = 0x6c7967656e657261ULL
    cdef uint64_t v3 = 0x7465646279746573ULL
    cdef uint64_t b
    cdef uint64_t k0 = u8to64_le(key)
    cdef uint64_t k1 = u8to64_le(key + 8)
    cdef uint64_t m
    cdef int i
    cdef uint8_t* end = data + datalen - (datalen % sizeof(uint64_t))
    cdef int left = datalen & 7
    cdef int left_byte

    b = (<uint64_t>datalen) << 56
    v3 ^= k1
    v2 ^= k0
    v1 ^= k1
    v0 ^= k0

    while (data != end):
        m = u8to64_le(data)
        v3 ^= m
        for i in range(cROUNDS):
            _sipround(&v0, &v1, &v2, &v3)
        v0 ^= m

        data += sizeof(uint64_t)

    for i in range(left-1, -1, -1):
        b |= (<uint64_t>data[i]) << (i * 8)

    v3 ^= b

    for i in range(cROUNDS):
        _sipround(&v0, &v1, &v2, &v3)

    v0 ^= b
    v2 ^= 0xff

    for i in range(dROUNDS):
        _sipround(&v0, &v1, &v2, &v3)

    b = v0 ^ v1 ^ v2 ^ v3

    return b
