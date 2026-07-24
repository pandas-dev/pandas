"""
Native accelerators for Parquet encoding and decoding.
"""
# cython: profile=False
# cython: linetrace=False
# cython: binding=False
# cython: language_level=3
# cython: initializedcheck=False
# cython: boundscheck=False
# cython: wraparound=False
# cython: overflowcheck=False
# cython: initializedcheck=False
# cython: cdivision=True
# cython: always_allow_keywords=False

from libc.string cimport memcpy

from cpython cimport (PyUnicode_AsUTF8String, PyUnicode_DecodeUTF8,
                      PyBytes_CheckExact, PyBytes_FromStringAndSize,
                      PyBytes_GET_SIZE, PyBytes_AS_STRING)
from cpython.unicode cimport PyUnicode_DecodeUTF8

import numpy as np
cimport numpy as np
import cython


_obj_dtype = np.dtype('object')


def array_encode_utf8(inp):
    """
    utf-8 encode all elements of a 1d ndarray of "object" dtype.
    A new ndarray of bytes objects is returned.

    Accepts both numpy object arrays and pandas Series (including
    ArrowStringArray-backed StringDtype series in pandas 3).
    """
    # TODO: combine with pack_byte_array as is done for unpack
    cdef:
        Py_ssize_t i, n
        np.ndarray[object, ndim=1] arr
        np.ndarray[object] result

    # Use np.asarray with dtype=object to handle both plain numpy object arrays
    # (zero-copy) and pandas ArrowStringArray-backed series (one copy, unavoidable
    # since arrow uses a different memory layout incompatible with numpy views).
    arr = np.asarray(inp, dtype=object)

    n = arr.shape[0]
    # TODO: why not inplace?
    result = np.empty(n, dtype=object)
    for i in range(n):
        # Fast utf-8 encoding, avoiding method call and codec lookup indirection
        result[i] = PyUnicode_AsUTF8String(arr[i])

    return result


def pack_byte_array(list items):
    """
    Pack a variable length byte array column.
    A bytes object is returned.
    """
    cdef:
        Py_ssize_t i, n, itemlen, total_size
        unsigned char *start
        unsigned char *data
        object val, out

    # Strategy: compute the total output size and allocate it in one go.
    n = len(items)
    total_size = 0
    for i in range(n):
        val = items[i]
        if not PyBytes_CheckExact(val):
            raise TypeError("expected list of bytes")
        total_size += 4 + PyBytes_GET_SIZE(val)

    out = PyBytes_FromStringAndSize(NULL, total_size)
    start = data = <unsigned char *> PyBytes_AS_STRING(out)

    # Copy data to output.
    for i in range(n):
        val = items[i]
        # `itemlen` should be >= 0, so no signed extension issues
        itemlen = PyBytes_GET_SIZE(val)
        (<int*> data)[0] = itemlen
        data += 4
        memcpy(data, PyBytes_AS_STRING(val), itemlen)
        data += itemlen

    assert (data - start) == total_size
    return out


@cython.boundscheck(False)
def unpack_byte_array(const unsigned char[::1] raw_bytes, Py_ssize_t n, const char utf=False):
    """
    Unpack a variable length byte array column.
    An array of bytes objects is returned.
    """
    cdef:
        Py_ssize_t i = 0
        char* ptr = <char*>&raw_bytes[0]
        int itemlen, bytecount
        np.ndarray[object, ndim=1, mode="c"] out = np.empty(n, dtype="object")

    assert out is not None
    bytecount = raw_bytes.shape[0]
    while i < n and bytecount > 0:

        itemlen = (<int*> ptr)[0]
        ptr += 4
        if utf:
            out[i] = PyUnicode_DecodeUTF8(ptr, itemlen, "ignore")
        else:
            out[i] = PyBytes_FromStringAndSize(ptr, itemlen)
        ptr += itemlen
        bytecount -= 4 + itemlen
        i += 1

    return out


def unpack_byte_array_arrow(const unsigned char[::1] raw_bytes, Py_ssize_t n,
                            validity=None, Py_ssize_t n_total=-1):
    """
    Unpack a variable-length BYTE_ARRAY parquet column directly into a
    ``pd.arrays.ArrowStringArray``, without creating intermediate Python
    ``str`` or ``bytes`` objects in the inner loop.

    Parameters
    ----------
    raw_bytes : 1-D uint8 array
        Packed parquet BYTE_ARRAY data for the *valid* (non-null) values only.
        Each value is encoded as a little-endian int32 length followed by the
        raw UTF-8 bytes.
    n : int
        Number of valid (non-null) string values in ``raw_bytes``.
    validity : 1-D uint8 array or None
        Optional definition-level array of length ``n_total``.
        A value equal to the maximum (non-zero) level means the row is valid;
        any other value means null.  When provided, ``n_total`` must also be
        given.  The convention matches the ``defi`` arrays produced by
        ``read_def``/``read_data_page_v2``: the caller passes the raw defi
        array; this function treats ``defi[i] == max(defi)`` as valid.
        Pass ``None`` (the default) for columns with no nulls.
    n_total : int
        Total number of rows including nulls.  Ignored when ``validity`` is
        ``None`` (n_total defaults to n in that case).

    Returns
    -------
    pd.arrays.ArrowStringArray
    """
    import pandas as pd
    import pyarrow as pa

    cdef:
        Py_ssize_t i, j, src_pos, dst_pos, total_bytes, length
        char* src_ptr
        unsigned char* dst_ptr
        unsigned char* val_ptr
        np.ndarray[np.int64_t, ndim=1] offsets
        np.ndarray[np.uint8_t, ndim=1] data_buf
        np.ndarray[np.uint8_t, ndim=1] null_bitmap
        np.ndarray[np.uint8_t, ndim=1] val_arr
        unsigned char max_defi
        Py_ssize_t n_out
        int has_nulls

    has_nulls = validity is not None
    if has_nulls:
        n_out = n_total if n_total >= 0 else n
        # Find the max definition level — it's the value that marks "valid"
        val_arr = np.asarray(validity, dtype=np.uint8)
        max_defi = val_arr.max() if len(val_arr) > 0 else 0
    else:
        n_out = n

    # Handle the empty/all-null case: n==0 means no valid (non-null) strings.
    if n == 0:
        if has_nulls and n_out > 0:
            # All rows are null: build a validity bitmap of all zeros.
            n_bitmap_bytes = (n_out + 7) // 8
            _null_bitmap = np.zeros(n_bitmap_bytes, dtype=np.uint8)
            _offsets = np.zeros(n_out + 1, dtype=np.int64)
            null_arr = pa.Array.from_buffers(
                pa.large_utf8(), n_out,
                [pa.py_buffer(_null_bitmap.tobytes()),
                 pa.py_buffer(_offsets.tobytes()),
                 pa.py_buffer(b"")]
            )
        else:
            # Truly empty column (0 rows).
            null_arr = pa.Array.from_buffers(
                pa.large_utf8(), 0,
                [None,
                 pa.py_buffer(np.zeros(1, dtype=np.int64).tobytes()),
                 pa.py_buffer(b"")]
            )
        return pd.arrays.ArrowStringArray(pa.chunked_array([null_arr]))

    # --- First pass: compute string lengths and total byte count ---
    # We need offsets of length n_out+1, where null positions have 0-length.
    offsets = np.empty(n_out + 1, dtype=np.int64)
    offsets[0] = 0

    src_ptr = <char*>&raw_bytes[0]
    src_pos = 0
    j = 0  # index into valid strings
    total_bytes = 0

    if has_nulls:
        for i in range(n_out):
            if val_arr[i] == max_defi:
                # Valid: read the length from raw_bytes
                length = (<int*>(src_ptr + src_pos))[0]
                src_pos += 4 + length
                total_bytes += length
                offsets[i + 1] = offsets[i] + length
                j += 1
            else:
                # Null: zero-length entry in offsets
                offsets[i + 1] = offsets[i]
    else:
        for i in range(n_out):
            length = (<int*>(src_ptr + src_pos))[0]
            src_pos += 4 + length
            offsets[i + 1] = offsets[i] + length
        total_bytes = offsets[n_out]

    # --- Allocate string data buffer ---
    data_buf = np.empty(total_bytes, dtype=np.uint8)

    # --- Second pass: copy string bytes into data_buf ---
    src_ptr = <char*>&raw_bytes[0]
    src_pos = 0
    dst_ptr = <unsigned char*>data_buf.data
    dst_pos = 0

    if has_nulls:
        for i in range(n_out):
            if val_arr[i] == max_defi:
                length = (<int*>(src_ptr + src_pos))[0]
                src_pos += 4
                memcpy(dst_ptr + dst_pos, src_ptr + src_pos, length)
                src_pos += length
                dst_pos += length
    else:
        for i in range(n_out):
            length = (<int*>(src_ptr + src_pos))[0]
            src_pos += 4
            memcpy(dst_ptr + dst_pos, src_ptr + src_pos, length)
            src_pos += length
            dst_pos += length

    # --- Build null bitmap if needed ---
    # Arrow convention: bit=1 means valid, bit=0 means null (LSB first within byte)
    if has_nulls:
        n_bitmap_bytes = (n_out + 7) // 8
        null_bitmap = np.zeros(n_bitmap_bytes, dtype=np.uint8)
        for i in range(n_out):
            if val_arr[i] == max_defi:
                null_bitmap[i // 8] |= (<unsigned char>1 << (i % 8))
        validity_buf = pa.py_buffer(null_bitmap.tobytes())
    else:
        validity_buf = None

    # --- Build pyarrow LargeString array from buffers ---
    arr = pa.Array.from_buffers(
        pa.large_utf8(),
        n_out,
        [validity_buf,
         pa.py_buffer(offsets.tobytes()),
         pa.py_buffer(data_buf.tobytes())]
    )
    return pd.arrays.ArrowStringArray(pa.chunked_array([arr]))
