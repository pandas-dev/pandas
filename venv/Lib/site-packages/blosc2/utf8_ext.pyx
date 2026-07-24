#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################
# cython: boundscheck=False, wraparound=False, initializedcheck=False

"""Bulk StringDType construction and bulk UTF-8 encoding for utf8 columns.

NumPy provides no bulk constructor for a ``StringDType`` array from an
offsets+bytes buffer pair -- the only Python-level ways in are per-element
assignment or conversion from another array.  This kernel uses NumPy's C
API for ``StringDType`` (``NpyString_pack``) to fill every element of a
preallocated ``StringDType`` array directly from the raw offsets/bytes
representation, in a single C loop with no per-row Python object churn.

A second kernel (``encode_utf8_span``) goes the other direction: it turns a
list of ``str`` into the concatenated UTF-8 bytes plus per-row lengths
needed to write a utf8 column, using each string's cached UTF-8
representation instead of allocating one ``bytes`` object per row.
"""

import numpy as np
cimport numpy as cnp

from cpython.unicode cimport PyUnicode_AsUTF8AndSize
from libc.stdint cimport int64_t, uint8_t
from libc.stdlib cimport free, malloc
from libc.string cimport memcpy

cnp.import_array()


# Declared here instead of relying on `cimport numpy`: the NpyString C API
# has been part of the NumPy headers since 2.0, but its Cython declarations
# only appear in the numpy/__init__.pxd of newer NumPy versions, and some
# build environments (e.g. the Pyodide cross-build) pin an older one.  The
# functions resolve through the API table populated by cnp.import_array().
cdef extern from "numpy/ndarraytypes.h":
    ctypedef struct npy_string_allocator:
        pass
    ctypedef struct npy_packed_static_string:
        pass
    ctypedef struct PyArray_StringDTypeObject:
        pass

cdef extern from "numpy/arrayobject.h":
    npy_string_allocator* NpyString_acquire_allocator(
        const PyArray_StringDTypeObject* descr
    ) nogil
    void NpyString_release_allocator(npy_string_allocator* allocator) nogil
    int NpyString_pack(
        npy_string_allocator* allocator,
        npy_packed_static_string* packed_string,
        const char* buf,
        size_t size,
    ) nogil


def pack_utf8_span(cnp.ndarray rel not None, cnp.ndarray data not None, cnp.ndarray out not None):
    """Fill *out* in place with rows carved out of *data* using *rel*.

    *rel* is ``int64``, length ``len(out) + 1``, ``rel[0] == 0``: the
    relative byte offset of each row within *data*.  *data* is ``uint8``
    and holds valid UTF-8 bytes (this packs bytes, it does not validate
    the encoding).  *out* is a ``numpy.dtypes.StringDType`` array of length
    ``len(rel) - 1``, already allocated by the caller.
    """
    if rel.ndim != 1 or data.ndim != 1 or out.ndim != 1:
        raise ValueError("rel, data and out must be 1-D arrays")
    cdef Py_ssize_t n = out.shape[0]
    if rel.shape[0] != n + 1:
        raise ValueError("rel must have length len(out) + 1")
    if rel.dtype != np.dtype(np.int64):
        raise TypeError("rel must have dtype int64")
    if data.dtype != np.dtype(np.uint8):
        raise TypeError("data must have dtype uint8")
    if not (rel.flags["C_CONTIGUOUS"] and data.flags["C_CONTIGUOUS"] and out.flags["C_CONTIGUOUS"]):
        raise ValueError("rel, data and out must be C-contiguous")
    if n == 0:
        return
    # Cheap, vectorized well-formedness checks: a malformed rel (decreasing,
    # negative, or reaching past the end of data) would otherwise drive the
    # unchecked pointer arithmetic below out of bounds.
    if int(rel[0]) != 0:
        raise ValueError("rel[0] must be 0")
    if bool((np.diff(rel) < 0).any()):
        raise ValueError("rel must be non-decreasing")
    if int(rel[n]) > data.shape[0]:
        raise ValueError("rel values must not exceed len(data)")

    cdef const int64_t* rel_ptr = <const int64_t*>cnp.PyArray_DATA(rel)
    cdef const uint8_t* data_ptr = <const uint8_t*>cnp.PyArray_DATA(data)
    cdef char* out_data = <char*>cnp.PyArray_DATA(out)
    cdef cnp.npy_intp itemsize = cnp.PyArray_ITEMSIZE(out)
    cdef npy_string_allocator* allocator = NpyString_acquire_allocator(
        <const PyArray_StringDTypeObject*>cnp.PyArray_DESCR(out)
    )
    if allocator == NULL:
        raise TypeError("out must be a StringDType array")

    cdef Py_ssize_t i
    cdef int64_t start, length
    cdef int ret = 0
    try:
        for i in range(n):
            start = rel_ptr[i]
            length = rel_ptr[i + 1] - start
            ret = NpyString_pack(
                allocator,
                <npy_packed_static_string*>(out_data + i * itemsize),
                <const char*>(data_ptr + start),
                <size_t>length,
            )
            if ret == -1:
                break
    finally:
        NpyString_release_allocator(allocator)

    if ret == -1:
        raise MemoryError("Failed to pack a UTF-8 row into the StringDType array")


def encode_utf8_span(list values not None):
    """Return ``(data, lengths)`` for *values*, a list of ``str``.

    *data* is a ``uint8`` NDArray holding the concatenated UTF-8 encoding of
    every value, in order.  *lengths* is an ``int64`` NDArray of length
    ``len(values)`` giving each value's UTF-8 byte length.  Each value's
    UTF-8 bytes are taken from its cached representation
    (``PyUnicode_AsUTF8AndSize``) rather than allocating a new ``bytes``
    object per row.
    """
    cdef Py_ssize_t n = len(values)
    if n == 0:
        return np.empty(0, dtype=np.uint8), np.empty(0, dtype=np.int64)

    cdef cnp.ndarray lengths = np.empty(n, dtype=np.int64)
    cdef int64_t* lengths_ptr = <int64_t*>cnp.PyArray_DATA(lengths)
    cdef const char** ptrs = <const char**>malloc(<size_t>n * sizeof(char*))
    if ptrs == NULL:
        raise MemoryError("Failed to allocate temporary pointer buffer")

    cdef Py_ssize_t i
    cdef Py_ssize_t size
    cdef int64_t total = 0
    cdef int64_t offset
    cdef cnp.ndarray data
    cdef uint8_t* data_ptr
    try:
        for i in range(n):
            ptrs[i] = PyUnicode_AsUTF8AndSize(values[i], &size)
            lengths_ptr[i] = <int64_t>size
            total += size

        data = np.empty(total if total > 0 else 1, dtype=np.uint8)
        data_ptr = <uint8_t*>cnp.PyArray_DATA(data)
        offset = 0
        for i in range(n):
            size = <Py_ssize_t>lengths_ptr[i]
            if size:
                memcpy(data_ptr + offset, ptrs[i], <size_t>size)
            offset += size
    finally:
        free(ptrs)

    return data, lengths
