# -*- coding: utf-8 -*-

import cython
from cython import Py_ssize_t

from cpython cimport PyBytes_GET_SIZE, PyUnicode_GET_SIZE

try:
    from cpython cimport PyString_GET_SIZE
except ImportError:
    from cpython cimport PyUnicode_GET_SIZE as PyString_GET_SIZE

import numpy as np
from numpy cimport ndarray, uint8_t


ctypedef fused pandas_string:
    str
    unicode
    bytes


@cython.boundscheck(False)
@cython.wraparound(False)
def write_csv_rows(list data, ndarray data_index,
                   int nlevels, ndarray cols, object writer):
    """
    Write the given data to the writer object, pre-allocating where possible
    for performance improvements.

    Parameters
    ----------
    data : list
    data_index : ndarray
    nlevels : int
    cols : ndarray
    writer : object
    """
    cdef:
        int N, j, i, ncols
        list rows
        object val

    # In crude testing, N>100 yields little marginal improvement
    N = 100

    # pre-allocate rows
    ncols = len(cols)
    rows = [[None] * (nlevels + ncols) for x in range(N)]

    j = -1
    if nlevels == 1:
        for j in range(len(data_index)):
            row = rows[j % N]
            row[0] = data_index[j]
            for i in range(ncols):
                row[1 + i] = data[i][j]

            if j >= N - 1 and j % N == N - 1:
                writer.writerows(rows)
    elif nlevels > 1:
        for j in range(len(data_index)):
            row = rows[j % N]
            row[:nlevels] = list(data_index[j])
            for i in range(ncols):
                row[nlevels + i] = data[i][j]

            if j >= N - 1 and j % N == N - 1:
                writer.writerows(rows)
    else:
        for j in range(len(data_index)):
            row = rows[j % N]
            for i in range(ncols):
                row[i] = data[i][j]

            if j >= N - 1 and j % N == N - 1:
                writer.writerows(rows)

    if j >= 0 and (j < N - 1 or (j % N) != N - 1):
        writer.writerows(rows[:((j + 1) % N)])


@cython.boundscheck(False)
@cython.wraparound(False)
def convert_json_to_lines(object arr):
    """
    replace comma separated json with line feeds, paying special attention
    to quotes & brackets
    """
    cdef:
        Py_ssize_t i = 0, num_open_brackets_seen = 0, length
        bint in_quotes = 0, is_escaping = 0
        ndarray[uint8_t] narr
        unsigned char v, comma, left_bracket, right_brack, newline

    newline = ord('\n')
    comma = ord(',')
    left_bracket = ord('{')
    right_bracket = ord('}')
    quote = ord('"')
    backslash = ord('\\')

    narr = np.frombuffer(arr.encode('utf-8'), dtype='u1').copy()
    length = narr.shape[0]
    for i in range(length):
        v = narr[i]
        if v == quote and i > 0 and not is_escaping:
            in_quotes = ~in_quotes
        if v == backslash or is_escaping:
            is_escaping = ~is_escaping
        if v == comma:  # commas that should be \n
            if num_open_brackets_seen == 0 and not in_quotes:
                narr[i] = newline
        elif v == left_bracket:
            if not in_quotes:
                num_open_brackets_seen += 1
        elif v == right_bracket:
            if not in_quotes:
                num_open_brackets_seen -= 1

    return narr.tostring().decode('utf-8')


# stata, pytables
@cython.boundscheck(False)
@cython.wraparound(False)
def max_len_string_array(pandas_string[:] arr) -> Py_ssize_t:
    """ return the maximum size of elements in a 1-dim string array """
    cdef:
        Py_ssize_t i, m = 0, l = 0, length = arr.shape[0]
        pandas_string v

    for i in range(length):
        v = arr[i]
        if isinstance(v, str):
            l = PyString_GET_SIZE(v)
        elif isinstance(v, bytes):
            l = PyBytes_GET_SIZE(v)
        elif isinstance(v, unicode):
            l = PyUnicode_GET_SIZE(v)

        if l > m:
            m = l

    return m


# ------------------------------------------------------------------
# PyTables Helpers


@cython.boundscheck(False)
@cython.wraparound(False)
def string_array_replace_from_nan_rep(
        ndarray[object, ndim=1] arr, object nan_rep,
        object replace=None):
    """
    Replace the values in the array with 'replacement' if
    they are 'nan_rep'. Return the same array.
    """
    cdef:
        int length = arr.shape[0], i = 0

    if replace is None:
        replace = np.nan

    for i in range(length):
        if arr[i] == nan_rep:
            arr[i] = replace

    return arr
