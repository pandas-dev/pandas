# -*- coding: utf-8 -*-

cimport cython
from cython cimport Py_ssize_t

from cpython cimport (PyString_Check, PyBytes_Check, PyUnicode_Check,
                      PyBytes_GET_SIZE, PyUnicode_GET_SIZE)

try:
    from cpython cimport PyString_GET_SIZE
except ImportError:
    from cpython cimport PyUnicode_GET_SIZE as PyString_GET_SIZE

import numpy as np
cimport numpy as cnp
from numpy cimport ndarray, uint8_t
cnp.import_array()

cimport util


ctypedef fused pandas_string:
    str
    unicode
    bytes


def sanitize_objects(ndarray[object] values, set na_values,
                     convert_empty=True):
    cdef:
        Py_ssize_t i, n
        object val, onan
        Py_ssize_t na_count = 0
        dict memo = {}

    n = len(values)
    onan = np.nan

    for i from 0 <= i < n:
        val = values[i]
        if (convert_empty and val == '') or (val in na_values):
            values[i] = onan
            na_count += 1
        elif val in memo:
            values[i] = memo[val]
        else:
            memo[val] = val

    return na_count


@cython.boundscheck(False)
@cython.wraparound(False)
def write_csv_rows(list data, ndarray data_index,
                   int nlevels, ndarray cols, object writer):

    cdef int N, j, i, ncols
    cdef list rows
    cdef object val

    # In crude testing, N>100 yields little marginal improvement
    N=100

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
cpdef Py_ssize_t max_len_string_array(pandas_string[:] arr):
    """ return the maximum size of elements in a 1-dim string array """
    cdef:
        Py_ssize_t i, m = 0, l = 0, length = arr.shape[0]
        pandas_string v

    for i in range(length):
        v = arr[i]
        if PyString_Check(v):
            l = PyString_GET_SIZE(v)
        elif PyBytes_Check(v):
            l = PyBytes_GET_SIZE(v)
        elif PyUnicode_Check(v):
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

    cdef int length = arr.shape[0], i = 0
    if replace is None:
        replace = np.nan

    for i from 0 <= i < length:
        if arr[i] == nan_rep:
            arr[i] = replace

    return arr


def convert_timestamps(ndarray values):
    cdef:
        object val, f, result
        dict cache = {}
        Py_ssize_t i, n = len(values)
        ndarray[object] out

    # for HDFStore, a bit temporary but...

    from datetime import datetime
    f = datetime.fromtimestamp

    out = np.empty(n, dtype='O')

    for i in range(n):
        val = util.get_value_1d(values, i)
        if val in cache:
            out[i] = cache[val]
        else:
            cache[val] = out[i] = f(val)

    return out


@cython.wraparound(False)
@cython.boundscheck(False)
def fast_unique(ndarray[object] values):
    cdef:
        Py_ssize_t i, n = len(values)
        list uniques = []
        dict table = {}
        object val, stub = 0

    for i from 0 <= i < n:
        val = values[i]
        if val not in table:
            table[val] = stub
            uniques.append(val)
    try:
        uniques.sort()
    except Exception:
        pass

    return uniques
