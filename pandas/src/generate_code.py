import os
from cStringIO import StringIO

header = """
cimport numpy as np
cimport cython

from libc.string cimport memmove

from numpy cimport *

from cpython cimport (PyDict_New, PyDict_GetItem, PyDict_SetItem,
                      PyDict_Contains, PyDict_Keys,
                      Py_INCREF, PyTuple_SET_ITEM,
                      PyTuple_SetItem,
                      PyTuple_New)
from cpython cimport PyFloat_Check
cimport cpython

import numpy as np
isnan = np.isnan

from datetime import datetime as pydatetime

# this is our datetime.pxd
from datetime cimport *

from khash cimport *

ctypedef unsigned char UChar

cimport util
from util cimport is_array, _checknull, _checknan

# import datetime C API
PyDateTime_IMPORT

# initialize numpy
import_array()
import_ufunc()

cdef int PLATFORM_INT = (<ndarray> np.arange(0, dtype=np.int_)).descr.type_num

cpdef ensure_platform_int(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == PLATFORM_INT:
            return arr
        else:
            return arr.astype(np.int_)
    else:
        return np.array(arr, dtype=np.int_)

"""


take_1d_template = """@cython.wraparound(False)
def take_1d_%(name)s_%(dest)s(ndarray[%(c_type_in)s] values,
                              ndarray[int64_t] indexer,
                              ndarray[%(c_type_out)s] out,
                              fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        %(c_type_out)s fv

    n = len(indexer)

    fv = fill_value
    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            out[i] = fv
        else:
            out[i] = %(preval)svalues[idx]%(postval)s

"""

take_2d_axis0_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_%(name)s_%(dest)s(ndarray[%(c_type_in)s, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[%(c_type_out)s, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        %(c_type_out)s fv

    n = len(indexer)
    k = values.shape[1]

    fv = fill_value

    IF %(can_copy)s:
        cdef:
            %(c_type_out)s *v, *o

        #GH3130
        if (values.strides[1] == out.strides[1] and
            values.strides[1] == sizeof(%(c_type_out)s) and
            sizeof(%(c_type_out)s) * n >= 256):

            for i from 0 <= i < n:
                idx = indexer[i]
                if idx == -1:
                    for j from 0 <= j < k:
                        out[i, j] = fv
                else:
                    v = &values[idx, 0]
                    o = &out[i, 0]
                    memmove(o, v, <size_t>(sizeof(%(c_type_out)s) * k))
            return

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                out[i, j] = %(preval)svalues[idx, j]%(postval)s

"""

take_2d_axis1_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_%(name)s_%(dest)s(ndarray[%(c_type_in)s, ndim=2] values,
                                    ndarray[int64_t] indexer,
                                    ndarray[%(c_type_out)s, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        %(c_type_out)s fv

    n = len(values)
    k = len(indexer)

    if n == 0 or k == 0:
        return

    fv = fill_value

    IF %(can_copy)s:
        cdef:
            %(c_type_out)s *v, *o

        #GH3130
        if (values.strides[0] == out.strides[0] and
            values.strides[0] == sizeof(%(c_type_out)s) and
            sizeof(%(c_type_out)s) * n >= 256):

            for j from 0 <= j < k:
                idx = indexer[j]
                if idx == -1:
                    for i from 0 <= i < n:
                        out[i, j] = fv
                else:
                    v = &values[0, idx]
                    o = &out[0, j]
                    memmove(o, v, <size_t>(sizeof(%(c_type_out)s) * n))
            return

    for j from 0 <= j < k:
        idx = indexer[j]
        if idx == -1:
            for i from 0 <= i < n:
                out[i, j] = fv
        else:
            for i from 0 <= i < n:
                out[i, j] = %(preval)svalues[i, idx]%(postval)s

"""

take_2d_multi_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_multi_%(name)s_%(dest)s(ndarray[%(c_type_in)s, ndim=2] values,
                                    indexer,
                                    ndarray[%(c_type_out)s, ndim=2] out,
                                    fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t] idx0 = indexer[0]
        ndarray[int64_t] idx1 = indexer[1]
        %(c_type_out)s fv

    n = len(idx0)
    k = len(idx1)

    fv = fill_value
    for i from 0 <= i < n:
        idx = idx0[i]
        if idx == -1:
            for j from 0 <= j < k:
                out[i, j] = fv
        else:
            for j from 0 <= j < k:
                if idx1[j] == -1:
                    out[i, j] = fv
                else:
                    out[i, j] = %(preval)svalues[idx, idx1[j]]%(postval)s

"""



'''
Backfilling logic for generating fill vector

Diagram of what's going on

Old      New    Fill vector    Mask
         .        0               1
         .        0               1
         .        0               1
A        A        0               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
B        B        1               1
         .        2               1
         .        2               1
         .        2               1
C        C        2               1
         .                        0
         .                        0
D
'''

backfill_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_%(name)s(ndarray[%(c_type)s] old, ndarray[%(c_type)s] new,
                      limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef %(c_type)s cur, prev
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[0] > old[nleft - 1]:
        return indexer

    i = nleft - 1
    j = nright - 1

    cur = old[nleft - 1]

    while j >= 0 and new[j] > cur:
        j -= 1

    while True:
        if j < 0:
            break

        if i == 0:
            while j >= 0:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] < cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j -= 1
            break

        prev = old[i - 1]

        while j >= 0 and prev < new[j] <= cur:
            if new[j] == cur:
                indexer[j] = i
            elif new[j] < cur and fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j -= 1

        fill_count = 0
        i -= 1
        cur = prev

    return indexer

"""


pad_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def pad_%(name)s(ndarray[%(c_type)s] old, ndarray[%(c_type)s] new,
                   limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int64_t, ndim=1] indexer
    cdef %(c_type)s cur, next
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int64)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[nright - 1] < old[0]:
        return indexer

    i = j = 0

    cur = old[0]

    while j <= nright - 1 and new[j] < cur:
        j += 1

    while True:
        if j == nright:
            break

        if i == nleft - 1:
            while j < nright:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] > cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j += 1
            break

        next = old[i + 1]

        while j < nright and cur <= new[j] < next:
            if new[j] == cur:
                indexer[j] = i
            elif fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j += 1

        fill_count = 0
        i += 1
        cur = next

    return indexer

"""

pad_1d_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def pad_inplace_%(name)s(ndarray[%(c_type)s] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef %(c_type)s val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[0]
    for i in range(N):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]

"""

pad_2d_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_%(name)s(ndarray[%(c_type)s, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef %(c_type)s val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, 0]
        for i in range(N):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]
"""

backfill_2d_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_2d_inplace_%(name)s(ndarray[%(c_type)s, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef %(c_type)s val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1 , -1):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]
"""

backfill_1d_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_inplace_%(name)s(ndarray[%(c_type)s] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef %(c_type)s val
    cdef int lim, fill_count = 0

    N = len(values)

    # GH 2778
    if N == 0:
        return

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1 , -1):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]
"""


diff_2d_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def diff_2d_%(name)s(ndarray[%(c_type)s, ndim=2] arr,
                     ndarray[%(dest_type2)s, ndim=2] out,
                    Py_ssize_t periods, int axis):
    cdef:
        Py_ssize_t i, j, sx, sy

    sx, sy = (<object> arr).shape
    if arr.flags.f_contiguous:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for j in range(sy):
                for i in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for j in range(start, stop):
                for i in range(sx):
                    out[i, j] = arr[i, j] - arr[i, j - periods]
    else:
        if axis == 0:
            if periods >= 0:
                start, stop = periods, sx
            else:
                start, stop = 0, sx + periods
            for i in range(start, stop):
                for j in range(sy):
                    out[i, j] = arr[i, j] - arr[i - periods, j]
        else:
            if periods >= 0:
                start, stop = periods, sy
            else:
                start, stop = 0, sy + periods
            for i in range(sx):
                for j in range(start, stop):
                    out[i, j] = arr[i, j] - arr[i, j - periods]
"""

is_monotonic_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_%(name)s(ndarray[%(c_type)s] arr):
    '''
    Returns
    -------
    is_monotonic, is_unique
    '''
    cdef:
        Py_ssize_t i, n
        %(c_type)s prev, cur
        bint is_unique = 1

    n = len(arr)

    if n < 2:
        return True, True

    prev = arr[0]
    for i in range(1, n):
        cur = arr[i]
        if cur < prev:
            return False, None
        elif cur == prev:
            is_unique = 0
        prev = cur
    return True, is_unique
"""

map_indices_template = """@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_%(name)s(ndarray[%(c_type)s] index):
    '''
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    '''
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result

"""

groupby_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_%(name)s(ndarray[%(c_type)s] index, ndarray labels):
    cdef dict result = {}
    cdef Py_ssize_t i, length
    cdef list members
    cdef object idx, key

    length = len(index)

    if not length == len(labels):
       raise AssertionError("len(index) != len(labels)")

    for i in range(length):
        key = util.get_value_1d(labels, i)

        if _checknull(key):
            continue

        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result

"""

group_last_template = """@cython.wraparound(False)
@cython.wraparound(False)
def group_last_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[%(c_type)s, ndim=2] values,
               ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    for i in range(N):
        lab = labels[i]
        if lab < 0:
            continue

        counts[lab] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[lab, j] += 1
                resx[lab, j] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]
"""

group_last_bin_template = """@cython.wraparound(False)
@cython.wraparound(False)
def group_last_bin_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[%(c_type)s, ndim=2] values,
                   ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] resx, nobs

    nobs = np.zeros_like(out)
    resx = np.empty_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    for i in range(N):
        while b < ngroups - 1 and i >= bins[b]:
            b += 1

        counts[b] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[b, j] += 1
                resx[b, j] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]
"""

group_nth_bin_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth_bin_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[%(c_type)s, ndim=2] values,
                  ndarray[int64_t] bins, int64_t rank):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] resx, nobs

    nobs = np.zeros_like(out)
    resx = np.empty_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    for i in range(N):
        while b < ngroups - 1 and i >= bins[b]:
            b += 1

        counts[b] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[b, j] += 1
                if nobs[b, j] == rank:
                    resx[b, j] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]
"""

group_nth_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[%(c_type)s, ndim=2] values,
              ndarray[int64_t] labels, int64_t rank):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    for i in range(N):
        lab = labels[i]
        if lab < 0:
            continue

        counts[lab] += 1
        for j in range(K):
            val = values[i, j]

            # not nan
            if val == val:
                nobs[lab, j] += 1
                if nobs[lab, j] == rank:
                    resx[lab, j] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = resx[i, j]
"""

group_add_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def group_add_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[%(c_type)s, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] sumx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    sumx[lab, j] += val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                sumx[lab, 0] += val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j]
"""

group_add_bin_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def group_add_bin_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[%(dest_type2)s, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b, nbins
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] sumx, nobs

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    sumx[b, j] += val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                sumx[b, 0] += val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j]
"""

group_prod_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def group_prod_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[%(c_type)s, ndim=2] values,
               ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] prodx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    prodx = np.ones_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    prodx[lab, j] *= val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                prodx[lab, 0] *= val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = prodx[i, j]
"""

group_prod_bin_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def group_prod_bin_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[%(dest_type2)s, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] prodx, nobs

    nobs = np.zeros_like(out)
    prodx = np.ones_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    prodx[b, j] *= val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                prodx[b, 0] *= val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = prodx[i, j]
"""

group_var_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def group_var_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[%(dest_type2)s, ndim=2] values,
              ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab
        %(dest_type2)s val, ct
        ndarray[%(dest_type2)s, ndim=2] nobs, sumx, sumxx

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)
    sumxx = np.zeros_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):

            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1

            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    sumx[lab, j] += val
                    sumxx[lab, j] += val * val
    else:
        for i in range(N):

            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]
            # not nan
            if val == val:
                nobs[lab, 0] += 1
                sumx[lab, 0] += val
                sumxx[lab, 0] += val * val


    for i in range(len(counts)):
        for j in range(K):
            ct = nobs[i, j]
            if ct < 2:
                out[i, j] = nan
            else:
                out[i, j] = ((ct * sumxx[i, j] - sumx[i, j] * sumx[i, j]) /
                             (ct * ct - ct))
"""

group_var_bin_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def group_var_bin_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[%(dest_type2)s, ndim=2] values,
                  ndarray[int64_t] bins):

    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        %(dest_type2)s val, ct
        ndarray[%(dest_type2)s, ndim=2] nobs, sumx, sumxx

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)
    sumxx = np.zeros_like(out)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1

            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    sumx[b, j] += val
                    sumxx[b, j] += val * val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                sumx[b, 0] += val
                sumxx[b, 0] += val * val

    for i in range(ngroups):
        for j in range(K):
            ct = nobs[i, j]
            if ct < 2:
                out[i, j] = nan
            else:
                out[i, j] = ((ct * sumxx[i, j] - sumx[i, j] * sumx[i, j]) /
                             (ct * ct - ct))
"""

# add passing bin edges, instead of labels


#----------------------------------------------------------------------
# group_min, group_max

group_min_bin_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def group_min_bin_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[%(dest_type2)s, ndim=2] values,
                   ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] minx, nobs

    nobs = np.zeros_like(out)

    minx = np.empty_like(out)
    minx.fill(np.inf)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    if val < minx[b, j]:
                        minx[b, j] = val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                if val < minx[b, 0]:
                    minx[b, 0] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = minx[i, j]
"""

group_max_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def group_max_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[%(dest_type2)s, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] maxx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    maxx = np.empty_like(out)
    maxx.fill(-np.inf)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    if val > maxx[lab, j]:
                        maxx[lab, j] = val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                if val > maxx[lab, 0]:
                    maxx[lab, 0] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = maxx[i, j]
"""

group_max_bin_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def group_max_bin_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[%(dest_type2)s, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] maxx, nobs

    nobs = np.zeros_like(out)
    maxx = np.empty_like(out)
    maxx.fill(-np.inf)

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    if val > maxx[b, j]:
                        maxx[b, j] = val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                if val > maxx[b, 0]:
                    maxx[b, 0] = val

    for i in range(ngroups):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = maxx[i, j]
"""


group_min_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def group_min_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[%(dest_type2)s, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] minx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    minx = np.empty_like(out)
    minx.fill(np.inf)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    if val < minx[lab, j]:
                        minx[lab, j] = val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[lab, 0] += 1
                if val < minx[lab, 0]:
                    minx[lab, 0] = val

    for i in range(len(counts)):
        for j in range(K):
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = minx[i, j]
"""


group_mean_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def group_mean_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[%(dest_type2)s, ndim=2] values,
               ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] sumx, nobs

    if not len(values) == len(labels):
       raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape

    if K > 1:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]
                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    sumx[lab, j] += val
    else:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            val = values[i, 0]
            # not nan
            if val == val:
                nobs[lab, 0] += 1
                sumx[lab, 0] += val

    for i in range(len(counts)):
        for j in range(K):
            count = nobs[i, j]
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j] / count
"""

group_mean_bin_template = """
def group_mean_bin_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[%(dest_type2)s, ndim=2] values,
                   ndarray[int64_t] bins):
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        %(dest_type2)s val, count
        ndarray[%(dest_type2)s, ndim=2] sumx, nobs

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape
    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    b = 0
    if K > 1:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[b, j] += 1
                    sumx[b, j] += val
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                b += 1

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                nobs[b, 0] += 1
                sumx[b, 0] += val

    for i in range(ngroups):
        for j in range(K):
            count = nobs[i, j]
            if nobs[i, j] == 0:
                out[i, j] = nan
            else:
                out[i, j] = sumx[i, j] / count
"""

group_ohlc_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def group_ohlc_%(name)s(ndarray[%(dest_type2)s, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[%(dest_type2)s, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        %(dest_type2)s val, count
        %(dest_type2)s vopen, vhigh, vlow, vclose, NA
        bint got_first = 0

    if bins[len(bins) - 1] == len(values):
        ngroups = len(bins)
    else:
        ngroups = len(bins) + 1

    N, K = (<object> values).shape

    if out.shape[1] != 4:
        raise ValueError('Output array must have 4 columns')

    NA = np.nan

    b = 0
    if K > 1:
        raise NotImplementedError
    else:
        for i in range(N):
            while b < ngroups - 1 and i >= bins[b]:
                if not got_first:
                    out[b, 0] = NA
                    out[b, 1] = NA
                    out[b, 2] = NA
                    out[b, 3] = NA
                else:
                    out[b, 0] = vopen
                    out[b, 1] = vhigh
                    out[b, 2] = vlow
                    out[b, 3] = vclose
                b += 1
                got_first = 0

            counts[b] += 1
            val = values[i, 0]

            # not nan
            if val == val:
                if not got_first:
                    got_first = 1
                    vopen = val
                    vlow = val
                    vhigh = val
                else:
                    if val < vlow:
                        vlow = val
                    if val > vhigh:
                        vhigh = val
                vclose = val

        if not got_first:
            out[b, 0] = NA
            out[b, 1] = NA
            out[b, 2] = NA
            out[b, 3] = NA
        else:
            out[b, 0] = vopen
            out[b, 1] = vhigh
            out[b, 2] = vlow
            out[b, 3] = vclose
"""

arrmap_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_%(name)s(ndarray[%(c_type)s] index, object func):
    cdef Py_ssize_t length = index.shape[0]
    cdef Py_ssize_t i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    from pandas.lib import maybe_convert_objects

    for i in range(length):
        result[i] = func(index[i])

    return maybe_convert_objects(result)

"""

#----------------------------------------------------------------------
# Joins on ordered, unique indices

# right might contain non-unique values

left_join_unique_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_indexer_unique_%(name)s(ndarray[%(c_type)s] left,
                                      ndarray[%(c_type)s] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int64_t] indexer
        %(c_type)s lval, rval

    i = 0
    j = 0
    nleft = len(left)
    nright = len(right)

    indexer = np.empty(nleft, dtype=np.int64)
    while True:
        if i == nleft:
            break

        if j == nright:
            indexer[i] = -1
            i += 1
            continue

        rval = right[j]

        while i < nleft - 1 and left[i] == rval:
            indexer[i] = j
            i += 1

        if left[i] == right[j]:
            indexer[i] = j
            i += 1
            while i < nleft - 1 and left[i] == rval:
                indexer[i] = j
                i += 1
            j += 1
        elif left[i] > rval:
            indexer[i] = -1
            j += 1
        else:
            indexer[i] = -1
            i += 1
    return indexer

"""

# @cython.wraparound(False)
# @cython.boundscheck(False)

left_join_template = """
def left_join_indexer_%(name)s(ndarray[%(c_type)s] left,
                              ndarray[%(c_type)s] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        %(c_type)s lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[%(c_type)s] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft > 0:
        while i < nleft:
            if j == nright:
                count += nleft - i
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                count += 1
                i += 1
            else:
                j += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=%(dtype)s)

    i = 0
    j = 0
    count = 0
    if nleft > 0:
        while i < nleft:
            if j == nright:
                while i < nleft:
                    lindexer[count] = i
                    rindexer[count] = -1
                    result[count] = left[i]
                    i += 1
                    count += 1
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = left[i]
                count += 1
                i += 1
            else:
                j += 1

    return result, lindexer, rindexer

"""


inner_join_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer_%(name)s(ndarray[%(c_type)s] left,
                              ndarray[%(c_type)s] right):
    '''
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        %(c_type)s lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[%(c_type)s] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft > 0 and nright > 0:
        while True:
            if i == nleft:
                break
            if j == nright:
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                i += 1
            else:
                j += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=%(dtype)s)

    i = 0
    j = 0
    count = 0
    if nleft > 0 and nright > 0:
        while True:
            if i == nleft:
                break
            if j == nright:
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = rval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                i += 1
            else:
                j += 1

    return result, lindexer, rindexer

"""


outer_join_template2 = """@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer_%(name)s(ndarray[%(c_type)s] left,
                                ndarray[%(c_type)s] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        %(c_type)s lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[%(c_type)s] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    if nleft == 0:
        count = nright
    elif nright == 0:
        count = nleft
    else:
        while True:
            if i == nleft:
                count += nright - j
                break
            if j == nright:
                count += nleft - i
                break

            lval = left[i]
            rval = right[j]
            if lval == rval:
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                count += 1
                i += 1
            else:
                count += 1
                j += 1

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=%(dtype)s)

    # do it again, but populate the indexers / result

    i = 0
    j = 0
    count = 0
    if nleft == 0:
        for j in range(nright):
            lindexer[j] = -1
            rindexer[j] = j
            result[j] = right[j]
    elif nright == 0:
        for i in range(nright):
            lindexer[i] = i
            rindexer[i] = -1
            result[i] = left[i]
    else:
        while True:
            if i == nleft:
                while j < nright:
                    lindexer[count] = -1
                    rindexer[count] = j
                    result[count] = right[j]
                    count += 1
                    j += 1
                break
            if j == nright:
                while i < nleft:
                    lindexer[count] = i
                    rindexer[count] = -1
                    result[count] = left[i]
                    count += 1
                    i += 1
                break

            lval = left[i]
            rval = right[j]

            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                count += 1
                if i < nleft - 1:
                    if j < nright - 1 and right[j + 1] == rval:
                        j += 1
                    else:
                        i += 1
                        if left[i] != rval:
                            j += 1
                elif j < nright - 1:
                    j += 1
                    if lval != right[j]:
                        i += 1
                else:
                    # end of the road
                    break
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = lval
                count += 1
                i += 1
            else:
                lindexer[count] = -1
                rindexer[count] = j
                result[count] = rval
                count += 1
                j += 1

    return result, lindexer, rindexer

"""

outer_join_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer_%(name)s(ndarray[%(c_type)s] left,
                                ndarray[%(c_type)s] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        %(c_type)s lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[%(c_type)s] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    while True:
        if i == nleft:
            if j == nright:
                # we are done
                break
            else:
                while j < nright:
                    j += 1
                    count += 1
                break
        elif j == nright:
            while i < nleft:
                i += 1
                count += 1
            break
        else:
            if left[i] == right[j]:
                i += 1
                j += 1
            elif left[i] < right[j]:
                i += 1
            else:
                j += 1

            count += 1

    lindexer = np.empty(count, dtype=np.int64)
    rindexer = np.empty(count, dtype=np.int64)
    result = np.empty(count, dtype=%(dtype)s)

    # do it again, but populate the indexers / result

    i = 0
    j = 0
    count = 0
    while True:
        if i == nleft:
            if j == nright:
                # we are done
                break
            else:
                while j < nright:
                    lindexer[count] = -1
                    rindexer[count] = j
                    result[count] = right[j]
                    j += 1
                    count += 1
                break
        elif j == nright:
            while i < nleft:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = left[i]
                i += 1
                count += 1
            break
        else:
            lval = left[i]
            rval = right[j]
            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                i += 1
                j += 1
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = lval
                i += 1
            else:
                lindexer[count] = -1
                rindexer[count] = j
                result[count] = rval
                j += 1

            count += 1

    return result, lindexer, rindexer

"""

# ensure_dtype functions

ensure_dtype_template = """
cpdef ensure_%(name)s(object arr):
    if util.is_array(arr):
        if (<ndarray> arr).descr.type_num == NPY_%(ctype)s:
            return arr
        else:
            return arr.astype(np.%(dtype)s)
    else:
        return np.array(arr, dtype=np.%(dtype)s)

"""

ensure_functions = [
    ('float64', 'FLOAT64', 'float64'),
    ('float32', 'FLOAT32', 'float32'),
    ('int8', 'INT8', 'int8'),
    ('int16', 'INT16', 'int16'),
    ('int32', 'INT32', 'int32'),
    ('int64', 'INT64', 'int64'),
    # ('platform_int', 'INT', 'int_'),
    ('object', 'OBJECT', 'object_'),
]

def generate_ensure_dtypes():
    output = StringIO()
    for name, ctype, dtype in ensure_functions:
        filled = ensure_dtype_template % locals()
        output.write(filled)
    return output.getvalue()

#----------------------------------------------------------------------
# Fast "put" logic for speeding up interleaving logic

put2d_template = """
def put2d_%(name)s_%(dest_type)s(ndarray[%(c_type)s, ndim=2, cast=True] values,
                              ndarray[int64_t] indexer, Py_ssize_t loc,
                              ndarray[%(dest_type2)s] out):
    cdef:
        Py_ssize_t i, j, k

    k = len(values)
    for j from 0 <= j < k:
        i = indexer[j]
        out[i] = values[j, loc]
"""


#-------------------------------------------------------------------------
# Generators

def generate_put_template(template, use_ints = True, use_floats = True):
    floats_list = [
        ('float64', 'float64_t', 'float64_t', 'np.float64'),
        ('float32', 'float32_t', 'float32_t', 'np.float32'),
        ]
    ints_list = [
        ('int8',  'int8_t',  'float32_t', 'np.float32'),
        ('int16', 'int16_t', 'float32_t', 'np.float32'),
        ('int32', 'int32_t', 'float64_t', 'np.float64'),
        ('int64', 'int64_t', 'float64_t', 'np.float64'),
        ]
    function_list = []
    if use_floats:
        function_list.extend(floats_list)
    if use_ints:
        function_list.extend(ints_list)

    output = StringIO()
    for name, c_type, dest_type, dest_dtype in function_list:
        func = template % {'name' : name,
                           'c_type' : c_type,
                           'dest_type' : dest_type.replace('_t', ''),
                           'dest_type2' : dest_type,
                           'dest_dtype' : dest_dtype}
        output.write(func)
    return output.getvalue()

def generate_take_template(template, exclude=None):
    # name, dest, ctypein, ctypeout, preval, postval, cancopy
    function_list = [
        ('bool', 'bool', 'uint8_t', 'uint8_t', '', '', True),
        ('bool', 'object', 'uint8_t', 'object',
         'True if ', ' > 0 else False', False),
        ('int8', 'int8', 'int8_t', 'int8_t', '', '', True),
        ('int8', 'int32', 'int8_t', 'int32_t', '', '', False),
        ('int8', 'int64', 'int8_t', 'int64_t', '', '', False),
        ('int8', 'float64', 'int8_t', 'float64_t', '', '', False),
        ('int16', 'int16', 'int16_t', 'int16_t', '', '', True),
        ('int16', 'int32', 'int16_t', 'int32_t', '', '', False),
        ('int16', 'int64', 'int16_t', 'int64_t', '', '', False),
        ('int16', 'float64', 'int16_t', 'float64_t', '', '', False),
        ('int32', 'int32', 'int32_t', 'int32_t', '', '', True),
        ('int32', 'int64', 'int32_t', 'int64_t', '', '', False),
        ('int32', 'float64', 'int32_t', 'float64_t', '', '', False),
        ('int64', 'int64', 'int64_t', 'int64_t', '', '', True),
        ('int64', 'float64', 'int64_t', 'float64_t', '', '', False),
        ('float32', 'float32', 'float32_t', 'float32_t', '', '', True),
        ('float32', 'float64', 'float32_t', 'float64_t', '', '', False),
        ('float64', 'float64', 'float64_t', 'float64_t', '', '', True),
        ('object', 'object', 'object', 'object', '', '', False)
    ]

    output = StringIO()
    for (name, dest, c_type_in, c_type_out,
         preval, postval, can_copy) in function_list:
        if exclude is not None and name in exclude:
            continue

        func = template % {'name': name, 'dest': dest,
                           'c_type_in': c_type_in, 'c_type_out': c_type_out,
                           'preval': preval, 'postval': postval,
                           'can_copy': 'True' if can_copy else 'False'}
        output.write(func)
    return output.getvalue()

def generate_from_template(template, exclude=None):
    # name, ctype, capable of holding NA
    function_list = [
        ('float64', 'float64_t', 'np.float64', True),
        ('float32', 'float32_t', 'np.float32', True),
        ('object', 'object', 'object', True),
        ('int32', 'int32_t', 'np.int32', False),
        ('int64', 'int64_t', 'np.int64', False),
        ('bool', 'uint8_t', 'np.bool', False)
    ]

    output = StringIO()
    for name, c_type, dtype, can_hold_na in function_list:
        if exclude is not None and name in exclude:
            continue

        func = template % {'name': name, 'c_type': c_type,
                           'dtype': dtype,
                           'raise_on_na': 'False' if can_hold_na else 'True'}
        output.write(func)
    return output.getvalue()

put_2d = [diff_2d_template]
groupbys = [group_last_template,
            group_last_bin_template,
            group_nth_template,
            group_nth_bin_template,
            group_add_template,
            group_add_bin_template,
            group_prod_template,
            group_prod_bin_template,
            group_var_template,
            group_var_bin_template,
            group_mean_template,
            group_mean_bin_template,
            group_min_template,
            group_min_bin_template,
            group_max_template,
            group_max_bin_template,
            group_ohlc_template]

templates_1d = [map_indices_template,
                pad_template,
                backfill_template,
                pad_1d_template,
                backfill_1d_template,
                pad_2d_template,
                backfill_2d_template,
                is_monotonic_template,
                groupby_template,
                arrmap_template]

nobool_1d_templates = [left_join_unique_template,
                       left_join_template,
                       outer_join_template2,
                       inner_join_template]

take_templates = [take_1d_template,
                  take_2d_axis0_template,
                  take_2d_axis1_template,
                  take_2d_multi_template]

def generate_take_cython_file(path='generated.pyx'):
    with open(path, 'w') as f:
        print >> f, header

        print >> f, generate_ensure_dtypes()

        for template in templates_1d:
            print >> f, generate_from_template(template)

        for template in take_templates:
            print >> f, generate_take_template(template)

        for template in put_2d:
            print >> f, generate_put_template(template)

        for template in groupbys:
            print >> f, generate_put_template(template, use_ints = False)

        # for template in templates_1d_datetime:
        #     print >> f, generate_from_template_datetime(template)

        # for template in templates_2d_datetime:
        #     print >> f, generate_from_template_datetime(template, ndim=2)

        for template in nobool_1d_templates:
            print >> f, generate_from_template(template, exclude=['bool'])

if __name__ == '__main__':
    generate_take_cython_file()
