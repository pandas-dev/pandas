# -*- coding: utf-8 -*-

cimport cython
from cython cimport Py_ssize_t

import numpy as np
cimport numpy as cnp
from numpy cimport (ndarray,
                    int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t,
                    uint32_t, uint64_t, float32_t, float64_t)
cnp.import_array()


cdef double NaN = <double> np.NaN
cdef double nan = NaN

from pandas._libs.algos import groupsort_indexer, ensure_platform_int
from pandas.core.algorithms import take_nd

include "join_func_helper.pxi"


def inner_join(ndarray[int64_t] left, ndarray[int64_t] right,
               Py_ssize_t max_groups):
    cdef:
        Py_ssize_t i, j, k, count = 0
        ndarray[int64_t] left_count, right_count, left_sorter, right_sorter
        ndarray[int64_t] left_indexer, right_indexer
        int64_t lc, rc

    # NA group in location 0

    left_sorter, left_count = groupsort_indexer(left, max_groups)
    right_sorter, right_count = groupsort_indexer(right, max_groups)

    # First pass, determine size of result set, do not use the NA group
    for i in range(1, max_groups + 1):
        lc = left_count[i]
        rc = right_count[i]

        if rc > 0 and lc > 0:
            count += lc * rc

    # group 0 is the NA group
    cdef:
        Py_ssize_t loc, left_pos = 0, right_pos = 0, position = 0
        Py_ssize_t offset

    # exclude the NA group
    left_pos = left_count[0]
    right_pos = right_count[0]

    left_indexer = np.empty(count, dtype=np.int64)
    right_indexer = np.empty(count, dtype=np.int64)

    for i in range(1, max_groups + 1):
        lc = left_count[i]
        rc = right_count[i]

        if rc > 0 and lc > 0:
            for j in range(lc):
                offset = position + j * rc
                for k in range(rc):
                    left_indexer[offset + k] = left_pos + j
                    right_indexer[offset + k] = right_pos + k
            position += lc * rc
        left_pos += lc
        right_pos += rc

    return (_get_result_indexer(left_sorter, left_indexer),
            _get_result_indexer(right_sorter, right_indexer))


def left_outer_join(ndarray[int64_t] left, ndarray[int64_t] right,
                    Py_ssize_t max_groups, sort=True):
    cdef:
        Py_ssize_t i, j, k, count = 0
        ndarray[int64_t] left_count, right_count
        ndarray left_sorter, right_sorter, rev
        ndarray[int64_t] left_indexer, right_indexer
        int64_t lc, rc

    # NA group in location 0

    left_sorter, left_count = groupsort_indexer(left, max_groups)
    right_sorter, right_count = groupsort_indexer(right, max_groups)

    # First pass, determine size of result set, do not use the NA group
    for i in range(1, max_groups + 1):
        if right_count[i] > 0:
            count += left_count[i] * right_count[i]
        else:
            count += left_count[i]

    # group 0 is the NA group
    cdef:
        Py_ssize_t loc, left_pos = 0, right_pos = 0, position = 0
        Py_ssize_t offset

    # exclude the NA group
    left_pos = left_count[0]
    right_pos = right_count[0]

    left_indexer = np.empty(count, dtype=np.int64)
    right_indexer = np.empty(count, dtype=np.int64)

    for i in range(1, max_groups + 1):
        lc = left_count[i]
        rc = right_count[i]

        if rc == 0:
            for j in range(lc):
                left_indexer[position + j] = left_pos + j
                right_indexer[position + j] = -1
            position += lc
        else:
            for j in range(lc):
                offset = position + j * rc
                for k in range(rc):
                    left_indexer[offset + k] = left_pos + j
                    right_indexer[offset + k] = right_pos + k
            position += lc * rc
        left_pos += lc
        right_pos += rc

    left_indexer = _get_result_indexer(left_sorter, left_indexer)
    right_indexer = _get_result_indexer(right_sorter, right_indexer)

    if not sort:  # if not asked to sort, revert to original order
        if len(left) == len(left_indexer):
            # no multiple matches for any row on the left
            # this is a short-cut to avoid groupsort_indexer
            # otherwise, the `else` path also works in this case
            left_sorter = ensure_platform_int(left_sorter)

            rev = np.empty(len(left), dtype=np.intp)
            rev.put(left_sorter, np.arange(len(left)))
        else:
            rev, _ = groupsort_indexer(left_indexer, len(left))

        rev = ensure_platform_int(rev)
        right_indexer = right_indexer.take(rev)
        left_indexer = left_indexer.take(rev)

    return left_indexer, right_indexer


def full_outer_join(ndarray[int64_t] left, ndarray[int64_t] right,
                    Py_ssize_t max_groups):
    cdef:
        Py_ssize_t i, j, k, count = 0
        ndarray[int64_t] left_count, right_count, left_sorter, right_sorter
        ndarray[int64_t] left_indexer, right_indexer
        int64_t lc, rc

    # NA group in location 0

    left_sorter, left_count = groupsort_indexer(left, max_groups)
    right_sorter, right_count = groupsort_indexer(right, max_groups)

    # First pass, determine size of result set, do not use the NA group
    for i in range(1, max_groups + 1):
        lc = left_count[i]
        rc = right_count[i]

        if rc > 0 and lc > 0:
            count += lc * rc
        else:
            count += lc + rc

    # group 0 is the NA group
    cdef:
        int64_t left_pos = 0, right_pos = 0
        Py_ssize_t offset, position = 0

    # exclude the NA group
    left_pos = left_count[0]
    right_pos = right_count[0]

    left_indexer = np.empty(count, dtype=np.int64)
    right_indexer = np.empty(count, dtype=np.int64)

    for i in range(1, max_groups + 1):
        lc = left_count[i]
        rc = right_count[i]

        if rc == 0:
            for j in range(lc):
                left_indexer[position + j] = left_pos + j
                right_indexer[position + j] = -1
            position += lc
        elif lc == 0:
            for j in range(rc):
                left_indexer[position + j] = -1
                right_indexer[position + j] = right_pos + j
            position += rc
        else:
            for j in range(lc):
                offset = position + j * rc
                for k in range(rc):
                    left_indexer[offset + k] = left_pos + j
                    right_indexer[offset + k] = right_pos + k
            position += lc * rc
        left_pos += lc
        right_pos += rc

    return (_get_result_indexer(left_sorter, left_indexer),
            _get_result_indexer(right_sorter, right_indexer))


def _get_result_indexer(sorter, indexer):
    if len(sorter) > 0:
        res = take_nd(sorter, indexer, fill_value=-1)
    else:
        # length-0 case
        res = np.empty(len(indexer), dtype=np.int64)
        res.fill(-1)

    return res


def ffill_indexer(ndarray[int64_t] indexer):
    cdef:
        Py_ssize_t i, n = len(indexer)
        ndarray[int64_t] result
        int64_t val, last_obs

    result = np.empty(n, dtype=np.int64)
    last_obs = -1

    for i in range(n):
        val = indexer[i]
        if val == -1:
            result[i] = last_obs
        else:
            result[i] = val
            last_obs = val

    return result


# ----------------------------------------------------------------------
# left_join_indexer, inner_join_indexer, outer_join_indexer
# ----------------------------------------------------------------------

ctypedef fused join_t:
    float64_t
    float32_t
    object
    int32_t
    int64_t
    uint64_t


# Joins on ordered, unique indices

# right might contain non-unique values

@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_indexer_unique(ndarray[join_t] left, ndarray[join_t] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int64_t] indexer
        join_t lval, rval

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


left_join_indexer_unique_float64 = left_join_indexer_unique["float64_t"]
left_join_indexer_unique_float32 = left_join_indexer_unique["float32_t"]
left_join_indexer_unique_object = left_join_indexer_unique["object"]
left_join_indexer_unique_int32 = left_join_indexer_unique["int32_t"]
left_join_indexer_unique_int64 = left_join_indexer_unique["int64_t"]
left_join_indexer_unique_uint64 = left_join_indexer_unique["uint64_t"]


# @cython.wraparound(False)
# @cython.boundscheck(False)
def left_join_indexer(ndarray[join_t] left, ndarray[join_t] right):
    """
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    """
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        join_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[join_t] result

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
    result = np.empty(count, dtype=left.dtype)

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


left_join_indexer_float64 = left_join_indexer["float64_t"]
left_join_indexer_float32 = left_join_indexer["float32_t"]
left_join_indexer_object = left_join_indexer["object"]
left_join_indexer_int32 = left_join_indexer["int32_t"]
left_join_indexer_int64 = left_join_indexer["int64_t"]
left_join_indexer_uint64 = left_join_indexer["uint64_t"]


@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer(ndarray[join_t] left, ndarray[join_t] right):
    """
    Two-pass algorithm for monotonic indexes. Handles many-to-one merges
    """
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        join_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[join_t] result

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
    result = np.empty(count, dtype=left.dtype)

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


inner_join_indexer_float64 = inner_join_indexer["float64_t"]
inner_join_indexer_float32 = inner_join_indexer["float32_t"]
inner_join_indexer_object = inner_join_indexer["object"]
inner_join_indexer_int32 = inner_join_indexer["int32_t"]
inner_join_indexer_int64 = inner_join_indexer["int64_t"]
inner_join_indexer_uint64 = inner_join_indexer["uint64_t"]


@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer(ndarray[join_t] left, ndarray[join_t] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        join_t lval, rval
        ndarray[int64_t] lindexer, rindexer
        ndarray[join_t] result

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
    result = np.empty(count, dtype=left.dtype)

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
        for i in range(nleft):
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


outer_join_indexer_float64 = outer_join_indexer["float64_t"]
outer_join_indexer_float32 = outer_join_indexer["float32_t"]
outer_join_indexer_object = outer_join_indexer["object"]
outer_join_indexer_int32 = outer_join_indexer["int32_t"]
outer_join_indexer_int64 = outer_join_indexer["int64_t"]
outer_join_indexer_uint64 = outer_join_indexer["uint64_t"]
