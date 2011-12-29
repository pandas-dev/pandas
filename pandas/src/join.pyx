
def inner_join(ndarray[int32_t] left, ndarray[int32_t] right):
    cdef:
        Py_ssize_t count

def left_outer_join(ndarray[int32_t] left, ndarray[int32_t] right,
                    Py_ssize_t max_groups):
    cdef:
        Py_ssize_t i, j, k, count = 0
        ndarray[int32_t] left_count, right_count, left_sorter, right_sorter
        ndarray[int32_t] left_indexer, right_indexer
        int32_t lc, rc

    # NA group in location 0

    left_sorter, left_count = groupsort_indexer(left, max_groups)
    right_sorter, right_count = groupsort_indexer(right, max_groups)

    # First pass, determine size of result set, do not use the NA group
    for i in range(1, max_groups + 1):
        if right_count[i] > 0:
            count += left_count[i] * right_count[i]
        else:
            count += left_count[i]

    left_indexer = np.empty(count, dtype='i4')
    right_indexer = np.empty(count, dtype='i4')

    left = left.take(left_sorter)
    right = right.take(right_sorter)

    # group 0 is the NA group
    cdef:
        Py_ssize_t loc, left_pos = 0, right_pos = 0, position = 0
        Py_ssize_t offset

    # exclude the NA group
    left_pos = left_count[0]
    right_pos = right_count[0]

    for i in range(1, max_groups + 1):
        lc = left_count[i]
        rc = right_count[i]

        if rc == 0:
            for j in range(lc):
                left_indexer[position + j] = left_pos + j
                right_indexer[position + j] = -1
            left_pos += lc
            position += lc
        else:
            for j in range(lc):
                offset = position + j * rc
                for k in range(rc):
                    left_indexer[offset + k] = left_pos + j
                    right_indexer[offset + k] = right_pos + k
            left_pos += lc
            right_pos += rc
            position += lc * rc

    return left_sorter, left_indexer, right_sorter, right_indexer


def full_outer_join(ndarray[int32_t] left, ndarray[int32_t] right):
    pass


@cython.boundscheck(False)
@cython.wraparound(False)
def groupsort_indexer(ndarray[int32_t] index, Py_ssize_t ngroups):
    cdef:
        Py_ssize_t i, loc, label, n
        ndarray[int32_t] counts, where, result

    # count group sizes, location 0 for NA
    counts = np.zeros(ngroups + 1, dtype='i4')
    n = len(index)
    for i from 0 <= i < n:
        counts[index[i] + 1] += 1

    # mark the start of each contiguous group of like-indexed data
    where = np.zeros(ngroups + 1, dtype='i4')
    for i from 1 <= i < ngroups + 1:
        where[i] = where[i - 1] + counts[i - 1]

    # this is our indexer
    result = np.zeros(n, dtype='i4')
    for i from 0 <= i < n:
        label = index[i] + 1
        result[where[label]] = i
        where[label] += 1

    return result, counts
