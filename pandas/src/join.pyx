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

    if not sort:
        if left_sorter.dtype != np.int_:
            left_sorter = left_sorter.astype(np.int_)

        rev = np.empty(len(left), dtype=np.int_)
        rev.put(left_sorter, np.arange(len(left)))

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
    if indexer.dtype != np.int_:
        indexer = indexer.astype(np.int_)
    if len(sorter) > 0:
        res = sorter.take(indexer)
        np.putmask(res, indexer == -1, -1)
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


def ffill_by_group(ndarray[int64_t] indexer, ndarray[int64_t] group_ids,
                   int64_t max_group):
    cdef:
        Py_ssize_t i, n = len(indexer)
        ndarray[int64_t] result, last_obs
        int64_t gid, val

    result = np.empty(n, dtype=np.int64)

    last_obs = np.empty(max_group, dtype=np.int64)
    last_obs.fill(-1)

    for i in range(n):
        gid = group_ids[i]
        val = indexer[i]
        if val == -1:
            result[i] = last_obs[gid]
        else:
            result[i] = val
            last_obs[gid] = val

    return result

