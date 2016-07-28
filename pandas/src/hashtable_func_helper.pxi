"""
Template for each `dtype` helper function for hashtable

WARNING: DO NOT edit .pxi FILE directly, .pxi is generated from .pxi.in
"""

#----------------------------------------------------------------------
# VectorData
#----------------------------------------------------------------------


@cython.wraparound(False)
@cython.boundscheck(False)
cdef build_count_table_float64(float64_t[:] values,
                                 kh_float64_t *table, bint dropna):
    cdef:
        khiter_t k
        Py_ssize_t i, n = len(values)
        float64_t val
        int ret = 0

    with nogil:
        kh_resize_float64(table, n)

        for i in range(n):
            val = values[i]
            if val == val or not dropna:
                k = kh_get_float64(table, val)
                if k != table.n_buckets:
                    table.vals[k] += 1
                else:
                    k = kh_put_float64(table, val, &ret)
                    table.vals[k] = 1


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef value_count_float64(float64_t[:] values, bint dropna):
    cdef:
        Py_ssize_t i=0
        kh_float64_t *table
        float64_t[:] result_keys
        int64_t[:] result_counts
        int k

    table = kh_init_float64()
    build_count_table_float64(values, table, dropna)

    result_keys = np.empty(table.n_occupied, dtype=np.float64)
    result_counts = np.zeros(table.n_occupied, dtype=np.int64)

    with nogil:
        for k in range(table.n_buckets):
            if kh_exist_float64(table, k):
                result_keys[i] = table.keys[k]
                result_counts[i] = table.vals[k]
                i += 1
    kh_destroy_float64(table)

    return np.asarray(result_keys), np.asarray(result_counts)


@cython.wraparound(False)
@cython.boundscheck(False)
def duplicated_float64(float64_t[:] values,
                         object keep='first'):
    cdef:
        int ret = 0, k
        float64_t value
        Py_ssize_t i, n = len(values)
        kh_float64_t * table = kh_init_float64()
        ndarray[uint8_t, ndim=1, cast=True] out = np.empty(n, dtype='bool')

    kh_resize_float64(table, min(n, _SIZE_HINT_LIMIT))

    if keep not in ('last', 'first', False):
        raise ValueError('keep must be either "first", "last" or False')

    if keep == 'last':
        with nogil:
            for i from n > i >=0:
                kh_put_float64(table, values[i], &ret)
                out[i] = ret == 0
    elif keep == 'first':
        with nogil:
            for i from 0 <= i < n:
                kh_put_float64(table, values[i], &ret)
                out[i] = ret == 0
    else:
        with nogil:
            for i from 0 <= i < n:
                value = values[i]
                k = kh_get_float64(table, value)
                if k != table.n_buckets:
                    out[table.vals[k]] = 1
                    out[i] = 1
                else:
                    k = kh_put_float64(table, value, &ret)
                    table.keys[k] = value
                    table.vals[k] = i
                    out[i] = 0
    kh_destroy_float64(table)
    return out


@cython.wraparound(False)
@cython.boundscheck(False)
cdef build_count_table_int64(int64_t[:] values,
                                 kh_int64_t *table, bint dropna):
    cdef:
        khiter_t k
        Py_ssize_t i, n = len(values)
        int64_t val
        int ret = 0

    with nogil:
        kh_resize_int64(table, n)

        for i in range(n):
            val = values[i]
            if val == val or not dropna:
                k = kh_get_int64(table, val)
                if k != table.n_buckets:
                    table.vals[k] += 1
                else:
                    k = kh_put_int64(table, val, &ret)
                    table.vals[k] = 1


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef value_count_int64(int64_t[:] values, bint dropna):
    cdef:
        Py_ssize_t i=0
        kh_int64_t *table
        int64_t[:] result_keys
        int64_t[:] result_counts
        int k

    table = kh_init_int64()
    build_count_table_int64(values, table, dropna)

    result_keys = np.empty(table.n_occupied, dtype=np.int64)
    result_counts = np.zeros(table.n_occupied, dtype=np.int64)

    with nogil:
        for k in range(table.n_buckets):
            if kh_exist_int64(table, k):
                result_keys[i] = table.keys[k]
                result_counts[i] = table.vals[k]
                i += 1
    kh_destroy_int64(table)

    return np.asarray(result_keys), np.asarray(result_counts)


@cython.wraparound(False)
@cython.boundscheck(False)
def duplicated_int64(int64_t[:] values,
                         object keep='first'):
    cdef:
        int ret = 0, k
        int64_t value
        Py_ssize_t i, n = len(values)
        kh_int64_t * table = kh_init_int64()
        ndarray[uint8_t, ndim=1, cast=True] out = np.empty(n, dtype='bool')

    kh_resize_int64(table, min(n, _SIZE_HINT_LIMIT))

    if keep not in ('last', 'first', False):
        raise ValueError('keep must be either "first", "last" or False')

    if keep == 'last':
        with nogil:
            for i from n > i >=0:
                kh_put_int64(table, values[i], &ret)
                out[i] = ret == 0
    elif keep == 'first':
        with nogil:
            for i from 0 <= i < n:
                kh_put_int64(table, values[i], &ret)
                out[i] = ret == 0
    else:
        with nogil:
            for i from 0 <= i < n:
                value = values[i]
                k = kh_get_int64(table, value)
                if k != table.n_buckets:
                    out[table.vals[k]] = 1
                    out[i] = 1
                else:
                    k = kh_put_int64(table, value, &ret)
                    table.keys[k] = value
                    table.vals[k] = i
                    out[i] = 0
    kh_destroy_int64(table)
    return out
