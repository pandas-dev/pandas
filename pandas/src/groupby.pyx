#-------------------------------------------------------------------------------
# Groupby-related functions

@cython.boundscheck(False)
def arrmap(ndarray[object] index, object func):
    cdef int length = index.shape[0]
    cdef int i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    for i from 0 <= i < length:
        result[i] = func(index[i])

    return result

@cython.boundscheck(False)
def groupby_func(object index, object mapper):
    cdef dict result = {}
    cdef ndarray[object] mapped_index
    cdef ndarray[object] index_buf
    cdef ndarray[int8_t] mask
    cdef int i, length
    cdef list members
    cdef object idx, key

    length = len(index)

    index_buf = np.asarray(index)
    mapped_index = arrmap(index_buf, mapper)
    mask = isnullobj(mapped_index)

    for i from 0 <= i < length:
        if mask[i]:
            continue

        key = mapped_index[i]
        idx = index_buf[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result


def func_groupby_indices(object index, object mapper):
    return groupby_indices_naive(arrmap(index, mapper))

@cython.boundscheck(False)
cpdef groupby_indices_naive(ndarray[object] values):
    cdef dict result
    cdef ndarray[int8_t] mask
    cdef Py_ssize_t i, length = len(values)
    cdef object key

    result = {}
    mask = isnullobj(values)
    for i from 0 <= i < length:
        if mask[i]:
            continue

        key = values[i]
        if key in result:
            (<list> result[key]).append(i)
        else:
            result[key] = [i]

    return result

@cython.boundscheck(False)
def groupby_indices(ndarray values):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int64_t] labels, counts, arr, seen
        int64_t loc
        dict ids = {}
        object val
        int64_t k

    ids, labels, counts = group_labels(values)
    seen = np.zeros_like(counts)

    # try not to get in trouble here...
    cdef int64_t **vecs = <int64_t **> malloc(len(ids) * sizeof(int64_t*))
    result = {}
    for i from 0 <= i < len(counts):
        arr = np.empty(counts[i], dtype=np.int64)
        result[ids[i]] = arr
        vecs[i] = <int64_t *> arr.data

    for i from 0 <= i < n:
        k = labels[i]

        # was NaN
        if k == -1:
            continue

        loc = seen[k]
        vecs[k][loc] = i
        seen[k] = loc + 1

    free(vecs)

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def is_lexsorted(list list_of_arrays):
    cdef:
        int i
        Py_ssize_t n, nlevels
        int64_t k, cur, pre
        ndarray arr

    nlevels = len(list_of_arrays)
    n = len(list_of_arrays[0])

    cdef int64_t **vecs = <int64_t**> malloc(nlevels * sizeof(int64_t*))
    for i from 0 <= i < nlevels:
        # vecs[i] = <int64_t *> (<ndarray> list_of_arrays[i]).data

        arr = list_of_arrays[i]
        vecs[i] = <int64_t *> arr.data
    # assume uniqueness??

    for i from 1 <= i < n:
        for k from 0 <= k < nlevels:
            cur = vecs[k][i]
            pre = vecs[k][i-1]
            if cur == pre:
                continue
            elif cur > pre:
                break
            else:
                return False
    free(vecs)
    return True

@cython.wraparound(False)
@cython.boundscheck(False)
def group_labels(ndarray[object] values):
    '''
    Compute label vector from input values and associated useful data

    Returns
    -------
    '''
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int64_t] labels = np.empty(n, dtype=np.int64)
        ndarray[int64_t] counts = np.empty(n, dtype=np.int64)
        dict ids = {}, reverse = {}
        int64_t idx
        object val
        int64_t count = 0

    for i from 0 <= i < n:
        val = values[i]

        # is NaN
        if val != val:
            labels[i] = -1
            continue

        # for large number of groups, not doing try: except: makes a big
        # difference
        if val in ids:
            idx = ids[val]
            labels[i] = idx
            counts[idx] = counts[idx] + 1
        else:
            ids[val] = count
            reverse[count] = val
            labels[i] = count
            counts[count] = 1
            count += 1

    return reverse, labels, counts[:count].copy()


@cython.wraparound(False)
@cython.boundscheck(False)
def get_unique_labels(ndarray[object] values, dict idMap):
    cdef int i, length
    cdef object idx
    cdef ndarray[int64_t] fillVec
    length = len(values)
    fillVec = np.empty(length, dtype=np.int64)
    for i from 0 <= i < length:
        idx = values[i]
        fillVec[i] = idMap[idx]

    return fillVec

@cython.boundscheck(False)
@cython.wraparound(False)
def groupsort_indexer(ndarray[int64_t] index, Py_ssize_t ngroups):
    cdef:
        Py_ssize_t i, loc, label, n
        ndarray[int64_t] counts, where, result

    # count group sizes, location 0 for NA
    counts = np.zeros(ngroups + 1, dtype=np.int64)
    n = len(index)
    for i from 0 <= i < n:
        counts[index[i] + 1] += 1

    # mark the start of each contiguous group of like-indexed data
    where = np.zeros(ngroups + 1, dtype=np.int64)
    for i from 1 <= i < ngroups + 1:
        where[i] = where[i - 1] + counts[i - 1]

    # this is our indexer
    result = np.zeros(n, dtype=np.int64)
    for i from 0 <= i < n:
        label = index[i] + 1
        result[where[label]] = i
        where[label] += 1

    return result, counts

# TODO: aggregate multiple columns in single pass

@cython.boundscheck(False)
@cython.wraparound(False)
def group_add(ndarray[float64_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

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


@cython.boundscheck(False)
@cython.wraparound(False)
def group_prod(ndarray[float64_t, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[float64_t, ndim=2] values,
               ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] prodx, nobs

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

#----------------------------------------------------------------------
# first, nth, last

@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth(ndarray[float64_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int64_t] labels, int64_t rank):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

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

@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth_object(ndarray[object, ndim=2] out,
                     ndarray[int64_t] counts,
                     ndarray[object, ndim=2] values,
                     ndarray[int64_t] labels,
                     int64_t rank):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        object val
        float64_t count
        ndarray[int64_t, ndim=2] nobs
        ndarray[object, ndim=2] resx

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty((<object> out).shape, dtype=object)

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
                out[i, j] = <object> nan
            else:
                out[i, j] = resx[i, j]

@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth_bin(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins, int64_t rank):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] resx, nobs

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

@cython.boundscheck(False)
@cython.wraparound(False)
def group_nth_bin_object(ndarray[object, ndim=2] out,
                         ndarray[int64_t] counts,
                         ndarray[object, ndim=2] values,
                         ndarray[int64_t] bins, int64_t rank):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        object val
        float64_t count
        ndarray[object, ndim=2] resx
        ndarray[float64_t, ndim=2] nobs

    nobs = np.zeros((<object> out).shape, dtype=np.float64)
    resx = np.empty((<object> out).shape, dtype=object)

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

@cython.boundscheck(False)
@cython.wraparound(False)
def group_last(ndarray[float64_t, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[float64_t, ndim=2] values,
               ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

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

@cython.boundscheck(False)
@cython.wraparound(False)
def group_last_object(ndarray[object, ndim=2] out,
                      ndarray[int64_t] counts,
                      ndarray[object, ndim=2] values,
                      ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        object val
        float64_t count
        ndarray[object, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty((<object> out).shape, dtype=object)

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


@cython.boundscheck(False)
@cython.wraparound(False)
def group_last_bin(ndarray[float64_t, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[float64_t, ndim=2] values,
                   ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] resx, nobs

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

@cython.boundscheck(False)
@cython.wraparound(False)
def group_last_bin_object(ndarray[object, ndim=2] out,
                          ndarray[int64_t] counts,
                          ndarray[object, ndim=2] values,
                          ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        object val
        float64_t count
        ndarray[object, ndim=2] resx
        ndarray[float64_t, ndim=2] nobs

    nobs = np.zeros((<object> out).shape, dtype=np.float64)
    resx = np.empty((<object> out).shape, dtype=object)

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

#----------------------------------------------------------------------
# group_min, group_max


@cython.boundscheck(False)
@cython.wraparound(False)
def group_min(ndarray[float64_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] minx, nobs

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


@cython.boundscheck(False)
@cython.wraparound(False)
def group_max(ndarray[float64_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] maxx, nobs

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


@cython.boundscheck(False)
@cython.wraparound(False)
def group_mean(ndarray[float64_t, ndim=2] out,
               ndarray[int64_t] counts,
               ndarray[float64_t, ndim=2] values,
               ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

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


def group_median(ndarray[float64_t, ndim=2] out,
                 ndarray[int64_t] counts,
                 ndarray[float64_t, ndim=2] values,
                 ndarray[int64_t] labels):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, size
        ndarray[int64_t] _counts
        ndarray data
        float64_t* ptr

    from pandas._algos import take_2d_axis1_float64

    ngroups = len(counts)
    N, K = (<object> values).shape

    indexer, _counts = groupsort_indexer(labels, ngroups)
    counts[:] = _counts[1:]

    data = np.empty((K, N), dtype=np.float64)
    ptr = <float64_t*> data.data

    take_2d_axis1_float64(values.T, indexer, out=data)

    for i in range(K):
        # exclude NA group
        ptr += _counts[0]
        for j in range(ngroups):
            size = _counts[j + 1]
            out[j, i] = _median_linear(ptr, size)
            ptr += size


cdef inline float64_t _median_linear(float64_t* a, int n):
    cdef int i, j, na_count = 0
    cdef float64_t result
    cdef float64_t* tmp

    # count NAs
    for i in range(n):
        if a[i] != a[i]:
            na_count += 1

    if na_count:
        if na_count == n:
            return NaN

        tmp = <float64_t*> malloc((n - na_count) * sizeof(float64_t))

        j = 0
        for i in range(n):
            if a[i] == a[i]:
                tmp[j] = a[i]
                j += 1

        a = tmp
        n -= na_count


    if n % 2:
        result = kth_smallest_c(a, n / 2, n)
    else:
        result = (kth_smallest_c(a, n / 2, n) +
                  kth_smallest_c(a, n / 2 - 1, n)) / 2

    if na_count:
        free(a)

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def group_var(ndarray[float64_t, ndim=2] out,
              ndarray[int64_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, ct
        ndarray[float64_t, ndim=2] nobs, sumx, sumxx

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

# TODO: could do even better if we know something about the data. eg, index has
# 1-min data, binner has 5-min data, then  bins are just strides in index. This
# is a general, O(max(len(values), len(binner))) method.

@cython.boundscheck(False)
@cython.wraparound(False)
def generate_bins_dt64(ndarray[int64_t] values, ndarray[int64_t] binner,
                       object closed='left'):
    """
    Int64 (datetime64) version of generic python version in groupby.py
    """
    cdef:
        Py_ssize_t lenidx, lenbin, i, j, bc, vc
        ndarray[int64_t] bins
        int64_t l_bin, r_bin
        bint right_closed = closed == 'right'

    lenidx = len(values)
    lenbin = len(binner)

    if lenidx <= 0 or lenbin <= 0:
        raise ValueError("Invalid length for values or for binner")

    # check binner fits data
    if values[0] < binner[0]:
        raise ValueError("Values falls before first bin")

    if values[lenidx-1] > binner[lenbin-1]:
        raise ValueError("Values falls after last bin")

    bins   = np.empty(lenbin - 1, dtype=np.int64)

    j  = 0 # index into values
    bc = 0 # bin count

    # linear scan
    for i in range(0, lenbin - 1):
        l_bin = binner[i]
        r_bin = binner[i+1]

        # count values in current bin, advance to next bin
        while j < lenidx and (values[j] < r_bin or
                              (right_closed and values[j] == r_bin)):
            j += 1

        bins[bc] = j
        bc += 1

    return bins

# add passing bin edges, instead of labels

@cython.boundscheck(False)
@cython.wraparound(False)
def group_add_bin(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b, nbins
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

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

@cython.boundscheck(False)
@cython.wraparound(False)
def group_prod_bin(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] prodx, nobs

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

@cython.boundscheck(False)
@cython.wraparound(False)
def group_min_bin(ndarray[float64_t, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[float64_t, ndim=2] values,
                   ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] minx, nobs

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


@cython.boundscheck(False)
@cython.wraparound(False)
def group_max_bin(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] maxx, nobs

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


@cython.boundscheck(False)
@cython.wraparound(False)
def group_ohlc(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        float64_t vopen, vhigh, vlow, vclose, NA
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


# @cython.boundscheck(False)
# @cython.wraparound(False)
def group_mean_bin(ndarray[float64_t, ndim=2] out,
                   ndarray[int64_t] counts,
                   ndarray[float64_t, ndim=2] values,
                   ndarray[int64_t] bins):
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

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

@cython.boundscheck(False)
@cython.wraparound(False)
def group_var_bin(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] bins):

    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, ct
        ndarray[float64_t, ndim=2] nobs, sumx, sumxx

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



@cython.boundscheck(False)
@cython.wraparound(False)
def row_bool_subset(ndarray[float64_t, ndim=2] values,
                    ndarray[uint8_t, cast=True] mask):
    cdef:
        Py_ssize_t i, j, n, k, pos = 0
        ndarray[float64_t, ndim=2] out

    n, k = (<object> values).shape
    assert(n == len(mask))

    out = np.empty((mask.sum(), k), dtype=np.float64)

    for i in range(n):
        if mask[i]:
            for j in range(k):
                out[pos, j] = values[i, j]
            pos += 1

    return out

@cython.boundscheck(False)
@cython.wraparound(False)
def row_bool_subset_object(ndarray[object, ndim=2] values,
                           ndarray[uint8_t, cast=True] mask):
    cdef:
        Py_ssize_t i, j, n, k, pos = 0
        ndarray[object, ndim=2] out

    n, k = (<object> values).shape
    assert(n == len(mask))

    out = np.empty((mask.sum(), k), dtype=object)

    for i in range(n):
        if mask[i]:
            for j in range(k):
                out[pos, j] = values[i, j]
            pos += 1

    return out


def group_count(ndarray[int64_t] values, Py_ssize_t size):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int64_t] counts

    counts = np.zeros(size, dtype=np.int64)
    for i in range(n):
        counts[values[i]] += 1
    return counts

def lookup_values(ndarray[object] values, dict mapping):
    cdef:
        Py_ssize_t i, n = len(values)

    result = np.empty(n, dtype='O')
    for i in range(n):
        result[i] = mapping[values[i]]
    return maybe_convert_objects(result)


def count_level_1d(ndarray[uint8_t, cast=True] mask,
                   ndarray[int64_t] labels, Py_ssize_t max_bin):
    cdef:
        Py_ssize_t i, n
        ndarray[int64_t] counts

    counts = np.zeros(max_bin, dtype='i8')

    n = len(mask)

    for i from 0 <= i < n:
        if mask[i]:
            counts[labels[i]] += 1

    return counts


def count_level_2d(ndarray[uint8_t, ndim=2, cast=True] mask,
                   ndarray[int64_t] labels, Py_ssize_t max_bin):
    cdef:
        Py_ssize_t i, j, k, n
        ndarray[int64_t, ndim=2] counts

    n, k = (<object> mask).shape
    counts = np.zeros((max_bin, k), dtype='i8')

    for i from 0 <= i < n:
        for j from 0 <= j < k:
            if mask[i, j]:
                counts[labels[i], j] += 1

    return counts

cdef class _PandasNull:

    def __richcmp__(_PandasNull self, object other, int op):
        if op == 2: # ==
            return isinstance(other, _PandasNull)
        elif op == 3: # !=
            return not isinstance(other, _PandasNull)
        else:
            return False

    def __hash__(self):
        return 0

pandas_null = _PandasNull()

def fast_zip_fillna(list ndarrays, fill_value=pandas_null):
    '''
    For zipping multiple ndarrays into an ndarray of tuples
    '''
    cdef:
        Py_ssize_t i, j, k, n
        ndarray[object] result
        flatiter it
        object val, tup

    k = len(ndarrays)
    n = len(ndarrays[0])

    result = np.empty(n, dtype=object)

    # initialize tuples on first pass
    arr = ndarrays[0]
    it = <flatiter> PyArray_IterNew(arr)
    for i in range(n):
        val = PyArray_GETITEM(arr, PyArray_ITER_DATA(it))
        tup = PyTuple_New(k)

        if val != val:
            val = fill_value

        PyTuple_SET_ITEM(tup, 0, val)
        Py_INCREF(val)
        result[i] = tup
        PyArray_ITER_NEXT(it)

    for j in range(1, k):
        arr = ndarrays[j]
        it = <flatiter> PyArray_IterNew(arr)
        if len(arr) != n:
            raise ValueError('all arrays must be same length')

        for i in range(n):
            val = PyArray_GETITEM(arr, PyArray_ITER_DATA(it))
            if val != val:
                val = fill_value

            PyTuple_SET_ITEM(result[i], j, val)
            Py_INCREF(val)
            PyArray_ITER_NEXT(it)

    return result

def duplicated(ndarray[object] values, take_last=False):
    cdef:
        Py_ssize_t i, n
        dict seen = {}
        object row

    n = len(values)
    cdef ndarray[uint8_t] result = np.zeros(n, dtype=np.uint8)

    if take_last:
        for i from n > i >= 0:
            row = values[i]

            if row in seen:
                result[i] = 1
            else:
                seen[row] = None
                result[i] = 0
    else:
        for i from 0 <= i < n:
            row = values[i]
            if row in seen:
                result[i] = 1
            else:
                seen[row] = None
                result[i] = 0

    return result.view(np.bool_)

def generate_slices(ndarray[int64_t] labels, Py_ssize_t ngroups):
    cdef:
        Py_ssize_t i, group_size, n, lab, start
        object slobj
        ndarray[int64_t] starts

    n = len(labels)

    starts = np.zeros(ngroups, dtype=np.int64)
    ends = np.zeros(ngroups, dtype=np.int64)

    start = 0
    group_size = 0
    for i in range(n):
        group_size += 1
        lab = labels[i]
        if i == n - 1 or lab != labels[i + 1]:
            starts[lab] = start
            ends[lab] = start + group_size
            start += group_size
            group_size = 0

    return starts, ends


def groupby_arrays(ndarray index, ndarray[int64_t] labels, sort=True):
    cdef:
        Py_ssize_t i, lab, cur, start, n = len(index)
        dict result = {}

    index = np.asarray(index)

    # this is N log N. If this is a bottleneck may we worth fixing someday
    if sort:
        indexer = labels.argsort(kind='mergesort')

        labels = labels.take(indexer)
        index = index.take(indexer)

    if n == 0:
        return result

    start = 0
    cur = labels[0]
    for i in range(1, n):
        lab = labels[i]

        if lab != cur:
            if lab != -1:
                result[cur] = index[start:i]
            start = i
        cur = lab

    result[cur] = index[start:]
    return result

def indices_fast(object index, ndarray[int64_t] labels, list keys,
                 list sorted_labels):
    cdef:
        Py_ssize_t i, j, k, lab, cur, start, n = len(labels)
        dict result = {}
        object tup

    k = len(keys)

    if n == 0:
        return result

    start = 0
    cur = labels[0]
    for i in range(1, n):
        lab = labels[i]

        if lab != cur:
            if lab != -1:
                tup = PyTuple_New(k)
                for j in range(k):
                    val = util.get_value_at(keys[j],
                                            sorted_labels[j][i-1])
                    PyTuple_SET_ITEM(tup, j, val)
                    Py_INCREF(val)

                result[tup] = index[start:i]
            start = i
        cur = lab

    tup = PyTuple_New(k)
    for j in range(k):
        val = util.get_value_at(keys[j],
                                sorted_labels[j][n - 1])
        PyTuple_SET_ITEM(tup, j, val)
        Py_INCREF(val)
    result[tup] = index[start:]

    return result
