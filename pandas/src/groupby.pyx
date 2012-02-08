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
        ndarray[int32_t] labels, counts, arr, seen
        int32_t loc
        dict ids = {}
        object val
        int32_t k

    ids, labels, counts = group_labels(values)
    seen = np.zeros_like(counts)

    # try not to get in trouble here...
    cdef int32_t **vecs = <int32_t **> malloc(len(ids) * sizeof(int32_t*))
    result = {}
    for i from 0 <= i < len(counts):
        arr = np.empty(counts[i], dtype=np.int32)
        result[ids[i]] = arr
        vecs[i] = <int32_t *> arr.data

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
        int32_t k, cur, pre

    nlevels = len(list_of_arrays)
    n = len(list_of_arrays[0])

    cdef int32_t **vecs = <int32_t**> malloc(nlevels * sizeof(int32_t*))
    for i from 0 <= i < nlevels:
        vecs[i] = <int32_t *> (<ndarray> list_of_arrays[i]).data
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
        ndarray[int32_t] labels = np.empty(n, dtype=np.int32)
        ndarray[int32_t] counts = np.empty(n, dtype=np.int32)
        dict ids = {}, reverse = {}
        int32_t idx
        object val
        int32_t count = 0

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
    cdef ndarray[int32_t] fillVec
    length = len(values)
    fillVec = np.empty(length, dtype=np.int32)
    for i from 0 <= i < length:
        idx = values[i]
        fillVec[i] = idMap[idx]

    return fillVec

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


# TODO: aggregate multiple columns in single pass

@cython.boundscheck(False)
@cython.wraparound(False)
def group_add(ndarray[float64_t, ndim=2] out,
              ndarray[int32_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int32_t] labels):
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
def group_mean(ndarray[float64_t, ndim=2] out,
               ndarray[int32_t] counts,
               ndarray[float64_t, ndim=2] values,
               ndarray[int32_t] labels):
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

@cython.boundscheck(False)
@cython.wraparound(False)
def group_var(ndarray[float64_t, ndim=2] out,
              ndarray[int32_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int32_t] labels):
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

def group_count(ndarray[int32_t] values, Py_ssize_t size):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int32_t] counts

    counts = np.zeros(size, dtype='i4')
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

def reduce_mean(ndarray[object] indices,
                ndarray[object] buckets,
                ndarray[float64_t] values,
                inclusive=False):
    cdef:
        Py_ssize_t i, j, nbuckets, nvalues
        ndarray[float64_t] output
        float64_t the_sum, val, nobs



    nbuckets = len(buckets)
    nvalues = len(indices)

    assert(len(values) == len(indices))

    output = np.empty(nbuckets, dtype=float)
    output.fill(np.NaN)

    j = 0
    for i from 0 <= i < nbuckets:
        next_bound = buckets[i]
        the_sum = 0
        nobs = 0
        if inclusive:
            while j < nvalues and indices[j] <= next_bound:
                val = values[j]
                # not NaN
                if val == val:
                    the_sum += val
                    nobs += 1
                j += 1
        else:
            while j < nvalues and indices[j] < next_bound:
                val = values[j]
                # not NaN
                if val == val:
                    the_sum += val
                    nobs += 1
                j += 1

        if nobs > 0:
            output[i] = the_sum / nobs

        if j >= nvalues:
            break

    return output

def _bucket_locs(index, buckets, inclusive=False):
    if inclusive:
        locs = index.searchsorted(buckets, side='left')
    else:
        locs = index.searchsorted(buckets, side='right')

    return locs

def count_level_1d(ndarray[uint8_t, cast=True] mask,
                   ndarray[int32_t] labels, Py_ssize_t max_bin):
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
                   ndarray[int32_t] labels, Py_ssize_t max_bin):
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

def duplicated(list values, take_last=False):
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


def generate_slices(ndarray[int32_t] labels, Py_ssize_t ngroups):
    cdef:
        Py_ssize_t i, group_size, n, lab, start
        object slobj
        ndarray[int32_t] starts

    n = len(labels)

    starts = np.zeros(ngroups, dtype='i4')
    ends = np.zeros(ngroups, dtype='i4')

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

'''

def ts_upsample_mean(ndarray[object] indices,
                     ndarray[object] buckets,
                     ndarray[float64_t] values,
                     inclusive=False):
    cdef:
        Py_ssize_t i, j, nbuckets, nvalues
        ndarray[float64_t] output
        object next_bound
        float64_t the_sum, val, nobs

    nbuckets = len(buckets)
    nvalues = len(indices)

    assert(len(values) == len(indices))

    output = np.empty(nbuckets, dtype=float)
    output.fill(np.NaN)

    j = 0
    for i from 0 <= i < nbuckets:
        next_bound = buckets[i]
        the_sum = 0
        nobs = 0
        if inclusive:
            while j < nvalues and indices[j] <= next_bound:
                val = values[j]
                # not NaN
                if val == val:
                    the_sum += val
                    nobs += 1
                j += 1
        else:
            while j < nvalues and indices[j] < next_bound:

    cdef:
        Py_ssize_t i, j, nbuckets, nvalues
        ndarray[float64_t] output
        object next_bound
        float64_t the_sum, val, nobs

    nbuckets = len(buckets)
    nvalues = len(indices)

    assert(len(values) == len(indices))

    output = np.empty(nbuckets, dtype=float)
    output.fill(np.NaN)

    j = 0
    for i from 0 <= i < nbuckets:
        next_bound = buckets[i]
        the_sum = 0
        nobs = 0
        if inclusive:
            while j < nvalues and indices[j] <= next_bound:
                val = values[j]
                # not NaN
                if val == val:
                    the_sum += val
                    nobs += 1
                j += 1
        else:
            while j < nvalues and indices[j] < next_bound:
                val = values[j]
                # not NaN
                if val == val:
                    the_sum += val
                    nobs += 1
                j += 1

        if nobs > 0:
            output[i] = the_sum / nobs

        if j >= nvalues:
            break

    return output
                val = values[j]
                # not NaN
                if val == val:
                    the_sum += val
                    nobs += 1
                j += 1

        if nobs > 0:
            output[i] = the_sum / nobs

        if j >= nvalues:
            break

    return output
'''

def ts_upsample_generic(ndarray[object] indices,
                        ndarray[object] buckets,
                        ndarray[float64_t] values,
                        object aggfunc,
                        inclusive=False):
    '''
    put something here
    '''
    cdef:
        Py_ssize_t i, j, jstart, nbuckets, nvalues
        ndarray[float64_t] output
        object next_bound
        float64_t the_sum, val, nobs

    nbuckets = len(buckets)
    nvalues = len(indices)

    assert(len(values) == len(indices))

    output = np.empty(nbuckets, dtype=float)
    output.fill(np.NaN)

    j = 0
    for i from 0 <= i < nbuckets:
        next_bound = buckets[i]
        the_sum = 0
        nobs = 0

        jstart = j
        if inclusive:
            while j < nvalues and indices[j] <= next_bound:
                j += 1
        else:
            while j < nvalues and indices[j] < next_bound:
                j += 1

        if nobs > 0:
            output[i] = aggfunc(values[jstart:j])

        if j >= nvalues:
            break

    return output

