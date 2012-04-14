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
        ndarray arr

    nlevels = len(list_of_arrays)
    n = len(list_of_arrays[0])

    cdef int32_t **vecs = <int32_t**> malloc(nlevels * sizeof(int32_t*))
    for i from 0 <= i < nlevels:
        # vecs[i] = <int32_t *> (<ndarray> list_of_arrays[i]).data

        arr = list_of_arrays[i]
        vecs[i] = <int32_t *> arr.data
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
def group_prod(ndarray[float64_t, ndim=2] out,
               ndarray[int32_t] counts,
               ndarray[float64_t, ndim=2] values,
               ndarray[int32_t] labels):
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


@cython.boundscheck(False)
@cython.wraparound(False)
def group_min(ndarray[float64_t, ndim=2] out,
              ndarray[int32_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int32_t] labels):
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
              ndarray[int32_t] counts,
              ndarray[float64_t, ndim=2] values,
              ndarray[int32_t] labels):
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

# TODO: could do even better if we know something about the data. eg, index has
# 1-min data, binner has 5-min data, then  bins are just strides in index. This
# is a general, O(max(len(values), len(binner))) method.

@cython.boundscheck(False)
@cython.wraparound(False)
def generate_bins_dt64(ndarray[int64_t] values, ndarray[int64_t] binner,
                       object closed='left', object label='left'):
    """
    Int64 (datetime64) version of generic python version in groupby.py
    """
    cdef:
        Py_ssize_t lenidx, lenbin, i, j, bc, vc
        ndarray[int64_t] labels
        ndarray[int32_t] bins
        int64_t l_bin, r_bin

    lenidx = len(values)
    lenbin = len(binner)

    if lenidx <= 0 or lenbin <= 0:
        raise ValueError("Invalid length for values or for binner")

    # check binner fits data
    if values[0] < binner[0]:
        raise ValueError("Values falls before first bin")

    if values[lenidx-1] > binner[lenbin-1]:
        raise ValueError("Values falls after last bin")

    labels = np.empty(lenbin, dtype=np.int64)
    bins   = np.empty(lenbin, dtype=np.int32)

    j  = 0 # index into values
    bc = 0 # bin count
    vc = 0 # value count

    # linear scan, presume nothing about values/binner except that it
    # fits ok
    for i in range(0, lenbin-1):
        l_bin = binner[i]
        r_bin = binner[i+1]

        # set label of bin
        if label == 'left':
            labels[bc] = l_bin
        else:
            labels[bc] = r_bin

        # count values in current bin, advance to next bin
        while values[j] < r_bin or closed == 'right' and values[j] == r_bin:
            j += 1
            vc += 1
            if j >= lenidx:
                break

        # check we have data left to scan
        if j >= lenidx:
            break

        # if we've seen some values, mark bin
        if vc != 0:
            bins[bc] = j
            bc += 1
            vc = 0

    labels = np.resize(labels, bc + 1)
    bins = np.resize(bins, bc)

    return bins, labels

# add passing bin edges, instead of labels

@cython.boundscheck(False)
@cython.wraparound(False)
def group_add_bin(ndarray[float64_t, ndim=2] out,
                  ndarray[int32_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int32_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            if b < ngroups - 1 and i >= bins[b]:
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
            if b < ngroups - 1 and i >= bins[b]:
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
                  ndarray[int32_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int32_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] prodx, nobs

    nobs = np.zeros_like(out)
    prodx = np.ones_like(out)

    ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            if b < ngroups - 1 and i >= bins[b]:
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
            if b < ngroups - 1 and i >= bins[b]:
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
                   ndarray[int32_t] counts,
                   ndarray[float64_t, ndim=2] values,
                   ndarray[int32_t] bins):
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


    ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            if b < ngroups - 1 and i >= bins[b]:
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
            if b < ngroups - 1 and i >= bins[b]:
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
                  ndarray[int32_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int32_t] bins):
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

    ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            if b < ngroups - 1 and i >= bins[b]:
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
            if b < ngroups - 1 and i >= bins[b]:
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
                  ndarray[int32_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int32_t] bins):
    '''
    Only aggregates on axis=0
    '''
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        float64_t vopen, vhigh, vlow, vclose, NA
        bint got_first = 0

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
            if b < ngroups - 1 and i >= bins[b]:
                if not got_first:
                    out[b, 0] = NA
                    out[b, 1] = NA
                    out[b, 2] = NA
                    out[b, 3] = NA
                else:
                    out[b, 0] = vopen
                    out[b, 1] = vlow
                    out[b, 2] = vhigh
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
            out[b, 1] = vlow
            out[b, 2] = vhigh
            out[b, 3] = vclose


@cython.boundscheck(False)
@cython.wraparound(False)
def group_mean_bin(ndarray[float64_t, ndim=2] out,
                   ndarray[int32_t] counts,
                   ndarray[float64_t, ndim=2] values,
                   ndarray[int32_t] bins):
    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            if b < ngroups - 1 and i >= bins[b]:
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
            if b < ngroups - 1 and i >= bins[b]:
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
                  ndarray[int32_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int32_t] bins):

    cdef:
        Py_ssize_t i, j, N, K, ngroups, b
        float64_t val, ct
        ndarray[float64_t, ndim=2] nobs, sumx, sumxx

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)
    sumxx = np.zeros_like(out)

    ngroups = len(bins) + 1
    N, K = (<object> values).shape

    b = 0
    if K > 1:
        for i in range(N):
            if b < ngroups - 1 and i >= bins[b]:
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
            if b < ngroups - 1 and i >= bins[b]:
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


def groupby_arrays(ndarray index, ndarray _labels):
    cdef:
        Py_ssize_t i, lab, cur, start, n = len(index)
        ndarray[int32_t] labels
        dict result = {}

    if _labels.dtype == np.int32:
        labels = _labels
    else:
        labels = _labels.astype(np.int32)

    index = np.asarray(index)

    # this is N log N. If this is a bottleneck may we worth fixing someday
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

    return result
