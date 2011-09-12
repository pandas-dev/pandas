#-------------------------------------------------------------------------------
# Groupby-related functions

cdef inline _isnan(object o):
    return o != o

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

@cython.boundscheck(False)
def groupby(ndarray[object] index, ndarray[object] labels):
    cdef dict result = {}
    cdef ndarray[int8_t] mask
    cdef int i, length
    cdef list members
    cdef object idx, key

    length = len(index)
    mask = isnullobj(labels)

    for i from 0 <= i < length:
        if mask[i]:
            continue

        key = labels[i]
        idx = index[i]
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

    cdef int32_t **vecs = <int32_t **> malloc(nlevels * sizeof(int32_t*))
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
def group_labels2(ndarray[object] values):
    '''
    Compute label vector from input values and associated useful data

    Returns
    -------
    '''
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int32_t] labels = np.empty(n, dtype=np.int32)
        dict ids = {}
        dict reverse = {}
        int32_t idx
        object val
        int32_t count = 0

    for i from 0 <= i < n:
        val = values[i]

        # is NaN
        if val != val:
            labels[i] = -1
            continue

        if val in ids:
            idx = ids[val]
            labels[i] = idx
        else:
            ids[val] = count
            reverse[count] = val
            labels[i] = count
            count += 1

    return reverse, labels

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

@cython.wraparound(False)
@cython.boundscheck(False)
def fast_unique_multiple(list arrays):
    cdef:
        ndarray[object] buf
        Py_ssize_t k = len(arrays)
        Py_ssize_t i, j, n
        list uniques = []
        dict table = {}
        object val, stub = 0

    for i from 0 <= i < k:
        buf = arrays[i]
        n = len(buf)
        for j from 0 <= j < n:
            val = buf[j]
            if val not in table:
                table[val] = stub
                uniques.append(val)
    try:
        uniques.sort()
    except Exception:
        pass

    return uniques

# from libcpp.set cimport set as stlset

# cdef fast_unique_int32(ndarray arr):
#     cdef:
#         cdef stlset[int] table

#         Py_ssize_t i, n = len(arr)
#         int32_t* values
#         list uniques = []
#         int32_t val

#     values = <int32_t*> arr.data

#     for i from 0 <= i < n:
#         val = values[i]
#         if table.count(val) == 0:
#             table.insert(val)
#             uniques.append(val)
#     return np.asarray(sorted(uniques), dtype=object)

ctypedef double_t (* agg_func)(double_t *out, int32_t *counts, double_t *values,
                               int32_t *labels, int start, int end,
                               Py_ssize_t offset)

cdef agg_func get_agg_func(object how):
    if how == 'add':
        return _group_add
    elif how == 'mean':
        return _group_mean

@cython.boundscheck(False)
def group_aggregate(ndarray[double_t] values, list label_list,
                    object shape, how='add'):
    cdef:
        list sorted_labels
        ndarray result, counts
        agg_func func

    func = get_agg_func(how)

    values, sorted_labels = _group_reorder(values, label_list)
    result = np.empty(shape, dtype=np.float64)
    result.fill(nan)

    counts = np.zeros(shape, dtype=np.int32)

    if not values.flags.c_contiguous:
        values = values.copy()

    _aggregate_group(<double_t*> result.data,
                     <int32_t*> counts.data,
                     <double_t*> values.data,
                     sorted_labels, 0, len(values),
                     shape, 0, 0, func)

    return result, counts

def _group_reorder(values, label_list):
    indexer = np.lexsort(label_list[::-1])
    sorted_labels = [labels.take(indexer) for labels in label_list]
    sorted_values = values.take(indexer)
    return sorted_values, sorted_labels

cdef void _aggregate_group(double_t *out, int32_t *counts, double_t *values,
                           list labels, int start, int end, tuple shape,
                           Py_ssize_t which, Py_ssize_t offset,
                           agg_func func):
    cdef:
        ndarray[int32_t] axis
        cdef Py_ssize_t stride

    # time to actually aggregate
    if which == len(labels) - 1:
        axis = labels[which]

        while axis[start] == -1 and start < end:
            start += 1
        func(out, counts, values, <int32_t*> axis.data, start, end, offset)
    else:
        axis = labels[which][start:end]
        stride = np.prod(shape[which+1:])
        # get group counts on axisp
        edges = axis.searchsorted(np.arange(1, shape[which] + 1), side='left')
        # print edges, axis

        left = axis.searchsorted(0) # ignore NA values coded as -1

        # aggregate each subgroup
        for right in edges:
            _aggregate_group(out, counts, values, labels, start + left,
                             start + right, shape, which + 1, offset, func)
            offset += stride
            left = right

# TODO: aggregate multiple columns in single pass

cdef double_t _group_add(double_t *out, int32_t *counts, double_t *values,
                         int32_t *labels, int start, int end,
                         Py_ssize_t offset):
    cdef:
        Py_ssize_t i, it = start
        int32_t lab
        int32_t count = 0, tot = 0
        double_t val, cum = 0

    while it < end:
        i = labels[it]
        val = values[it]
        tot += 1

        # not nan
        if val == val:
            count += 1
            cum += val

        if it == end - 1 or labels[it + 1] > i:
            if count == 0:
                out[offset + i] = nan
            else:
                out[offset + i] = cum

            counts[offset + i] = tot

            count = 0
            cum = 0
            tot = 0

        it += 1

cdef double_t _group_mean(double_t *out, int32_t *counts, double_t *values,
                          int32_t *labels, int start, int end,
                          Py_ssize_t offset):
    cdef:
        Py_ssize_t i, it = start
        int32_t lab
        int32_t count = 0, tot = 0
        double_t val, cum = 0

    while it < end:
        i = labels[it]
        val = values[it]
        tot += 1

        # not nan
        if val == val:
            count += 1
            cum += val

        if it == end - 1 or labels[it + 1] > i:
            if count == 0:
                out[offset + i] = nan
            else:
                out[offset + i] = cum / count

            counts[offset + i] = tot

            count = 0
            cum = 0
            tot = 0

        it += 1

def _result_shape(label_list):
    # assumed sorted
    shape = []
    for labels in label_list:
        shape.append(1 + labels[-1])
    return tuple(shape)

def reduce_mean(ndarray[object] indices,
                ndarray[object] buckets,
                ndarray[double_t] values,
                inclusive=False):
    cdef:
        Py_ssize_t i, j, nbuckets, nvalues
        ndarray[double_t] output
        double_t the_sum, val, nobs



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

'''

def ts_upsample_mean(ndarray[object] indices,
                     ndarray[object] buckets,
                     ndarray[double_t] values,
                     inclusive=False):
    cdef:
        Py_ssize_t i, j, nbuckets, nvalues
        ndarray[double_t] output
        object next_bound
        double_t the_sum, val, nobs

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
        ndarray[double_t] output
        object next_bound
        double_t the_sum, val, nobs

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
                        ndarray[double_t] values,
                        object aggfunc,
                        inclusive=False):
    '''
    put something here
    '''
    cdef:
        Py_ssize_t i, j, jstart, nbuckets, nvalues
        ndarray[double_t] output
        object next_bound
        double_t the_sum, val, nobs

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

