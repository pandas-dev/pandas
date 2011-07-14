
#-------------------------------------------------------------------------------
# Groupby-related functions

cdef inline _isnan(object o):
    return o != o

@cython.boundscheck(False)
def arrmap(ndarray[object, ndim=1] index, object func):
    cdef int length = index.shape[0]
    cdef int i = 0

    cdef ndarray[object, ndim=1] result = np.empty(length, dtype=np.object_)

    for i from 0 <= i < length:
        result[i] = func(index[i])

    return result

@cython.boundscheck(False)
def groupby(object index, object mapper, output=None):
    cdef dict result
    cdef ndarray[object, ndim=1] mapped_index
    cdef ndarray[object, ndim=1] index_buf
    cdef ndarray[int8_t, ndim=1] mask
    cdef int i, length
    cdef list members
    cdef object idx, key

    length = len(index)

    if output is None:
        result = {}
    else:
        result = <dict> output

    index_buf = np.asarray(index)
    mapped_index = arrmap(index_buf, mapper)
    mask = isnullobj(mapped_index)

    # nullkeys = index_buf[mask.astype(bool)]
    # if len(nullkeys) > 0:
    #     result[np.NaN] = nullkeys

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
def groupby_indices(object index, object mapper):
    cdef dict result
    cdef ndarray[object, ndim=1] mapped_index
    cdef ndarray[int8_t, ndim=1] mask
    cdef int i, length
    cdef list members, null_list
    cdef object key

    length = len(index)

    result = {}
    index = np.asarray(index)
    mapped_index = arrmap(index, mapper)
    mask = isnullobj(mapped_index)
    for i from 0 <= i < length:
        if mask[i]:
            continue

        key = mapped_index[i]
        if key in result:
            (<list> result[key]).append(i)
        else:
            result[key] = [i]

    return result

def group_labels(ndarray[object] values):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int32_t] labels = np.empty(n, dtype=np.int32)
        dict ids = {}
        object val
        int32_t count = 0

    for i from 0 <= i < n:
        val = values[i]
        try:
            labels[i] = ids[val]
        except KeyError:
            ids[val] = count
            labels[i] = count
            count += 1

    return count, labels

def labelize(*key_arrays):
    shape = []
    labels = []
    for key_arr in key_arrays:
        ct, lab = group_labels(key_arrays)
        shape.append(ct)
        labels.append(lab)

    return tuple(shape), labels

ctypedef double_t (* agg_func)(double_t *out, int32_t *counts,
                               double_t *values, int32_t *labels,
                               int start, int end)

cdef agg_func get_agg_func(object how):
    if how == 'add':
        return _group_add

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
    counts = = np.zeros(shape, dtype=np.int32)

    _aggregate_group(<double_t*> result.data,
                     <double_t*> values.data,
                     <int32_t*> counts.data
                      sorted_labels, 0, len(values),
                     shape, 0, func)

    return result, sorted_labels

cdef _aggregate_group(double_t *out, int32_t *counts, double_t *values,
                      list labels, int start, int end, tuple shape,
                      Py_ssize_t which, agg_func func):
    cdef:
        ndarray[int32_t] axis0
        int32_t label_end

    axis0 = labels[which][start:end]
    label_end = shape[which]

    # time to actually aggregate
    if which == len(labels) - 1:
        func(out, values, axis0, start, end)
    else:
        # get group counts on axis
        edges = axis0.searchsorted(np.arange(1, label_end + 1), side='left')
        start = 0
        # aggregate each subgroup
        for end in edges:
            _aggregate_group(out, counts, values, sorted_labels[1:],
                             start, end, shape, which + 1, func)
            start = end

cdef _group_add(double_t *out, int32_t *counts, double_t *values,
                int32_t *labels, int start, int end, int rng):
    cdef:
        Py_ssize_t i, it = start
        int32_t lab
        int32_t count = 0
        double_t val, cum = 0

    for i in range(rng):
        while it < end:
            if labels[it] > i:
                counts[i] = count
                out[i] = cum
                break

            val = values[it]
            # not nan
            if val == val:
                count += 1
                cum += val

        count = 0
        cum = 0

def _group_reorder(values, label_list)
    indexer = np.lexsort(label_list[::-1])
    sorted_labels = [labels.take(indexer) for labels in label_list]
    sorted_values = values.take(indexer)
    return sorted_values, sorted_labels

def _result_shape(label_list):
    # assumed sorted
    shape = []
    for labels in label_list:
        size.append(1 + labels[-1])
    return tuple(shape)

def reduce_mean(ndarray[object, ndim=1] indices,
                ndarray[object, ndim=1] buckets,
                ndarray[double_t, ndim=1] values,
                inclusive=False):
    cdef:
        Py_ssize_t i, j, nbuckets, nvalues
        ndarray[double_t, ndim=1] output
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

def ts_upsample_mean(ndarray[object, ndim=1] indices,
                     ndarray[object, ndim=1] buckets,
                     ndarray[double_t, ndim=1] values,
                     inclusive=False):
    cdef:
        Py_ssize_t i, j, nbuckets, nvalues
        ndarray[double_t, ndim=1] output
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
        ndarray[double_t, ndim=1] output
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

def ts_upsample_generic(ndarray[object, ndim=1] indices,
                        ndarray[object, ndim=1] buckets,
                        ndarray[double_t, ndim=1] values,
                        object aggfunc,
                        inclusive=False):
    '''
    put something here
    '''
    cdef:
        Py_ssize_t i, j, jstart, nbuckets, nvalues
        ndarray[double_t, ndim=1] output
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

