
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

    return labels

ctypedef double_t (* agg_func)(double_t *out, double_t *values,
                               int32_t *labels, int start, int end)

cdef agg_func get_agg_func(object how):
    pass

def group_aggregate(ndarray[double_t] values, list label_list,
                    how='add'):
    cdef:
        list sorted_labels
        ndarray[double_t] result
        double_t *resbuf
        double_t *valbuf
        Py_ssize_t i
        agg_func func

    values, sorted_labels = _group_reorder(values, label_list)
    result = np.empty(_result_shape(sorted_labels), dtype=np.float64)

    resbuf = <double_t*> result.data
    valbuf = <double_t*> values.data

    _aggregate_group(resbuf, valbuf, sorted_labels, 0, len(values),
                     func)

    return result, sorted_labels

cdef _aggregate_group(double_t *out, double_t *values, list labels,
                      int start, int end, agg_func func):
    cdef:
        ndarray[int32_t] axis0
        int32_t label_end

    axis0 = sorted_labels[0]
    label_end = axis0[0]
    if len(labels) == 1:
        func(out, values, axis0, start, end)
    else:
        # get group counts on axis
        edges = axis0.searchsorted(np.arange(1, label_end + 1), side='left')
        start = 0
        for end in edges:
            _aggregate_group(resbuf, valbuf, sorted_labels[1:],
                             start, end, func)
            start = end

cdef _group_add(double_t *out, double_t *values, int32_t *labels,
                int start, int end):
    cdef:
        int i
        int32_t lab
        double_t val

    for i from start <= i < end:


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

