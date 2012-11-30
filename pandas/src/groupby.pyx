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
