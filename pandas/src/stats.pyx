cdef float64_t FP_ERR = 1e-13

cimport util

def rank_1d_float64(object in_arr):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, n, dups = 0
        ndarray[float64_t] sorted_data, ranks, values
        ndarray[int64_t] argsorted
        int32_t idx
        float64_t val, nan_value
        float64_t sum_ranks = 0

    values = np.asarray(in_arr).copy()

    nan_value = np.inf
    mask = np.isnan(values)
    np.putmask(values, mask, nan_value)

    n = len(values)
    ranks = np.empty(n, dtype='f8')

    # py2.5/win32 hack, can't pass i8
    _as = values.argsort()
    sorted_data = values.take(_as)

    argsorted = _as.astype('i8')

    for i in range(n):
        sum_ranks += i + 1
        dups += 1
        val = sorted_data[i]
        if val == nan_value:
            ranks[argsorted[i]] = nan
            continue
        if i == n - 1 or fabs(sorted_data[i + 1] - val) > FP_ERR:
            for j in range(i - dups + 1, i + 1):
                ranks[argsorted[j]] = sum_ranks / dups
            sum_ranks = dups = 0
    return ranks

def rank_2d_float64(object in_arr, axis=0):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, z, k, n, dups = 0
        ndarray[float64_t, ndim=2] ranks, values
        ndarray[int64_t, ndim=2] argsorted
        int32_t idx
        float64_t val, nan_value
        float64_t sum_ranks = 0

    in_arr = np.asarray(in_arr)

    if axis == 0:
        values = in_arr.T.copy()
    else:
        values = in_arr.copy()

    nan_value = np.inf
    np.putmask(values, np.isnan(values), nan_value)

    n, k = (<object> values).shape
    ranks = np.empty((n, k), dtype='f8')
    argsorted = values.argsort(1).astype('i8')
    values.sort(axis=1)

    for i in range(n):
        dups = sum_ranks = 0
        for j in range(k):
            sum_ranks += j + 1
            dups += 1
            val = values[i, j]
            if val == nan_value:
                ranks[i, argsorted[i, j]] = nan
                continue
            if j == k - 1 or fabs(values[i, j + 1] - val) > FP_ERR:
                for z in range(j - dups + 1, j + 1):
                    ranks[i, argsorted[i, z]] = sum_ranks / dups
                sum_ranks = dups = 0

    if axis == 0:
        return ranks.T
    else:
        return ranks

def rank_1d_generic(object in_arr, bint retry=1):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, n, dups = 0
        ndarray[float64_t] ranks
        ndarray sorted_data, values
        ndarray[int64_t] argsorted
        int32_t idx
        object val, nan_value
        float64_t sum_ranks = 0

    values = np.array(in_arr, copy=True)

    if values.dtype != np.object_:
        values = values.astype('O')

    nan_value = Infinity() # always greater than everything
    mask = isnullobj(values)
    np.putmask(values, mask, nan_value)

    n = len(values)
    ranks = np.empty(n, dtype='f8')

    # py2.5/win32 hack, can't pass i8
    try:
        _as = values.argsort()
    except TypeError:
        if not retry:
            raise

        valid_locs = (-mask).nonzero()[0]
        ranks.put(valid_locs, rank_1d_generic(values.take(valid_locs), 0))
        np.putmask(ranks, mask, np.nan)
        return ranks

    sorted_data = values.take(_as)
    argsorted = _as.astype('i8')

    for i in range(n):
        sum_ranks += i + 1
        dups += 1
        val = util.get_value_at(sorted_data, i)
        if val is nan_value:
            ranks[argsorted[i]] = nan
            continue
        if (i == n - 1 or
            are_diff(util.get_value_at(sorted_data, i + 1), val)):
            for j in range(i - dups + 1, i + 1):
                ranks[argsorted[j]] = sum_ranks / dups
            sum_ranks = dups = 0
    return ranks

cdef inline are_diff(object left, object right):
    try:
        return fabs(left - right) > FP_ERR
    except TypeError:
        return left != right

class Infinity(object):

    return_false = lambda self, other: False
    return_true = lambda self, other: True
    __lt__ = return_false
    __le__ = return_false
    __eq__ = return_false
    __ne__ = return_true
    __gt__ = return_true
    __ge__ = return_true
    __cmp__ = return_false

def rank_2d_generic(object in_arr, axis=0):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, z, k, n, infs, dups = 0
        ndarray[float64_t, ndim=2] ranks
        ndarray[object, ndim=2] values
        ndarray[int64_t, ndim=2] argsorted
        int32_t idx
        object val, nan_value
        float64_t sum_ranks = 0

    in_arr = np.asarray(in_arr)

    if axis == 0:
        values = in_arr.T.copy()
    else:
        values = in_arr.copy()

    if values.dtype != np.object_:
        values = values.astype('O')

    nan_value = Infinity() # always greater than everything
    mask = isnullobj2d(values)
    np.putmask(values, mask, nan_value)

    n, k = (<object> values).shape
    ranks = np.empty((n, k), dtype='f8')

    try:
        argsorted = values.argsort(1).astype('i8')
    except TypeError:
        values = in_arr
        for i in range(len(values)):
            ranks[i] = rank_1d_generic(in_arr[i])
        if axis == 0:
            return ranks.T
        else:
            return ranks

    values.sort(axis=1)

    for i in range(n):
        dups = sum_ranks = infs = 0
        for j in range(k):
            val = values[i, j]
            if val is nan_value:
                ranks[i, argsorted[i, j]] = nan
                infs += 1
                continue
            sum_ranks += (j - infs) + 1
            dups += 1
            if j == k - 1 or are_diff(values[i, j + 1], val):
                for z in range(j - dups + 1, j + 1):
                    ranks[i, argsorted[i, z]] = sum_ranks / dups
                sum_ranks = dups = 0

    if axis == 0:
        return ranks.T
    else:
        return ranks
