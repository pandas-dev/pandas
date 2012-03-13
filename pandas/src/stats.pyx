cdef float64_t FP_ERR = 1e-13

cimport util

cdef:
    int TIEBREAK_AVERAGE = 0
    int TIEBREAK_MIN = 1
    int TIEBREAK_MAX = 2
    int TIEBREAK_FIRST = 3
    int TIEBREAK_FIRST_DESCENDING = 4

tiebreakers = {
    'average' : TIEBREAK_AVERAGE,
    'min' : TIEBREAK_MIN,
    'max' : TIEBREAK_MAX,
    'first' : TIEBREAK_FIRST
}


# ctypedef fused pvalue_t:
#     float64_t
#     int64_t
#     object

# from cython cimport floating, integral

cdef _take_2d_float64(ndarray[float64_t, ndim=2] values,
                      object idx):
    cdef:
        Py_ssize_t i, j, N, K
        ndarray[Py_ssize_t, ndim=2, cast=True] indexer = idx
        ndarray[float64_t, ndim=2] result
        object val

    N, K = (<object> values).shape
    result = np.empty_like(values)
    for i in range(N):
        for j in range(K):
            result[i, j] = values[i, indexer[i, j]]
    return result

cdef _take_2d_int64(ndarray[int64_t, ndim=2] values,
                      object idx):
    cdef:
        Py_ssize_t i, j, N, K
        ndarray[Py_ssize_t, ndim=2, cast=True] indexer = idx
        ndarray[int64_t, ndim=2] result
        object val

    N, K = (<object> values).shape
    result = np.empty_like(values)
    for i in range(N):
        for j in range(K):
            result[i, j] = values[i, indexer[i, j]]
    return result

cdef _take_2d_object(ndarray[object, ndim=2] values,
                     object idx):
    cdef:
        Py_ssize_t i, j, N, K
        ndarray[Py_ssize_t, ndim=2, cast=True] indexer = idx
        ndarray[object, ndim=2] result
        object val

    N, K = (<object> values).shape
    result = values.copy()
    for i in range(N):
        for j in range(K):
            result[i, j] = values[i, indexer[i, j]]
    return result


def rank_1d_float64(object in_arr, ties_method='average', ascending=True):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, n, dups = 0
        ndarray[float64_t] sorted_data, ranks, values
        ndarray[int64_t] argsorted
        float64_t val, nan_value
        float64_t sum_ranks = 0
        int tiebreak = 0
    tiebreak = tiebreakers[ties_method]

    values = np.asarray(in_arr).copy()

    if ascending:
        nan_value = np.inf
    else:
        nan_value = -np.inf
    mask = np.isnan(values)
    np.putmask(values, mask, nan_value)

    n = len(values)
    ranks = np.empty(n, dtype='f8')

    # py2.5/win32 hack, can't pass i8
    if tiebreak == TIEBREAK_FIRST:
        # need to use a stable sort here
        _as = values.argsort(kind='mergesort')
        if not ascending:
            tiebreak = TIEBREAK_FIRST_DESCENDING
    else:
        _as = values.argsort()

    if not ascending:
        _as = _as[::-1]

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
            if tiebreak == TIEBREAK_AVERAGE:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = sum_ranks / dups
            elif tiebreak == TIEBREAK_MIN:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i - dups + 2
            elif tiebreak == TIEBREAK_MAX:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i + 1
            elif tiebreak == TIEBREAK_FIRST:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = j + 1
            elif tiebreak == TIEBREAK_FIRST_DESCENDING:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = 2 * i - j - dups + 2
            sum_ranks = dups = 0
    return ranks


def rank_1d_int64(object in_arr, ties_method='average', ascending=True):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, n, dups = 0
        ndarray[int64_t] sorted_data, values
        ndarray[float64_t] ranks
        ndarray[int64_t] argsorted
        int64_t val
        float64_t sum_ranks = 0
        int tiebreak = 0
    tiebreak = tiebreakers[ties_method]

    values = np.asarray(in_arr)

    n = len(values)
    ranks = np.empty(n, dtype='f8')

    # py2.5/win32 hack, can't pass i8
    if tiebreak == TIEBREAK_FIRST:
        # need to use a stable sort here
        _as = values.argsort(kind='mergesort')
        if not ascending:
            tiebreak = TIEBREAK_FIRST_DESCENDING
    else:
        _as = values.argsort()

    if not ascending:
        _as = _as[::-1]

    sorted_data = values.take(_as)
    argsorted = _as.astype('i8')

    for i in range(n):
        sum_ranks += i + 1
        dups += 1
        val = sorted_data[i]
        if i == n - 1 or fabs(sorted_data[i + 1] - val) > 0:
            if tiebreak == TIEBREAK_AVERAGE:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = sum_ranks / dups
            elif tiebreak == TIEBREAK_MIN:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i - dups + 2
            elif tiebreak == TIEBREAK_MAX:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i + 1
            elif tiebreak == TIEBREAK_FIRST:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = j + 1
            elif tiebreak == TIEBREAK_FIRST_DESCENDING:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = 2 * i - j - dups + 2
            sum_ranks = dups = 0
    return ranks


def rank_2d_float64(object in_arr, axis=0, ties_method='average',
                    ascending=True):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, z, k, n, dups = 0
        ndarray[float64_t, ndim=2] ranks, values
        ndarray[int64_t, ndim=2] argsorted
        float64_t val, nan_value
        float64_t sum_ranks = 0
        int tiebreak = 0
    tiebreak = tiebreakers[ties_method]

    in_arr = np.asarray(in_arr)

    if axis == 0:
        values = in_arr.T.copy()
    else:
        values = in_arr.copy()

    if ascending:
        nan_value = np.inf
    else:
        nan_value = -np.inf

    np.putmask(values, np.isnan(values), nan_value)

    n, k = (<object> values).shape
    ranks = np.empty((n, k), dtype='f8')

    if tiebreak == TIEBREAK_FIRST:
        # need to use a stable sort here
        _as = values.argsort(axis=1, kind='mergesort')
        if not ascending:
            tiebreak = TIEBREAK_FIRST_DESCENDING
    else:
        _as = values.argsort(1)

    if not ascending:
        _as = _as[:, ::-1]

    values = _take_2d_float64(values, _as)
    argsorted = _as.astype('i8')

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
                if tiebreak == TIEBREAK_AVERAGE:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = sum_ranks / dups
                elif tiebreak == TIEBREAK_MIN:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j - dups + 2
                elif tiebreak == TIEBREAK_MAX:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j + 1
                elif tiebreak == TIEBREAK_FIRST:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = z + 1
                elif tiebreak == TIEBREAK_FIRST_DESCENDING:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = 2 * j - z - dups + 2
                sum_ranks = dups = 0

    if axis == 0:
        return ranks.T
    else:
        return ranks


def rank_2d_int64(object in_arr, axis=0, ties_method='average',
                    ascending=True):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, z, k, n, dups = 0
        ndarray[float64_t, ndim=2] ranks
        ndarray[int64_t, ndim=2] argsorted
        ndarray[int64_t, ndim=2, cast=True] values
        int64_t val
        float64_t sum_ranks = 0
        int tiebreak = 0
    tiebreak = tiebreakers[ties_method]

    if axis == 0:
        values = np.asarray(in_arr).T
    else:
        values = np.asarray(in_arr)

    n, k = (<object> values).shape
    ranks = np.empty((n, k), dtype='f8')

    if tiebreak == TIEBREAK_FIRST:
        # need to use a stable sort here
        _as = values.argsort(axis=1, kind='mergesort')
        if not ascending:
            tiebreak = TIEBREAK_FIRST_DESCENDING
    else:
        _as = values.argsort(1)

    if not ascending:
        _as = _as[:, ::-1]

    values = _take_2d_int64(values, _as)
    argsorted = _as.astype('i8')

    for i in range(n):
        dups = sum_ranks = 0
        for j in range(k):
            sum_ranks += j + 1
            dups += 1
            val = values[i, j]
            if j == k - 1 or fabs(values[i, j + 1] - val) > FP_ERR:
                if tiebreak == TIEBREAK_AVERAGE:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = sum_ranks / dups
                elif tiebreak == TIEBREAK_MIN:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j - dups + 2
                elif tiebreak == TIEBREAK_MAX:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j + 1
                elif tiebreak == TIEBREAK_FIRST:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = z + 1
                elif tiebreak == TIEBREAK_FIRST_DESCENDING:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = 2 * j - z - dups + 2
                sum_ranks = dups = 0

    if axis == 0:
        return ranks.T
    else:
        return ranks


def rank_1d_generic(object in_arr, bint retry=1, ties_method='average',
                    ascending=True):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, n, dups = 0
        ndarray[float64_t] ranks
        ndarray sorted_data, values
        ndarray[int64_t] argsorted
        object val, nan_value
        float64_t sum_ranks = 0
        int tiebreak = 0
    tiebreak = tiebreakers[ties_method]

    values = np.array(in_arr, copy=True)

    if values.dtype != np.object_:
        values = values.astype('O')

    if ascending:
        # always greater than everything
        nan_value = Infinity()
    else:
        nan_value = NegInfinity()

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
        ranks.put(valid_locs, rank_1d_generic(values.take(valid_locs), 0,
                                              ties_method=ties_method,
                                              ascending=ascending))
        np.putmask(ranks, mask, np.nan)
        return ranks

    if not ascending:
        _as = _as[::-1]

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
            if tiebreak == TIEBREAK_AVERAGE:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = sum_ranks / dups
            elif tiebreak == TIEBREAK_MIN:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i - dups + 2
            elif tiebreak == TIEBREAK_MAX:
                for j in range(i - dups + 1, i + 1):
                    ranks[argsorted[j]] = i + 1
            elif tiebreak == TIEBREAK_FIRST:
                raise ValueError('first not supported for non-numeric data')
            sum_ranks = dups = 0
    return ranks

cdef inline are_diff(object left, object right):
    try:
        return fabs(left - right) > FP_ERR
    except TypeError:
        return left != right

_return_false = lambda self, other: False
_return_true = lambda self, other: True

class Infinity(object):

    __lt__ = _return_false
    __le__ = _return_false
    __eq__ = _return_false
    __ne__ = _return_true
    __gt__ = _return_true
    __ge__ = _return_true
    __cmp__ = _return_false

class NegInfinity(object):

    __lt__ = _return_true
    __le__ = _return_true
    __eq__ = _return_false
    __ne__ = _return_true
    __gt__ = _return_false
    __ge__ = _return_false
    __cmp__ = _return_true

def rank_2d_generic(object in_arr, axis=0, ties_method='average',
                    ascending=True):
    """
    Fast NaN-friendly version of scipy.stats.rankdata
    """

    cdef:
        Py_ssize_t i, j, z, k, n, infs, dups = 0
        ndarray[float64_t, ndim=2] ranks
        ndarray[object, ndim=2] values
        ndarray[int64_t, ndim=2] argsorted
        object val, nan_value
        float64_t sum_ranks = 0
        int tiebreak = 0
    tiebreak = tiebreakers[ties_method]

    in_arr = np.asarray(in_arr)

    if axis == 0:
        values = in_arr.T.copy()
    else:
        values = in_arr.copy()

    if values.dtype != np.object_:
        values = values.astype('O')

    if ascending:
        # always greater than everything
        nan_value = Infinity()
    else:
        nan_value = NegInfinity()

    mask = isnullobj2d(values)
    np.putmask(values, mask, nan_value)

    n, k = (<object> values).shape
    ranks = np.empty((n, k), dtype='f8')

    try:
        _as = values.argsort(1)
    except TypeError:
        values = in_arr
        for i in range(len(values)):
            ranks[i] = rank_1d_generic(in_arr[i],
                                       ties_method=ties_method,
                                       ascending=ascending)
        if axis == 0:
            return ranks.T
        else:
            return ranks

    if not ascending:
        _as = _as[:, ::-1]

    values = _take_2d_object(values, _as)
    argsorted = _as.astype('i8')

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
                if tiebreak == TIEBREAK_AVERAGE:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = sum_ranks / dups
                elif tiebreak == TIEBREAK_MIN:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j - dups + 2
                elif tiebreak == TIEBREAK_MAX:
                    for z in range(j - dups + 1, j + 1):
                        ranks[i, argsorted[i, z]] = j + 1
                elif tiebreak == TIEBREAK_FIRST:
                    raise ValueError('first not supported for '
                                     'non-numeric data')
                sum_ranks = dups = 0

    if axis == 0:
        return ranks.T
    else:
        return ranks

# def _take_indexer_2d(ndarray[float64_t, ndim=2] values,
#                      ndarray[Py_ssize_t, ndim=2, cast=True] indexer):
#     cdef:
#         Py_ssize_t i, j, N, K
#         ndarray[float64_t, ndim=2] result

#     N, K = (<object> values).shape
#     result = np.empty_like(values)
#     for i in range(N):
#         for j in range(K):
#             result[i, j] = values[i, indexer[i, j]]
#     return result
