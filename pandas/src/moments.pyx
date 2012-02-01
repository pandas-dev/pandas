# Cython implementations of rolling sum, mean, variance, skewness,
# other statistical moment functions
#
# Misc implementation notes
# -------------------------
#
# - In Cython x * x is faster than x ** 2 for C types, this should be
#   periodically revisited to see if it's still true.
#
# -

# original C implementation by N. Devillard.
# This code in public domain.
# Function :   kth_smallest()
# In       :   array of elements, # of elements in the array, rank k
# Out      :   one element
# Job      :   find the kth smallest element in the array

#             Reference:

#               Author: Wirth, Niklaus
#                Title: Algorithms + data structures = programs
#            Publisher: Englewood Cliffs: Prentice-Hall, 1976
# Physical description: 366 p.
#               Series: Prentice-Hall Series in Automatic Computation


def kth_smallest(ndarray[double_t] a, Py_ssize_t k):
    cdef:
        Py_ssize_t i,j,l,m,n
        double_t x, t

    n = len(a)

    l = 0
    m = n-1
    while (l<m):
        x = a[k]
        i = l
        j = m

        while 1:
            while a[i] < x: i += 1
            while x < a[j]: j -= 1
            if i <= j:
                t = a[i]
                a[i] = a[j]
                a[j] = t
                i += 1; j -= 1

            if i > j: break

        if j < k: l = i
        if k < i: m = j
    return a[k]


def median(ndarray arr):
    '''
    A faster median
    '''
    cdef int n = len(arr)

    if len(arr) == 0:
        return np.NaN

    arr = arr.copy()

    if n % 2:
        return kth_smallest(arr, n / 2)
    else:
        return (kth_smallest(arr, n / 2) +
                kth_smallest(arr, n / 2 - 1)) / 2

#-------------------------------------------------------------------------------
# Rolling sum

def roll_sum(ndarray[double_t] input, int win, int minp):
    cdef double val, prev, sum_x = 0
    cdef int nobs = 0, i
    cdef int N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    minp = _check_minp(minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            sum_x += val

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if i > win - 1:
            prev = input[i - win]
            if prev == prev:
                sum_x -= prev
                nobs -= 1

        if val == val:
            nobs += 1
            sum_x += val

        if nobs >= minp:
            output[i] = sum_x
        else:
            output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling mean

def roll_mean(ndarray[double_t] input,
               int win, int minp):
    cdef double val, prev, sum_x = 0
    cdef Py_ssize_t nobs = 0, i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    minp = _check_minp(minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            sum_x += val

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if i > win - 1:
            prev = input[i - win]
            if prev == prev:
                sum_x -= prev
                nobs -= 1

        if val == val:
            nobs += 1
            sum_x += val

        if nobs >= minp:
            output[i] = sum_x / nobs
        else:
            output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Exponentially weighted moving average

def ewma(ndarray[double_t] input, double_t com):
    '''
    Compute exponentially-weighted moving average using center-of-mass.

    Parameters
    ----------
    input : ndarray (float64 type)
    com : float64

    Returns
    -------
    y : ndarray
    '''

    cdef double cur, prev, neww, oldw, adj
    cdef Py_ssize_t i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)


    neww = 1. / (1. + com)
    oldw = 1. - neww
    adj = oldw

    output[0] = neww * input[0]

    for i from 1 <= i < N:
        cur = input[i]
        prev = output[i - 1]

        if cur == cur:
            if prev == prev:
                output[i] = oldw * prev + neww * cur
            else:
                output[i] = neww * cur
        else:
            output[i] = prev

    for i from 0 <= i < N:
        cur = input[i]
        output[i] = output[i] / (1. - adj)

        if cur == cur:
            adj *= oldw

    return output

#-------------------------------------------------------------------------------
# Rolling variance

def _check_minp(minp, N):
    if minp > N:
        minp = N + 1
    elif minp == 0:
        minp = 1
    elif minp < 0:
        raise ValueError('min_periods must be >= 0')
    return minp

def roll_var(ndarray[double_t] input, int win, int minp):
    cdef double val, prev, sum_x = 0, sum_xx = 0, nobs = 0
    cdef Py_ssize_t i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    minp = _check_minp(minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            sum_x += val
            sum_xx += val * val

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if i > win - 1:
            prev = input[i - win]
            if prev == prev:
                sum_x -= prev
                sum_xx -= prev * prev
                nobs -= 1

        if val == val:
            nobs += 1
            sum_x += val
            sum_xx += val * val

        if nobs >= minp:
            output[i] = (nobs * sum_xx - sum_x * sum_x) / (nobs * nobs - nobs)
        else:
            output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling skewness

def roll_skew(ndarray[double_t] input, int win, int minp):
    cdef double val, prev
    cdef double x = 0, xx = 0, xxx = 0
    cdef Py_ssize_t nobs = 0, i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    # 3 components of the skewness equation
    cdef double A, B, C, R

    minp = _check_minp(minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            x += val
            xx += val * val
            xxx += val * val * val

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if i > win - 1:
            prev = input[i - win]
            if prev == prev:
                x -= prev
                xx -= prev * prev
                xxx -= prev * prev * prev

                nobs -= 1

        if val == val:
            nobs += 1
            x += val
            xx += val * val
            xxx += val * val * val

        if nobs >= minp:
            A = x / nobs
            B = xx / nobs - A * A
            C = xxx / nobs - A * A * A - 3 * A * B

            R = sqrt(B)

            output[i] = ((sqrt(nobs * (nobs - 1.)) * C) /
                         ((nobs-2) * R * R * R))
        else:
            output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling kurtosis


def roll_kurt(ndarray[double_t] input,
               int win, int minp):
    cdef double val, prev
    cdef double x = 0, xx = 0, xxx = 0, xxxx = 0
    cdef Py_ssize_t nobs = 0, i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    # 5 components of the kurtosis equation
    cdef double A, B, C, D, R, K

    minp = _check_minp(minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1

            # seriously don't ask me why this is faster
            x += val
            xx += val * val
            xxx += val * val * val
            xxxx += val * val * val * val

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if i > win - 1:
            prev = input[i - win]
            if prev == prev:
                x -= prev
                xx -= prev * prev
                xxx -= prev * prev * prev
                xxxx -= prev * prev * prev * prev

                nobs -= 1

        if val == val:
            nobs += 1
            x += val
            xx += val * val
            xxx += val * val * val
            xxxx += val * val * val * val

        if nobs >= minp:
            A = x / nobs
            R = A * A
            B = xx / nobs - R
            R = R * A
            C = xxx / nobs - R - 3 * A * B
            R = R * A
            D = xxxx / nobs - R - 6*B*A*A - 4*C*A

            K = (nobs * nobs - 1.)*D/(B*B) - 3*((nobs-1.)**2)
            K = K / ((nobs - 2.)*(nobs-3.))

            output[i] = K
        else:
            output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling median, min, max

ctypedef double_t (* skiplist_f)(object sl, int n, int p)

cdef _roll_skiplist_op(ndarray arg, int win, int minp, skiplist_f op):
    cdef ndarray[double_t] input = arg
    cdef double val, prev, midpoint
    cdef IndexableSkiplist skiplist
    cdef Py_ssize_t nobs = 0, i

    cdef Py_ssize_t N = len(input)
    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    skiplist = IndexableSkiplist(win)

    minp = _check_minp(minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            skiplist.insert(val)

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if i > win - 1:
            prev = input[i - win]

            if prev == prev:
                skiplist.remove(prev)
                nobs -= 1

        if val == val:
            nobs += 1
            skiplist.insert(val)

        output[i] = op(skiplist, nobs, minp)

    return output

from skiplist cimport *

def roll_median_c(ndarray[float64_t] arg, int win, int minp):
    cdef double val, res, prev
    cdef:
        int ret=0
        skiplist_t *sl
        Py_ssize_t midpoint, nobs = 0, i


    cdef Py_ssize_t N = len(arg)
    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    sl = skiplist_init(win)

    minp = _check_minp(minp, N)

    for i from 0 <= i < minp - 1:
        val = arg[i]

        # Not NaN
        if val == val:
            nobs += 1
            skiplist_insert(sl, val)

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = arg[i]

        if i > win - 1:
            prev = arg[i - win]

            if prev == prev:
                skiplist_remove(sl, prev)
                nobs -= 1

        if val == val:
            nobs += 1
            skiplist_insert(sl, val)

        if nobs >= minp:
            midpoint = nobs / 2
            if nobs % 2:
                res = skiplist_get(sl, midpoint, &ret)
            else:
                res = (skiplist_get(sl, midpoint, &ret) +
                       skiplist_get(sl, (midpoint - 1), &ret)) / 2
        else:
            res = NaN

        output[i] = res

    skiplist_destroy(sl)

    return output

def roll_median_cython(ndarray input, int win, int minp):
    '''
    O(N log(window)) implementation using skip list
    '''
    return _roll_skiplist_op(input, win, minp, _get_median)

# Unfortunately had to resort to some hackery here, would like for
# Cython to be able to get this right.

cdef double_t _get_median(object sl, int nobs, int minp):
    cdef Py_ssize_t midpoint
    cdef IndexableSkiplist skiplist = <IndexableSkiplist> sl
    if nobs >= minp:
        midpoint = nobs / 2
        if nobs % 2:
            return skiplist.get(midpoint)
        else:
            return (skiplist.get(midpoint) +
                    skiplist.get(midpoint - 1)) / 2
    else:
        return NaN

def roll_max(ndarray input, int win, int minp):
    '''
    O(N log(window)) implementation using skip list
    '''
    return _roll_skiplist_op(input, win, minp, _get_max)

cdef double_t _get_max(object skiplist, int nobs, int minp):
    if nobs >= minp:
        return <IndexableSkiplist> skiplist.get(nobs - 1)
    else:
        return NaN

def roll_min(ndarray input, int win, int minp):
    '''
    O(N log(window)) implementation using skip list
    '''
    return _roll_skiplist_op(input, win, minp, _get_min)

cdef double_t _get_min(object skiplist, int nobs, int minp):
    if nobs >= minp:
        return <IndexableSkiplist> skiplist.get(0)
    else:
        return NaN

def roll_quantile(ndarray[float64_t, cast=True] input, int win,
                  int minp, double quantile):
    '''
    O(N log(window)) implementation using skip list
    '''
    cdef double val, prev, midpoint
    cdef IndexableSkiplist skiplist
    cdef Py_ssize_t nobs = 0, i
    cdef Py_ssize_t N = len(input)
    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    skiplist = IndexableSkiplist(win)

    minp = _check_minp(minp, N)

    for i from 0 <= i < minp - 1:
        val = input[i]

        # Not NaN
        if val == val:
            nobs += 1
            skiplist.insert(val)

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = input[i]

        if i > win - 1:
            prev = input[i - win]

            if prev == prev:
                skiplist.remove(prev)
                nobs -= 1

        if val == val:
            nobs += 1
            skiplist.insert(val)

        if nobs >= minp:
            idx = int((quantile / 1.) * (nobs - 1))
            output[i] = skiplist.get(idx)
        else:
            output[i] = NaN

    return output

def roll_generic(ndarray[float64_t, cast=True] input, int win,
                 int minp, object func):
    cdef ndarray[double_t] output, counts, bufarr
    cdef Py_ssize_t i, n
    cdef float64_t *buf, *oldbuf

    if not input.flags.c_contiguous:
        input = input.copy('C')

    buf = <float64_t*> input.data

    n = len(input)
    minp = _check_minp(minp, n)
    output = np.empty(n, dtype=float)
    counts = roll_sum(np.isfinite(input).astype(float), win, minp)

    bufarr = np.empty(win, dtype=float)
    oldbuf = <float64_t*> bufarr.data

    n = len(input)
    for i from 0 <= i < win:
        if counts[i] >= minp:
            output[i] = func(input[int_max(i - win + 1, 0) : i + 1])
        else:
            output[i] = NaN

    for i from win <= i < n:
        buf = buf + 1
        bufarr.data = <char*> buf
        if counts[i] >= minp:
            output[i] = func(bufarr)
        else:
            output[i] = NaN

    bufarr.data = <char*> oldbuf

    return output
