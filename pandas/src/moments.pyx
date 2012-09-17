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

def _check_minp(win, minp, N):
    if minp > win:
        raise ValueError('min_periods (%d) must be <= window (%d)'
                        % (minp, win))
    elif minp > N:
        minp = N + 1
    elif minp == 0:
        minp = 1
    elif minp < 0:
        raise ValueError('min_periods must be >= 0')
    return minp

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

cdef inline kth_smallest_c(float64_t* a, Py_ssize_t k, Py_ssize_t n):
    cdef:
        Py_ssize_t i,j,l,m
        double_t x, t

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


# -------------- Min, Max subsequence

def max_subseq(ndarray[double_t] arr):
    cdef:
        Py_ssize_t i=0,s=0,e=0,T,n
        double m, S

    n = len(arr)

    if len(arr) == 0:
        return (-1,-1,None)

    m = arr[0]
    S = m
    T = 0

    for i in range(1, n):
        # S = max { S + A[i], A[i] )
        if (S > 0):
            S = S + arr[i]
        else:
            S = arr[i]
            T = i
        if S > m:
            s = T
            e = i
            m = S

    return (s, e, m)

def min_subseq(ndarray[double_t] arr):
    cdef:
        Py_ssize_t s, e
        double m

    (s, e, m) = max_subseq(-arr)

    return (s, e, -m)

#-------------------------------------------------------------------------------
# Rolling sum

def roll_sum(ndarray[double_t] input, int win, int minp):
    cdef double val, prev, sum_x = 0
    cdef int nobs = 0, i
    cdef int N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    minp = _check_minp(win, minp, N)

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

    minp = _check_minp(win, minp, N)

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

def ewma(ndarray[double_t] input, double_t com, int adjust):
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

    if N == 0:
        return output

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

    if adjust:
        for i from 0 <= i < N:
            cur = input[i]
            output[i] = output[i] / (1. - adj)

            if cur == cur:
                adj *= oldw

    return output

#----------------------------------------------------------------------
# Pairwise correlation/covariance

@cython.boundscheck(False)
@cython.wraparound(False)
def nancorr(ndarray[float64_t, ndim=2] mat, cov=False):
    cdef:
        Py_ssize_t i, j, xi, yi, N, K
        ndarray[float64_t, ndim=2] result
        ndarray[uint8_t, ndim=2] mask
        int64_t nobs = 0
        float64_t vx, vy, sumx, sumy, sumxx, sumyy, meanx, meany, divisor

    N, K = (<object> mat).shape

    result = np.empty((K, K), dtype=np.float64)
    mask = np.isfinite(mat).view(np.uint8)

    for xi in range(K):
        for yi in range(xi + 1):
            nobs = sumxx = sumyy = sumx = sumy = 0
            for i in range(N):
                if mask[i, xi] and mask[i, yi]:
                    vx = mat[i, xi]
                    vy = mat[i, yi]
                    nobs += 1
                    sumx += vx
                    sumy += vy

            if nobs == 0:
                result[xi, yi] = result[yi, xi] = np.NaN
            else:
                meanx = sumx / nobs
                meany = sumy / nobs

                # now the cov numerator
                sumx = 0

                for i in range(N):
                    if mask[i, xi] and mask[i, yi]:
                        vx = mat[i, xi] - meanx
                        vy = mat[i, yi] - meany

                        sumx += vx * vy
                        sumxx += vx * vx
                        sumyy += vy * vy

                divisor = (nobs - 1.0) if cov else sqrt(sumxx * sumyy)

                if divisor != 0:
                    result[xi, yi] = result[yi, xi] = sumx / divisor
                else:
                    result[xi, yi] = result[yi, xi] = np.NaN

    return result

#----------------------------------------------------------------------
# Rolling variance

def roll_var(ndarray[double_t] input, int win, int minp, int ddof=1):
    cdef double val, prev, sum_x = 0, sum_xx = 0, nobs = 0
    cdef Py_ssize_t i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    minp = _check_minp(win, minp, N)

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
            # pathological case
            if nobs == 1:
                output[i] = 0
                continue

            output[i] = (nobs * sum_xx - sum_x * sum_x) / (nobs * (nobs - ddof))
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

    minp = _check_minp(win, minp, N)

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

    minp = _check_minp(win, minp, N)

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

    minp = _check_minp(win, minp, N)

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

    minp = _check_minp(win, minp, N)

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

#----------------------------------------------------------------------

# Moving maximum / minimum code taken from Bottleneck under the terms
# of its Simplified BSD license
# https://github.com/kwgoodman/bottleneck

cdef struct pairs:
    double value
    int death

from libc cimport stdlib

@cython.boundscheck(False)
@cython.wraparound(False)
def roll_max2(ndarray[float64_t] a, int window, int minp):
    "Moving max of 1d array of dtype=float64 along axis=0 ignoring NaNs."
    cdef np.float64_t ai, aold
    cdef Py_ssize_t count
    cdef pairs* ring
    cdef pairs* minpair
    cdef pairs* end
    cdef pairs* last
    cdef Py_ssize_t i0
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef np.npy_intp *dims = [n0]
    cdef np.ndarray[np.float64_t, ndim=1] y = PyArray_EMPTY(1, dims,
		NPY_float64, 0)

    if window < 1:
        raise ValueError('Invalid window size %d'
                         % (window))

    if minp > window:
        raise ValueError('Invalid min_periods size %d greater than window %d'
                        % (minp, window))

    minp = _check_minp(window, minp, n0)

    window = min(window, n0)

    ring = <pairs*>stdlib.malloc(window * sizeof(pairs))
    end = ring + window
    last = ring

    minpair = ring
    ai = a[0]
    if ai == ai:
        minpair.value = ai
    else:
        minpair.value = MINfloat64
    minpair.death = window

    count = 0
    for i0 in range(n0):
        ai = a[i0]
        if ai == ai:
            count += 1
        else:
            ai = MINfloat64
        if i0 >= window:
            aold = a[i0 - window]
            if aold == aold:
                count -= 1
        if minpair.death == i0:
            minpair += 1
            if minpair >= end:
                minpair = ring
        if ai >= minpair.value:
            minpair.value = ai
            minpair.death = i0 + window
            last = minpair
        else:
            while last.value <= ai:
                if last == ring:
                    last = end
                last -= 1
            last += 1
            if last == end:
                last = ring
            last.value = ai
            last.death = i0 + window
        if count >= minp:
            y[i0] = minpair.value
        else:
            y[i0] = NaN

    for i0 in range(minp - 1):
        y[i0] = NaN

    stdlib.free(ring)
    return y

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

@cython.boundscheck(False)
@cython.wraparound(False)
def roll_min2(np.ndarray[np.float64_t, ndim=1] a, int window, int minp):
    "Moving min of 1d array of dtype=float64 along axis=0 ignoring NaNs."
    cdef np.float64_t ai, aold
    cdef Py_ssize_t count
    cdef pairs* ring
    cdef pairs* minpair
    cdef pairs* end
    cdef pairs* last
    cdef Py_ssize_t i0
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef np.npy_intp *dims = [n0]
    cdef np.ndarray[np.float64_t, ndim=1] y = PyArray_EMPTY(1, dims,
		NPY_float64, 0)

    if window < 1:
        raise ValueError('Invalid window size %d'
                         % (window))

    if minp > window:
        raise ValueError('Invalid min_periods size %d greater than window %d'
                        % (minp, window))

    window = min(window, n0)

    minp = _check_minp(window, minp, n0)

    ring = <pairs*>stdlib.malloc(window * sizeof(pairs))
    end = ring + window
    last = ring

    minpair = ring
    ai = a[0]
    if ai == ai:
        minpair.value = ai
    else:
        minpair.value = MAXfloat64
    minpair.death = window

    count = 0
    for i0 in range(n0):
        ai = a[i0]
        if ai == ai:
            count += 1
        else:
            ai = MAXfloat64
        if i0 >= window:
            aold = a[i0 - window]
            if aold == aold:
                count -= 1
        if minpair.death == i0:
            minpair += 1
            if minpair >= end:
                minpair = ring
        if ai <= minpair.value:
            minpair.value = ai
            minpair.death = i0 + window
            last = minpair
        else:
            while last.value >= ai:
                if last == ring:
                    last = end
                last -= 1
            last += 1
            if last == end:
                last = ring
            last.value = ai
            last.death = i0 + window
        if count >= minp:
            y[i0] = minpair.value
        else:
            y[i0] = NaN

    for i0 in range(minp - 1):
        y[i0] = NaN

    stdlib.free(ring)
    return y

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

    minp = _check_minp(win, minp, N)

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
    if n == 0:
        return input

    minp = _check_minp(win, minp, n)
    output = np.empty(n, dtype=float)
    counts = roll_sum(np.isfinite(input).astype(float), win, minp)

    bufarr = np.empty(win, dtype=float)
    oldbuf = <float64_t*> bufarr.data

    n = len(input)
    for i from 0 <= i < int_min(win, n):
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
