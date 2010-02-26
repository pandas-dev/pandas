# Cython implementations of rolling sum, mean, variance, skewness,
# other statistical moment functions

cdef extern from "wirth.h":
    double kth_smallest(double *a, int n, int k)

def median(ndarray arr):
    '''
    A faster median
    '''
    cdef double *values
    cdef int n = len(arr)

    if len(arr) == 0:
        return np.NaN

    if not np.PyArray_CHKFLAGS(arr, np.NPY_C_CONTIGUOUS):
        arr = np.array(arr)

    values = <double *> arr.data

    if n % 2:
        return kth_smallest(values, n, n / 2)
    else:
        return (kth_smallest(values, n, n / 2) +
                kth_smallest(values, n, n / 2 - 1)) / 2

#-------------------------------------------------------------------------------
# Rolling sum

def _roll_sum(ndarray[double_t, ndim=1] input,
              int win, int minp):
    cdef double val, prev, sum_x = 0
    cdef int nobs = 0, i
    cdef int N = len(input)

    cdef ndarray[double_t, ndim=1] output = np.empty(N, dtype=float)

    if minp > N:
        minp = N + 1

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

def _roll_mean(ndarray[double_t, ndim=1] input,
               int win, int minp):
    cdef double val, prev, sum_x = 0
    cdef int nobs = 0, i
    cdef int N = len(input)

    cdef ndarray[double_t, ndim=1] output = np.empty(N, dtype=float)

    if minp > N:
        minp = N + 1

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

def _ewma(ndarray[double_t, ndim=1] input,
          int com):
    cdef double cur, prev, neww, oldw, adj
    cdef int i
    cdef int N = len(input)

    cdef ndarray[double_t, ndim=1] output = np.empty(N, dtype=float)


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

def _roll_var(ndarray[double_t, ndim=1] input,
              int win, int minp):
    cdef double val, prev, sum_x = 0, sum_xx = 0
    cdef int nobs = 0, i
    cdef int N = len(input)

    cdef ndarray[double_t, ndim=1] output = np.empty(N, dtype=float)

    if minp > N:
        minp = N + 1

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

def _roll_skew(ndarray[double_t, ndim=1] input,
               int win, int minp):
    cdef double val, prev
    cdef double x = 0, xx = 0, xxx = 0
    cdef int nobs = 0, i
    cdef int N = len(input)

    cdef ndarray[double_t, ndim=1] output = np.empty(N, dtype=float)

    # 3 components of the skewness equation
    cdef double A, B, C, R

    if minp > N:
        minp = N + 1

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


def _roll_kurt(ndarray[double_t, ndim=1] input,
               int win, int minp):
    cdef double val, prev
    cdef double x = 0, xx = 0, xxx = 0, xxxx = 0
    cdef int nobs = 0, i
    cdef int N = len(input)

    cdef ndarray[double_t, ndim=1] output = np.empty(N, dtype=float)

    # 5 components of the kurtosis equation
    cdef double A, B, C, D, R, K

    if minp > N:
        minp = N + 1

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
# Rolling median

def _roll_median(ndarray[double_t, ndim=1] input,
                 int win, int minp):

    cdef double val, prev, midpoint
    cdef IndexableSkiplist skiplist
    cdef int nobs = 0, i
    cdef int N = len(input)

    cdef ndarray[double_t, ndim=1] output = np.empty(N, dtype=float)

    skiplist = IndexableSkiplist(win)

    if minp > N:
        minp = N + 1

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
            midpoint = nobs / 2
            if nobs % 2:
                output[i] = skiplist.get(midpoint)
            else:
                output[i] = (skiplist.get(midpoint) +
                             skiplist.get(midpoint - 1)) / 2
        else:
            output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling min / max

def _roll_max(ndarray[double_t, ndim=1] input,
              int win):
    cdef double val
    cdef int i, j
    cdef int N = len(input)

    cdef ndarray[double_t, ndim=1] output = np.empty(N, dtype=float)

    for i from 0 <= i < win - 1:
        output[i] = NaN

    for i from win - 1 <= i < N:
        val = input[i - win + 1]

        for j from (i - win + 1) <= j <= i:
            val = max(val, input[j])

        output[j] = val

    return output

def _roll_min(ndarray[double_t, ndim=1] input,
              int win):
    cdef double val
    cdef int i, j
    cdef int N = len(input)

    cdef ndarray[double_t, ndim=1] output = np.empty(N, dtype=float)

    for i from 0 <= i < win - 1:
        output[i] = NaN

    for i from win - 1 <= i < N:
        val = input[i - win + 1]

        for j from (i - win + 1) <= j <= i:
            val = min(val, input[j])

        output[j] = val

    return output

#-------------------------------------------------------------------------------
# Python interface

def rolling_sum(ndarray input, window, minp=None):
    '''
    Compute rolling sum of input array

    Parameters
    ----------
    input: ndarray
    window: int
    minp: int

    Returns
    -------
    ndarray with length of input
    '''
    if minp is None:
        minp = window

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    return _roll_sum(input, window, minp)

def rolling_mean(ndarray input, window, minp=None):
    if minp is None:
        minp = window

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    return _roll_mean(input, window, minp)

def rolling_median(ndarray input, window, minp=None):
    if minp is None:
        minp = window

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    return _roll_median(input, window, minp)

def rolling_var(ndarray input, window, minp=None):
    if minp is None:
        minp = window
    else:
        minp = max(2, minp)

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    return _roll_var(input, window, minp)

def rolling_std(ndarray input, window, minp=None):
    output = rolling_var(input, window, minp=minp)
    return np.sqrt(output)

def rolling_skew(ndarray input, window, minp=None):
    if minp is None:
        minp = window
    else:
        minp = max(2, minp)

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    return _roll_skew(input, window, minp)

def rolling_kurt(ndarray input, window, minp=None):
    if minp is None:
        minp = window
    else:
        minp = max(2, minp)

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    return _roll_kurt(input, window, minp)

def ewma(ndarray input, com):
    N = len(input)

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    return _ewma(input, com, N)
