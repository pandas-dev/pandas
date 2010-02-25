# Cython implementations of rolling sum, mean, variance, skewness,
# other statistical moment functions

include "skiplist.pyx"

cdef extern from "wirth.h":
    double kth_smallest(double *a, int n, int k)

#-------------------------------------------------------------------------------
# Rolling sum

cdef void _roll_sum(double_t *input, double_t *output,
                    int win, int minp, int N):
    cdef double val, prev, sum_x = 0
    cdef int nobs = 0, i

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

def _roll_sum_noncontig(ndarray[double_t, ndim=1] input,
                        ndarray[double_t, ndim=1] output,
                        int win, int minp, int N):
    cdef double val, prev, sum_x = 0
    cdef int nobs = 0, i

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

#-------------------------------------------------------------------------------
# Rolling mean

cdef void _roll_mean(double_t *input, double_t *output,
                    int win, int minp, int N):
    cdef double val, prev, sum_x = 0
    cdef int nobs = 0, i

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

def _roll_mean_noncontig(ndarray[double_t, ndim=1] input,
                         ndarray[double_t, ndim=1] output,
                         int win, int minp, int N):
    cdef double val, prev, sum_x = 0
    cdef int nobs = 0, i

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

#-------------------------------------------------------------------------------
# Exponentially weighted moving average

cdef void _ewma(double_t *input, double_t *output,
                      int com, int N):
    cdef double cur, prev, neww, oldw, adj
    cdef int i

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

def _ewma_noncontig(ndarray[double_t, ndim=1] input,
                          ndarray[double_t, ndim=1] output,
                          int com, int N):
    cdef double cur, prev, neww, oldw, adj
    cdef int i

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

#-------------------------------------------------------------------------------
# Rolling variance

cdef void _roll_var(double_t *input, double_t *output,
                    int win, int minp, int N):
    cdef double val, prev, sum_x = 0, sum_xx = 0
    cdef int nobs = 0, i

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


def _roll_var_noncontig(ndarray[double_t, ndim=1] input,
                        ndarray[double_t, ndim=1] output,
                        int win, int minp, int N):
    cdef double val, prev, sum_x = 0, sum_xx = 0
    cdef int nobs = 0, i

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

#-------------------------------------------------------------------------------
# Rolling skewness

cdef void _roll_skew(double_t *input, double_t *output,
                    int win, int minp, int N):
    cdef double val, prev
    cdef double x = 0, xx = 0, xxx = 0
    cdef int nobs = 0, i

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


def _roll_skew_noncontig(ndarray[double_t, ndim=1] input,
                         ndarray[double_t, ndim=1] output,
                         int win, int minp, int N):
    cdef double val, prev
    cdef double x = 0, xx = 0, xxx = 0
    cdef int nobs = 0, i

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

#-------------------------------------------------------------------------------
# Rolling kurtosis


cdef void _roll_kurt(double_t *input, double_t *output,
                     int win, int minp, int N):
    cdef double val, prev
    cdef double x = 0, xx = 0, xxx = 0, xxxx = 0
    cdef int nobs = 0, i

    # 5 components of the kurtosis equation
    cdef double A, B, C, D, R, K

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

def _roll_kurt_noncontig(ndarray[double_t, ndim=1] input,
                         ndarray[double_t, ndim=1] output,
                         int win, int minp, int N):
    cdef double val, prev
    cdef double x = 0, xx = 0, xxx = 0, xxxx = 0
    cdef int nobs = 0, i

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

#-------------------------------------------------------------------------------
# Rolling median

def median(ndarray arr):
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
# Rolling min / max

cdef void _roll_max(double_t *input, double_t *output,
                    int win, int N):
    cdef double val
    cdef int i, j

    for i from 0 <= i < win - 1:
        output[i] = NaN

    for i from win - 1 <= i < N:
        val = input[i - win + 1]

        for j from (i - win + 1) <= j <= i:
            val = max(val, input[j])

        output[j] = val

def _roll_max_noncontig(ndarray[double_t, ndim=1] input,
                        ndarray[double_t, ndim=1] output,
                        int win, int N):
    cdef double val
    cdef int i, j

    for i from 0 <= i < win - 1:
        output[i] = NaN

    for i from win - 1 <= i < N:
        val = input[i - win + 1]

        for j from (i - win + 1) <= j <= i:
            val = max(val, input[j])

        output[j] = val

cdef void _roll_min(double_t *input, double_t *output,
                    int win, int N):
    cdef double val
    cdef int i, j

    for i from 0 <= i < win - 1:
        output[i] = NaN

    for i from win - 1 <= i < N:
        val = input[i - win + 1]

        for j from (i - win + 1) <= j <= i:
            val = min(val, input[j])

        output[j] = val

def _roll_min_noncontig(ndarray[double_t, ndim=1] input,
                        ndarray[double_t, ndim=1] output,
                        int win, int N):
    cdef double val
    cdef int i, j

    for i from 0 <= i < win - 1:
        output[i] = NaN

    for i from win - 1 <= i < N:
        val = input[i - win + 1]

        for j from (i - win + 1) <= j <= i:
            val = min(val, input[j])

        output[j] = val

def roll_median(ndarray[npy_float64, ndim=1] arr, int window, int minp):
    cdef char *mask_data
    cdef ndarray mask
    cdef ndarray[npy_float64, ndim=1] result
    cdef int i, n

    n = len(arr)
    arr = arr.copy()

    mask = <ndarray> np.isfinite(arr)
    mask_data = <char *> mask.data

    result = np.empty(len(arr), dtype=float)

    for i from minp <= i <= n:
        pass

def rolling_median(ndarray[double, ndim=1] arr, int window):
    cdef int i, n, midpoint
    cdef IndexableSkiplist skiplist
    cdef ndarray[double, ndim=1] result

    n = len(arr)
    skiplist = IndexableSkiplist(window)
    result = np.empty(n, dtype=float)
    result[:window] = NaN

    for i from 0 <= i < window:
        skiplist.insert(arr[i])

    midpoint = window / 2

    cdef int flag = window % 2

    for i from window <= i < n:
        skiplist.remove(arr[i - window])
        skiplist.insert(arr[i])

        if flag:
            result[i] = skiplist.get(midpoint)
        else:
            result[i] = (skiplist.get(midpoint) +
                         skiplist.get(midpoint - 1)) / 2

    return result

#-------------------------------------------------------------------------------
# Python interface

# I could use function pointers here and trim down the code duplication

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
    N = len(input)

    if minp is None:
        minp = window

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    cdef ndarray output = np.empty(N, dtype=float)

    if is_contiguous(input):
        _roll_sum(get_double_ptr(input), get_double_ptr(output),
                   window, minp, N)
    else:
        _roll_sum_noncontig(input, output, window, minp, N)

    return output

def rolling_mean(ndarray input, window, minp=None):
    N = len(input)

    if minp is None:
        minp = window

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    cdef ndarray output = np.empty(N, dtype=float)

    if is_contiguous(input):
        _roll_mean(get_double_ptr(input), get_double_ptr(output),
                   window, minp, N)
    else:
        _roll_mean_noncontig(input, output, window, minp, N)

    return output


def ewma(ndarray input, com):
    N = len(input)

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    cdef ndarray output = np.empty(N, dtype=float)

    if is_contiguous(input):
        _ewma(get_double_ptr(input), get_double_ptr(output),
              com, N)
    else:
        _ewma_noncontig(input, output, com, N)

    return output

def rolling_var(ndarray input, window, minp=None):
    N = len(input)

    if minp is None:
        minp = window
    else:
        minp = max(2, minp)

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    cdef ndarray output = np.empty(N, dtype=float)

    if is_contiguous(input):
        _roll_var(get_double_ptr(input), get_double_ptr(output),
                   window, minp, N)
    else:
        _roll_var_noncontig(input, output, window, minp, N)

    return output

def rolling_std(ndarray input, window, minp=None):
    output = rolling_var(input, window, minp=minp)

    return np.sqrt(output)

def rolling_skew(ndarray input, window, minp=None):
    N = len(input)

    if minp is None:
        minp = window
    else:
        minp = max(2, minp)

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    cdef ndarray output = np.empty(N, dtype=float)

    if is_contiguous(input):
        _roll_skew(get_double_ptr(input), get_double_ptr(output),
                   window, minp, N)
    else:
        _roll_skew_noncontig(input, output, window, minp, N)

    return output

def rolling_kurt(ndarray input, window, minp=None):
    N = len(input)

    if minp is None:
        minp = window
    else:
        minp = max(2, minp)

    if not issubclass(input.dtype.type, float):
        input = input.astype(float)

    cdef ndarray output = np.empty(N, dtype=float)

    if is_contiguous(input):
        _roll_kurt(get_double_ptr(input), get_double_ptr(output),
                   window, minp, N)
    else:
        _roll_kurt_noncontig(input, output, window, minp, N)

    return output

