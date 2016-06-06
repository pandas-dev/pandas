from numpy cimport *
cimport numpy as np
import numpy as np

cimport cython

import_array()

cimport util

from libc.stdlib cimport malloc, free

from numpy cimport NPY_INT8 as NPY_int8
from numpy cimport NPY_INT16 as NPY_int16
from numpy cimport NPY_INT32 as NPY_int32
from numpy cimport NPY_INT64 as NPY_int64
from numpy cimport NPY_FLOAT16 as NPY_float16
from numpy cimport NPY_FLOAT32 as NPY_float32
from numpy cimport NPY_FLOAT64 as NPY_float64

from numpy cimport (int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t,
                    uint32_t, uint64_t, float16_t, float32_t, float64_t)

int8 = np.dtype(np.int8)
int16 = np.dtype(np.int16)
int32 = np.dtype(np.int32)
int64 = np.dtype(np.int64)
float16 = np.dtype(np.float16)
float32 = np.dtype(np.float32)
float64 = np.dtype(np.float64)

cdef np.int8_t MINint8 = np.iinfo(np.int8).min
cdef np.int16_t MINint16 = np.iinfo(np.int16).min
cdef np.int32_t MINint32 = np.iinfo(np.int32).min
cdef np.int64_t MINint64 = np.iinfo(np.int64).min
cdef np.float16_t MINfloat16 = np.NINF
cdef np.float32_t MINfloat32 = np.NINF
cdef np.float64_t MINfloat64 = np.NINF

cdef np.int8_t MAXint8 = np.iinfo(np.int8).max
cdef np.int16_t MAXint16 = np.iinfo(np.int16).max
cdef np.int32_t MAXint32 = np.iinfo(np.int32).max
cdef np.int64_t MAXint64 = np.iinfo(np.int64).max
cdef np.float16_t MAXfloat16 = np.inf
cdef np.float32_t MAXfloat32 = np.inf
cdef np.float64_t MAXfloat64 = np.inf

cdef double NaN = <double> np.NaN
cdef double nan = NaN

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

# this is our util.pxd
from util cimport numeric

cdef extern from "src/headers/math.h":
    double sqrt(double x) nogil
    int signbit(double) nogil

include "skiplist.pyx"

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

def _check_minp(win, minp, N, floor=1):
    if minp > win:
        raise ValueError('min_periods (%d) must be <= window (%d)'
                        % (minp, win))
    elif minp > N:
        minp = N + 1
    elif minp < 0:
        raise ValueError('min_periods must be >= 0')
    return max(minp, floor)

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

#-------------------------------------------------------------------------------
# Rolling sum
@cython.boundscheck(False)
@cython.wraparound(False)
def roll_sum(ndarray[double_t] input, int win, int minp):
    cdef double val, prev, sum_x = 0
    cdef int nobs = 0, i
    cdef int N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    minp = _check_minp(win, minp, N)
    with nogil:
        for i from 0 <= i < minp - 1:
            val = input[i]

            # Not NaN
            if val == val:
                nobs += 1
                sum_x += val

            output[i] = NaN

        for i from minp - 1 <= i < N:
            val = input[i]

            if val == val:
                nobs += 1
                sum_x += val

            if i > win - 1:
                prev = input[i - win]
                if prev == prev:
                    sum_x -= prev
                    nobs -= 1

            if nobs >= minp:
                output[i] = sum_x
            else:
                output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling mean
@cython.boundscheck(False)
@cython.wraparound(False)
def roll_mean(ndarray[double_t] input,
               int win, int minp):
    cdef:
        double val, prev, result, sum_x = 0
        Py_ssize_t nobs = 0, i, neg_ct = 0
        Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)
    minp = _check_minp(win, minp, N)
    with nogil:
        for i from 0 <= i < minp - 1:
            val = input[i]

            # Not NaN
            if val == val:
                nobs += 1
                sum_x += val
                if signbit(val):
                    neg_ct += 1

            output[i] = NaN

        for i from minp - 1 <= i < N:
            val = input[i]

            if val == val:
                nobs += 1
                sum_x += val
                if signbit(val):
                    neg_ct += 1

            if i > win - 1:
                prev = input[i - win]
                if prev == prev:
                    sum_x -= prev
                    nobs -= 1
                    if signbit(prev):
                        neg_ct -= 1

            if nobs >= minp:
                result = sum_x / nobs
                if neg_ct == 0 and result < 0:
                    # all positive
                    output[i] = 0
                elif neg_ct == nobs and result > 0:
                    # all negative
                    output[i] = 0
                else:
                    output[i] = result
            else:
                output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Exponentially weighted moving average

def ewma(ndarray[double_t] input, double_t com, int adjust, int ignore_na, int minp):
    """
    Compute exponentially-weighted moving average using center-of-mass.

    Parameters
    ----------
    input : ndarray (float64 type)
    com : float64
    adjust: int
    ignore_na: int
    minp: int

    Returns
    -------
    y : ndarray
    """

    cdef Py_ssize_t N = len(input)
    cdef ndarray[double_t] output = np.empty(N, dtype=float)
    if N == 0:
        return output

    minp = max(minp, 1)

    cdef double alpha, old_wt_factor, new_wt, weighted_avg, old_wt, cur
    cdef Py_ssize_t i, nobs

    alpha = 1. / (1. + com)
    old_wt_factor = 1. - alpha
    new_wt = 1. if adjust else alpha

    weighted_avg = input[0]
    is_observation = (weighted_avg == weighted_avg)
    nobs = int(is_observation)
    output[0] = weighted_avg if (nobs >= minp) else NaN
    old_wt = 1.

    for i from 1 <= i < N:
        cur = input[i]
        is_observation = (cur == cur)
        nobs += int(is_observation)
        if weighted_avg == weighted_avg:
            if is_observation or (not ignore_na):
                old_wt *= old_wt_factor
                if is_observation:
                    if weighted_avg != cur:  # avoid numerical errors on constant series
                        weighted_avg = ((old_wt * weighted_avg) + (new_wt * cur)) / (old_wt + new_wt)
                    if adjust:
                        old_wt += new_wt
                    else:
                        old_wt = 1.
        elif is_observation:
            weighted_avg = cur

        output[i] = weighted_avg if (nobs >= minp) else NaN

    return output

#-------------------------------------------------------------------------------
# Exponentially weighted moving covariance

def ewmcov(ndarray[double_t] input_x, ndarray[double_t] input_y,
           double_t com, int adjust, int ignore_na, int minp, int bias):
    """
    Compute exponentially-weighted moving variance using center-of-mass.

    Parameters
    ----------
    input_x : ndarray (float64 type)
    input_y : ndarray (float64 type)
    com : float64
    adjust: int
    ignore_na: int
    minp: int
    bias: int

    Returns
    -------
    y : ndarray
    """

    cdef Py_ssize_t N = len(input_x)
    if len(input_y) != N:
        raise ValueError('arrays are of different lengths (%d and %d)' % (N, len(input_y)))
    cdef ndarray[double_t] output = np.empty(N, dtype=float)
    if N == 0:
        return output

    minp = max(minp, 1)

    cdef double alpha, old_wt_factor, new_wt, mean_x, mean_y, cov
    cdef double sum_wt, sum_wt2, old_wt, cur_x, cur_y, old_mean_x, old_mean_y
    cdef Py_ssize_t i, nobs

    alpha = 1. / (1. + com)
    old_wt_factor = 1. - alpha
    new_wt = 1. if adjust else alpha

    mean_x = input_x[0]
    mean_y = input_y[0]
    is_observation = ((mean_x == mean_x) and (mean_y == mean_y))
    nobs = int(is_observation)
    if not is_observation:
        mean_x = NaN
        mean_y = NaN
    output[0] = (0. if bias else NaN) if (nobs >= minp) else NaN
    cov = 0.
    sum_wt = 1.
    sum_wt2 = 1.
    old_wt = 1.

    for i from 1 <= i < N:
        cur_x = input_x[i]
        cur_y = input_y[i]
        is_observation = ((cur_x == cur_x) and (cur_y == cur_y))
        nobs += int(is_observation)
        if mean_x == mean_x:
            if is_observation or (not ignore_na):
                sum_wt *= old_wt_factor
                sum_wt2 *= (old_wt_factor * old_wt_factor)
                old_wt *= old_wt_factor
                if is_observation:
                    old_mean_x = mean_x
                    old_mean_y = mean_y
                    if mean_x != cur_x:  # avoid numerical errors on constant series
                        mean_x = ((old_wt * old_mean_x) + (new_wt * cur_x)) / (old_wt + new_wt)
                    if mean_y != cur_y:  # avoid numerical errors on constant series
                        mean_y = ((old_wt * old_mean_y) + (new_wt * cur_y)) / (old_wt + new_wt)
                    cov = ((old_wt * (cov + ((old_mean_x - mean_x) * (old_mean_y - mean_y)))) +
                           (new_wt * ((cur_x - mean_x) * (cur_y - mean_y)))) / (old_wt + new_wt)
                    sum_wt += new_wt
                    sum_wt2 += (new_wt * new_wt)
                    old_wt += new_wt
                    if not adjust:
                        sum_wt /= old_wt
                        sum_wt2 /= (old_wt * old_wt)
                        old_wt = 1.
        elif is_observation:
            mean_x = cur_x
            mean_y = cur_y

        if nobs >= minp:
            if not bias:
                numerator = sum_wt * sum_wt
                denominator = numerator - sum_wt2
                output[i] = ((numerator / denominator) * cov) if (denominator > 0.) else NaN
            else:
                output[i] = cov
        else:
            output[i] = NaN

    return output

#----------------------------------------------------------------------
# Rolling variance

@cython.boundscheck(False)
@cython.wraparound(False)
def roll_var(ndarray[double_t] input, int win, int minp, int ddof=1):
    """
    Numerically stable implementation using Welford's method.
    """
    cdef double val, prev, mean_x = 0, ssqdm_x = 0, nobs = 0, delta
    cdef Py_ssize_t i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    minp = _check_minp(win, minp, N)

    # Check for windows larger than array, addresses #7297
    win = min(win, N)

    with nogil:
        # Over the first window, observations can only be added, never removed
        for i from 0 <= i < win:
            val = input[i]

            # Not NaN
            if val == val:
                nobs += 1
                delta = (val - mean_x)
                mean_x += delta / nobs
                ssqdm_x += delta * (val - mean_x)

            if (nobs >= minp) and (nobs > ddof):
                #pathological case
                if nobs == 1:
                    val = 0
                else:
                    val = ssqdm_x / (nobs - ddof)
                    if val < 0:
                        val = 0
            else:
                val = NaN

            output[i] = val

        # After the first window, observations can both be added and removed
        for i from win <= i < N:
            val = input[i]
            prev = input[i - win]

            if val == val:
                if prev == prev:
                    # Adding one observation and removing another one
                    delta = val - prev
                    prev -= mean_x
                    mean_x += delta / nobs
                    val -= mean_x
                    ssqdm_x += (val + prev) * delta
                else:
                    # Adding one observation and not removing any
                    nobs += 1
                    delta = (val - mean_x)
                    mean_x += delta / nobs
                    ssqdm_x += delta * (val - mean_x)
            elif prev == prev:
                # Adding no new observation, but removing one
                nobs -= 1
                if nobs:
                    delta = (prev - mean_x)
                    mean_x -= delta  / nobs
                    ssqdm_x -= delta * (prev - mean_x)
                else:
                    mean_x = 0
                    ssqdm_x = 0
            # Variance is unchanged if no observation is added or removed

            if (nobs >= minp) and (nobs > ddof):
                #pathological case
                if nobs == 1:
                    val = 0
                else:
                    val = ssqdm_x / (nobs - ddof)
                    if val < 0:
                        val = 0
            else:
                val = NaN

            output[i] = val

    return output


#-------------------------------------------------------------------------------
# Rolling skewness
@cython.boundscheck(False)
@cython.wraparound(False)
def roll_skew(ndarray[double_t] input, int win, int minp):
    cdef double val, prev
    cdef double x = 0, xx = 0, xxx = 0
    cdef Py_ssize_t nobs = 0, i
    cdef Py_ssize_t N = len(input)

    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    # 3 components of the skewness equation
    cdef double A, B, C, R

    minp = _check_minp(win, minp, N)
    with nogil:
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

            if val == val:
                nobs += 1
                x += val
                xx += val * val
                xxx += val * val * val

            if i > win - 1:
                prev = input[i - win]
                if prev == prev:
                    x -= prev
                    xx -= prev * prev
                    xxx -= prev * prev * prev

                    nobs -= 1
            if nobs >= minp:
                A = x / nobs
                B = xx / nobs - A * A
                C = xxx / nobs - A * A * A - 3 * A * B
                if B <= 0 or nobs < 3:
                    output[i] = NaN
                else:
                    R = sqrt(B)
                    output[i] = ((sqrt(nobs * (nobs - 1.)) * C) /
                                 ((nobs-2) * R * R * R))
            else:
                output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling kurtosis
@cython.boundscheck(False)
@cython.wraparound(False)
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
    with nogil:
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

            if val == val:
                nobs += 1
                x += val
                xx += val * val
                xxx += val * val * val
                xxxx += val * val * val * val

            if i > win - 1:
                prev = input[i - win]
                if prev == prev:
                    x -= prev
                    xx -= prev * prev
                    xxx -= prev * prev * prev
                    xxxx -= prev * prev * prev * prev

                    nobs -= 1

            if nobs >= minp:
                A = x / nobs
                R = A * A
                B = xx / nobs - R
                R = R * A
                C = xxx / nobs - R - 3 * A * B
                R = R * A
                D = xxxx / nobs - R - 6*B*A*A - 4*C*A

                if B == 0 or nobs < 4:
                    output[i] = NaN

                else:
                    K = (nobs * nobs - 1.)*D/(B*B) - 3*((nobs-1.)**2)
                    K = K / ((nobs - 2.)*(nobs-3.))

                    output[i] = K

            else:
                output[i] = NaN

    return output

#-------------------------------------------------------------------------------
# Rolling median, min, max

from skiplist cimport *

@cython.boundscheck(False)
@cython.wraparound(False)
def roll_median_c(ndarray[float64_t] arg, int win, int minp):
    cdef:
        double val, res, prev
        bint err=0
        int ret=0
        skiplist_t *sl
        Py_ssize_t midpoint, nobs = 0, i


    cdef Py_ssize_t N = len(arg)
    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    sl = skiplist_init(win)
    if sl == NULL:
        raise MemoryError("skiplist_init failed")

    minp = _check_minp(win, minp, N)

    with nogil:
        for i from 0 <= i < minp - 1:
            val = arg[i]

            # Not NaN
            if val == val:
                nobs += 1
                err = skiplist_insert(sl, val) != 1
                if err:
                    break
            output[i] = NaN

    with nogil:
        if not err:
            for i from minp - 1 <= i < N:

                val = arg[i]

                if i > win - 1:
                    prev = arg[i - win]

                    if prev == prev:
                        skiplist_remove(sl, prev)
                        nobs -= 1

                if val == val:
                    nobs += 1
                    err = skiplist_insert(sl, val) != 1
                    if err:
                        break

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
    if err:
        raise MemoryError("skiplist_insert failed")
    return output

#----------------------------------------------------------------------

# Moving maximum / minimum code taken from Bottleneck under the terms
# of its Simplified BSD license
# https://github.com/kwgoodman/bottleneck

@cython.boundscheck(False)
@cython.wraparound(False)
def roll_max(ndarray[numeric] a, int window, int minp):
    """
    Moving max of 1d array of any numeric type along axis=0 ignoring NaNs.

    Parameters
    ----------
    a: numpy array
    window: int, size of rolling window
    minp: if number of observations in window
          is below this, output a NaN
    """
    return _roll_min_max(a, window, minp, 1)

@cython.boundscheck(False)
@cython.wraparound(False)
def roll_min(ndarray[numeric] a, int window, int minp):
    """
    Moving max of 1d array of any numeric type along axis=0 ignoring NaNs.

    Parameters
    ----------
    a: numpy array
    window: int, size of rolling window
    minp: if number of observations in window
          is below this, output a NaN
    """
    return _roll_min_max(a, window, minp, 0)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef _roll_min_max(ndarray[numeric] a, int window, int minp, bint is_max):
    "Moving min/max of 1d array of any numeric type along axis=0 ignoring NaNs."
    cdef numeric ai, aold
    cdef Py_ssize_t count
    cdef Py_ssize_t* death
    cdef numeric* ring
    cdef numeric* minvalue
    cdef numeric* end
    cdef numeric* last
    cdef Py_ssize_t i0
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef np.npy_intp *dims = [n0]
    cdef bint should_replace
    cdef np.ndarray[numeric, ndim=1] y = PyArray_EMPTY(1, dims, PyArray_TYPE(a), 0)

    if window < 1:
        raise ValueError('Invalid window size %d'
                         % (window))

    if minp > window:
        raise ValueError('Invalid min_periods size %d greater than window %d'
                        % (minp, window))

    minp = _check_minp(window, minp, n0)
    with nogil:
        ring = <numeric*>malloc(window * sizeof(numeric))
        death = <Py_ssize_t*>malloc(window * sizeof(Py_ssize_t))
        end = ring + window
        last = ring

        minvalue = ring
        ai = a[0]
        if numeric in cython.floating:
            if ai == ai:
                minvalue[0] = ai
            elif is_max:
                minvalue[0] = MINfloat64
            else:
                minvalue[0] = MAXfloat64
        else:
            minvalue[0] = ai
        death[0] = window

        count = 0
        for i0 in range(n0):
            ai = a[i0]
            if numeric in cython.floating:
                if ai == ai:
                    count += 1
                elif is_max:
                    ai = MINfloat64
                else:
                    ai = MAXfloat64
            else:
                count += 1
            if i0 >= window:
                aold = a[i0 - window]
                if aold == aold:
                    count -= 1
            if death[minvalue-ring] == i0:
                minvalue += 1
                if minvalue >= end:
                    minvalue = ring
            should_replace = ai >= minvalue[0] if is_max else ai <= minvalue[0]
            if should_replace:
                minvalue[0] = ai
                death[minvalue-ring] = i0 + window
                last = minvalue
            else:
                should_replace = last[0] <= ai if is_max else last[0] >= ai
                while should_replace:
                    if last == ring:
                        last = end
                    last -= 1
                    should_replace = last[0] <= ai if is_max else last[0] >= ai
                last += 1
                if last == end:
                    last = ring
                last[0] = ai
                death[last - ring] = i0 + window
            if numeric in cython.floating:
                if count >= minp:
                    y[i0] = minvalue[0]
                else:
                    y[i0] = NaN
            else:
                y[i0] = minvalue[0]

        for i0 in range(minp - 1):
            if numeric in cython.floating:
                y[i0] = NaN
            else:
                y[i0] = 0

        free(ring)
        free(death)
    return y

def roll_quantile(ndarray[float64_t, cast=True] input, int win,
                  int minp, double quantile):
    """
    O(N log(window)) implementation using skip list
    """
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

def roll_generic(ndarray[float64_t, cast=True] input,
                 int win, int minp, int offset,
                 object func, object args, object kwargs):
    cdef ndarray[double_t] output, counts, bufarr
    cdef Py_ssize_t i, n
    cdef float64_t *buf
    cdef float64_t *oldbuf

    if not input.flags.c_contiguous:
        input = input.copy('C')

    n = len(input)
    if n == 0:
        return input

    minp = _check_minp(win, minp, n, floor=0)
    output = np.empty(n, dtype=float)
    counts = roll_sum(np.concatenate((np.isfinite(input).astype(float), np.array([0.] * offset))), win, minp)[offset:]

    # truncated windows at the beginning, through first full-length window
    for i from 0 <= i < (int_min(win, n) - offset):
        if counts[i] >= minp:
            output[i] = func(input[0 : (i + offset + 1)], *args, **kwargs)
        else:
            output[i] = NaN

    # remaining full-length windows
    buf = <float64_t*> input.data
    bufarr = np.empty(win, dtype=float)
    oldbuf = <float64_t*> bufarr.data
    for i from (win - offset) <= i < (n - offset):
        buf = buf + 1
        bufarr.data = <char*> buf
        if counts[i] >= minp:
            output[i] = func(bufarr, *args, **kwargs)
        else:
            output[i] = NaN
    bufarr.data = <char*> oldbuf

    # truncated windows at the end
    for i from int_max(n - offset, 0) <= i < n:
        if counts[i] >= minp:
            output[i] = func(input[int_max(i + offset - win + 1, 0) : n], *args, **kwargs)
        else:
            output[i] = NaN

    return output


def roll_window(ndarray[float64_t, ndim=1, cast=True] input,
                ndarray[float64_t, ndim=1, cast=True] weights,
                int minp, bint avg=True):
    """
    Assume len(weights) << len(input)
    """
    cdef:
        ndarray[double_t] output, tot_wgt, counts
        Py_ssize_t in_i, win_i, win_n, win_k, in_n, in_k
        float64_t val_in, val_win, c, w

    in_n = len(input)
    win_n = len(weights)
    output = np.zeros(in_n, dtype=float)
    counts = np.zeros(in_n, dtype=float)
    if avg:
        tot_wgt = np.zeros(in_n, dtype=float)

    minp = _check_minp(len(weights), minp, in_n)

    if avg:
        for win_i from 0 <= win_i < win_n:
            val_win = weights[win_i]
            if val_win != val_win:
                continue

            for in_i from 0 <= in_i < in_n - (win_n - win_i) + 1:
                val_in = input[in_i]
                if val_in == val_in:
                    output[in_i + (win_n - win_i) - 1] += val_in * val_win
                    counts[in_i + (win_n - win_i) - 1] += 1
                    tot_wgt[in_i + (win_n - win_i) - 1] += val_win

        for in_i from 0 <= in_i < in_n:
            c = counts[in_i]
            if c < minp:
                output[in_i] = NaN
            else:
                w = tot_wgt[in_i]
                if w == 0:
                    output[in_i] = NaN
                else:
                    output[in_i] /= tot_wgt[in_i]

    else:
        for win_i from 0 <= win_i < win_n:
            val_win = weights[win_i]
            if val_win != val_win:
                continue

            for in_i from 0 <= in_i < in_n - (win_n - win_i) + 1:
                val_in = input[in_i]

                if val_in == val_in:
                    output[in_i + (win_n - win_i) - 1] += val_in * val_win
                    counts[in_i + (win_n - win_i) - 1] += 1

        for in_i from 0 <= in_i < in_n:
            c = counts[in_i]
            if c < minp:
                output[in_i] = NaN

    return output
