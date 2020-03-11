# cython: boundscheck=False, wraparound=False, cdivision=True

import cython
from cython import Py_ssize_t
from libcpp.deque cimport deque

from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as cnp
from numpy cimport ndarray, int64_t, float64_t, float32_t
cnp.import_array()


cdef extern from "src/headers/cmath" namespace "std":
    bint isnan(float64_t) nogil
    bint notnan(float64_t) nogil
    int signbit(float64_t) nogil
    float64_t sqrt(float64_t x) nogil

from pandas._libs.algos import is_monotonic

from pandas._libs.util cimport numeric

cdef extern from "../src/skiplist.h":
    ctypedef struct node_t:
        node_t **next
        int *width
        double value
        int is_nil
        int levels
        int ref_count

    ctypedef struct skiplist_t:
        node_t *head
        node_t **tmp_chain
        int *tmp_steps
        int size
        int maxlevels

    skiplist_t* skiplist_init(int) nogil
    void skiplist_destroy(skiplist_t*) nogil
    double skiplist_get(skiplist_t*, int, int*) nogil
    int skiplist_insert(skiplist_t*, double) nogil
    int skiplist_remove(skiplist_t*, double) nogil

cdef:
    float32_t MINfloat32 = np.NINF
    float64_t MINfloat64 = np.NINF

    float32_t MAXfloat32 = np.inf
    float64_t MAXfloat64 = np.inf

    float64_t NaN = <float64_t>np.NaN

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

cdef inline bint is_monotonic_start_end_bounds(ndarray[int64_t, ndim=1] start,
                                               ndarray[int64_t, ndim=1] end):
    return is_monotonic(start, False)[0] and is_monotonic(end, False)[0]

# Cython implementations of rolling sum, mean, variance, skewness,
# other statistical moment functions
#
# Misc implementation notes
# -------------------------
#
# - In Cython x * x is faster than x ** 2 for C types, this should be
#   periodically revisited to see if it's still true.
#

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

# ----------------------------------------------------------------------
# Rolling count
# this is only an impl for index not None, IOW, freq aware


def roll_count(ndarray[float64_t] values, ndarray[int64_t] start, ndarray[int64_t] end,
               int64_t minp):
    cdef:
        float64_t val, count_x = 0.0
        int64_t s, e, nobs, N = len(values)
        Py_ssize_t i, j
        ndarray[float64_t] output

    output = np.empty(N, dtype=float)

    with nogil:

        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0:

                # setup
                count_x = 0.0
                for j in range(s, e):
                    val = values[j]
                    if notnan(val):
                        count_x += 1.0

            else:

                # calculate deletes
                for j in range(start[i - 1], s):
                    val = values[j]
                    if notnan(val):
                        count_x -= 1.0

                # calculate adds
                for j in range(end[i - 1], e):
                    val = values[j]
                    if notnan(val):
                        count_x += 1.0

            if count_x >= minp:
                output[i] = count_x
            else:
                output[i] = NaN

    return output


# ----------------------------------------------------------------------
# Rolling sum


cdef inline float64_t calc_sum(int64_t minp, int64_t nobs, float64_t sum_x) nogil:
    cdef:
        float64_t result

    if nobs >= minp:
        result = sum_x
    else:
        result = NaN

    return result


cdef inline void add_sum(float64_t val, int64_t *nobs, float64_t *sum_x) nogil:
    """ add a value from the sum calc """

    # Not NaN
    if notnan(val):
        nobs[0] = nobs[0] + 1
        sum_x[0] = sum_x[0] + val


cdef inline void remove_sum(float64_t val, int64_t *nobs, float64_t *sum_x) nogil:
    """ remove a value from the sum calc """

    if notnan(val):
        nobs[0] = nobs[0] - 1
        sum_x[0] = sum_x[0] - val


def roll_sum_variable(ndarray[float64_t] values, ndarray[int64_t] start,
                      ndarray[int64_t] end, int64_t minp):
    cdef:
        float64_t sum_x = 0
        int64_t s, e
        int64_t nobs = 0, i, j, N = len(values)
        ndarray[float64_t] output
        bint is_monotonic_bounds

    is_monotonic_bounds = is_monotonic_start_end_bounds(start, end)
    output = np.empty(N, dtype=float)

    with nogil:

        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0 or not is_monotonic_bounds:

                # setup

                for j in range(s, e):
                    add_sum(values[j], &nobs, &sum_x)

            else:

                # calculate deletes
                for j in range(start[i - 1], s):
                    remove_sum(values[j], &nobs, &sum_x)

                # calculate adds
                for j in range(end[i - 1], e):
                    add_sum(values[j], &nobs, &sum_x)

            output[i] = calc_sum(minp, nobs, sum_x)

            if not is_monotonic_bounds:
                for j in range(s, e):
                    remove_sum(values[j], &nobs, &sum_x)

    return output


def roll_sum_fixed(ndarray[float64_t] values, ndarray[int64_t] start,
                   ndarray[int64_t] end, int64_t minp, int64_t win):
    cdef:
        float64_t val, prev_x, sum_x = 0
        int64_t range_endpoint
        int64_t nobs = 0, i, N = len(values)
        ndarray[float64_t] output

    output = np.empty(N, dtype=float)

    range_endpoint = int_max(minp, 1) - 1

    with nogil:

        for i in range(0, range_endpoint):
            add_sum(values[i], &nobs, &sum_x)
            output[i] = NaN

        for i in range(range_endpoint, N):
            val = values[i]
            add_sum(val, &nobs, &sum_x)

            if i > win - 1:
                prev_x = values[i - win]
                remove_sum(prev_x, &nobs, &sum_x)

            output[i] = calc_sum(minp, nobs, sum_x)

    return output

# ----------------------------------------------------------------------
# Rolling mean


cdef inline float64_t calc_mean(int64_t minp, Py_ssize_t nobs,
                                Py_ssize_t neg_ct, float64_t sum_x) nogil:
    cdef:
        float64_t result

    if nobs >= minp:
        result = sum_x / <float64_t>nobs
        if neg_ct == 0 and result < 0:
            # all positive
            result = 0
        elif neg_ct == nobs and result > 0:
            # all negative
            result = 0
        else:
            pass
    else:
        result = NaN
    return result


cdef inline void add_mean(float64_t val, Py_ssize_t *nobs, float64_t *sum_x,
                          Py_ssize_t *neg_ct) nogil:
    """ add a value from the mean calc """

    # Not NaN
    if notnan(val):
        nobs[0] = nobs[0] + 1
        sum_x[0] = sum_x[0] + val
        if signbit(val):
            neg_ct[0] = neg_ct[0] + 1


cdef inline void remove_mean(float64_t val, Py_ssize_t *nobs, float64_t *sum_x,
                             Py_ssize_t *neg_ct) nogil:
    """ remove a value from the mean calc """

    if notnan(val):
        nobs[0] = nobs[0] - 1
        sum_x[0] = sum_x[0] - val
        if signbit(val):
            neg_ct[0] = neg_ct[0] - 1


def roll_mean_fixed(ndarray[float64_t] values, ndarray[int64_t] start,
                    ndarray[int64_t] end, int64_t minp, int64_t win):
    cdef:
        float64_t val, prev_x, sum_x = 0
        Py_ssize_t nobs = 0, i, neg_ct = 0, N = len(values)
        ndarray[float64_t] output

    output = np.empty(N, dtype=float)

    with nogil:
        for i in range(minp - 1):
            val = values[i]
            add_mean(val, &nobs, &sum_x, &neg_ct)
            output[i] = NaN

        for i in range(minp - 1, N):
            val = values[i]
            add_mean(val, &nobs, &sum_x, &neg_ct)

            if i > win - 1:
                prev_x = values[i - win]
                remove_mean(prev_x, &nobs, &sum_x, &neg_ct)

            output[i] = calc_mean(minp, nobs, neg_ct, sum_x)

    return output


def roll_mean_variable(ndarray[float64_t] values, ndarray[int64_t] start,
                       ndarray[int64_t] end, int64_t minp):
    cdef:
        float64_t val, sum_x = 0
        int64_t s, e
        Py_ssize_t nobs = 0, i, j, neg_ct = 0, N = len(values)
        ndarray[float64_t] output
        bint is_monotonic_bounds

    is_monotonic_bounds = is_monotonic_start_end_bounds(start, end)
    output = np.empty(N, dtype=float)

    with nogil:

        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0 or not is_monotonic_bounds:

                # setup
                for j in range(s, e):
                    val = values[j]
                    add_mean(val, &nobs, &sum_x, &neg_ct)

            else:

                # calculate deletes
                for j in range(start[i - 1], s):
                    val = values[j]
                    remove_mean(val, &nobs, &sum_x, &neg_ct)

                # calculate adds
                for j in range(end[i - 1], e):
                    val = values[j]
                    add_mean(val, &nobs, &sum_x, &neg_ct)

            output[i] = calc_mean(minp, nobs, neg_ct, sum_x)

            if not is_monotonic_bounds:
                for j in range(s, e):
                    val = values[j]
                    remove_mean(val, &nobs, &sum_x, &neg_ct)
    return output

# ----------------------------------------------------------------------
# Rolling variance


cdef inline float64_t calc_var(int64_t minp, int ddof, float64_t nobs,
                               float64_t ssqdm_x) nogil:
    cdef:
        float64_t result

    # Variance is unchanged if no observation is added or removed
    if (nobs >= minp) and (nobs > ddof):

        # pathological case
        if nobs == 1:
            result = 0
        else:
            result = ssqdm_x / (nobs - <float64_t>ddof)
            if result < 0:
                result = 0
    else:
        result = NaN

    return result


cdef inline void add_var(float64_t val, float64_t *nobs, float64_t *mean_x,
                         float64_t *ssqdm_x) nogil:
    """ add a value from the var calc """
    cdef:
        float64_t delta

    # `isnan` instead of equality as fix for GH-21813, msvc 2017 bug
    if isnan(val):
        return

    nobs[0] = nobs[0] + 1
    # a part of Welford's method for the online variance-calculation
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    delta = val - mean_x[0]
    mean_x[0] = mean_x[0] + delta / nobs[0]
    ssqdm_x[0] = ssqdm_x[0] + ((nobs[0] - 1) * delta ** 2) / nobs[0]


cdef inline void remove_var(float64_t val, float64_t *nobs, float64_t *mean_x,
                            float64_t *ssqdm_x) nogil:
    """ remove a value from the var calc """
    cdef:
        float64_t delta

    if notnan(val):
        nobs[0] = nobs[0] - 1
        if nobs[0]:
            # a part of Welford's method for the online variance-calculation
            # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
            delta = val - mean_x[0]
            mean_x[0] = mean_x[0] - delta / nobs[0]
            ssqdm_x[0] = ssqdm_x[0] - ((nobs[0] + 1) * delta ** 2) / nobs[0]
        else:
            mean_x[0] = 0
            ssqdm_x[0] = 0


def roll_var_fixed(ndarray[float64_t] values, ndarray[int64_t] start,
                   ndarray[int64_t] end, int64_t minp, int64_t win, int ddof=1):
    """
    Numerically stable implementation using Welford's method.
    """
    cdef:
        float64_t mean_x = 0, ssqdm_x = 0, nobs = 0,
        float64_t val, prev, delta, mean_x_old
        int64_t s, e
        Py_ssize_t i, j, N = len(values)
        ndarray[float64_t] output

    output = np.empty(N, dtype=float)

    # Check for windows larger than array, addresses #7297
    win = min(win, N)

    with nogil:

        # Over the first window, observations can only be added, never
        # removed
        for i in range(win):
            add_var(values[i], &nobs, &mean_x, &ssqdm_x)
            output[i] = calc_var(minp, ddof, nobs, ssqdm_x)

        # a part of Welford's method for the online variance-calculation
        # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

        # After the first window, observations can both be added and
        # removed
        for i in range(win, N):
            val = values[i]
            prev = values[i - win]

            if notnan(val):
                if prev == prev:

                    # Adding one observation and removing another one
                    delta = val - prev
                    mean_x_old = mean_x

                    mean_x += delta / nobs
                    ssqdm_x += ((nobs - 1) * val
                                + (nobs + 1) * prev
                                - 2 * nobs * mean_x_old) * delta / nobs

                else:
                    add_var(val, &nobs, &mean_x, &ssqdm_x)
            elif prev == prev:
                remove_var(prev, &nobs, &mean_x, &ssqdm_x)

            output[i] = calc_var(minp, ddof, nobs, ssqdm_x)

    return output


def roll_var_variable(ndarray[float64_t] values, ndarray[int64_t] start,
                      ndarray[int64_t] end, int64_t minp, int ddof=1):
    """
    Numerically stable implementation using Welford's method.
    """
    cdef:
        float64_t mean_x = 0, ssqdm_x = 0, nobs = 0,
        float64_t val, prev, delta, mean_x_old
        int64_t s, e
        Py_ssize_t i, j, N = len(values)
        ndarray[float64_t] output
        bint is_monotonic_bounds

    is_monotonic_bounds = is_monotonic_start_end_bounds(start, end)
    output = np.empty(N, dtype=float)

    with nogil:

        for i in range(0, N):

            s = start[i]
            e = end[i]

            # Over the first window, observations can only be added
            # never removed
            if i == 0 or not is_monotonic_bounds:

                for j in range(s, e):
                    add_var(values[j], &nobs, &mean_x, &ssqdm_x)

            else:

                # After the first window, observations can both be added
                # and removed

                # calculate adds
                for j in range(end[i - 1], e):
                    add_var(values[j], &nobs, &mean_x, &ssqdm_x)

                # calculate deletes
                for j in range(start[i - 1], s):
                    remove_var(values[j], &nobs, &mean_x, &ssqdm_x)

            output[i] = calc_var(minp, ddof, nobs, ssqdm_x)

            if not is_monotonic_bounds:
                for j in range(s, e):
                    remove_var(values[j], &nobs, &mean_x, &ssqdm_x)

    return output

# ----------------------------------------------------------------------
# Rolling skewness


cdef inline float64_t calc_skew(int64_t minp, int64_t nobs,
                                float64_t x, float64_t xx,
                                float64_t xxx) nogil:
    cdef:
        float64_t result, dnobs
        float64_t A, B, C, R

    if nobs >= minp:
        dnobs = <float64_t>nobs
        A = x / dnobs
        B = xx / dnobs - A * A
        C = xxx / dnobs - A * A * A - 3 * A * B

        # #18044: with uniform distribution, floating issue will
        #         cause B != 0. and cause the result is a very
        #         large number.
        #
        #         in core/nanops.py nanskew/nankurt call the function
        #         _zero_out_fperr(m2) to fix floating error.
        #         if the variance is less than 1e-14, it could be
        #         treat as zero, here we follow the original
        #         skew/kurt behaviour to check B <= 1e-14
        if B <= 1e-14 or nobs < 3:
            result = NaN
        else:
            R = sqrt(B)
            result = ((sqrt(dnobs * (dnobs - 1.)) * C) /
                      ((dnobs - 2) * R * R * R))
    else:
        result = NaN

    return result


cdef inline void add_skew(float64_t val, int64_t *nobs,
                          float64_t *x, float64_t *xx,
                          float64_t *xxx) nogil:
    """ add a value from the skew calc """

    # Not NaN
    if notnan(val):
        nobs[0] = nobs[0] + 1

        # seriously don't ask me why this is faster
        x[0] = x[0] + val
        xx[0] = xx[0] + val * val
        xxx[0] = xxx[0] + val * val * val


cdef inline void remove_skew(float64_t val, int64_t *nobs,
                             float64_t *x, float64_t *xx,
                             float64_t *xxx) nogil:
    """ remove a value from the skew calc """

    # Not NaN
    if notnan(val):
        nobs[0] = nobs[0] - 1

        # seriously don't ask me why this is faster
        x[0] = x[0] - val
        xx[0] = xx[0] - val * val
        xxx[0] = xxx[0] - val * val * val


def roll_skew_fixed(ndarray[float64_t] values, ndarray[int64_t] start,
                    ndarray[int64_t] end, int64_t minp, int64_t win):
    cdef:
        float64_t val, prev
        float64_t x = 0, xx = 0, xxx = 0
        int64_t nobs = 0, i, j, N = len(values)
        int64_t s, e
        ndarray[float64_t] output

    output = np.empty(N, dtype=float)

    with nogil:
        for i in range(minp - 1):
            val = values[i]
            add_skew(val, &nobs, &x, &xx, &xxx)
            output[i] = NaN

        for i in range(minp - 1, N):
            val = values[i]
            add_skew(val, &nobs, &x, &xx, &xxx)

            if i > win - 1:
                prev = values[i - win]
                remove_skew(prev, &nobs, &x, &xx, &xxx)

            output[i] = calc_skew(minp, nobs, x, xx, xxx)

    return output


def roll_skew_variable(ndarray[float64_t] values, ndarray[int64_t] start,
                       ndarray[int64_t] end, int64_t minp):
    cdef:
        float64_t val, prev
        float64_t x = 0, xx = 0, xxx = 0
        int64_t nobs = 0, i, j, N = len(values)
        int64_t s, e
        ndarray[float64_t] output
        bint is_monotonic_bounds

    is_monotonic_bounds = is_monotonic_start_end_bounds(start, end)
    output = np.empty(N, dtype=float)

    with nogil:

        for i in range(0, N):

            s = start[i]
            e = end[i]

            # Over the first window, observations can only be added
            # never removed
            if i == 0 or not is_monotonic_bounds:

                for j in range(s, e):
                    val = values[j]
                    add_skew(val, &nobs, &x, &xx, &xxx)

            else:

                # After the first window, observations can both be added
                # and removed

                # calculate adds
                for j in range(end[i - 1], e):
                    val = values[j]
                    add_skew(val, &nobs, &x, &xx, &xxx)

                # calculate deletes
                for j in range(start[i - 1], s):
                    val = values[j]
                    remove_skew(val, &nobs, &x, &xx, &xxx)

            output[i] = calc_skew(minp, nobs, x, xx, xxx)

            if not is_monotonic_bounds:
                for j in range(s, e):
                    val = values[j]
                    remove_skew(val, &nobs, &x, &xx, &xxx)

    return output

# ----------------------------------------------------------------------
# Rolling kurtosis


cdef inline float64_t calc_kurt(int64_t minp, int64_t nobs,
                                float64_t x, float64_t xx,
                                float64_t xxx, float64_t xxxx) nogil:
    cdef:
        float64_t result, dnobs
        float64_t A, B, C, D, R, K

    if nobs >= minp:
        dnobs = <float64_t>nobs
        A = x / dnobs
        R = A * A
        B = xx / dnobs - R
        R = R * A
        C = xxx / dnobs - R - 3 * A * B
        R = R * A
        D = xxxx / dnobs - R - 6 * B * A * A - 4 * C * A

        # #18044: with uniform distribution, floating issue will
        #         cause B != 0. and cause the result is a very
        #         large number.
        #
        #         in core/nanops.py nanskew/nankurt call the function
        #         _zero_out_fperr(m2) to fix floating error.
        #         if the variance is less than 1e-14, it could be
        #         treat as zero, here we follow the original
        #         skew/kurt behaviour to check B <= 1e-14
        if B <= 1e-14 or nobs < 4:
            result = NaN
        else:
            K = (dnobs * dnobs - 1.) * D / (B * B) - 3 * ((dnobs - 1.) ** 2)
            result = K / ((dnobs - 2.) * (dnobs - 3.))
    else:
        result = NaN

    return result


cdef inline void add_kurt(float64_t val, int64_t *nobs,
                          float64_t *x, float64_t *xx,
                          float64_t *xxx, float64_t *xxxx) nogil:
    """ add a value from the kurotic calc """

    # Not NaN
    if notnan(val):
        nobs[0] = nobs[0] + 1

        # seriously don't ask me why this is faster
        x[0] = x[0] + val
        xx[0] = xx[0] + val * val
        xxx[0] = xxx[0] + val * val * val
        xxxx[0] = xxxx[0] + val * val * val * val


cdef inline void remove_kurt(float64_t val, int64_t *nobs,
                             float64_t *x, float64_t *xx,
                             float64_t *xxx, float64_t *xxxx) nogil:
    """ remove a value from the kurotic calc """

    # Not NaN
    if notnan(val):
        nobs[0] = nobs[0] - 1

        # seriously don't ask me why this is faster
        x[0] = x[0] - val
        xx[0] = xx[0] - val * val
        xxx[0] = xxx[0] - val * val * val
        xxxx[0] = xxxx[0] - val * val * val * val


def roll_kurt_fixed(ndarray[float64_t] values, ndarray[int64_t] start,
                    ndarray[int64_t] end, int64_t minp, int64_t win):
    cdef:
        float64_t val, prev
        float64_t x = 0, xx = 0, xxx = 0, xxxx = 0
        int64_t nobs = 0, i, j, N = len(values)
        int64_t s, e
        ndarray[float64_t] output

    output = np.empty(N, dtype=float)

    with nogil:

        for i in range(minp - 1):
            add_kurt(values[i], &nobs, &x, &xx, &xxx, &xxxx)
            output[i] = NaN

        for i in range(minp - 1, N):
            add_kurt(values[i], &nobs, &x, &xx, &xxx, &xxxx)

            if i > win - 1:
                prev = values[i - win]
                remove_kurt(prev, &nobs, &x, &xx, &xxx, &xxxx)

            output[i] = calc_kurt(minp, nobs, x, xx, xxx, xxxx)

    return output


def roll_kurt_variable(ndarray[float64_t] values, ndarray[int64_t] start,
                       ndarray[int64_t] end, int64_t minp):
    cdef:
        float64_t val, prev
        float64_t x = 0, xx = 0, xxx = 0, xxxx = 0
        int64_t nobs = 0, i, j, s, e, N = len(values)
        ndarray[float64_t] output
        bint is_monotonic_bounds

    is_monotonic_bounds = is_monotonic_start_end_bounds(start, end)
    output = np.empty(N, dtype=float)

    with nogil:

        for i in range(0, N):

            s = start[i]
            e = end[i]

            # Over the first window, observations can only be added
            # never removed
            if i == 0 or not is_monotonic_bounds:

                for j in range(s, e):
                    add_kurt(values[j], &nobs, &x, &xx, &xxx, &xxxx)

            else:

                # After the first window, observations can both be added
                # and removed

                # calculate adds
                for j in range(end[i - 1], e):
                    add_kurt(values[j], &nobs, &x, &xx, &xxx, &xxxx)

                # calculate deletes
                for j in range(start[i - 1], s):
                    remove_kurt(values[j], &nobs, &x, &xx, &xxx, &xxxx)

            output[i] = calc_kurt(minp, nobs, x, xx, xxx, xxxx)

            if not is_monotonic_bounds:
                for j in range(s, e):
                    remove_kurt(values[j], &nobs, &x, &xx, &xxx, &xxxx)

    return output


# ----------------------------------------------------------------------
# Rolling median, min, max


def roll_median_c(ndarray[float64_t] values, ndarray[int64_t] start,
                  ndarray[int64_t] end, int64_t minp, int64_t win):
    cdef:
        float64_t val, res, prev
        bint err = 0
        int ret = 0
        skiplist_t *sl
        Py_ssize_t i, j
        int64_t nobs = 0, N = len(values), s, e
        int midpoint
        ndarray[float64_t] output

    # we use the Fixed/Variable Indexer here as the
    # actual skiplist ops outweigh any window computation costs
    output = np.empty(N, dtype=float)

    if win == 0 or (end - start).max() == 0:
        output[:] = NaN
        return output
    win = (end - start).max()
    sl = skiplist_init(<int>win)
    if sl == NULL:
        raise MemoryError("skiplist_init failed")

    with nogil:

        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0:

                # setup
                for j in range(s, e):
                    val = values[j]
                    if notnan(val):
                        nobs += 1
                        err = skiplist_insert(sl, val) != 1
                        if err:
                            break

            else:

                # calculate adds
                for j in range(end[i - 1], e):
                    val = values[j]
                    if notnan(val):
                        nobs += 1
                        err = skiplist_insert(sl, val) != 1
                        if err:
                            break

                # calculate deletes
                for j in range(start[i - 1], s):
                    val = values[j]
                    if notnan(val):
                        skiplist_remove(sl, val)
                        nobs -= 1

            if nobs >= minp:
                midpoint = <int>(nobs / 2)
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


# ----------------------------------------------------------------------

# Moving maximum / minimum code taken from Bottleneck under the terms
# of its Simplified BSD license
# https://github.com/pydata/bottleneck


cdef inline numeric init_mm(numeric ai, Py_ssize_t *nobs, bint is_max) nogil:

    if numeric in cython.floating:
        if ai == ai:
            nobs[0] = nobs[0] + 1
        elif is_max:
            if numeric == cython.float:
                ai = MINfloat32
            else:
                ai = MINfloat64
        else:
            if numeric == cython.float:
                ai = MAXfloat32
            else:
                ai = MAXfloat64

    else:
        nobs[0] = nobs[0] + 1

    return ai


cdef inline void remove_mm(numeric aold, Py_ssize_t *nobs) nogil:
    """ remove a value from the mm calc """
    if numeric in cython.floating and aold == aold:
        nobs[0] = nobs[0] - 1


cdef inline numeric calc_mm(int64_t minp, Py_ssize_t nobs,
                            numeric value) nogil:
    cdef:
        numeric result

    if numeric in cython.floating:
        if nobs >= minp:
            result = value
        else:
            result = NaN
    else:
        result = value

    return result


def roll_max_fixed(ndarray[float64_t] values, ndarray[int64_t] start,
                   ndarray[int64_t] end, int64_t minp, int64_t win):
    """
    Moving max of 1d array of any numeric type along axis=0 ignoring NaNs.

    Parameters
    ----------
    values : np.ndarray[np.float64]
    window : int, size of rolling window
    minp : if number of observations in window
          is below this, output a NaN
    index : ndarray, optional
       index for window computation
    closed : 'right', 'left', 'both', 'neither'
            make the interval closed on the right, left,
            both or neither endpoints
    """
    return _roll_min_max_fixed(values, start, end, minp, win, is_max=1)


def roll_max_variable(ndarray[float64_t] values, ndarray[int64_t] start,
                      ndarray[int64_t] end, int64_t minp):
    """
    Moving max of 1d array of any numeric type along axis=0 ignoring NaNs.

    Parameters
    ----------
    values : np.ndarray[np.float64]
    window : int, size of rolling window
    minp : if number of observations in window
          is below this, output a NaN
    index : ndarray, optional
       index for window computation
    closed : 'right', 'left', 'both', 'neither'
            make the interval closed on the right, left,
            both or neither endpoints
    """
    return _roll_min_max_variable(values, start, end, minp, is_max=1)


def roll_min_fixed(ndarray[float64_t] values, ndarray[int64_t] start,
                   ndarray[int64_t] end, int64_t minp, int64_t win):
    """
    Moving min of 1d array of any numeric type along axis=0 ignoring NaNs.

    Parameters
    ----------
    values : np.ndarray[np.float64]
    window : int, size of rolling window
    minp : if number of observations in window
          is below this, output a NaN
    index : ndarray, optional
       index for window computation
    """
    return _roll_min_max_fixed(values, start, end, minp, win, is_max=0)


def roll_min_variable(ndarray[float64_t] values, ndarray[int64_t] start,
                      ndarray[int64_t] end, int64_t minp):
    """
    Moving min of 1d array of any numeric type along axis=0 ignoring NaNs.

    Parameters
    ----------
    values : np.ndarray[np.float64]
    window : int, size of rolling window
    minp : if number of observations in window
          is below this, output a NaN
    index : ndarray, optional
       index for window computation
    """
    return _roll_min_max_variable(values, start, end, minp, is_max=0)


cdef _roll_min_max_variable(ndarray[numeric] values,
                            ndarray[int64_t] starti,
                            ndarray[int64_t] endi,
                            int64_t minp,
                            bint is_max):
    cdef:
        numeric ai
        int64_t i, close_offset, curr_win_size
        Py_ssize_t nobs = 0, N = len(values)
        deque Q[int64_t]  # min/max always the front
        deque W[int64_t]  # track the whole window for nobs compute
        ndarray[float64_t, ndim=1] output

    output = np.empty(N, dtype=float)
    Q = deque[int64_t]()
    W = deque[int64_t]()

    with nogil:

        # This is using a modified version of the C++ code in this
        # SO post: http://bit.ly/2nOoHlY
        # The original impl didn't deal with variable window sizes
        # So the code was optimized for that

        for i in range(starti[0], endi[0]):
            ai = init_mm(values[i], &nobs, is_max)

            # Discard previous entries if we find new min or max
            if is_max:
                while not Q.empty() and ((ai >= values[Q.back()]) or
                                         values[Q.back()] != values[Q.back()]):
                    Q.pop_back()
            else:
                while not Q.empty() and ((ai <= values[Q.back()]) or
                                         values[Q.back()] != values[Q.back()]):
                    Q.pop_back()
            Q.push_back(i)
            W.push_back(i)

        # if right is open then the first window is empty
        close_offset = 0 if endi[0] > starti[0] else 1
        # first window's size
        curr_win_size = endi[0] - starti[0]

        for i in range(endi[0], endi[N-1]):
            if not Q.empty() and curr_win_size > 0:
                output[i-1+close_offset] = calc_mm(
                    minp, nobs, values[Q.front()])
            else:
                output[i-1+close_offset] = NaN

            ai = init_mm(values[i], &nobs, is_max)

            # Discard previous entries if we find new min or max
            if is_max:
                while not Q.empty() and ((ai >= values[Q.back()]) or
                                         values[Q.back()] != values[Q.back()]):
                    Q.pop_back()
            else:
                while not Q.empty() and ((ai <= values[Q.back()]) or
                                         values[Q.back()] != values[Q.back()]):
                    Q.pop_back()

            # Maintain window/nobs retention
            curr_win_size = endi[i + close_offset] - starti[i + close_offset]
            while not Q.empty() and Q.front() <= i - curr_win_size:
                Q.pop_front()
            while not W.empty() and W.front() <= i - curr_win_size:
                remove_mm(values[W.front()], &nobs)
                W.pop_front()

            Q.push_back(i)
            W.push_back(i)

        if not Q.empty() and curr_win_size > 0:
            output[N-1] = calc_mm(minp, nobs, values[Q.front()])
        else:
            output[N-1] = NaN

    return output


cdef _roll_min_max_fixed(ndarray[numeric] values,
                         ndarray[int64_t] starti,
                         ndarray[int64_t] endi,
                         int64_t minp,
                         int64_t win,
                         bint is_max):
    cdef:
        numeric ai
        bint should_replace
        int64_t i, removed, window_i,
        Py_ssize_t nobs = 0, N = len(values)
        int64_t* death
        numeric* ring
        numeric* minvalue
        numeric* end
        numeric* last
        ndarray[float64_t, ndim=1] output

    output = np.empty(N, dtype=float)
    # setup the rings of death!
    ring = <numeric *>malloc(win * sizeof(numeric))
    death = <int64_t *>malloc(win * sizeof(int64_t))

    end = ring + win
    last = ring
    minvalue = ring
    ai = values[0]
    minvalue[0] = init_mm(values[0], &nobs, is_max)
    death[0] = win
    nobs = 0

    with nogil:

        for i in range(N):
            ai = init_mm(values[i], &nobs, is_max)

            if i >= win:
                remove_mm(values[i - win], &nobs)

            if death[minvalue - ring] == i:
                minvalue = minvalue + 1
                if minvalue >= end:
                    minvalue = ring

            if is_max:
                should_replace = ai >= minvalue[0]
            else:
                should_replace = ai <= minvalue[0]
            if should_replace:

                minvalue[0] = ai
                death[minvalue - ring] = i + win
                last = minvalue

            else:

                if is_max:
                    should_replace = last[0] <= ai
                else:
                    should_replace = last[0] >= ai
                while should_replace:
                    if last == ring:
                        last = end
                    last -= 1
                    if is_max:
                        should_replace = last[0] <= ai
                    else:
                        should_replace = last[0] >= ai

                last += 1
                if last == end:
                    last = ring
                last[0] = ai
                death[last - ring] = i + win

            output[i] = calc_mm(minp, nobs, minvalue[0])

        for i in range(minp - 1):
            if numeric in cython.floating:
                output[i] = NaN
            else:
                output[i] = 0

        free(ring)
        free(death)

    return output


cdef enum InterpolationType:
    LINEAR,
    LOWER,
    HIGHER,
    NEAREST,
    MIDPOINT


interpolation_types = {
    'linear': LINEAR,
    'lower': LOWER,
    'higher': HIGHER,
    'nearest': NEAREST,
    'midpoint': MIDPOINT,
}


def roll_quantile(ndarray[float64_t, cast=True] values, ndarray[int64_t] start,
                  ndarray[int64_t] end, int64_t minp, int64_t win,
                  float64_t quantile, str interpolation):
    """
    O(N log(window)) implementation using skip list
    """
    cdef:
        float64_t val, prev, midpoint, idx_with_fraction
        skiplist_t *skiplist
        int64_t nobs = 0, i, j, s, e, N = len(values)
        Py_ssize_t idx
        ndarray[float64_t] output
        float64_t vlow, vhigh
        InterpolationType interpolation_type
        int ret = 0

    if quantile <= 0.0 or quantile >= 1.0:
        raise ValueError(f"quantile value {quantile} not in [0, 1]")

    try:
        interpolation_type = interpolation_types[interpolation]
    except KeyError:
        raise ValueError(f"Interpolation '{interpolation}' is not supported")

    # we use the Fixed/Variable Indexer here as the
    # actual skiplist ops outweigh any window computation costs
    output = np.empty(N, dtype=float)

    if win == 0 or (end - start).max() == 0:
        output[:] = NaN
        return output
    win = (end - start).max()
    skiplist = skiplist_init(<int>win)
    if skiplist == NULL:
        raise MemoryError("skiplist_init failed")

    with nogil:
        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0:

                # setup
                for j in range(s, e):
                    val = values[j]
                    if notnan(val):
                        nobs += 1
                        skiplist_insert(skiplist, val)

            else:

                # calculate adds
                for j in range(end[i - 1], e):
                    val = values[j]
                    if notnan(val):
                        nobs += 1
                        skiplist_insert(skiplist, val)

                # calculate deletes
                for j in range(start[i - 1], s):
                    val = values[j]
                    if notnan(val):
                        skiplist_remove(skiplist, val)
                        nobs -= 1

            if nobs >= minp:
                if nobs == 1:
                    # Single value in skip list
                    output[i] = skiplist_get(skiplist, 0, &ret)
                else:
                    idx_with_fraction = quantile * (nobs - 1)
                    idx = <int>idx_with_fraction

                    if idx_with_fraction == idx:
                        # no need to interpolate
                        output[i] = skiplist_get(skiplist, idx, &ret)
                        continue

                    if interpolation_type == LINEAR:
                        vlow = skiplist_get(skiplist, idx, &ret)
                        vhigh = skiplist_get(skiplist, idx + 1, &ret)
                        output[i] = ((vlow + (vhigh - vlow) *
                                      (idx_with_fraction - idx)))
                    elif interpolation_type == LOWER:
                        output[i] = skiplist_get(skiplist, idx, &ret)
                    elif interpolation_type == HIGHER:
                        output[i] = skiplist_get(skiplist, idx + 1, &ret)
                    elif interpolation_type == NEAREST:
                        # the same behaviour as round()
                        if idx_with_fraction - idx == 0.5:
                            if idx % 2 == 0:
                                output[i] = skiplist_get(skiplist, idx, &ret)
                            else:
                                output[i] = skiplist_get(
                                    skiplist, idx + 1, &ret)
                        elif idx_with_fraction - idx < 0.5:
                            output[i] = skiplist_get(skiplist, idx, &ret)
                        else:
                            output[i] = skiplist_get(skiplist, idx + 1, &ret)
                    elif interpolation_type == MIDPOINT:
                        vlow = skiplist_get(skiplist, idx, &ret)
                        vhigh = skiplist_get(skiplist, idx + 1, &ret)
                        output[i] = <float64_t>(vlow + vhigh) / 2
            else:
                output[i] = NaN

    skiplist_destroy(skiplist)

    return output


def roll_generic_fixed(object obj,
                       ndarray[int64_t] start, ndarray[int64_t] end,
                       int64_t minp, int64_t win,
                       int offset, object func, bint raw,
                       object args, object kwargs):
    cdef:
        ndarray[float64_t] output, counts, bufarr
        ndarray[float64_t, cast=True] arr
        float64_t *buf
        float64_t *oldbuf
        int64_t nobs = 0, i, j, s, e, N = len(start)

    n = len(obj)
    if n == 0:
        return obj

    arr = np.asarray(obj)

    # ndarray input
    if raw:
        if not arr.flags.c_contiguous:
            arr = arr.copy('C')

    counts = roll_sum_fixed(np.concatenate([np.isfinite(arr).astype(float),
                                            np.array([0.] * offset)]),
                            start, end, minp, win)[offset:]

    output = np.empty(N, dtype=float)

    if not raw:
        # series
        for i in range(N):
            if counts[i] >= minp:
                sl = slice(int_max(i + offset - win + 1, 0),
                           int_min(i + offset + 1, N))
                output[i] = func(obj.iloc[sl], *args, **kwargs)
            else:
                output[i] = NaN

    else:

        # truncated windows at the beginning, through first full-length window
        for i in range((int_min(win, N) - offset)):
            if counts[i] >= minp:
                output[i] = func(arr[0: (i + offset + 1)], *args, **kwargs)
            else:
                output[i] = NaN

        # remaining full-length windows
        buf = <float64_t *>arr.data
        bufarr = np.empty(win, dtype=float)
        oldbuf = <float64_t *>bufarr.data
        for i in range((win - offset), (N - offset)):
            buf = buf + 1
            bufarr.data = <char *>buf
            if counts[i] >= minp:
                output[i] = func(bufarr, *args, **kwargs)
            else:
                output[i] = NaN
        bufarr.data = <char *>oldbuf

        # truncated windows at the end
        for i in range(int_max(N - offset, 0), N):
            if counts[i] >= minp:
                output[i] = func(arr[int_max(i + offset - win + 1, 0): N],
                                 *args,
                                 **kwargs)
            else:
                output[i] = NaN

    return output


def roll_generic_variable(object obj,
                          ndarray[int64_t] start, ndarray[int64_t] end,
                          int64_t minp,
                          int offset, object func, bint raw,
                          object args, object kwargs):
    cdef:
        ndarray[float64_t] output, counts, bufarr
        ndarray[float64_t, cast=True] arr
        float64_t *buf
        float64_t *oldbuf
        int64_t nobs = 0, i, j, s, e, N = len(start)

    n = len(obj)
    if n == 0:
        return obj

    arr = np.asarray(obj)

    # ndarray input
    if raw:
        if not arr.flags.c_contiguous:
            arr = arr.copy('C')

    counts = roll_sum_variable(np.concatenate([np.isfinite(arr).astype(float),
                                               np.array([0.] * offset)]),
                               start, end, minp)[offset:]

    output = np.empty(N, dtype=float)

    if offset != 0:
        raise ValueError("unable to roll_generic with a non-zero offset")

    for i in range(0, N):
        s = start[i]
        e = end[i]

        if counts[i] >= minp:
            if raw:
                output[i] = func(arr[s:e], *args, **kwargs)
            else:
                output[i] = func(obj.iloc[s:e], *args, **kwargs)
        else:
            output[i] = NaN

    return output


# ----------------------------------------------------------------------
# Rolling sum and mean for weighted window


def roll_weighted_sum(float64_t[:] values, float64_t[:] weights, int minp):
    return _roll_weighted_sum_mean(values, weights, minp, avg=0)


def roll_weighted_mean(float64_t[:] values, float64_t[:] weights, int minp):
    return _roll_weighted_sum_mean(values, weights, minp, avg=1)


cdef ndarray[float64_t] _roll_weighted_sum_mean(float64_t[:] values,
                                                float64_t[:] weights,
                                                int minp, bint avg):
    """
    Assume len(weights) << len(values)
    """
    cdef:
        float64_t[:] output, tot_wgt, counts
        Py_ssize_t in_i, win_i, win_n, in_n
        float64_t val_in, val_win, c, w

    in_n = len(values)
    win_n = len(weights)

    output = np.zeros(in_n, dtype=np.float64)
    counts = np.zeros(in_n, dtype=np.float64)
    if avg:
        tot_wgt = np.zeros(in_n, dtype=np.float64)

    if minp > win_n:
        raise ValueError(f"min_periods (minp) must be <= "
                         f"window (win)")
    elif minp > in_n:
        minp = in_n + 1
    elif minp < 0:
        raise ValueError('min_periods must be >= 0')

    minp = max(minp, 1)

    with nogil:
        if avg:
            for win_i in range(win_n):
                val_win = weights[win_i]
                if val_win != val_win:
                    continue

                for in_i in range(in_n - (win_n - win_i) + 1):
                    val_in = values[in_i]
                    if val_in == val_in:
                        output[in_i + (win_n - win_i) - 1] += val_in * val_win
                        counts[in_i + (win_n - win_i) - 1] += 1
                        tot_wgt[in_i + (win_n - win_i) - 1] += val_win

            for in_i in range(in_n):
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
            for win_i in range(win_n):
                val_win = weights[win_i]
                if val_win != val_win:
                    continue

                for in_i in range(in_n - (win_n - win_i) + 1):
                    val_in = values[in_i]

                    if val_in == val_in:
                        output[in_i + (win_n - win_i) - 1] += val_in * val_win
                        counts[in_i + (win_n - win_i) - 1] += 1

            for in_i in range(in_n):
                c = counts[in_i]
                if c < minp:
                    output[in_i] = NaN

    return np.asarray(output)


# ----------------------------------------------------------------------
# Rolling var for weighted window


cdef inline float64_t calc_weighted_var(float64_t t,
                                        float64_t sum_w,
                                        Py_ssize_t win_n,
                                        unsigned int ddof,
                                        float64_t nobs,
                                        int64_t minp) nogil:
    """
    Calculate weighted variance for a window using West's method.

    Paper: https://dl.acm.org/citation.cfm?id=359153

    Parameters
    ----------
    t: float64_t
        sum of weighted squared differences
    sum_w: float64_t
        sum of weights
    win_n: Py_ssize_t
        window size
    ddof: unsigned int
        delta degrees of freedom
    nobs: float64_t
        number of observations
    minp: int64_t
        minimum number of observations

    Returns
    -------
    result : float64_t
        weighted variance of the window
    """

    cdef:
        float64_t result

    # Variance is unchanged if no observation is added or removed
    if (nobs >= minp) and (nobs > ddof):

        # pathological case
        if nobs == 1:
            result = 0
        else:
            result = t * win_n / ((win_n - ddof) * sum_w)
            if result < 0:
                result = 0
    else:
        result = NaN

    return result


cdef inline void add_weighted_var(float64_t val,
                                  float64_t w,
                                  float64_t *t,
                                  float64_t *sum_w,
                                  float64_t *mean,
                                  float64_t *nobs) nogil:
    """
    Update weighted mean, sum of weights and sum of weighted squared
    differences to include value and weight pair in weighted variance
    calculation using West's method.

    Paper: https://dl.acm.org/citation.cfm?id=359153

    Parameters
    ----------
    val: float64_t
        window values
    w: float64_t
        window weights
    t: float64_t
        sum of weighted squared differences
    sum_w: float64_t
        sum of weights
    mean: float64_t
        weighted mean
    nobs: float64_t
        number of observations
    """

    cdef:
        float64_t temp, q, r

    if isnan(val):
        return

    nobs[0] = nobs[0] + 1

    q = val - mean[0]
    temp = sum_w[0] + w
    r = q * w / temp

    mean[0] = mean[0] + r
    t[0] = t[0] + r * sum_w[0] * q
    sum_w[0] = temp


cdef inline void remove_weighted_var(float64_t val,
                                     float64_t w,
                                     float64_t *t,
                                     float64_t *sum_w,
                                     float64_t *mean,
                                     float64_t *nobs) nogil:
    """
    Update weighted mean, sum of weights and sum of weighted squared
    differences to remove value and weight pair from weighted variance
    calculation using West's method.

    Paper: https://dl.acm.org/citation.cfm?id=359153

    Parameters
    ----------
    val: float64_t
        window values
    w: float64_t
        window weights
    t: float64_t
        sum of weighted squared differences
    sum_w: float64_t
        sum of weights
    mean: float64_t
        weighted mean
    nobs: float64_t
        number of observations
    """

    cdef:
        float64_t temp, q, r

    if notnan(val):
        nobs[0] = nobs[0] - 1

        if nobs[0]:
            q = val - mean[0]
            temp = sum_w[0] - w
            r = q * w / temp

            mean[0] = mean[0] - r
            t[0] = t[0] - r * sum_w[0] * q
            sum_w[0] = temp

        else:
            t[0] = 0
            sum_w[0] = 0
            mean[0] = 0


def roll_weighted_var(float64_t[:] values, float64_t[:] weights,
                      int64_t minp, unsigned int ddof):
    """
    Calculates weighted rolling variance using West's online algorithm.

    Paper: https://dl.acm.org/citation.cfm?id=359153

    Parameters
    ----------
    values: float64_t[:]
        values to roll window over
    weights: float64_t[:]
        array of weights whose length is window size
    minp: int64_t
        minimum number of observations to calculate
        variance of a window
    ddof: unsigned int
         the divisor used in variance calculations
         is the window size - ddof

    Returns
    -------
    output: float64_t[:]
        weighted variances of windows
    """

    cdef:
        float64_t t = 0, sum_w = 0, mean = 0, nobs = 0
        float64_t val, pre_val, w, pre_w
        Py_ssize_t i, n, win_n
        float64_t[:] output

    n = len(values)
    win_n = len(weights)
    output = np.empty(n, dtype=float)

    with nogil:

        for i in range(win_n):
            add_weighted_var(values[i], weights[i], &t,
                             &sum_w, &mean, &nobs)

            output[i] = calc_weighted_var(t, sum_w, win_n,
                                          ddof, nobs, minp)

        for i in range(win_n, n):
            val = values[i]
            pre_val = values[i - win_n]

            w = weights[i % win_n]
            pre_w = weights[(i - win_n) % win_n]

            if notnan(val):
                if pre_val == pre_val:
                    remove_weighted_var(pre_val, pre_w, &t,
                                        &sum_w, &mean, &nobs)

                add_weighted_var(val, w, &t, &sum_w, &mean, &nobs)

            elif pre_val == pre_val:
                remove_weighted_var(pre_val, pre_w, &t,
                                    &sum_w, &mean, &nobs)

            output[i] = calc_weighted_var(t, sum_w, win_n,
                                          ddof, nobs, minp)

    return output


# ----------------------------------------------------------------------
# Exponentially weighted moving average


def ewma(float64_t[:] vals, float64_t com, int adjust, bint ignore_na, int minp):
    """
    Compute exponentially-weighted moving average using center-of-mass.

    Parameters
    ----------
    vals : ndarray (float64 type)
    com : float64
    adjust: int
    ignore_na: bool
    minp: int

    Returns
    -------
    ndarray
    """

    cdef:
        Py_ssize_t N = len(vals)
        ndarray[float64_t] output = np.empty(N, dtype=float)
        float64_t alpha, old_wt_factor, new_wt, weighted_avg, old_wt, cur
        Py_ssize_t i, nobs
        bint is_observation

    if N == 0:
        return output

    minp = max(minp, 1)

    alpha = 1. / (1. + com)
    old_wt_factor = 1. - alpha
    new_wt = 1. if adjust else alpha

    weighted_avg = vals[0]
    is_observation = (weighted_avg == weighted_avg)
    nobs = int(is_observation)
    output[0] = weighted_avg if (nobs >= minp) else NaN
    old_wt = 1.

    with nogil:
        for i in range(1, N):
            cur = vals[i]
            is_observation = (cur == cur)
            nobs += is_observation
            if weighted_avg == weighted_avg:

                if is_observation or (not ignore_na):

                    old_wt *= old_wt_factor
                    if is_observation:

                        # avoid numerical errors on constant series
                        if weighted_avg != cur:
                            weighted_avg = ((old_wt * weighted_avg) +
                                            (new_wt * cur)) / (old_wt + new_wt)
                        if adjust:
                            old_wt += new_wt
                        else:
                            old_wt = 1.
            elif is_observation:
                weighted_avg = cur

            output[i] = weighted_avg if (nobs >= minp) else NaN

    return output


# ----------------------------------------------------------------------
# Exponentially weighted moving covariance


def ewmcov(float64_t[:] input_x, float64_t[:] input_y,
           float64_t com, int adjust, bint ignore_na, int minp, int bias):
    """
    Compute exponentially-weighted moving variance using center-of-mass.

    Parameters
    ----------
    input_x : ndarray (float64 type)
    input_y : ndarray (float64 type)
    com : float64
    adjust: int
    ignore_na: bool
    minp: int
    bias: int

    Returns
    -------
    ndarray
    """

    cdef:
        Py_ssize_t N = len(input_x)
        float64_t alpha, old_wt_factor, new_wt, mean_x, mean_y, cov
        float64_t sum_wt, sum_wt2, old_wt, cur_x, cur_y, old_mean_x, old_mean_y
        float64_t numerator, denominator
        Py_ssize_t i, nobs
        ndarray[float64_t] output
        bint is_observation

    if <Py_ssize_t>len(input_y) != N:
        raise ValueError(f"arrays are of different lengths "
                         f"({N} and {len(input_y)})")

    output = np.empty(N, dtype=float)
    if N == 0:
        return output

    minp = max(minp, 1)

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

    with nogil:

        for i in range(1, N):
            cur_x = input_x[i]
            cur_y = input_y[i]
            is_observation = ((cur_x == cur_x) and (cur_y == cur_y))
            nobs += is_observation
            if mean_x == mean_x:
                if is_observation or (not ignore_na):
                    sum_wt *= old_wt_factor
                    sum_wt2 *= (old_wt_factor * old_wt_factor)
                    old_wt *= old_wt_factor
                    if is_observation:
                        old_mean_x = mean_x
                        old_mean_y = mean_y

                        # avoid numerical errors on constant series
                        if mean_x != cur_x:
                            mean_x = ((old_wt * old_mean_x) +
                                      (new_wt * cur_x)) / (old_wt + new_wt)

                        # avoid numerical errors on constant series
                        if mean_y != cur_y:
                            mean_y = ((old_wt * old_mean_y) +
                                      (new_wt * cur_y)) / (old_wt + new_wt)
                        cov = ((old_wt * (cov + ((old_mean_x - mean_x) *
                                                 (old_mean_y - mean_y)))) +
                               (new_wt * ((cur_x - mean_x) *
                                          (cur_y - mean_y)))) / (old_wt + new_wt)
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
                    if (denominator > 0.):
                        output[i] = ((numerator / denominator) * cov)
                    else:
                        output[i] = NaN
                else:
                    output[i] = cov
            else:
                output[i] = NaN

    return output
