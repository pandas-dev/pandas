# -*- coding: utf-8 -*-
# cython: boundscheck=False, wraparound=False, cdivision=True

cimport cython
from cython cimport Py_ssize_t
from libcpp.deque cimport deque

from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as cnp
from numpy cimport ndarray, double_t, int64_t, float64_t, float32_t
cnp.import_array()


cdef extern from "src/headers/cmath" namespace "std":
    bint isnan(double) nogil
    int signbit(double) nogil
    double sqrt(double x) nogil

cimport util
from util cimport numeric

from skiplist cimport (skiplist_t,
                       skiplist_init, skiplist_destroy,
                       skiplist_get, skiplist_insert, skiplist_remove)

cdef float32_t MINfloat32 = np.NINF
cdef float64_t MINfloat64 = np.NINF

cdef float32_t MAXfloat32 = np.inf
cdef float64_t MAXfloat64 = np.inf

cdef double NaN = <double> np.NaN

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b


# Cython implementations of rolling sum, mean, variance, skewness,
# other statistical moment functions
#
# Misc implementation notes
# -------------------------
#
# - In Cython x * x is faster than x ** 2 for C types, this should be
#   periodically revisited to see if it's still true.
#


def _check_minp(win, minp, N, floor=None):
    """
    Parameters
    ----------
    win: int
    minp: int or None
    N: len of window
    floor: int, optional
        default 1

    Returns
    -------
    minimum period
    """

    if minp is None:
        minp = 1
    if not util.is_integer_object(minp):
        raise ValueError("min_periods must be an integer")
    if minp > win:
        raise ValueError("min_periods (%d) must be <= "
                         "window (%d)" % (minp, win))
    elif minp > N:
        minp = N + 1
    elif minp < 0:
        raise ValueError('min_periods must be >= 0')
    if floor is None:
        floor = 1

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

# ----------------------------------------------------------------------
# The indexer objects for rolling
# These define start/end indexers to compute offsets


cdef class WindowIndexer:

    cdef:
        ndarray start, end
        int64_t N, minp, win
        bint is_variable

    def get_data(self):
        return (self.start, self.end, <int64_t>self.N,
                <int64_t>self.win, <int64_t>self.minp,
                self.is_variable)


cdef class MockFixedWindowIndexer(WindowIndexer):
    """

    We are just checking parameters of the indexer,
    and returning a consistent API with fixed/variable
    indexers.

    Parameters
    ----------
    input: ndarray
        input data array
    win: int64_t
        window size
    minp: int64_t
        min number of obs in a window to consider non-NaN
    index: object
        index of the input
    floor: optional
        unit for flooring
    left_closed: bint
        left endpoint closedness
    right_closed: bint
        right endpoint closedness

    """
    def __init__(self, ndarray input, int64_t win, int64_t minp,
                 bint left_closed, bint right_closed,
                 object index=None, object floor=None):

        assert index is None
        self.is_variable = 0
        self.N = len(input)
        self.minp = _check_minp(win, minp, self.N, floor=floor)
        self.start = np.empty(0, dtype='int64')
        self.end = np.empty(0, dtype='int64')
        self.win = win


cdef class FixedWindowIndexer(WindowIndexer):
    """
    create a fixed length window indexer object
    that has start & end, that point to offsets in
    the index object; these are defined based on the win
    arguments

    Parameters
    ----------
    input: ndarray
        input data array
    win: int64_t
        window size
    minp: int64_t
        min number of obs in a window to consider non-NaN
    index: object
        index of the input
    floor: optional
        unit for flooring the unit
    left_closed: bint
        left endpoint closedness
    right_closed: bint
        right endpoint closedness

    """
    def __init__(self, ndarray input, int64_t win, int64_t minp,
                 bint left_closed, bint right_closed,
                 object index=None, object floor=None):
        cdef ndarray start_s, start_e, end_s, end_e

        assert index is None
        self.is_variable = 0
        self.N = len(input)
        self.minp = _check_minp(win, minp, self.N, floor=floor)

        start_s = np.zeros(win, dtype='int64')
        start_e = np.arange(win, self.N, dtype='int64') - win + 1
        self.start = np.concatenate([start_s, start_e])

        end_s = np.arange(win, dtype='int64') + 1
        end_e = start_e + win
        self.end = np.concatenate([end_s, end_e])
        self.win = win


cdef class VariableWindowIndexer(WindowIndexer):
    """
    create a variable length window indexer object
    that has start & end, that point to offsets in
    the index object; these are defined based on the win
    arguments

    Parameters
    ----------
    input: ndarray
        input data array
    win: int64_t
        window size
    minp: int64_t
        min number of obs in a window to consider non-NaN
    index: ndarray
        index of the input
    left_closed: bint
        left endpoint closedness
        True if the left endpoint is closed, False if open
    right_closed: bint
        right endpoint closedness
        True if the right endpoint is closed, False if open
    floor: optional
        unit for flooring the unit
    """
    def __init__(self, ndarray input, int64_t win, int64_t minp,
                 bint left_closed, bint right_closed, ndarray index,
                 object floor=None):

        self.is_variable = 1
        self.N = len(index)
        self.minp = _check_minp(win, minp, self.N, floor=floor)

        self.start = np.empty(self.N, dtype='int64')
        self.start.fill(-1)

        self.end = np.empty(self.N, dtype='int64')
        self.end.fill(-1)

        self.build(index, win, left_closed, right_closed)

        # max window size
        self.win = (self.end - self.start).max()

    def build(self, ndarray[int64_t] index, int64_t win, bint left_closed,
              bint right_closed):

        cdef:
            ndarray[int64_t] start, end
            int64_t start_bound, end_bound, N
            Py_ssize_t i, j

        start = self.start
        end = self.end
        N = self.N

        start[0] = 0

        # right endpoint is closed
        if right_closed:
            end[0] = 1
        # right endpoint is open
        else:
            end[0] = 0

        with nogil:

            # start is start of slice interval (including)
            # end is end of slice interval (not including)
            for i in range(1, N):
                end_bound = index[i]
                start_bound = index[i] - win

                # left endpoint is closed
                if left_closed:
                    start_bound -= 1

                # advance the start bound until we are
                # within the constraint
                start[i] = i
                for j in range(start[i - 1], i):
                    if index[j] > start_bound:
                        start[i] = j
                        break

                # end bound is previous end
                # or current index
                if index[end[i - 1]] <= end_bound:
                    end[i] = i + 1
                else:
                    end[i] = end[i - 1]

                # right endpoint is open
                if not right_closed:
                    end[i] -= 1


def get_window_indexer(input, win, minp, index, closed,
                       floor=None, use_mock=True):
    """
    return the correct window indexer for the computation

    Parameters
    ----------
    input: 1d ndarray
    win: integer, window size
    minp: integer, minimum periods
    index: 1d ndarray, optional
        index to the input array
    closed: string, default None
        {'right', 'left', 'both', 'neither'}
        window endpoint closedness. Defaults to 'right' in
        VariableWindowIndexer and to 'both' in FixedWindowIndexer
    floor: optional
        unit for flooring the unit
    use_mock: boolean, default True
        if we are a fixed indexer, return a mock indexer
        instead of the FixedWindow Indexer. This is a type
        compat Indexer that allows us to use a standard
        code path with all of the indexers.


    Returns
    -------
    tuple of 1d int64 ndarrays of the offsets & data about the window

    """

    cdef:
        bint left_closed = False
        bint right_closed = False

    assert closed is None or closed in ['right', 'left', 'both', 'neither']

    # if windows is variable, default is 'right', otherwise default is 'both'
    if closed is None:
        closed = 'right' if index is not None else 'both'

    if closed in ['right', 'both']:
        right_closed = True

    if closed in ['left', 'both']:
        left_closed = True

    if index is not None:
        indexer = VariableWindowIndexer(input, win, minp, left_closed,
                                        right_closed, index, floor)
    elif use_mock:
        indexer = MockFixedWindowIndexer(input, win, minp, left_closed,
                                         right_closed, index, floor)
    else:
        indexer = FixedWindowIndexer(input, win, minp, left_closed,
                                     right_closed, index, floor)
    return indexer.get_data()

# ----------------------------------------------------------------------
# Rolling count
# this is only an impl for index not None, IOW, freq aware


def roll_count(ndarray[double_t] input, int64_t win, int64_t minp,
               object index, object closed):
    cdef:
        double val, count_x = 0.0
        int64_t s, e, nobs, N
        Py_ssize_t i, j
        ndarray[int64_t] start, end
        ndarray[double_t] output

    start, end, N, win, minp, _ = get_window_indexer(input, win,
                                                     minp, index, closed)
    output = np.empty(N, dtype=float)

    with nogil:

        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0:

                # setup
                count_x = 0.0
                for j in range(s, e):
                    val = input[j]
                    if val == val:
                        count_x += 1.0

            else:

                # calculate deletes
                for j in range(start[i - 1], s):
                    val = input[j]
                    if val == val:
                        count_x -= 1.0

                # calculate adds
                for j in range(end[i - 1], e):
                    val = input[j]
                    if val == val:
                        count_x += 1.0

            if count_x >= minp:
                output[i] = count_x
            else:
                output[i] = NaN

    return output

# ----------------------------------------------------------------------
# Rolling sum


cdef inline double calc_sum(int64_t minp, int64_t nobs, double sum_x) nogil:
    cdef double result

    if nobs >= minp:
        result = sum_x
    else:
        result = NaN

    return result


cdef inline void add_sum(double val, int64_t *nobs, double *sum_x) nogil:
    """ add a value from the sum calc """

    # Not NaN
    if val == val:
        nobs[0] = nobs[0] + 1
        sum_x[0] = sum_x[0] + val


cdef inline void remove_sum(double val, int64_t *nobs, double *sum_x) nogil:
    """ remove a value from the sum calc """

    if val == val:
        nobs[0] = nobs[0] - 1
        sum_x[0] = sum_x[0] - val


def roll_sum(ndarray[double_t] input, int64_t win, int64_t minp,
             object index, object closed):
    cdef:
        double val, prev_x, sum_x = 0
        int64_t s, e, range_endpoint
        int64_t nobs = 0, i, j, N
        bint is_variable
        ndarray[int64_t] start, end
        ndarray[double_t] output

    start, end, N, win, minp, is_variable = get_window_indexer(input, win,
                                                               minp, index,
                                                               closed,
                                                               floor=0)
    output = np.empty(N, dtype=float)

    # for performance we are going to iterate
    # fixed windows separately, makes the code more complex as we have 2 paths
    # but is faster

    if is_variable:

        # variable window
        with nogil:

            for i in range(0, N):
                s = start[i]
                e = end[i]

                if i == 0:

                    # setup
                    sum_x = 0.0
                    nobs = 0
                    for j in range(s, e):
                        add_sum(input[j], &nobs, &sum_x)

                else:

                    # calculate deletes
                    for j in range(start[i - 1], s):
                        remove_sum(input[j], &nobs, &sum_x)

                    # calculate adds
                    for j in range(end[i - 1], e):
                        add_sum(input[j], &nobs, &sum_x)

                output[i] = calc_sum(minp, nobs, sum_x)

    else:

        # fixed window

        range_endpoint = int_max(minp, 1) - 1

        with nogil:

            for i in range(0, range_endpoint):
                add_sum(input[i], &nobs, &sum_x)
                output[i] = NaN

            for i in range(range_endpoint, N):
                val = input[i]
                add_sum(val, &nobs, &sum_x)

                if i > win - 1:
                    prev_x = input[i - win]
                    remove_sum(prev_x, &nobs, &sum_x)

                output[i] = calc_sum(minp, nobs, sum_x)

    return output

# ----------------------------------------------------------------------
# Rolling mean


cdef inline double calc_mean(int64_t minp, Py_ssize_t nobs,
                             Py_ssize_t neg_ct, double sum_x) nogil:
    cdef double result

    if nobs >= minp:
        result = sum_x / <double>nobs
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


cdef inline void add_mean(double val, Py_ssize_t *nobs, double *sum_x,
                          Py_ssize_t *neg_ct) nogil:
    """ add a value from the mean calc """

    # Not NaN
    if val == val:
        nobs[0] = nobs[0] + 1
        sum_x[0] = sum_x[0] + val
        if signbit(val):
            neg_ct[0] = neg_ct[0] + 1


cdef inline void remove_mean(double val, Py_ssize_t *nobs, double *sum_x,
                             Py_ssize_t *neg_ct) nogil:
    """ remove a value from the mean calc """

    if val == val:
        nobs[0] = nobs[0] - 1
        sum_x[0] = sum_x[0] - val
        if signbit(val):
            neg_ct[0] = neg_ct[0] - 1


def roll_mean(ndarray[double_t] input, int64_t win, int64_t minp,
              object index, object closed):
    cdef:
        double val, prev_x, result, sum_x = 0
        int64_t s, e
        bint is_variable
        Py_ssize_t nobs = 0, i, j, neg_ct = 0, N
        ndarray[int64_t] start, end
        ndarray[double_t] output

    start, end, N, win, minp, is_variable = get_window_indexer(input, win,
                                                               minp, index,
                                                               closed)
    output = np.empty(N, dtype=float)

    # for performance we are going to iterate
    # fixed windows separately, makes the code more complex as we have 2 paths
    # but is faster

    if is_variable:

        with nogil:

            for i in range(0, N):
                s = start[i]
                e = end[i]

                if i == 0:

                    # setup
                    sum_x = 0.0
                    nobs = 0
                    for j in range(s, e):
                        val = input[j]
                        add_mean(val, &nobs, &sum_x, &neg_ct)

                else:

                    # calculate deletes
                    for j in range(start[i - 1], s):
                        val = input[j]
                        remove_mean(val, &nobs, &sum_x, &neg_ct)

                    # calculate adds
                    for j in range(end[i - 1], e):
                        val = input[j]
                        add_mean(val, &nobs, &sum_x, &neg_ct)

                output[i] = calc_mean(minp, nobs, neg_ct, sum_x)

    else:

        with nogil:
            for i in range(minp - 1):
                val = input[i]
                add_mean(val, &nobs, &sum_x, &neg_ct)
                output[i] = NaN

            for i in range(minp - 1, N):
                val = input[i]
                add_mean(val, &nobs, &sum_x, &neg_ct)

                if i > win - 1:
                    prev_x = input[i - win]
                    remove_mean(prev_x, &nobs, &sum_x, &neg_ct)

                output[i] = calc_mean(minp, nobs, neg_ct, sum_x)

    return output

# ----------------------------------------------------------------------
# Rolling variance


cdef inline double calc_var(int64_t minp, int ddof, double nobs,
                            double ssqdm_x) nogil:
    cdef double result

    # Variance is unchanged if no observation is added or removed
    if (nobs >= minp) and (nobs > ddof):

        # pathological case
        if nobs == 1:
            result = 0
        else:
            result = ssqdm_x / (nobs - <double>ddof)
            if result < 0:
                result = 0
    else:
        result = NaN

    return result


cdef inline void add_var(double val, double *nobs, double *mean_x,
                         double *ssqdm_x) nogil:
    """ add a value from the var calc """
    cdef double delta
    # `isnan` instead of equality as fix for GH-21813, msvc 2017 bug
    if isnan(val):
        return

    nobs[0] = nobs[0] + 1
    # a part of Welford's method for the online variance-calculation
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    delta = val - mean_x[0]
    mean_x[0] = mean_x[0] + delta / nobs[0]
    ssqdm_x[0] = ssqdm_x[0] + ((nobs[0] - 1) * delta ** 2) / nobs[0]


cdef inline void remove_var(double val, double *nobs, double *mean_x,
                            double *ssqdm_x) nogil:
    """ remove a value from the var calc """
    cdef double delta

    # Not NaN
    if val == val:
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


def roll_var(ndarray[double_t] input, int64_t win, int64_t minp,
             object index, object closed, int ddof=1):
    """
    Numerically stable implementation using Welford's method.
    """
    cdef:
        double val, prev, mean_x = 0, ssqdm_x = 0, nobs = 0, delta, mean_x_old
        int64_t s, e
        bint is_variable
        Py_ssize_t i, j, N
        ndarray[int64_t] start, end
        ndarray[double_t] output

    start, end, N, win, minp, is_variable = get_window_indexer(input, win,
                                                               minp, index,
                                                               closed)
    output = np.empty(N, dtype=float)

    # Check for windows larger than array, addresses #7297
    win = min(win, N)

    # for performance we are going to iterate
    # fixed windows separately, makes the code more complex as we
    # have 2 paths but is faster

    if is_variable:

        with nogil:

            for i in range(0, N):

                s = start[i]
                e = end[i]

                # Over the first window, observations can only be added
                # never removed
                if i == 0:

                    for j in range(s, e):
                        add_var(input[j], &nobs, &mean_x, &ssqdm_x)

                else:

                    # After the first window, observations can both be added
                    # and removed

                    # calculate adds
                    for j in range(end[i - 1], e):
                        add_var(input[j], &nobs, &mean_x, &ssqdm_x)

                    # calculate deletes
                    for j in range(start[i - 1], s):
                        remove_var(input[j], &nobs, &mean_x, &ssqdm_x)

                output[i] = calc_var(minp, ddof, nobs, ssqdm_x)

    else:

        with nogil:

            # Over the first window, observations can only be added, never
            # removed
            for i in range(win):
                add_var(input[i], &nobs, &mean_x, &ssqdm_x)
                output[i] = calc_var(minp, ddof, nobs, ssqdm_x)

            # a part of Welford's method for the online variance-calculation
            # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

            # After the first window, observations can both be added and
            # removed
            for i in range(win, N):
                val = input[i]
                prev = input[i - win]

                if val == val:
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


# ----------------------------------------------------------------------
# Rolling skewness

cdef inline double calc_skew(int64_t minp, int64_t nobs, double x, double xx,
                             double xxx) nogil:
    cdef double result, dnobs
    cdef double A, B, C, R

    if nobs >= minp:
        dnobs = <double>nobs
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


cdef inline void add_skew(double val, int64_t *nobs, double *x, double *xx,
                          double *xxx) nogil:
    """ add a value from the skew calc """

    # Not NaN
    if val == val:
        nobs[0] = nobs[0] + 1

        # seriously don't ask me why this is faster
        x[0] = x[0] + val
        xx[0] = xx[0] + val * val
        xxx[0] = xxx[0] + val * val * val


cdef inline void remove_skew(double val, int64_t *nobs, double *x, double *xx,
                             double *xxx) nogil:
    """ remove a value from the skew calc """

    # Not NaN
    if val == val:
        nobs[0] = nobs[0] - 1

        # seriously don't ask me why this is faster
        x[0] = x[0] - val
        xx[0] = xx[0] - val * val
        xxx[0] = xxx[0] - val * val * val


def roll_skew(ndarray[double_t] input, int64_t win, int64_t minp,
              object index, object closed):
    cdef:
        double val, prev
        double x = 0, xx = 0, xxx = 0
        int64_t nobs = 0, i, j, N
        int64_t s, e
        bint is_variable
        ndarray[int64_t] start, end
        ndarray[double_t] output

    start, end, N, win, minp, is_variable = get_window_indexer(input, win,
                                                               minp, index,
                                                               closed)
    output = np.empty(N, dtype=float)

    if is_variable:

        with nogil:

            for i in range(0, N):

                s = start[i]
                e = end[i]

                # Over the first window, observations can only be added
                # never removed
                if i == 0:

                    for j in range(s, e):
                        val = input[j]
                        add_skew(val, &nobs, &x, &xx, &xxx)

                else:

                    # After the first window, observations can both be added
                    # and removed

                    # calculate adds
                    for j in range(end[i - 1], e):
                        val = input[j]
                        add_skew(val, &nobs, &x, &xx, &xxx)

                    # calculate deletes
                    for j in range(start[i - 1], s):
                        val = input[j]
                        remove_skew(val, &nobs, &x, &xx, &xxx)

                output[i] = calc_skew(minp, nobs, x, xx, xxx)

    else:

        with nogil:
            for i in range(minp - 1):
                val = input[i]
                add_skew(val, &nobs, &x, &xx, &xxx)
                output[i] = NaN

            for i in range(minp - 1, N):
                val = input[i]
                add_skew(val, &nobs, &x, &xx, &xxx)

                if i > win - 1:
                    prev = input[i - win]
                    remove_skew(prev, &nobs, &x, &xx, &xxx)

                output[i] = calc_skew(minp, nobs, x, xx, xxx)

    return output

# ----------------------------------------------------------------------
# Rolling kurtosis


cdef inline double calc_kurt(int64_t minp, int64_t nobs, double x, double xx,
                             double xxx, double xxxx) nogil:
    cdef double result, dnobs
    cdef double A, B, C, D, R, K

    if nobs >= minp:
        dnobs = <double>nobs
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


cdef inline void add_kurt(double val, int64_t *nobs, double *x, double *xx,
                          double *xxx, double *xxxx) nogil:
    """ add a value from the kurotic calc """

    # Not NaN
    if val == val:
        nobs[0] = nobs[0] + 1

        # seriously don't ask me why this is faster
        x[0] = x[0] + val
        xx[0] = xx[0] + val * val
        xxx[0] = xxx[0] + val * val * val
        xxxx[0] = xxxx[0] + val * val * val * val


cdef inline void remove_kurt(double val, int64_t *nobs, double *x, double *xx,
                             double *xxx, double *xxxx) nogil:
    """ remove a value from the kurotic calc """

    # Not NaN
    if val == val:
        nobs[0] = nobs[0] - 1

        # seriously don't ask me why this is faster
        x[0] = x[0] - val
        xx[0] = xx[0] - val * val
        xxx[0] = xxx[0] - val * val * val
        xxxx[0] = xxxx[0] - val * val * val * val


def roll_kurt(ndarray[double_t] input, int64_t win, int64_t minp,
              object index, object closed):
    cdef:
        double val, prev
        double x = 0, xx = 0, xxx = 0, xxxx = 0
        int64_t nobs = 0, i, j, N
        int64_t s, e
        bint is_variable
        ndarray[int64_t] start, end
        ndarray[double_t] output

    start, end, N, win, minp, is_variable = get_window_indexer(input, win,
                                                               minp, index,
                                                               closed)
    output = np.empty(N, dtype=float)

    if is_variable:

        with nogil:

            for i in range(0, N):

                s = start[i]
                e = end[i]

                # Over the first window, observations can only be added
                # never removed
                if i == 0:

                    for j in range(s, e):
                        add_kurt(input[j], &nobs, &x, &xx, &xxx, &xxxx)

                else:

                    # After the first window, observations can both be added
                    # and removed

                    # calculate adds
                    for j in range(end[i - 1], e):
                        add_kurt(input[j], &nobs, &x, &xx, &xxx, &xxxx)

                    # calculate deletes
                    for j in range(start[i - 1], s):
                        remove_kurt(input[j], &nobs, &x, &xx, &xxx, &xxxx)

                output[i] = calc_kurt(minp, nobs, x, xx, xxx, xxxx)

    else:

        with nogil:

            for i in range(minp - 1):
                add_kurt(input[i], &nobs, &x, &xx, &xxx, &xxxx)
                output[i] = NaN

            for i in range(minp - 1, N):
                add_kurt(input[i], &nobs, &x, &xx, &xxx, &xxxx)

                if i > win - 1:
                    prev = input[i - win]
                    remove_kurt(prev, &nobs, &x, &xx, &xxx, &xxxx)

                output[i] = calc_kurt(minp, nobs, x, xx, xxx, xxxx)

    return output

# ----------------------------------------------------------------------
# Rolling median, min, max


def roll_median_c(ndarray[float64_t] input, int64_t win, int64_t minp,
                  object index, object closed):
    cdef:
        double val, res, prev
        bint err = 0, is_variable
        int ret = 0
        skiplist_t *sl
        Py_ssize_t i, j
        int64_t nobs = 0, N, s, e
        int midpoint
        ndarray[int64_t] start, end
        ndarray[double_t] output

    # we use the Fixed/Variable Indexer here as the
    # actual skiplist ops outweigh any window computation costs
    start, end, N, win, minp, is_variable = get_window_indexer(
        input, win,
        minp, index, closed,
        use_mock=False)
    output = np.empty(N, dtype=float)

    sl = skiplist_init(<int>win)
    if sl == NULL:
        raise MemoryError("skiplist_init failed")

    with nogil:

        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0:

                # setup
                val = input[i]
                if val == val:
                    nobs += 1
                    err = skiplist_insert(sl, val) != 1
                    if err:
                        break

            else:

                # calculate deletes
                for j in range(start[i - 1], s):
                    val = input[j]
                    if val == val:
                        skiplist_remove(sl, val)
                        nobs -= 1

                # calculate adds
                for j in range(end[i - 1], e):
                    val = input[j]
                    if val == val:
                        nobs += 1
                        err = skiplist_insert(sl, val) != 1
                        if err:
                            break

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
# https://github.com/kwgoodman/bottleneck


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
    cdef numeric result

    if numeric in cython.floating:
        if nobs >= minp:
            result = value
        else:
            result = NaN
    else:
        result = value

    return result


def roll_max(ndarray[numeric] input, int64_t win, int64_t minp,
             object index, object closed):
    """
    Moving max of 1d array of any numeric type along axis=0 ignoring NaNs.

    Parameters
    ----------
    input: numpy array
    window: int, size of rolling window
    minp: if number of observations in window
          is below this, output a NaN
    index: ndarray, optional
       index for window computation
    closed: 'right', 'left', 'both', 'neither'
            make the interval closed on the right, left,
            both or neither endpoints
    """
    return _roll_min_max(input, win, minp, index, closed=closed, is_max=1)


def roll_min(ndarray[numeric] input, int64_t win, int64_t minp,
             object index, object closed):
    """
    Moving max of 1d array of any numeric type along axis=0 ignoring NaNs.

    Parameters
    ----------
    input: numpy array
    window: int, size of rolling window
    minp: if number of observations in window
          is below this, output a NaN
    index: ndarray, optional
       index for window computation
    """
    return _roll_min_max(input, win, minp, index, is_max=0, closed=closed)


cdef _roll_min_max(ndarray[numeric] input, int64_t win, int64_t minp,
                   object index, object closed, bint is_max):
    """
    Moving min/max of 1d array of any numeric type along axis=0
    ignoring NaNs.
    """
    cdef:
        ndarray[int64_t] starti, endi
        int64_t N
        bint is_variable

    starti, endi, N, win, minp, is_variable = get_window_indexer(
        input, win,
        minp, index, closed)

    if is_variable:
        return _roll_min_max_variable(input, starti, endi, N, win, minp,
                                      is_max)
    else:
        return _roll_min_max_fixed(input, starti, endi, N, win, minp, is_max)


cdef _roll_min_max_variable(ndarray[numeric] input,
                            ndarray[int64_t] starti,
                            ndarray[int64_t] endi,
                            int64_t N,
                            int64_t win,
                            int64_t minp,
                            bint is_max):
    cdef:
        numeric ai
        int64_t i, close_offset, curr_win_size
        Py_ssize_t nobs = 0
        deque Q[int64_t]  # min/max always the front
        deque W[int64_t]  # track the whole window for nobs compute
        ndarray[double_t, ndim=1] output

    output = np.empty(N, dtype=float)
    Q = deque[int64_t]()
    W = deque[int64_t]()

    with nogil:

        # This is using a modified version of the C++ code in this
        # SO post: http://bit.ly/2nOoHlY
        # The original impl didn't deal with variable window sizes
        # So the code was optimized for that

        for i from starti[0] <= i < endi[0]:
            ai = init_mm(input[i], &nobs, is_max)

            # Discard previous entries if we find new min or max
            if is_max:
                while not Q.empty() and ((ai >= input[Q.back()]) or
                                         (input[Q.back()] != input[Q.back()])):
                    Q.pop_back()
            else:
                while not Q.empty() and ((ai <= input[Q.back()]) or
                                         (input[Q.back()] != input[Q.back()])):
                    Q.pop_back()
            Q.push_back(i)
            W.push_back(i)

        # if right is open then the first window is empty
        close_offset = 0 if endi[0] > starti[0] else 1

        for i in range(endi[0], endi[N-1]):
            if not Q.empty():
                output[i-1+close_offset] = calc_mm(
                    minp, nobs, input[Q.front()])
            else:
                output[i-1+close_offset] = NaN

            ai = init_mm(input[i], &nobs, is_max)

            # Discard previous entries if we find new min or max
            if is_max:
                while not Q.empty() and ((ai >= input[Q.back()]) or
                                         (input[Q.back()] != input[Q.back()])):
                    Q.pop_back()
            else:
                while not Q.empty() and ((ai <= input[Q.back()]) or
                                         (input[Q.back()] != input[Q.back()])):
                    Q.pop_back()

            # Maintain window/nobs retention
            curr_win_size = endi[i + close_offset] - starti[i + close_offset]
            while not Q.empty() and Q.front() <= i - curr_win_size:
                Q.pop_front()
            while not W.empty() and W.front() <= i - curr_win_size:
                remove_mm(input[W.front()], &nobs)
                W.pop_front()

            Q.push_back(i)
            W.push_back(i)

        output[N-1] = calc_mm(minp, nobs, input[Q.front()])

    return output


cdef _roll_min_max_fixed(ndarray[numeric] input,
                         ndarray[int64_t] starti,
                         ndarray[int64_t] endi,
                         int64_t N,
                         int64_t win,
                         int64_t minp,
                         bint is_max):
    cdef:
        numeric ai
        bint should_replace
        int64_t i, removed, window_i,
        Py_ssize_t nobs = 0
        int64_t* death
        numeric* ring
        numeric* minvalue
        numeric* end
        numeric* last
        ndarray[double_t, ndim=1] output

    output = np.empty(N, dtype=float)
    # setup the rings of death!
    ring = <numeric *>malloc(win * sizeof(numeric))
    death = <int64_t *>malloc(win * sizeof(int64_t))

    end = ring + win
    last = ring
    minvalue = ring
    ai = input[0]
    minvalue[0] = init_mm(input[0], &nobs, is_max)
    death[0] = win
    nobs = 0

    with nogil:

        for i in range(N):
            ai = init_mm(input[i], &nobs, is_max)

            if i >= win:
                remove_mm(input[i - win], &nobs)

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


def roll_quantile(ndarray[float64_t, cast=True] input, int64_t win,
                  int64_t minp, object index, object closed,
                  double quantile, str interpolation):
    """
    O(N log(window)) implementation using skip list
    """
    cdef:
        double val, prev, midpoint, idx_with_fraction
        skiplist_t *skiplist
        int64_t nobs = 0, i, j, s, e, N
        Py_ssize_t idx
        bint is_variable
        ndarray[int64_t] start, end
        ndarray[double_t] output
        double vlow, vhigh
        InterpolationType interpolation_type
        int ret = 0

    if quantile <= 0.0 or quantile >= 1.0:
        raise ValueError("quantile value {0} not in [0, 1]".format(quantile))

    try:
        interpolation_type = interpolation_types[interpolation]
    except KeyError:
        raise ValueError("Interpolation '{}' is not supported"
                         .format(interpolation))

    # we use the Fixed/Variable Indexer here as the
    # actual skiplist ops outweigh any window computation costs
    start, end, N, win, minp, is_variable = get_window_indexer(
        input, win,
        minp, index, closed,
        use_mock=False)
    output = np.empty(N, dtype=float)
    skiplist = skiplist_init(<int>win)
    if skiplist == NULL:
        raise MemoryError("skiplist_init failed")

    with nogil:
        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0:

                # setup
                val = input[i]
                if val == val:
                    nobs += 1
                    skiplist_insert(skiplist, val)

            else:

                # calculate deletes
                for j in range(start[i - 1], s):
                    val = input[j]
                    if val == val:
                        skiplist_remove(skiplist, val)
                        nobs -= 1

                # calculate adds
                for j in range(end[i - 1], e):
                    val = input[j]
                    if val == val:
                        nobs += 1
                        skiplist_insert(skiplist, val)

            if nobs >= minp:
                if nobs == 1:
                    # Single value in skip list
                    output[i] = skiplist_get(skiplist, 0, &ret)
                else:
                    idx_with_fraction = quantile * (nobs - 1)
                    idx = <int> idx_with_fraction

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
                        output[i] = <double> (vlow + vhigh) / 2
            else:
                output[i] = NaN

    skiplist_destroy(skiplist)

    return output


def roll_generic(object obj,
                 int64_t win, int64_t minp, object index, object closed,
                 int offset, object func, bint raw,
                 object args, object kwargs):
    cdef:
        ndarray[double_t] output, counts, bufarr
        ndarray[float64_t, cast=True] arr
        float64_t *buf
        float64_t *oldbuf
        int64_t nobs = 0, i, j, s, e, N
        bint is_variable
        ndarray[int64_t] start, end

    n = len(obj)
    if n == 0:
        return obj

    arr = np.asarray(obj)

    # ndarray input
    if raw:
        if not arr.flags.c_contiguous:
            arr = arr.copy('C')

    counts = roll_sum(np.concatenate([np.isfinite(arr).astype(float),
                                      np.array([0.] * offset)]),
                      win, minp, index, closed)[offset:]

    start, end, N, win, minp, is_variable = get_window_indexer(arr, win,
                                                               minp, index,
                                                               closed,
                                                               floor=0)

    output = np.empty(N, dtype=float)

    if is_variable:
        # variable window arr or series

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

    elif not raw:
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
        for i from 0 <= i < (int_min(win, N) - offset):
            if counts[i] >= minp:
                output[i] = func(arr[0: (i + offset + 1)], *args, **kwargs)
            else:
                output[i] = NaN

        # remaining full-length windows
        buf = <float64_t *> arr.data
        bufarr = np.empty(win, dtype=float)
        oldbuf = <float64_t *> bufarr.data
        for i from (win - offset) <= i < (N - offset):
            buf = buf + 1
            bufarr.data = <char *> buf
            if counts[i] >= minp:
                output[i] = func(bufarr, *args, **kwargs)
            else:
                output[i] = NaN
        bufarr.data = <char *> oldbuf

        # truncated windows at the end
        for i from int_max(N - offset, 0) <= i < N:
            if counts[i] >= minp:
                output[i] = func(arr[int_max(i + offset - win + 1, 0): N],
                                 *args,
                                 **kwargs)
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
        for win_i in range(win_n):
            val_win = weights[win_i]
            if val_win != val_win:
                continue

            for in_i from 0 <= in_i < in_n - (win_n - win_i) + 1:
                val_in = input[in_i]
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

            for in_i from 0 <= in_i < in_n - (win_n - win_i) + 1:
                val_in = input[in_i]

                if val_in == val_in:
                    output[in_i + (win_n - win_i) - 1] += val_in * val_win
                    counts[in_i + (win_n - win_i) - 1] += 1

        for in_i in range(in_n):
            c = counts[in_i]
            if c < minp:
                output[in_i] = NaN

    return output

# ----------------------------------------------------------------------
# Exponentially weighted moving average


def ewma(double_t[:] vals, double_t com, int adjust, int ignore_na, int minp):
    """
    Compute exponentially-weighted moving average using center-of-mass.

    Parameters
    ----------
    vals : ndarray (float64 type)
    com : float64
    adjust: int
    ignore_na: int
    minp: int

    Returns
    -------
    y : ndarray
    """

    cdef:
        Py_ssize_t N = len(vals)
        ndarray[double_t] output = np.empty(N, dtype=float)
        double alpha, old_wt_factor, new_wt, weighted_avg, old_wt, cur
        Py_ssize_t i, nobs

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

    for i in range(1, N):
        cur = vals[i]
        is_observation = (cur == cur)
        nobs += int(is_observation)
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


def ewmcov(double_t[:] input_x, double_t[:] input_y,
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

    cdef:
        Py_ssize_t N = len(input_x)
        double alpha, old_wt_factor, new_wt, mean_x, mean_y, cov
        double sum_wt, sum_wt2, old_wt, cur_x, cur_y, old_mean_x, old_mean_y
        Py_ssize_t i, nobs
        ndarray[double_t] output

    if len(input_y) != N:
        raise ValueError("arrays are of different lengths "
                         "({N} and {len_y})".format(N=N, len_y=len(input_y)))

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

    for i in range(1, N):
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
