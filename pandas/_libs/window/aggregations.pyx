# cython: boundscheck=False, wraparound=False, cdivision=True

from libc.math cimport (
    fabs,
    signbit,
    sqrt,
)
from libcpp.deque cimport deque
from libcpp.stack cimport stack
from libcpp.unordered_map cimport unordered_map

from pandas._libs.algos cimport TiebreakEnumType

import numpy as np

cimport numpy as cnp
from numpy cimport (
    float32_t,
    float64_t,
    int64_t,
    ndarray,
)

cnp.import_array()

import cython

from pandas._libs.algos import is_monotonic


cdef extern from "pandas/skiplist.h":
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
    int skiplist_rank(skiplist_t*, double) nogil
    int skiplist_min_rank(skiplist_t*, double) nogil

cdef:
    float32_t MINfloat32 = -np.inf
    float64_t MINfloat64 = -np.inf

    float32_t MAXfloat32 = np.inf
    float64_t MAXfloat64 = np.inf

    float64_t NaN = <float64_t>np.nan
    float64_t EpsF64 = np.finfo(np.float64).eps

    # Consider an operation ill-conditioned if
    # it will only have up to 3 significant digits in base 10 remaining.
    # https://en.wikipedia.org/wiki/Condition_number
    float64_t InvCondTol = EpsF64 * 1e3

cdef bint is_monotonic_increasing_start_end_bounds(
    ndarray[int64_t, ndim=1] start, ndarray[int64_t, ndim=1] end
):
    return is_monotonic(start, False)[0] and is_monotonic(end, False)[0]

# ----------------------------------------------------------------------
# Rolling sum


cdef float64_t calc_sum(int64_t minp, int64_t nobs, float64_t sum_x,
                        int64_t num_consecutive_same_value, float64_t prev_value
                        ) noexcept nogil:
    cdef:
        float64_t result

    if nobs == 0 == minp:
        result = 0
    elif nobs >= minp:
        if num_consecutive_same_value >= nobs:
            result = prev_value * nobs
        else:
            result = sum_x
    else:
        result = NaN

    return result


cdef void add_sum(float64_t val, int64_t *nobs, float64_t *sum_x,
                  float64_t *compensation, int64_t *num_consecutive_same_value,
                  float64_t *prev_value) noexcept nogil:
    """ add a value from the sum calc using Kahan summation """

    cdef:
        float64_t y, t

    # Not NaN
    if val == val:
        nobs[0] = nobs[0] + 1
        y = val - compensation[0]
        t = sum_x[0] + y
        compensation[0] = t - sum_x[0] - y
        sum_x[0] = t

        # GH#42064, record num of same values to remove floating point artifacts
        if val == prev_value[0]:
            num_consecutive_same_value[0] += 1
        else:
            # reset to 1 (include current value itself)
            num_consecutive_same_value[0] = 1
        prev_value[0] = val


cdef void remove_sum(float64_t val, int64_t *nobs, float64_t *sum_x,
                     float64_t *compensation) noexcept nogil:
    """ remove a value from the sum calc using Kahan summation """

    cdef:
        float64_t y, t

    # Not NaN
    if val == val:
        nobs[0] = nobs[0] - 1
        y = - val - compensation[0]
        t = sum_x[0] + y
        compensation[0] = t - sum_x[0] - y
        sum_x[0] = t


def roll_sum(const float64_t[:] values, ndarray[int64_t] start,
             ndarray[int64_t] end, int64_t minp) -> np.ndarray:
    cdef:
        Py_ssize_t i, j
        float64_t sum_x, compensation_add, compensation_remove, prev_value
        int64_t s, e, num_consecutive_same_value
        int64_t nobs = 0, N = len(start)
        ndarray[float64_t] output
        bint is_monotonic_increasing_bounds

    is_monotonic_increasing_bounds = is_monotonic_increasing_start_end_bounds(
        start, end
    )
    output = np.empty(N, dtype=np.float64)

    with nogil:

        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0 or not is_monotonic_increasing_bounds or s >= end[i - 1]:

                # setup
                prev_value = values[s]
                num_consecutive_same_value = 0
                sum_x = compensation_add = compensation_remove = 0
                nobs = 0
                for j in range(s, e):
                    add_sum(values[j], &nobs, &sum_x, &compensation_add,
                            &num_consecutive_same_value, &prev_value)

            else:

                # calculate deletes
                for j in range(start[i - 1], s):
                    remove_sum(values[j], &nobs, &sum_x, &compensation_remove)

                # calculate adds
                for j in range(end[i - 1], e):
                    add_sum(values[j], &nobs, &sum_x, &compensation_add,
                            &num_consecutive_same_value, &prev_value)

            output[i] = calc_sum(
                minp, nobs, sum_x, num_consecutive_same_value, prev_value
            )

            if not is_monotonic_increasing_bounds:
                nobs = 0
                sum_x = 0.0
                compensation_remove = 0.0

    return output


# ----------------------------------------------------------------------
# Rolling mean


cdef float64_t calc_mean(int64_t minp, Py_ssize_t nobs, Py_ssize_t neg_ct,
                         float64_t sum_x, int64_t num_consecutive_same_value,
                         float64_t prev_value) noexcept nogil:
    cdef:
        float64_t result

    if nobs >= minp and nobs > 0:
        result = sum_x / <float64_t>nobs
        if num_consecutive_same_value >= nobs:
            result = prev_value
        elif neg_ct == 0 and result < 0:
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


cdef void add_mean(
    float64_t val,
    Py_ssize_t *nobs,
    float64_t *sum_x,
    Py_ssize_t *neg_ct,
    float64_t *compensation,
    int64_t *num_consecutive_same_value,
    float64_t *prev_value
) noexcept nogil:
    """ add a value from the mean calc using Kahan summation """
    cdef:
        float64_t y, t

    # Not NaN
    if val == val:
        nobs[0] = nobs[0] + 1
        y = val - compensation[0]
        t = sum_x[0] + y
        compensation[0] = t - sum_x[0] - y
        sum_x[0] = t
        if signbit(val):
            neg_ct[0] = neg_ct[0] + 1

        # GH#42064, record num of same values to remove floating point artifacts
        if val == prev_value[0]:
            num_consecutive_same_value[0] += 1
        else:
            # reset to 1 (include current value itself)
            num_consecutive_same_value[0] = 1
        prev_value[0] = val


cdef void remove_mean(float64_t val, Py_ssize_t *nobs, float64_t *sum_x,
                      Py_ssize_t *neg_ct, float64_t *compensation) noexcept nogil:
    """ remove a value from the mean calc using Kahan summation """
    cdef:
        float64_t y, t

    if val == val:
        nobs[0] = nobs[0] - 1
        y = - val - compensation[0]
        t = sum_x[0] + y
        compensation[0] = t - sum_x[0] - y
        sum_x[0] = t
        if signbit(val):
            neg_ct[0] = neg_ct[0] - 1


def roll_mean(const float64_t[:] values, ndarray[int64_t] start,
              ndarray[int64_t] end, int64_t minp) -> np.ndarray:
    cdef:
        float64_t val, compensation_add, compensation_remove, sum_x, prev_value
        int64_t s, e, num_consecutive_same_value
        Py_ssize_t nobs, i, j, neg_ct, N = len(start)
        ndarray[float64_t] output
        bint is_monotonic_increasing_bounds

    is_monotonic_increasing_bounds = is_monotonic_increasing_start_end_bounds(
        start, end
    )
    output = np.empty(N, dtype=np.float64)

    with nogil:

        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0 or not is_monotonic_increasing_bounds or s >= end[i - 1]:

                # setup
                compensation_add = compensation_remove = sum_x = 0
                nobs = neg_ct = 0
                prev_value = values[s]
                num_consecutive_same_value = 0
                for j in range(s, e):
                    val = values[j]
                    add_mean(val, &nobs, &sum_x, &neg_ct, &compensation_add,
                             &num_consecutive_same_value, &prev_value)

            else:

                # calculate deletes
                for j in range(start[i - 1], s):
                    val = values[j]
                    remove_mean(val, &nobs, &sum_x, &neg_ct, &compensation_remove)

                # calculate adds
                for j in range(end[i - 1], e):
                    val = values[j]
                    add_mean(val, &nobs, &sum_x, &neg_ct, &compensation_add,
                             &num_consecutive_same_value, &prev_value)

            output[i] = calc_mean(
                minp, nobs, neg_ct, sum_x, num_consecutive_same_value, prev_value
            )

            if not is_monotonic_increasing_bounds:
                nobs = 0
                neg_ct = 0
                sum_x = 0.0
                compensation_remove = 0.0
    return output

# ----------------------------------------------------------------------
# Rolling variance


cdef float64_t calc_var(
    int64_t minp,
    int ddof,
    float64_t nobs,
    float64_t ssqdm_x,
) noexcept nogil:
    cdef:
        float64_t result

    # Variance is unchanged if no observation is added or removed
    if (nobs >= minp) and (nobs > ddof):
        result = ssqdm_x / (nobs - <float64_t>ddof)
    else:
        result = NaN

    return result


cdef void add_var(
    float64_t val,
    float64_t *nobs,
    float64_t *mean_x,
    float64_t *ssqdm_x,
    float64_t *compensation,
    bint *numerically_unstable,
) noexcept nogil:
    """ add a value from the var calc """
    cdef:
        float64_t delta, prev_mean, y, t
        float64_t prev_m2 = ssqdm_x[0]

    # GH#21813, if msvc 2017 bug is resolved, we should be OK with != instead of `isnan`
    if val != val:
        return

    nobs[0] = nobs[0] + 1

    # Welford's method for the online variance-calculation
    # using Kahan summation
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    prev_mean = mean_x[0] - compensation[0]
    y = val - compensation[0]
    t = y - mean_x[0]
    compensation[0] = t + mean_x[0] - y
    delta = t
    if nobs[0]:
        mean_x[0] = mean_x[0] + delta / nobs[0]
    else:
        mean_x[0] = 0
    ssqdm_x[0] = ssqdm_x[0] + (val - prev_mean) * (val - mean_x[0])

    if prev_m2 * InvCondTol > ssqdm_x[0]:
        # possible catastrophic cancellation
        numerically_unstable[0] = True


cdef void remove_var(
    float64_t val,
    float64_t *nobs,
    float64_t *mean_x,
    float64_t *ssqdm_x,
    float64_t *compensation,
    bint *numerically_unstable,
) noexcept nogil:
    """ remove a value from the var calc """
    cdef:
        float64_t delta, prev_mean, y, t
        float64_t prev_m2 = ssqdm_x[0]
    if val == val:
        nobs[0] = nobs[0] - 1
        if nobs[0]:
            # Welford's method for the online variance-calculation
            # using Kahan summation
            # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
            prev_mean = mean_x[0] - compensation[0]
            y = val - compensation[0]
            t = y - mean_x[0]
            compensation[0] = t + mean_x[0] - y
            delta = t
            mean_x[0] = mean_x[0] - delta / nobs[0]
            ssqdm_x[0] = ssqdm_x[0] - (val - prev_mean) * (val - mean_x[0])

            if prev_m2 * InvCondTol > ssqdm_x[0]:
                # possible catastrophic cancellation
                numerically_unstable[0] = True
        else:
            mean_x[0] = 0
            ssqdm_x[0] = 0
            numerically_unstable[0] = False


def roll_var(const float64_t[:] values, ndarray[int64_t] start,
             ndarray[int64_t] end, int64_t minp, int ddof=1) -> np.ndarray:
    """
    Numerically stable implementation using Welford's method.
    """
    cdef:
        float64_t mean_x, ssqdm_x, nobs, compensation_add,
        float64_t compensation_remove
        int64_t s, e
        Py_ssize_t i, j, N = len(start)
        ndarray[float64_t] output
        bint is_monotonic_increasing_bounds
        bint requires_recompute, numerically_unstable

    minp = max(minp, 1)
    is_monotonic_increasing_bounds = is_monotonic_increasing_start_end_bounds(
        start, end
    )
    output = np.empty(N, dtype=np.float64)

    with nogil:

        for i in range(0, N):

            s = start[i]
            e = end[i]

            # Over the first window, observations can only be added
            # never removed
            requires_recompute = (
                i == 0
                or not is_monotonic_increasing_bounds
                or s >= end[i - 1]
            )

            if not requires_recompute:
                # After the first window, observations can both be added
                # and removed

                # calculate deletes
                for j in range(start[i - 1], s):
                    remove_var(values[j], &nobs, &mean_x, &ssqdm_x,
                               &compensation_remove, &numerically_unstable)

                # calculate adds
                for j in range(end[i - 1], e):
                    add_var(values[j], &nobs, &mean_x, &ssqdm_x, &compensation_add,
                            &numerically_unstable)

            if requires_recompute or numerically_unstable:

                mean_x = ssqdm_x = nobs = compensation_add = compensation_remove = 0
                for j in range(s, e):
                    add_var(values[j], &nobs, &mean_x, &ssqdm_x, &compensation_add,
                            &numerically_unstable)
                numerically_unstable = False

            output[i] = calc_var(minp, ddof, nobs, ssqdm_x)

            if not is_monotonic_increasing_bounds:
                nobs = 0.0
                mean_x = 0.0
                ssqdm_x = 0.0
                compensation_remove = 0.0

    return output

# ----------------------------------------------------------------------
# Rolling skewness


cdef float64_t calc_skew(int64_t minp, int64_t nobs,
                         float64_t mean, float64_t m2, float64_t m3,
                         int64_t num_consecutive_same_value
                         ) noexcept nogil:
    cdef:
        float64_t result, dnobs
        float64_t moments_ratio, correction

    if nobs >= minp:
        dnobs = <float64_t>nobs

        if nobs < 3:
            result = NaN
        # GH 42064 46431
        # uniform case, force result to be 0
        elif num_consecutive_same_value >= nobs:
            result = 0.0
        # #18044: with degenerate distribution, floating issue will
        #         cause m2 != 0. and cause the result is a very
        #         large number.
        #
        #         in core/nanops.py nanskew/nankurt call the function
        #         _zero_out_fperr(m2) to fix floating error.
        #         if the variance is less than 1e-14, it could be
        #         treat as zero, here we follow the original
        #         skew/kurt behaviour to check m2 <= n * 1e-14
        elif m2 <= dnobs * 1e-14:
            result = NaN
        else:
            moments_ratio = m3 / (m2 * sqrt(m2))
            correction = dnobs * sqrt((dnobs - 1)) / (dnobs - 2)
            result = moments_ratio * correction
    else:
        result = NaN

    return result


cdef void add_skew(float64_t val, int64_t *nobs,
                   float64_t *mean, float64_t *m2,
                   float64_t *m3,
                   bint *numerically_unstable,
                   int64_t *num_consecutive_same_value,
                   float64_t *prev_value,
                   ) noexcept nogil:
    """ add a value from the skew calc """
    cdef:
        float64_t n, delta, delta_n, term1, m3_update, new_m3

    # Formulas adapted from
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics

    # Not NaN
    if val == val:
        nobs[0] += 1
        n = <float64_t>(nobs[0])
        delta = val - mean[0]
        delta_n = delta / n
        term1 = delta * delta_n * (n - 1.0)

        m3_update = delta_n * (term1 * (n - 2.0) - 3.0 * m2[0])
        new_m3 = m3[0] + m3_update
        if (fabs(m3_update) + fabs(m3[0])) * InvCondTol > fabs(new_m3):
            # possible catastrophic cancellation
            numerically_unstable[0] = True

        m3[0] = new_m3
        m2[0] += term1
        mean[0] += delta_n

        # GH#42064, record num of same values to remove floating point artifacts
        if val == prev_value[0]:
            num_consecutive_same_value[0] += 1
        else:
            # reset to 1 (include current value itself)
            num_consecutive_same_value[0] = 1
        prev_value[0] = val


cdef void remove_skew(float64_t val, int64_t *nobs,
                      float64_t *mean, float64_t *m2,
                      float64_t *m3,
                      bint *numerically_unstable) noexcept nogil:
    """ remove a value from the skew calc """
    cdef:
        float64_t n, delta, delta_n, term1, m3_update, new_m3

    # This is the online update for the central moments
    # when we remove an observation.
    #
    # δ = x - m_{n+1}
    # m_{n} = m_{n+1} - (δ / n)
    # m²_n = Σ_{i=1}^{n+1}(x_i - m_{n})² - (x - m_{n})² # uses new mean
    #      = m²_{n+1} - (δ²/n)*(n+1)
    # m³_n = Σ_{i=1}^{n+1}(x_i - m_{n})³ - (x - m_{n})³ # uses new mean
    #      = m³_{n+1} - (δ³/n²)*(n+1)*(n+2) + 3 * m²_{n+1}*(δ/n)

    # Not NaN
    if val == val:
        nobs[0] -= 1
        n = <float64_t>(nobs[0])
        delta = val - mean[0]
        delta_n = delta / n
        term1 = delta_n * delta * (n + 1.0)

        m3_update = delta_n * (term1 * (n + 2.0) - 3.0 * m2[0])
        new_m3 = m3[0] - m3_update

        if (fabs(m3_update) + fabs(m3[0])) * InvCondTol > fabs(new_m3):
            # possible catastrophic cancellation
            numerically_unstable[0] = True

        m3[0] = new_m3
        m2[0] -= term1
        mean[0] -= delta_n


def roll_skew(const float64_t[:] values, ndarray[int64_t] start,
              ndarray[int64_t] end, int64_t minp) -> np.ndarray:
    cdef:
        Py_ssize_t i, j
        float64_t val
        float64_t mean, m2, m3
        float64_t prev_value
        int64_t nobs = 0, N = len(start)
        int64_t s, e, num_consecutive_same_value
        ndarray[float64_t] output
        bint is_monotonic_increasing_bounds
        bint requires_recompute, numerically_unstable = False

    minp = max(minp, 3)
    is_monotonic_increasing_bounds = is_monotonic_increasing_start_end_bounds(
        start, end
    )
    output = np.empty(N, dtype=np.float64)

    with nogil:
        for i in range(0, N):

            s = start[i]
            e = end[i]

            # Over the first window, observations can only be added
            # never removed
            requires_recompute = (
                i == 0
                or not is_monotonic_increasing_bounds
                or s >= end[i - 1]
            )

            if not requires_recompute:
                # After the first window, observations can both be added
                # and removed
                # calculate deletes
                for j in range(start[i - 1], s):
                    val = values[j]
                    remove_skew(val, &nobs, &mean, &m2, &m3, &numerically_unstable)

                # calculate adds
                for j in range(end[i - 1], e):
                    val = values[j]
                    add_skew(val, &nobs, &mean, &m2, &m3, &numerically_unstable,
                             &num_consecutive_same_value, &prev_value)

            if requires_recompute or numerically_unstable:

                prev_value = values[s]
                num_consecutive_same_value = 0

                mean = m2 = m3 = 0.0
                nobs = 0

                for j in range(s, e):
                    val = values[j]
                    add_skew(val, &nobs, &mean, &m2, &m3, &numerically_unstable,
                             &num_consecutive_same_value, &prev_value)

                numerically_unstable = False

            output[i] = calc_skew(minp, nobs, mean, m2, m3, num_consecutive_same_value)

            if not is_monotonic_increasing_bounds:
                nobs = 0
                mean = 0.0
                m2 = 0.0
                m3 = 0.0

    return output

# ----------------------------------------------------------------------
# Rolling kurtosis


cdef float64_t calc_kurt(int64_t minp, int64_t nobs,
                         float64_t m2, float64_t m4,
                         int64_t num_consecutive_same_value,
                         ) noexcept nogil:
    cdef:
        float64_t result, dnobs, term1, term2, inner, correction
        float64_t moments_ratio

    if nobs >= minp:
        if nobs < 4:
            result = NaN
        # GH 42064 46431
        # degenerate case, force result to be -3.
        elif num_consecutive_same_value >= nobs:
            result = -3.
        else:
            dnobs = <float64_t>nobs

            # #18044: with uniform distribution, floating issue will
            #         cause B != 0. and cause the result is a very
            #         large number.
            #
            #         in core/nanops.py nanskew/nankurt call the function
            #         _zero_out_fperr(m2) to fix floating error.
            #         if the variance is less than 1e-14, it could be
            #         treat as zero, here we follow the original
            #         skew/kurt behaviour to check m2 <= n * 1e-14
            if m2 <= dnobs * 1e-14:
                result = NaN
            else:
                moments_ratio = m4 / (m2 * m2)
                term1 = dnobs * (dnobs + 1.0) * moments_ratio
                term2 = 3.0 * (dnobs - 1.0)
                inner = term1 - term2

                correction = (dnobs - 1.0) / ((dnobs - 2.0) * (dnobs - 3.0))
                result = correction * inner
    else:
        result = NaN

    return result


cdef void add_kurt(float64_t val, int64_t *nobs,
                   float64_t *mean, float64_t *m2,
                   float64_t *m3, float64_t *m4,
                   bint *numerically_unstable,
                   int64_t *num_consecutive_same_value,
                   float64_t *prev_value
                   ) noexcept nogil:
    """ add a value from the kurotic calc """
    cdef:
        float64_t n, delta, delta_n, term1, m4_update, new_m4

    # Formulas adapted from
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics

    # Not NaN
    if val == val:
        nobs[0] += 1
        n = <float64_t>(nobs[0])
        delta = val - mean[0]
        delta_n = delta / n
        term1 = delta * delta_n * (n - 1.0)

        m4_update = delta_n * (
                -4.0 * m3[0]
                + delta_n * (
                    6 * m2[0] + term1 * (n * n - 3.0 * n + 3.0)
                    )
                )
        new_m4 = m4[0] + m4_update

        if (fabs(m4_update) + fabs(m4[0])) * InvCondTol > fabs(new_m4):
            # possible catastrophic cancellation
            numerically_unstable[0] = True

        m4[0] = new_m4
        m3[0] += delta_n * (term1 * (n - 2.0) - 3.0 * m2[0])
        m2[0] += term1
        mean[0] += delta_n

        # GH#42064, record num of same values to remove floating point artifacts
        if val == prev_value[0]:
            num_consecutive_same_value[0] += 1
        else:
            # reset to 1 (include current value itself)
            num_consecutive_same_value[0] = 1
        prev_value[0] = val


cdef void remove_kurt(float64_t val, int64_t *nobs,
                      float64_t *mean, float64_t *m2,
                      float64_t *m3, float64_t *m4,
                      bint *numerically_unstable,
                      ) noexcept nogil:
    """ remove a value from the kurotic calc """
    cdef:
        float64_t n, delta, delta_n, term1, m4_update, new_m4

    # Not NaN
    if val == val:
        nobs[0] -= 1
        n = <float64_t>(nobs[0])
        delta = val - mean[0]
        delta_n = delta / n
        term1 = delta_n * delta * (n + 1.0)

        m4_update = delta_n * (
                4.0 * m3[0]
                + delta_n * (
                    6.0 * m2[0]
                    - term1 * (n * n + 3.0 * n + 3.0)
                    )
                )
        new_m4 = m4[0] + m4_update

        if (fabs(m4_update) + fabs(m4[0])) * InvCondTol > fabs(new_m4):
            # possible catastrophic cancellation
            numerically_unstable[0] = True

        m4[0] = new_m4
        m3[0] -= delta_n * (term1 * (n + 2.0) - 3.0 * m2[0])
        m2[0] -= term1
        mean[0] -= delta_n


def roll_kurt(const float64_t[:] values, ndarray[int64_t] start,
              ndarray[int64_t] end, int64_t minp) -> np.ndarray:
    cdef:
        Py_ssize_t i, j
        float64_t mean, m2, m3, m4
        float64_t prev_value
        int64_t nobs, s, e, num_consecutive_same_value
        int64_t N = len(start)
        ndarray[float64_t] output
        bint is_monotonic_increasing_bounds
        bint requires_recompute, numerically_unstable = False

    minp = max(minp, 4)
    is_monotonic_increasing_bounds = is_monotonic_increasing_start_end_bounds(
        start, end
    )
    output = np.empty(N, dtype=np.float64)

    with nogil:
        for i in range(0, N):

            s = start[i]
            e = end[i]

            # Over the first window, observations can only be added
            # never removed
            requires_recompute = (
                i == 0
                or not is_monotonic_increasing_bounds
                or s >= end[i - 1]
            )

            if not requires_recompute:

                # After the first window, observations can both be added
                # and removed
                # calculate deletes
                for j in range(start[i - 1], s):
                    remove_kurt(values[j], &nobs, &mean, &m2, &m3, &m4,
                                &numerically_unstable)

                # calculate adds
                for j in range(end[i - 1], e):
                    add_kurt(values[j], &nobs, &mean, &m2, &m3, &m4,
                             &numerically_unstable,
                             &num_consecutive_same_value, &prev_value)

            if requires_recompute or numerically_unstable:

                prev_value = values[s]
                num_consecutive_same_value = 0

                mean = m2 = m3 = m4 = 0.0
                nobs = 0
                for j in range(s, e):
                    add_kurt(values[j], &nobs, &mean, &m2, &m3, &m4,
                             &numerically_unstable,
                             &num_consecutive_same_value, &prev_value)

            output[i] = calc_kurt(minp, nobs, m2, m4,
                                  num_consecutive_same_value)

            if not is_monotonic_increasing_bounds:
                nobs = 0
                mean = 0.0
                m2 = 0.0
                m3 = 0.0
                m4 = 0.0

    return output


# ----------------------------------------------------------------------
# Rolling median, min, max


def roll_median_c(const float64_t[:] values, ndarray[int64_t] start,
                  ndarray[int64_t] end, int64_t minp) -> np.ndarray:
    cdef:
        Py_ssize_t i, j
        bint err = False, is_monotonic_increasing_bounds
        int midpoint, ret = 0
        int64_t nobs = 0, N = len(start), s, e, win
        float64_t val, res
        skiplist_t *sl
        ndarray[float64_t] output

    is_monotonic_increasing_bounds = is_monotonic_increasing_start_end_bounds(
        start, end
    )

    # we use the Fixed/Variable Indexer here as the
    # actual skiplist ops outweigh any window computation costs
    output = np.empty(N, dtype=np.float64)

    if (end - start).max() == 0:
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

            if i == 0 or not is_monotonic_increasing_bounds or s >= end[i - 1]:

                if i != 0:
                    skiplist_destroy(sl)
                    sl = skiplist_init(<int>win)
                    nobs = 0
                # setup
                for j in range(s, e):
                    val = values[j]
                    if val == val:
                        nobs += 1
                        err = skiplist_insert(sl, val) == -1
                        if err:
                            break

            else:

                # calculate adds
                for j in range(end[i - 1], e):
                    val = values[j]
                    if val == val:
                        nobs += 1
                        err = skiplist_insert(sl, val) == -1
                        if err:
                            break

                # calculate deletes
                for j in range(start[i - 1], s):
                    val = values[j]
                    if val == val:
                        skiplist_remove(sl, val)
                        nobs -= 1
            if nobs >= minp:
                midpoint = <int>(nobs / 2)
                if nobs % 2:
                    res = skiplist_get(sl, midpoint, &ret)
                else:
                    res = (skiplist_get(sl, midpoint, &ret) +
                           skiplist_get(sl, (midpoint - 1), &ret)) / 2
                if ret == 0:
                    res = NaN
            else:
                res = NaN

            output[i] = res

            if not is_monotonic_increasing_bounds:
                nobs = 0
                skiplist_destroy(sl)
                sl = skiplist_init(<int>win)

    skiplist_destroy(sl)
    if err:
        raise MemoryError("skiplist_insert failed")
    return output


# ----------------------------------------------------------------------

cdef int64_t bisect_left(
    deque[int64_t]& a,
    int64_t x,
    int64_t lo=0,
    int64_t hi=-1
) nogil:
    """Same as https://docs.python.org/3/library/bisect.html."""

    cdef int64_t mid
    if hi == -1:
        hi = a.size()
    while lo < hi:
        mid = (lo + hi) // 2
        if a.at(mid) < x:
            lo = mid + 1
        else:
            hi = mid
    return lo

from libc.math cimport isnan

# Prior version of moving maximum / minimum code taken from Bottleneck
# Licence at LICENSES/BOTTLENECK_LICENCE


def roll_max(ndarray[float64_t] values, ndarray[int64_t] start,
             ndarray[int64_t] end, int64_t minp) -> np.ndarray:
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

    Returns
    -------
    np.ndarray[float]
    """
    return _roll_min_max(values, start, end, minp, is_max=1)


def roll_min(ndarray[float64_t] values, ndarray[int64_t] start,
             ndarray[int64_t] end, int64_t minp) -> np.ndarray:
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

    Returns
    -------
    np.ndarray[float]
    """
    return _roll_min_max(values, start, end, minp, is_max=0)


def _roll_min_max(
    ndarray[float64_t] values,
    ndarray[int64_t] start,
    ndarray[int64_t] end,
    int64_t minp,
    bint is_max
):
    cdef:
        Py_ssize_t i, i_next, k, valid_start, last_end, last_start, N = len(start)
        # Indices of bounded extrema in `values`. `candidates[i]` is always increasing.
        # `values[candidates[i]]` is decreasing for max and increasing for min.
        deque candidates[int64_t]
        # Indices of largest windows that "cover" preceding windows.
        stack dominators[int64_t]
        ndarray[float64_t, ndim=1] output

        Py_ssize_t this_start, this_end, stash_start
        int64_t q_idx

    output = np.empty(N, dtype=np.float64)
    candidates = deque[int64_t]()
    dominators = stack[int64_t]()

    # This function was "ported" / translated from sliding_min_max()
    # in /pandas/core/_numba/kernels/min_max_.py.
    # (See there for credits and some comments.)
    # Code translation assumptions/rules:
    # - min_periods --> minp
    # - deque[0] --> front()
    # - deque[-1] --> back()
    # - stack[-1] --> top()
    # - bool(stack/deque) --> !empty()
    # - deque.append()    --> push_back()
    # - stack.append()    --> push()
    # - deque.popleft     --> pop_front()
    # - deque.pop()       --> pop_back()

    with nogil:
        if minp < 1:
            minp = 1

        if N>2:
            i_next = N - 1
            for i in range(N - 2, -1, -1):
                if start[i_next] < start[i] \
                    and (
                           dominators.empty()
                        or start[dominators.top()] > start[i_next]
                ):
                    dominators.push(i_next)
                i_next = i

        # NaN tracking to guarantee minp
        valid_start = -minp

        last_end = 0
        last_start = -1

        for i in range(N):
            this_start = start[i]
            this_end = end[i]

            if (not dominators.empty() and dominators.top() == i):
                dominators.pop()

            if not (this_end > last_end
                    or (this_end == last_end and this_start >= last_start)):
                raise ValueError(
                    "Start/End ordering requirement is violated at index {}".format(i))

            if dominators.empty():
                stash_start = this_start
            else:
                stash_start = min(this_start, start[dominators.top()])

            while not candidates.empty() and candidates.front() < stash_start:
                candidates.pop_front()

            for k in range(last_end, this_end):
                if not isnan(values[k]):
                    valid_start += 1
                    while valid_start >= 0 and isnan(values[valid_start]):
                        valid_start += 1

                    if is_max:
                        while (not candidates.empty()
                                and values[k] >= values[candidates.back()]):
                            candidates.pop_back()
                    else:
                        while (not candidates.empty()
                                and values[k] <= values[candidates.back()]):
                            candidates.pop_back()
                    candidates.push_back(k)

            if candidates.empty() or this_start > valid_start:
                output[i] = NaN
            elif candidates.front() >= this_start:
                # ^^ This is here to avoid costly bisection for fixed window sizes.
                output[i] = values[candidates.front()]
            else:
                q_idx = bisect_left(candidates, this_start, lo=1)
                output[i] = values[candidates[q_idx]]
            last_end = this_end
            last_start = this_start

    return output

# ----------------------------------------------------------------------
# Rolling first, last


def roll_first(const float64_t[:] values, ndarray[int64_t] start,
               ndarray[int64_t] end, int64_t minp) -> np.ndarray:
    return _roll_first_last(values, start, end, minp, is_first=1)


def roll_last(const float64_t[:] values, ndarray[int64_t] start,
              ndarray[int64_t] end, int64_t minp) -> np.ndarray:
    return _roll_first_last(values, start, end, minp, is_first=0)


cdef _roll_first_last(const float64_t[:] values, ndarray[int64_t] start,
                      ndarray[int64_t] end, int64_t minp, bint is_first):
    cdef:
        Py_ssize_t i, j, fl_idx
        bint is_monotonic_increasing_bounds
        int64_t nobs = 0, N = len(start), s, e
        float64_t val, res
        ndarray[float64_t] output

    is_monotonic_increasing_bounds = is_monotonic_increasing_start_end_bounds(
        start, end
    )

    output = np.empty(N, dtype=np.float64)

    if (end - start).max() == 0:
        output[:] = NaN
        return output

    with nogil:
        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0 or not is_monotonic_increasing_bounds or s >= end[i - 1]:
                fl_idx = -1
                nobs = 0
                for j in range(s, e):
                    val = values[j]
                    if val == val:
                        if not is_first or fl_idx < s:
                            fl_idx = j
                        nobs += 1
            else:
                # handle deletes
                for j in range(start[i - 1], s):
                    val = values[j]
                    if val == val:
                        nobs -= 1

                # update fl_idx if out of range, if first
                if is_first and fl_idx < s:
                    fl_idx = -1
                    for j in range(s, end[i - 1]):
                        val = values[j]
                        if val == val:
                            fl_idx = j
                            break

                # handle adds
                for j in range(end[i - 1], e):
                    val = values[j]
                    if val == val:
                        if not is_first or fl_idx < s:
                            fl_idx = j
                        nobs += 1

            if nobs >= minp and fl_idx >= s:
                res = values[fl_idx]
            else:
                res = NaN

            output[i] = res

            if not is_monotonic_increasing_bounds:
                nobs = 0

    return output


cdef enum InterpolationType:
    LINEAR,
    LOWER,
    HIGHER,
    NEAREST,
    MIDPOINT


interpolation_types = {
    "linear": LINEAR,
    "lower": LOWER,
    "higher": HIGHER,
    "nearest": NEAREST,
    "midpoint": MIDPOINT,
}


def roll_quantile(const float64_t[:] values, ndarray[int64_t] start,
                  ndarray[int64_t] end, int64_t minp,
                  float64_t quantile, str interpolation) -> np.ndarray:
    """
    O(N log(window)) implementation using skip list
    """
    cdef:
        Py_ssize_t i, j, s, e, N = len(start), idx
        int ret = 0
        int64_t nobs = 0, win
        float64_t val, idx_with_fraction
        float64_t vlow, vhigh
        skiplist_t *skiplist
        InterpolationType interpolation_type
        ndarray[float64_t] output

    if quantile <= 0.0 or quantile >= 1.0:
        raise ValueError(f"quantile value {quantile} not in [0, 1]")

    try:
        interpolation_type = interpolation_types[interpolation]
    except KeyError:
        raise ValueError(f"Interpolation '{interpolation}' is not supported")

    is_monotonic_increasing_bounds = is_monotonic_increasing_start_end_bounds(
        start, end
    )
    # we use the Fixed/Variable Indexer here as the
    # actual skiplist ops outweigh any window computation costs
    output = np.empty(N, dtype=np.float64)

    win = (end - start).max()
    if win == 0:
        output[:] = NaN
        return output
    skiplist = skiplist_init(<int>win)
    if skiplist == NULL:
        raise MemoryError("skiplist_init failed")

    with nogil:
        for i in range(0, N):
            s = start[i]
            e = end[i]

            if i == 0 or not is_monotonic_increasing_bounds or s >= end[i - 1]:
                if i != 0:
                    nobs = 0
                    skiplist_destroy(skiplist)
                    skiplist = skiplist_init(<int>win)

                # setup
                for j in range(s, e):
                    val = values[j]
                    if val == val:
                        nobs += 1
                        skiplist_insert(skiplist, val)

            else:
                # calculate adds
                for j in range(end[i - 1], e):
                    val = values[j]
                    if val == val:
                        nobs += 1
                        skiplist_insert(skiplist, val)

                # calculate deletes
                for j in range(start[i - 1], s):
                    val = values[j]
                    if val == val:
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
                        output[i] = (vlow + (vhigh - vlow) *
                                     (idx_with_fraction - idx))
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

                    if ret == 0:
                        output[i] = NaN
            else:
                output[i] = NaN

    skiplist_destroy(skiplist)

    return output


rolling_rank_tiebreakers = {
    "average": TiebreakEnumType.TIEBREAK_AVERAGE,
    "min": TiebreakEnumType.TIEBREAK_MIN,
    "max": TiebreakEnumType.TIEBREAK_MAX,
}


def roll_rank(const float64_t[:] values, ndarray[int64_t] start,
              ndarray[int64_t] end, int64_t minp, bint percentile,
              str method, bint ascending) -> np.ndarray:
    """
    O(N log(window)) implementation using skip list

    derived from roll_quantile
    """
    cdef:
        Py_ssize_t i, j, s, e, N = len(start)
        float64_t rank_min = 0, rank = 0
        int64_t nobs = 0, win
        float64_t val
        skiplist_t *skiplist
        float64_t[::1] output
        TiebreakEnumType rank_type

    try:
        rank_type = rolling_rank_tiebreakers[method]
    except KeyError:
        raise ValueError(f"Method '{method}' is not supported")

    is_monotonic_increasing_bounds = is_monotonic_increasing_start_end_bounds(
        start, end
    )
    # we use the Fixed/Variable Indexer here as the
    # actual skiplist ops outweigh any window computation costs
    output = np.empty(N, dtype=np.float64)

    win = (end - start).max()
    if win == 0:
        output[:] = NaN
        return np.asarray(output)
    skiplist = skiplist_init(<int>win)
    if skiplist == NULL:
        raise MemoryError("skiplist_init failed")

    with nogil:
        for i in range(N):
            s = start[i]
            e = end[i]

            if i == 0 or not is_monotonic_increasing_bounds or s >= end[i - 1]:
                if i != 0:
                    nobs = 0
                    skiplist_destroy(skiplist)
                    skiplist = skiplist_init(<int>win)

                # setup
                for j in range(s, e):
                    val = values[j] if ascending else -values[j]
                    if val == val:
                        nobs += 1
                        rank = skiplist_insert(skiplist, val)
                        if rank == -1:
                            raise MemoryError("skiplist_insert failed")
                        if rank_type == TiebreakEnumType.TIEBREAK_AVERAGE:
                            # The average rank of `val` is the sum of the ranks of all
                            # instances of `val` in the skip list divided by the number
                            # of instances. The sum of consecutive integers from 1 to N
                            # is N * (N + 1) / 2.
                            # The sum of the ranks is the sum of integers from the
                            # lowest rank to the highest rank, which is the sum of
                            # integers from 1 to the highest rank minus the sum of
                            # integers from 1 to one less than the lowest rank.
                            rank_min = skiplist_min_rank(skiplist, val)
                            rank = (((rank * (rank + 1) / 2)
                                    - ((rank_min - 1) * rank_min / 2))
                                    / (rank - rank_min + 1))
                        elif rank_type == TiebreakEnumType.TIEBREAK_MIN:
                            rank = skiplist_min_rank(skiplist, val)
                    else:
                        rank = NaN

            else:
                # calculate deletes
                for j in range(start[i - 1], s):
                    val = values[j] if ascending else -values[j]
                    if val == val:
                        skiplist_remove(skiplist, val)
                        nobs -= 1

                # calculate adds
                for j in range(end[i - 1], e):
                    val = values[j] if ascending else -values[j]
                    if val == val:
                        nobs += 1
                        rank = skiplist_insert(skiplist, val)
                        if rank == -1:
                            raise MemoryError("skiplist_insert failed")
                        if rank_type == TiebreakEnumType.TIEBREAK_AVERAGE:
                            rank_min = skiplist_min_rank(skiplist, val)
                            rank = (((rank * (rank + 1) / 2)
                                    - ((rank_min - 1) * rank_min / 2))
                                    / (rank - rank_min + 1))
                        elif rank_type == TiebreakEnumType.TIEBREAK_MIN:
                            rank = skiplist_min_rank(skiplist, val)
                    else:
                        rank = NaN
            if nobs >= minp:
                output[i] = rank / nobs if percentile else rank
            else:
                output[i] = NaN

    skiplist_destroy(skiplist)

    return np.asarray(output)


def roll_nunique(const float64_t[:] values, ndarray[int64_t] start,
                 ndarray[int64_t] end, int64_t minp) -> np.ndarray:
    """
    Rolling number of unique elements in the window
    """
    cdef:
        Py_ssize_t i, j, s, e, N = len(start)
        int64_t nobs = 0
        float64_t val
        float64_t[::1] output
        unordered_map[float64_t, int64_t] value_counts

    is_monotonic_increasing_bounds = is_monotonic_increasing_start_end_bounds(
        start, end
    )
    output = np.empty(N, dtype=np.float64)
    value_counts = unordered_map[float64_t, int64_t]()

    with nogil:
        for i in range(N):
            s = start[i]
            e = end[i]

            if i == 0 or not is_monotonic_increasing_bounds or s >= end[i - 1]:
                if i != 0:
                    nobs = 0
                    value_counts.clear()

                # setup
                for j in range(s, e):
                    val = values[j]
                    if val == val:
                        nobs += 1
                        value_counts[val] += 1

            else:
                # calculate deletes
                for j in range(start[i - 1], s):
                    val = values[j]
                    if val == val:
                        value_counts[val] -= 1
                        if value_counts[val] == 0:
                            value_counts.erase(val)
                        nobs -= 1

                # calculate adds
                for j in range(end[i - 1], e):
                    val = values[j]
                    if val == val:
                        nobs += 1
                        value_counts[val] += 1

            if nobs >= minp:
                output[i] = value_counts.size()
            else:
                output[i] = NaN

    return np.asarray(output)


def roll_apply(object obj,
               ndarray[int64_t] start, ndarray[int64_t] end,
               int64_t minp,
               object function, bint raw,
               tuple args, dict kwargs) -> np.ndarray:
    cdef:
        ndarray[float64_t] output, counts
        ndarray[float64_t, cast=True] arr
        Py_ssize_t i, s, e, N = len(start), n = len(obj)

    if n == 0:
        return np.array([], dtype=np.float64)

    arr = np.asarray(obj)

    # ndarray input
    if raw and not arr.flags.c_contiguous:
        arr = arr.copy("C")

    counts = roll_sum(np.isfinite(arr).astype(float), start, end, minp)

    output = np.empty(N, dtype=np.float64)

    for i in range(N):

        s = start[i]
        e = end[i]

        if counts[i] >= minp:
            if raw:
                output[i] = function(arr[s:e], *args, **kwargs)
            else:
                output[i] = function(obj.iloc[s:e], *args, **kwargs)
        else:
            output[i] = NaN

    return output


# ----------------------------------------------------------------------
# Rolling sum and mean for weighted window


def roll_weighted_sum(
    const float64_t[:] values, const float64_t[:] weights, int minp
) -> np.ndarray:
    return _roll_weighted_sum_mean(values, weights, minp, avg=0)


def roll_weighted_mean(
    const float64_t[:] values, const float64_t[:] weights, int minp
) -> np.ndarray:
    return _roll_weighted_sum_mean(values, weights, minp, avg=1)


cdef float64_t[:] _roll_weighted_sum_mean(const float64_t[:] values,
                                          const float64_t[:] weights,
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

    elif minp > in_n:
        minp = in_n + 1

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

    return output


# ----------------------------------------------------------------------
# Rolling var for weighted window


cdef float64_t calc_weighted_var(float64_t t,
                                 float64_t sum_w,
                                 Py_ssize_t win_n,
                                 unsigned int ddof,
                                 float64_t nobs,
                                 int64_t minp) noexcept nogil:
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


cdef void add_weighted_var(float64_t val,
                           float64_t w,
                           float64_t *t,
                           float64_t *sum_w,
                           float64_t *mean,
                           float64_t *nobs) noexcept nogil:
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

    if val != val:
        return

    nobs[0] = nobs[0] + 1

    q = val - mean[0]
    temp = sum_w[0] + w
    r = q * w / temp

    mean[0] = mean[0] + r
    t[0] = t[0] + r * sum_w[0] * q
    sum_w[0] = temp


cdef void remove_weighted_var(float64_t val,
                              float64_t w,
                              float64_t *t,
                              float64_t *sum_w,
                              float64_t *mean,
                              float64_t *nobs) noexcept nogil:
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

    if val == val:
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


def roll_weighted_var(const float64_t[:] values, const float64_t[:] weights,
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
    output = np.empty(n, dtype=np.float64)

    with nogil:

        for i in range(min(win_n, n)):
            add_weighted_var(values[i], weights[i], &t,
                             &sum_w, &mean, &nobs)

            output[i] = calc_weighted_var(t, sum_w, win_n,
                                          ddof, nobs, minp)

        for i in range(win_n, n):
            val = values[i]
            pre_val = values[i - win_n]

            w = weights[i % win_n]
            pre_w = weights[(i - win_n) % win_n]

            if val == val:
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
# Exponentially weighted moving
@cython.cpow(True)
def ewm(const float64_t[:] vals, const int64_t[:] start, const int64_t[:] end,
        int minp, float64_t com, bint adjust, bint ignore_na,
        const float64_t[:] deltas=None, bint normalize=True) -> np.ndarray:
    """
    Compute exponentially-weighted moving average or sum using center-of-mass.

    Parameters
    ----------
    vals : ndarray (float64 type)
    start: ndarray (int64 type)
    end: ndarray (int64 type)
    minp : int
    com : float64
    adjust : bool
    ignore_na : bool
    deltas : ndarray (float64 type), optional. If None, implicitly assumes equally
             spaced points (used when `times` is not passed)
    normalize : bool, optional.
                If True, calculate the mean. If False, calculate the sum.

    Returns
    -------
    np.ndarray[float64_t]
    """

    cdef:
        Py_ssize_t i, j, s, e, nobs, win_size, N = len(vals), M = len(start)
        const float64_t[:] sub_vals
        const float64_t[:] sub_deltas=None
        ndarray[float64_t] sub_output, output = np.empty(N, dtype=np.float64)
        float64_t alpha, old_wt_factor, new_wt, weighted, old_wt, cur
        bint is_observation, use_deltas

    if N == 0:
        return output

    use_deltas = deltas is not None

    alpha = 1. / (1. + com)
    old_wt_factor = 1. - alpha
    new_wt = 1. if adjust else alpha

    for j in range(M):
        s = start[j]
        e = end[j]
        sub_vals = vals[s:e]
        # note that len(deltas) = len(vals) - 1 and deltas[i] is to be used in
        # conjunction with vals[i+1]
        if use_deltas:
            sub_deltas = deltas[s:e - 1]
        win_size = len(sub_vals)
        sub_output = np.empty(win_size, dtype=np.float64)

        weighted = sub_vals[0]
        is_observation = weighted == weighted
        nobs = int(is_observation)
        sub_output[0] = weighted if nobs >= minp else NaN
        old_wt = 1.

        with nogil:
            for i in range(1, win_size):
                cur = sub_vals[i]
                is_observation = cur == cur
                nobs += is_observation
                if weighted == weighted:

                    if is_observation or not ignore_na:
                        if normalize:
                            if use_deltas:
                                old_wt *= old_wt_factor ** sub_deltas[i - 1]
                            else:
                                old_wt *= old_wt_factor
                        else:
                            weighted = old_wt_factor * weighted
                        if is_observation:
                            if normalize:
                                # avoid numerical errors on constant series
                                if weighted != cur:
                                    if not adjust and com == 1:
                                        # update in case of irregular-interval series
                                        new_wt = 1. - old_wt
                                    weighted = old_wt * weighted + new_wt * cur
                                    weighted /= (old_wt + new_wt)
                                if adjust:
                                    old_wt += new_wt
                                else:
                                    old_wt = 1.
                            else:
                                weighted += cur
                elif is_observation:
                    weighted = cur

                sub_output[i] = weighted if nobs >= minp else NaN

        output[s:e] = sub_output

    return output


def ewmcov(const float64_t[:] input_x, const int64_t[:] start, const int64_t[:] end,
           int minp, const float64_t[:] input_y, float64_t com, bint adjust,
           bint ignore_na, bint bias) -> np.ndarray:
    """
    Compute exponentially-weighted moving variance using center-of-mass.

    Parameters
    ----------
    input_x : ndarray (float64 type)
    start: ndarray (int64 type)
    end: ndarray (int64 type)
    minp : int
    input_y : ndarray (float64 type)
    com : float64
    adjust : bool
    ignore_na : bool
    bias : bool

    Returns
    -------
    np.ndarray[float64_t]
    """

    cdef:
        Py_ssize_t i, j, s, e, win_size, nobs
        Py_ssize_t N = len(input_x), M = len(input_y), L = len(start)
        float64_t alpha, old_wt_factor, new_wt, mean_x, mean_y, cov
        float64_t sum_wt, sum_wt2, old_wt, cur_x, cur_y, old_mean_x, old_mean_y
        float64_t numerator, denominator
        const float64_t[:] sub_x_vals, sub_y_vals
        ndarray[float64_t] sub_out, output = np.empty(N, dtype=np.float64)
        bint is_observation

    if M != N:
        raise ValueError(f"arrays are of different lengths ({N} and {M})")

    if N == 0:
        return output

    alpha = 1. / (1. + com)
    old_wt_factor = 1. - alpha
    new_wt = 1. if adjust else alpha

    for j in range(L):
        s = start[j]
        e = end[j]
        sub_x_vals = input_x[s:e]
        sub_y_vals = input_y[s:e]
        win_size = len(sub_x_vals)
        sub_out = np.empty(win_size, dtype=np.float64)

        mean_x = sub_x_vals[0]
        mean_y = sub_y_vals[0]
        is_observation = (mean_x == mean_x) and (mean_y == mean_y)
        nobs = int(is_observation)
        if not is_observation:
            mean_x = NaN
            mean_y = NaN
        sub_out[0] = (0. if bias else NaN) if nobs >= minp else NaN
        cov = 0.
        sum_wt = 1.
        sum_wt2 = 1.
        old_wt = 1.

        with nogil:
            for i in range(1, win_size):
                cur_x = sub_x_vals[i]
                cur_y = sub_y_vals[i]
                is_observation = (cur_x == cur_x) and (cur_y == cur_y)
                nobs += is_observation
                if mean_x == mean_x:
                    if is_observation or not ignore_na:
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
                        if denominator > 0:
                            sub_out[i] = (numerator / denominator) * cov
                        else:
                            sub_out[i] = NaN
                    else:
                        sub_out[i] = cov
                else:
                    sub_out[i] = NaN

        output[s:e] = sub_out

    return output
