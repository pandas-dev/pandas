cimport cython
from libc.math cimport (
    NAN,
    sqrt,
)
from numpy cimport (
    float64_t,
    int64_t,
)

from pandas._libs.dtypes cimport (
    numeric_object_t,
    numeric_t,
)


cdef numeric_t kth_smallest_c(numeric_t* arr, Py_ssize_t k, Py_ssize_t n) noexcept nogil


@cython.cdivision(True)
cdef inline void moments_add_value(
    float64_t val,
    int64_t* nobs,
    float64_t* mean,
    float64_t* m2,
    float64_t* m3,
    float64_t* m4,
    int max_moment,
) noexcept nogil:
    cdef:
        float64_t n
        float64_t delta = val - mean[0]
        float64_t delta_n, term1

    nobs[0] += 1
    n = <float64_t>nobs[0]
    delta_n = delta / n
    term1 = delta * delta_n * (n - 1.0)

    if max_moment >= 4:
        m4[0] += delta_n * (
                -4.0 * m3[0]
                + delta_n * (
                    6.0 * m2[0] + term1 * (n * n - 3.0 * n + 3.0)
                    )
                )
    if max_moment >= 3:
        m3[0] += delta_n * (term1 * (n - 2.0) - 3.0 * m2[0])
    m2[0] += term1
    mean[0] += delta_n


@cython.cdivision(True)
cdef inline float64_t calc_skew(
    int64_t nobs, float64_t m2, float64_t m3
) noexcept nogil:
    cdef:
        float64_t moments_ratio, correction, dnobs

    if nobs < 3:
        return NAN

    dnobs = <float64_t>nobs

    moments_ratio = m3 / (m2 * sqrt(m2))
    correction = (dnobs * sqrt(dnobs - 1.0)) / (dnobs - 2.0)
    return moments_ratio * correction


@cython.cdivision(True)
cdef inline float64_t calc_kurt(
    int64_t nobs, float64_t m2, float64_t m4
) noexcept nogil:
    cdef:
        float64_t result, dnobs, term1, term2, inner, correction
        float64_t moments_ratio

    if nobs < 4:
        return NAN
    dnobs = <float64_t>nobs
    moments_ratio = m4 / (m2 * m2)
    term1 = dnobs * (dnobs + 1.0) * moments_ratio
    term2 = 3.0 * (dnobs - 1.0)
    inner = term1 - term2

    correction = (dnobs - 1.0) / ((dnobs - 2.0) * (dnobs - 3.0))
    result = correction * inner

    return result


cdef enum TiebreakEnumType:
    TIEBREAK_AVERAGE
    TIEBREAK_MIN,
    TIEBREAK_MAX
    TIEBREAK_FIRST
    TIEBREAK_FIRST_DESCENDING
    TIEBREAK_DENSE


cdef numeric_object_t get_rank_nan_fill_val(
    bint rank_nans_highest,
    numeric_object_t val,
    bint is_datetimelike=*,
)
