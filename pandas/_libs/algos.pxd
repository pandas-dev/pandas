cimport cython
from libc.math cimport (
    NAN,
    sqrt,
)
from numpy cimport (
    float64_t,
    int64_t,
    uint8_t,
)

from pandas._libs.dtypes cimport (
    numeric_object_t,
    numeric_t,
)


cdef numeric_t kth_smallest_c(numeric_t* arr, Py_ssize_t k, Py_ssize_t n) noexcept nogil


cdef extern from "pandas/moments.h":
    ctypedef struct Moments:
        int64_t n
        float64_t mean
        float64_t m2
        float64_t m3
        float64_t m4

    void moments_add_value(
        Moments* moments,
        float64_t val,
        int max_moment,
    ) nogil

    Moments moments_merge(
        Moments a,
        Moments b,
        int max_moment,
    ) nogil

    Moments accumulate_moments_scalar(
        const float64_t* values,
        int64_t n,
        bint skipna,
        const uint8_t* mask,
        int max_moment,
    ) nogil

    float64_t calc_skew(Moments) nogil

    float64_t calc_kurt(Moments) nogil


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
