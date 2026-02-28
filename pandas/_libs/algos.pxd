from numpy cimport (
    float64_t,
    int64_t,
)

from pandas._libs.dtypes cimport (
    numeric_object_t,
    numeric_t,
)


cdef numeric_t kth_smallest_c(numeric_t* arr, Py_ssize_t k, Py_ssize_t n) noexcept nogil

cdef void moments_add_value(
    float64_t val,
    int64_t* nobs,
    float64_t* mean,
    float64_t* m2,
    float64_t* m3,
    float64_t* m4,
    int max_moment,
) noexcept nogil

cdef float64_t calc_skew(int64_t nobs, float64_t m2, float64_t m3) noexcept nogil

cdef float64_t calc_kurt(int64_t nobs, float64_t m2, float64_t m4) noexcept nogil

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
