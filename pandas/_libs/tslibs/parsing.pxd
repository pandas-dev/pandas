from cpython.datetime cimport datetime
from numpy cimport int64_t

from pandas._libs.tslibs.np_datetime cimport NPY_DATETIMEUNIT


cpdef str get_rule_month(str source)
cpdef quarter_to_myear(int year, int quarter, str freq)
cdef void warn_quarter_deprecated(str date_string, str freq)

# GH#65381 date strings parsed in a deprecated way
cdef enum DeprecatedParse:
    NO_DEPRECATION
    DATEUTIL_ONLY  # only the dateutil fallback parses it
    CONTRADICTS_DAYFIRST  # day/month order chosen against dayfirst
    YEAR_DAY_MONTH  # year-first string read as year-day-month by dayfirst

cdef void warn_parse_deprecated(
    str date_string, DeprecatedParse kind, bint dayfirst, bint* warned
)

cdef datetime parse_datetime_string(
    str date_string,
    bint dayfirst,
    bint yearfirst,
    NPY_DATETIMEUNIT* out_bestunit,
    int64_t* nanos,
    bint* out_is_quarter=*,
    DeprecatedParse* out_deprecated=*,
)
