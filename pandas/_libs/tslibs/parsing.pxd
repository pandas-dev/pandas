
cpdef str get_rule_month(str source)
cpdef quarter_to_myear(int year, int quarter, str freq)
cpdef ISOInfo null_iso_info()
from pandas._libs.tslibs.np_datetime cimport (
    NPY_DATETIMEUNIT,
    npy_datetimestruct,
)


cdef extern from "src/datetime/np_datetime_strings.h":
    ctypedef struct ISOInfo:
        const char *format
        int format_len
        const char *date_sep
        const char *time_sep
        const char *micro_or_tz
        int year
        int month
        int day
        int hour
        int minute
        int second
        int exact
