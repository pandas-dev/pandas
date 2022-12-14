from pandas._libs.tslibs.np_datetime cimport npy_datetimestruct


cdef strptime(
    val,
    str fmt,
    bint exact,
    format_regex,
    locale_time,
    npy_datetimestruct dts,
)
