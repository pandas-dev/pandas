from numpy cimport int64_t
from cpython cimport PyObject

cdef extern from "datetime.h":

    ctypedef class datetime.date [object PyDateTime_Date]:
        pass

    ctypedef class datetime.datetime [object PyDateTime_DateTime]:
        pass

    ctypedef class datetime.timedelta [object PyDateTime_Delta]:
        pass

    void PyDateTime_IMPORT()

    int PyDateTime_GET_YEAR(date)
    int PyDateTime_GET_MONTH(date)
    int PyDateTime_GET_DAY(date)
    int PyDateTime_DATE_GET_HOUR(object o)
    int PyDateTime_DATE_GET_MINUTE(object o)
    int PyDateTime_DATE_GET_SECOND(object o)
    int PyDateTime_DATE_GET_MICROSECOND(object o)
    int PyDateTime_TIME_GET_HOUR(object o)
    int PyDateTime_TIME_GET_MINUTE(object o)
    int PyDateTime_TIME_GET_SECOND(object o)
    int PyDateTime_TIME_GET_MICROSECOND(object o)
    bint PyDateTime_Check(object o)
    bint PyDate_Check(object o)
    object PyDateTime_FromDateAndTime(int year, int month, int day, int hour,
                                      int minute, int second, int us)

cdef extern from "numpy/ndarrayobject.h":

    ctypedef int64_t npy_timedelta
    ctypedef int64_t npy_datetime

    ctypedef struct npy_datetimestruct:
        int64_t year
        int month, day, hour, min, sec, us, ps, as

    ctypedef enum NPY_DATETIMEUNIT:
        #NPY_FR_Y
        #NPY_FR_M
        #NPY_FR_W
        #NPY_FR_B
        #NPY_FR_D
        #NPY_FR_h
        #NPY_FR_m
        #NPY_FR_s
        #NPY_FR_ms
        NPY_FR_us
        #NPY_FR_ns
        #NPY_FR_ps
        #NPY_FR_fs
        #NPY_FR_as

    ctypedef enum NPY_CASTING:
            NPY_NO_CASTING
            NPY_EQUIV_CASTING
            NPY_SAFE_CASTING
            NPY_SAME_KIND_CASTING
            NPY_UNSAFE_CASTING

    npy_datetime PyArray_DatetimeStructToDatetime(NPY_DATETIMEUNIT fr,
                                                  npy_datetimestruct *d)

    void PyArray_DatetimeToDatetimeStruct(npy_datetime val,
                                          NPY_DATETIMEUNIT fr,
                                          npy_datetimestruct *result)

cdef extern from "numpy/npy_common.h":

    ctypedef unsigned char npy_bool

cdef extern from "np_datetime.h":

    int convert_pydatetime_to_datetimestruct(PyObject *obj, npy_datetimestruct *out,
                                             NPY_DATETIMEUNIT *out_bestunit,
                                             int apply_tzinfo)
    int _days_per_month_table[2][12]

    int dayofweek(int y, int m, int d)
    int is_leapyear(int64_t year)

cdef extern from "np_datetime_strings.h":

    int parse_iso_8601_datetime(char *str, int len, NPY_DATETIMEUNIT unit,
                                NPY_CASTING casting, npy_datetimestruct *out,
                                npy_bool *out_local, NPY_DATETIMEUNIT *out_bestunit,
                                npy_bool *out_special)

    int make_iso_8601_datetime(npy_datetimestruct *dts, char *outstr, int outlen,
                               int local, NPY_DATETIMEUNIT base, int tzoffset,
                               NPY_CASTING casting)

    int get_datetime_iso_8601_strlen(int local, NPY_DATETIMEUNIT base)

cdef extern from "period.h":
    ctypedef struct date_info:
        long absdate
        double abstime
        double second
        int minute
        int hour
        int day
        int month
        int quarter
        int year
        int day_of_week
        int day_of_year
        int calendar

    long asfreq(long dtordinal, int freq1, int freq2, char relation) except -1

    long get_period_ordinal(int year, int month, int day,
                          int hour, int minute, int second,
                          int freq) except -1

    long get_python_ordinal(long period_ordinal, int freq) except -1

    char *skts_strftime(long value, int freq, PyObject *args)
    char *period_to_string(long value, int freq)
    char *period_to_string2(long value, int freq, char *fmt)

    int get_date_info(long ordinal, int freq, date_info *dinfo) except -1

    int pyear(long ordinal, int freq) except -1
    int pqyear(long ordinal, int freq) except -1
    int pquarter(long ordinal, int freq) except -1
    int pmonth(long ordinal, int freq) except -1
    int pday(long ordinal, int freq) except -1
    int pweekday(long ordinal, int freq) except -1
    int pday_of_week(long ordinal, int freq) except -1
    int pday_of_year(long ordinal, int freq) except -1
    int pweek(long ordinal, int freq) except -1
    int phour(long ordinal, int freq) except -1
    int pminute(long ordinal, int freq) except -1
    int psecond(long ordinal, int freq) except -1

cdef extern from "stdint.h":
    enum: INT64_MIN

