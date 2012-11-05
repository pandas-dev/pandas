from numpy cimport int64_t, int32_t, npy_int64, npy_int32
from cpython cimport PyObject


cdef extern from "stdint.h":
    enum: INT64_MIN
    enum: INT32_MIN



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
    bint PyTime_Check(object o)
    object PyDateTime_FromDateAndTime(int year, int month, int day, int hour,
                                      int minute, int second, int us)

cdef extern from "datetime_helper.h":
    void mangle_nat(object o)

cdef extern from "numpy/ndarrayobject.h":

    ctypedef int64_t npy_timedelta
    ctypedef int64_t npy_datetime

    ctypedef enum NPY_CASTING:
            NPY_NO_CASTING
            NPY_EQUIV_CASTING
            NPY_SAFE_CASTING
            NPY_SAME_KIND_CASTING
            NPY_UNSAFE_CASTING


cdef extern from "numpy_helper.h":
    npy_datetime get_datetime64_value(object o)

cdef extern from "numpy/npy_common.h":

    ctypedef unsigned char npy_bool

cdef extern from "datetime/np_datetime.h":

    ctypedef enum PANDAS_DATETIMEUNIT:
        PANDAS_FR_Y
        PANDAS_FR_M
        PANDAS_FR_W
        PANDAS_FR_D
        PANDAS_FR_B
        PANDAS_FR_h
        PANDAS_FR_m
        PANDAS_FR_s
        PANDAS_FR_ms
        PANDAS_FR_us
        PANDAS_FR_ns
        PANDAS_FR_ps
        PANDAS_FR_fs
        PANDAS_FR_as

    ctypedef struct pandas_datetimestruct:
        npy_int64 year
        npy_int32 month, day, hour, min, sec, us, ps, as

    int convert_pydatetime_to_datetimestruct(PyObject *obj,
                                             pandas_datetimestruct *out,
                                             PANDAS_DATETIMEUNIT *out_bestunit,
                                             int apply_tzinfo)

    npy_datetime pandas_datetimestruct_to_datetime(PANDAS_DATETIMEUNIT fr,
                                                   pandas_datetimestruct *d)
    void pandas_datetime_to_datetimestruct(npy_datetime val,
                                           PANDAS_DATETIMEUNIT fr,
                                           pandas_datetimestruct *result)
    int days_per_month_table[2][12]

    int dayofweek(int y, int m, int d)
    int is_leapyear(int64_t year)
    PANDAS_DATETIMEUNIT get_datetime64_unit(object o)

cdef extern from "datetime/np_datetime_strings.h":

    int parse_iso_8601_datetime(char *str, int len, PANDAS_DATETIMEUNIT unit,
                                NPY_CASTING casting, pandas_datetimestruct *out,
                                npy_bool *out_local, PANDAS_DATETIMEUNIT *out_bestunit,
                                npy_bool *out_special)

    int make_iso_8601_datetime(pandas_datetimestruct *dts, char *outstr, int outlen,
                               int local, PANDAS_DATETIMEUNIT base, int tzoffset,
                               NPY_CASTING casting)

    int get_datetime_iso_8601_strlen(int local, PANDAS_DATETIMEUNIT base)

    # int parse_python_string(object obj, pandas_datetimestruct *out) except -1
