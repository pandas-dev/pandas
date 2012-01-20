from numpy cimport int64_t

cdef extern from "datetime.h":

    ctypedef class datetime.datetime [object PyDateTime_DateTime]:
        # cdef int *data
        # cdef long hashcode
        # cdef char hastzinfo
        pass

    int PyDateTime_GET_YEAR(datetime o)
    int PyDateTime_GET_MONTH(datetime o)
    int PyDateTime_GET_DAY(datetime o)
    int PyDateTime_DATE_GET_HOUR(datetime o)
    int PyDateTime_DATE_GET_MINUTE(datetime o)
    int PyDateTime_DATE_GET_SECOND(datetime o)
    int PyDateTime_DATE_GET_MICROSECOND(datetime o)
    int PyDateTime_TIME_GET_HOUR(datetime o)
    int PyDateTime_TIME_GET_MINUTE(datetime o)
    int PyDateTime_TIME_GET_SECOND(datetime o)
    int PyDateTime_TIME_GET_MICROSECOND(datetime o)
    bint PyDateTime_Check(object o)
    void PyDateTime_IMPORT()

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

    npy_datetime PyArray_DatetimeStructToDatetime(NPY_DATETIMEUNIT fr,
                                                  npy_datetimestruct *d)

    void PyArray_DatetimeToDatetimeStruct(npy_datetime val,
                                          NPY_DATETIMEUNIT fr,
                                          npy_datetimestruct *result)
