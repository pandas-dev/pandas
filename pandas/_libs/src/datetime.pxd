# cython: profile=False
from numpy cimport int64_t, npy_int64, npy_int32

from cpython cimport PyUnicode_Check, PyUnicode_AsASCIIString


cdef extern from "numpy/ndarrayobject.h":
    ctypedef int64_t npy_datetime

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

    void pandas_datetime_to_datetimestruct(npy_datetime val,
                                           PANDAS_DATETIMEUNIT fr,
                                           pandas_datetimestruct *result) nogil


cdef extern from "datetime/np_datetime_strings.h":
    int parse_iso_8601_datetime(char *str, int len, PANDAS_DATETIMEUNIT unit,
                                pandas_datetimestruct *out,
                                int *out_local, int *out_tzoffset,
                                PANDAS_DATETIMEUNIT *out_bestunit,
                                npy_bool *out_special)

cdef inline int _string_to_dts(object val, pandas_datetimestruct* dts,
                           int* out_local, int* out_tzoffset) except? -1:
    cdef int result
    cdef char *tmp

    if PyUnicode_Check(val):
        val = PyUnicode_AsASCIIString(val);

    tmp = val
    result = _cstring_to_dts(tmp, len(val), dts, out_local, out_tzoffset)

    if result == -1:
        raise ValueError('Unable to parse %s' % str(val))
    return result

cdef inline int _cstring_to_dts(char *val, int length,
                                pandas_datetimestruct* dts,
                                int* out_local, int* out_tzoffset) except? -1:
    cdef:
        npy_bool special
        PANDAS_DATETIMEUNIT out_bestunit
        int result

    result = parse_iso_8601_datetime(val, length, PANDAS_FR_ns,
                                     dts, out_local, out_tzoffset,
                                     &out_bestunit, &special)
    return result
