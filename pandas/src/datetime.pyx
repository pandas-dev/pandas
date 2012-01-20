cimport numpy as cnp
cimport cython
cimport cpython
import numpy as np

from numpy cimport int64_t, import_array

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

# import datetime C API
PyDateTime_IMPORT

# initialize numpy
import_array()

cdef class Date:
    cdef:
        int64_t timestamp
        object freq
        object tzinfo

    def __init__(self, int64_t stamp, object freq = None,
                 object tzinfo = None):
        self.timestamp = stamp
        self.freq = <char *>freq
        self.tzinfo = tzinfo

_unbox_cache = dict()
def dt_unbox(object key):
    '''
    Unbox datetime to datetime64
    '''
    cdef int64_t y, M, d, h, m, s, u
    cdef int64_t val

    # NAT bit pattern
    val = 0x8000000000000000

    if PyDateTime_Check(key):
        u = PyDateTime_TIME_GET_MICROSECOND(key)
        s = PyDateTime_TIME_GET_SECOND(key)
        m = PyDateTime_TIME_GET_MINUTE(key)
        h = PyDateTime_TIME_GET_HOUR(key)
        d = PyDateTime_GET_DAY(key)
        M = PyDateTime_GET_MONTH(key)
        y = PyDateTime_GET_YEAR(key)
        val = y + M + d + h + m + s + u

    return np.datetime64(val)

_box_cache = dict()
def dt_box(object key):
    '''
    Box datetime64 to datetime
    '''
    return key.astype('O')

