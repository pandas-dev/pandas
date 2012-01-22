cimport numpy as cnp
cimport cython
cimport cpython
import numpy as np

from numpy cimport int64_t, import_array

# this is our datetime.pxd
from datetime cimport *

# import datetime C API
PyDateTime_IMPORT

# initialize numpy
import_array()

# in numpy 1.7, will prop need this
# numpy_pydatetime_import

cdef class Date:
    '''
    This is the custom pandas Date box for the numpy datetime64 dtype.
    '''
    cdef:
        int64_t timestamp
        object freq
        object tzinfo
        npy_datetimestruct dts

    def __init__(self, int64_t ts, object freq = None, object tzinfo = None):
        self.timestamp = ts
        self.freq = freq
        self.tzinfo = tzinfo

        # datetime64 decomposition to components
        PyArray_DatetimeToDatetimeStruct(self.timestamp, NPY_FR_us, &self.dts)

    # TODO: we'll probably need factory methods to construct this box from:
    #       -- datetime64 scalar python object
    #       -- datetime python object
    #       -- int64_t

    # --- the following properties to make it compatible with datetime

    property year:
        def __get__(self):
            return self.dts.year

    property month:
        def __get__(self):
            return self.dts.month

    property day:
        def __get__(self):
            return self.dts.day

    property hour:
        def __get__(self):
            return self.dts.hour

    property minute:
        def __get__(self):
            return self.dts.min

    property second:
        def __get__(self):
            return self.dts.sec

    property microsecond:
        def __get__(self):
            return self.dts.us

cdef:
    npy_datetimestruct g_dts
    NPY_DATETIMEUNIT g_out_bestunit

def pydt_to_dt64(object pydt):
    if PyDateTime_Check(pydt):
        convert_pydatetime_to_datetimestruct(<PyObject *>pydt, &g_dts, &g_out_bestunit, 1)
        return PyArray_DatetimeStructToDatetime(g_out_bestunit, &g_dts)

    raise ValueError("Expected a datetime, received a %s" % type(pydt))
