# cython: profile=False

cimport numpy as np
import numpy as np

cimport cython
from cython cimport Py_ssize_t

np.import_array()

from numpy cimport (ndarray,
                    int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t,
                    uint32_t, uint64_t, float32_t, float64_t)

cdef double NaN = <double> np.NaN
cdef double nan = NaN

include "reshape_helper.pxi"
