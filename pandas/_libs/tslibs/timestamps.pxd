# -*- coding: utf-8 -*-
# cython: profile=False

from numpy cimport int64_t
from np_datetime cimport npy_datetimestruct

cdef object create_timestamp_from_ts(int64_t value,
                                     npy_datetimestruct dts,
                                     object tz, object freq)

cdef int64_t _NS_UPPER_BOUND, _NS_LOWER_BOUND
