# -*- coding: utf-8 -*-

from numpy cimport int64_t
from pandas._libs.tslibs.np_datetime cimport npy_datetimestruct

cdef object create_timestamp_from_ts(int64_t value,
                                     npy_datetimestruct dts,
                                     object tz, object freq)
