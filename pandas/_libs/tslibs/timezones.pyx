# -*- coding: utf-8 -*-
# cython: profile=False

from pandas.compat import string_types, is_platform_windows

cdef extern from "Python.h":
    Py_ssize_t PY_SSIZE_T_MAX

cdef extern from "datetime.h":
    void PyDateTime_IMPORT()

# import datetime C API
PyDateTime_IMPORT

import numpy as np
cimport numpy as np
from numpy cimport ndarray, int64_t, float64_t
np.import_array()


# dateutil compat
from dateutil.tz import (tzoffset,
                         tzlocal as _dateutil_tzlocal,
                         tzfile as _dateutil_tzfile,
                         tzutc as _dateutil_tzutc,
                         tzstr as _dateutil_tzstr)

if is_platform_windows():
    from dateutil.zoneinfo import gettz as _dateutil_gettz
else:
    from dateutil.tz import gettz as _dateutil_gettz


from pytz.tzinfo import BaseTzInfo as _pytz_BaseTzInfo
import pytz
UTC = pytz.utc


cdef inline bint _is_utc(object tz):
    return tz is UTC or isinstance(tz, _dateutil_tzutc)
