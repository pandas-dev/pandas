# -*- coding: utf-8 -*-
# cython: profile=False

# dateutil compat
from dateutil.tz import tzutc as _dateutil_tzutc

import pytz
UTC = pytz.utc


cdef inline bint _is_utc(object tz):
    return tz is UTC or isinstance(tz, _dateutil_tzutc)
