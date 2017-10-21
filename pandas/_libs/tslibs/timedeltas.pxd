# -*- coding: utf-8 -*-
# cython: profile=False

from numpy cimport int64_t

# Exposed for tslib, not intended for outside use.
cdef parse_timedelta_string(object ts)
cpdef int64_t cast_from_unit(object ts, object unit) except? -1
