# -*- coding: utf-8 -*-
# cython: profile=False

cimport numpy as cnp
from numpy cimport int64_t, ndarray, float64_t
import numpy as np
cnp.import_array()


from cpython cimport PyFloat_Check, PyUnicode_Check

from util cimport (is_integer_object, is_float_object, is_string_object,
                   is_datetime64_object)

from cpython.datetime cimport (PyDateTime_Check, PyDate_Check,
                               PyDateTime_CheckExact,
                               PyDateTime_IMPORT,
                               timedelta, datetime, date, time)
# import datetime C API
PyDateTime_IMPORT


from tslibs.np_datetime cimport (check_dts_bounds,
                                 pandas_datetimestruct,
                                 _string_to_dts,
                                 dt64_to_dtstruct, dtstruct_to_dt64,
                                 pydatetime_to_dt64, pydate_to_dt64,
                                 get_datetime64_value)
from tslibs.np_datetime import OutOfBoundsDatetime

from tslibs.parsing import parse_datetime_string

cimport cython
from cython cimport Py_ssize_t


import pytz
UTC = pytz.utc

from tslibs.arraylike import (  # noqa:F841
    format_array_from_datetime, array_to_datetime, array_with_unit_to_datetime,
    ints_to_pytimedelta, ints_to_pydatetime)

from tslibs.timedeltas cimport cast_from_unit
from tslibs.timedeltas import Timedelta
from tslibs.timezones cimport (is_utc, is_tzlocal, is_fixed_offset,
                               treat_tz_as_pytz, get_dst_info)
from tslibs.conversion cimport (tz_convert_single, _TSObject,
                                convert_datetime_to_tsobject,
                                get_datetime64_nanos,
                                tz_convert_utc_to_tzlocal)
from tslibs.conversion import tz_convert_single

from tslibs.nattype import NaT, nat_strings, iNaT
from tslibs.nattype cimport checknull_with_nat, NPY_NAT

from tslibs.offsets cimport to_offset

from tslibs.timestamps cimport (create_timestamp_from_ts,
                                _NS_UPPER_BOUND, _NS_LOWER_BOUND)
from tslibs.timestamps import Timestamp


def _test_parse_iso8601(object ts):
    """
    TESTING ONLY: Parse string into Timestamp using iso8601 parser. Used
    only for testing, actual construction uses `convert_str_to_tsobject`
    """
    cdef:
        _TSObject obj
        int out_local = 0, out_tzoffset = 0

    obj = _TSObject()

    if ts == 'now':
        return Timestamp.utcnow()
    elif ts == 'today':
        return Timestamp.now().normalize()

    _string_to_dts(ts, &obj.dts, &out_local, &out_tzoffset)
    obj.value = dtstruct_to_dt64(&obj.dts)
    check_dts_bounds(&obj.dts)
    if out_local == 1:
        obj.tzinfo = pytz.FixedOffset(out_tzoffset)
        obj.value = tz_convert_single(obj.value, obj.tzinfo, 'UTC')
        return Timestamp(obj.value, tz=obj.tzinfo)
    else:
        return Timestamp(obj.value)


cpdef inline object _localize_pydatetime(object dt, object tz):
    """
    Take a datetime/Timestamp in UTC and localizes to timezone tz.
    """
    if tz is None:
        return dt
    elif isinstance(dt, Timestamp):
        return dt.tz_localize(tz)
    elif tz == 'UTC' or tz is UTC:
        return UTC.localize(dt)
    try:
        # datetime.replace with pytz may be incorrect result
        return tz.localize(dt)
    except AttributeError:
        return dt.replace(tzinfo=tz)


# ----------------------------------------------------------------------
# Some general helper functions


cpdef normalize_date(object dt):
    """
    Normalize datetime.datetime value to midnight. Returns datetime.date as a
    datetime.datetime at midnight

    Returns
    -------
    normalized : datetime.datetime or Timestamp
    """
    if PyDateTime_Check(dt):
        if not PyDateTime_CheckExact(dt):
            # i.e. a Timestamp object
            return dt.replace(hour=0, minute=0, second=0, microsecond=0,
                              nanosecond=0)
        else:
            # regular datetime object
            return dt.replace(hour=0, minute=0, second=0, microsecond=0)
            # TODO: Make sure DST crossing is handled correctly here
    elif PyDate_Check(dt):
        return datetime(dt.year, dt.month, dt.day)
    else:
        raise TypeError('Unrecognized type: %s' % type(dt))
