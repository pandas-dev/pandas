from cpython.datetime cimport tzinfo

import numpy as np
from numpy cimport ndarray, int64_t, int32_t

from pandas._libs.tslibs.util cimport get_nat

from pandas._libs.tslibs.dtypes import Resolution
from pandas._libs.tslibs.np_datetime cimport (
    npy_datetimestruct, dt64_to_dtstruct)
from pandas._libs.tslibs.timezones cimport get_tzconverter, TZConvertInfo
from pandas._libs.tslibs.ccalendar cimport get_days_in_month
from pandas._libs.tslibs.tzconversion cimport tz_convert_utc_to_tzlocal

# ----------------------------------------------------------------------
# Constants

cdef:
    int64_t NPY_NAT = get_nat()

    int RESO_NS = 0
    int RESO_US = 1
    int RESO_MS = 2
    int RESO_SEC = 3
    int RESO_MIN = 4
    int RESO_HR = 5
    int RESO_DAY = 6
    int RESO_MTH = 7
    int RESO_QTR = 8
    int RESO_YR = 9


# ----------------------------------------------------------------------


def get_resolution(const int64_t[:] stamps, tzinfo tz=None):
    cdef:
        Py_ssize_t i, n = len(stamps)
        npy_datetimestruct dts
        int reso = RESO_DAY, curr_reso
        int64_t local_val
        TZConvertInfo info

    info = get_tzconverter(tz, stamps)

    if info.use_fixed:
        assert info.delta != NPY_NAT
    elif not info.use_utc and not info.use_tzlocal:
        assert info.utcoffsets is not NULL
        assert info.positions is not NULL

    for i in range(n):
        if stamps[i] == NPY_NAT:
            continue

        if info.use_utc:
            local_val = stamps[i]
        elif info.use_tzlocal:
            local_val = tz_convert_utc_to_tzlocal(stamps[i], tz)
        elif info.use_fixed:
            local_val = stamps[i] + info.delta
        else:
            local_val = stamps[i] + info.utcoffsets[info.positions[i]]

        dt64_to_dtstruct(local_val, &dts)
        curr_reso = _reso_stamp(&dts)
        if curr_reso < reso:
            reso = curr_reso

    return Resolution(reso)


cdef inline int _reso_stamp(npy_datetimestruct *dts):
    if dts.us != 0:
        if dts.us % 1000 == 0:
            return RESO_MS
        return RESO_US
    elif dts.sec != 0:
        return RESO_SEC
    elif dts.min != 0:
        return RESO_MIN
    elif dts.hour != 0:
        return RESO_HR
    return RESO_DAY


# ----------------------------------------------------------------------
# Frequency Inference

def month_position_check(fields, weekdays):
    cdef:
        int32_t daysinmonth, y, m, d
        bint calendar_end = True
        bint business_end = True
        bint calendar_start = True
        bint business_start = True
        bint cal
        int32_t[:] years
        int32_t[:] months
        int32_t[:] days

    years = fields['Y']
    months = fields['M']
    days = fields['D']

    for y, m, d, wd in zip(years, months, days, weekdays):
        if calendar_start:
            calendar_start &= d == 1
        if business_start:
            business_start &= d == 1 or (d <= 3 and wd == 0)

        if calendar_end or business_end:
            daysinmonth = get_days_in_month(y, m)
            cal = d == daysinmonth
            if calendar_end:
                calendar_end &= cal
            if business_end:
                business_end &= cal or (daysinmonth - d < 3 and wd == 4)
        elif not calendar_start and not business_start:
            break

    if calendar_end:
        return 'ce'
    elif business_end:
        return 'be'
    elif calendar_start:
        return 'cs'
    elif business_start:
        return 'bs'
    else:
        return None
