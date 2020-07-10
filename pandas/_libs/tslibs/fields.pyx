"""
Functions for accessing attributes of Timestamp/datetime64/datetime-like
objects and arrays
"""
from locale import LC_TIME

import cython
from cython import Py_ssize_t

import numpy as np
cimport numpy as cnp
from numpy cimport ndarray, int64_t, int32_t, int8_t, uint32_t
cnp.import_array()

from pandas._config.localization import set_locale

from pandas._libs.tslibs.ccalendar import MONTHS_FULL, DAYS_FULL
from pandas._libs.tslibs.ccalendar cimport (
    get_days_in_month, is_leapyear, dayofweek, get_week_of_year,
    get_day_of_year, get_iso_calendar, iso_calendar_t,
    month_offset,
    get_firstbday,
    get_lastbday,
)
from pandas._libs.tslibs.np_datetime cimport (
    npy_datetimestruct, pandas_timedeltastruct, dt64_to_dtstruct,
    td64_to_tdstruct)
from pandas._libs.tslibs.nattype cimport NPY_NAT
from pandas._libs.tslibs.strptime import LocaleTime


@cython.wraparound(False)
@cython.boundscheck(False)
def build_field_sarray(const int64_t[:] dtindex):
    """
    Datetime as int64 representation to a structured array of fields
    """
    cdef:
        Py_ssize_t i, count = len(dtindex)
        npy_datetimestruct dts
        ndarray[int32_t] years, months, days, hours, minutes, seconds, mus

    sa_dtype = [
        ("Y", "i4"),  # year
        ("M", "i4"),  # month
        ("D", "i4"),  # day
        ("h", "i4"),  # hour
        ("m", "i4"),  # min
        ("s", "i4"),  # second
        ("u", "i4"),  # microsecond
    ]

    out = np.empty(count, dtype=sa_dtype)

    years = out['Y']
    months = out['M']
    days = out['D']
    hours = out['h']
    minutes = out['m']
    seconds = out['s']
    mus = out['u']

    for i in range(count):
        dt64_to_dtstruct(dtindex[i], &dts)
        years[i] = dts.year
        months[i] = dts.month
        days[i] = dts.day
        hours[i] = dts.hour
        minutes[i] = dts.min
        seconds[i] = dts.sec
        mus[i] = dts.us

    return out


def month_position_check(fields, weekdays):
    cdef:
        int32_t daysinmonth, y, m, d
        bint calendar_end = True
        bint business_end = True
        bint calendar_start = True
        bint business_start = True
        bint cal
        int32_t[:] years = fields["Y"]
        int32_t[:] months = fields["M"]
        int32_t[:] days = fields["D"]

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
        return "ce"
    elif business_end:
        return "be"
    elif calendar_start:
        return "cs"
    elif business_start:
        return "bs"
    else:
        return None


@cython.wraparound(False)
@cython.boundscheck(False)
def get_date_name_field(const int64_t[:] dtindex, str field, object locale=None):
    """
    Given a int64-based datetime index, return array of strings of date
    name based on requested field (e.g. day_name)
    """
    cdef:
        Py_ssize_t i, count = len(dtindex)
        ndarray[object] out, names
        npy_datetimestruct dts
        int dow

    out = np.empty(count, dtype=object)

    if field == 'day_name':
        if locale is None:
            names = np.array(DAYS_FULL, dtype=np.object_)
        else:
            names = np.array(get_locale_names('f_weekday', locale),
                             dtype=np.object_)
        for i in range(count):
            if dtindex[i] == NPY_NAT:
                out[i] = np.nan
                continue

            dt64_to_dtstruct(dtindex[i], &dts)
            dow = dayofweek(dts.year, dts.month, dts.day)
            out[i] = names[dow].capitalize()

    elif field == 'month_name':
        if locale is None:
            names = np.array(MONTHS_FULL, dtype=np.object_)
        else:
            names = np.array(get_locale_names('f_month', locale),
                             dtype=np.object_)
        for i in range(count):
            if dtindex[i] == NPY_NAT:
                out[i] = np.nan
                continue

            dt64_to_dtstruct(dtindex[i], &dts)
            out[i] = names[dts.month].capitalize()

    else:
        raise ValueError(f"Field {field} not supported")

    return out


@cython.wraparound(False)
@cython.boundscheck(False)
def get_start_end_field(const int64_t[:] dtindex, str field,
                        object freqstr=None, int month_kw=12):
    """
    Given an int64-based datetime index return array of indicators
    of whether timestamps are at the start/end of the month/quarter/year
    (defined by frequency).
    """
    cdef:
        Py_ssize_t i
        int count = len(dtindex)
        bint is_business = 0
        int end_month = 12
        int start_month = 1
        ndarray[int8_t] out
        npy_datetimestruct dts

    out = np.zeros(count, dtype='int8')

    if freqstr:
        if freqstr == 'C':
            raise ValueError(f"Custom business days is not supported by {field}")
        is_business = freqstr[0] == 'B'

        # YearBegin(), BYearBegin() use month = starting month of year.
        # QuarterBegin(), BQuarterBegin() use startingMonth = starting
        # month of year. Other offsets use month, startingMonth as ending
        # month of year.

        if (freqstr[0:2] in ['MS', 'QS', 'AS']) or (
                freqstr[1:3] in ['MS', 'QS', 'AS']):
            end_month = 12 if month_kw == 1 else month_kw - 1
            start_month = month_kw
        else:
            end_month = month_kw
            start_month = (end_month % 12) + 1
    else:
        end_month = 12
        start_month = 1

    if field == 'is_month_start':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if dts.day == get_firstbday(dts.year, dts.month):
                    out[i] = 1

        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if dts.day == 1:
                    out[i] = 1

    elif field == 'is_month_end':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if dts.day == get_lastbday(dts.year, dts.month):
                    out[i] = 1

        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if dts.day == get_days_in_month(dts.year, dts.month):
                    out[i] = 1

    elif field == 'is_quarter_start':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if ((dts.month - start_month) % 3 == 0) and (
                        dts.day == get_firstbday(dts.year, dts.month)):
                    out[i] = 1

        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if ((dts.month - start_month) % 3 == 0) and dts.day == 1:
                    out[i] = 1

    elif field == 'is_quarter_end':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if ((dts.month - end_month) % 3 == 0) and (
                        dts.day == get_lastbday(dts.year, dts.month)):
                    out[i] = 1

        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if ((dts.month - end_month) % 3 == 0) and (
                        dts.day == get_days_in_month(dts.year, dts.month)):
                    out[i] = 1

    elif field == 'is_year_start':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if (dts.month == start_month) and (
                        dts.day == get_firstbday(dts.year, dts.month)):
                    out[i] = 1

        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if (dts.month == start_month) and dts.day == 1:
                    out[i] = 1

    elif field == 'is_year_end':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if (dts.month == end_month) and (
                        dts.day == get_lastbday(dts.year, dts.month)):
                    out[i] = 1

        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = 0
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)

                if (dts.month == end_month) and (
                        dts.day == get_days_in_month(dts.year, dts.month)):
                    out[i] = 1

    else:
        raise ValueError(f"Field {field} not supported")

    return out.view(bool)


@cython.wraparound(False)
@cython.boundscheck(False)
def get_date_field(const int64_t[:] dtindex, str field):
    """
    Given a int64-based datetime index, extract the year, month, etc.,
    field and return an array of these values.
    """
    cdef:
        Py_ssize_t i, count = len(dtindex)
        ndarray[int32_t] out
        npy_datetimestruct dts

    out = np.empty(count, dtype='i4')

    if field == 'Y':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = dts.year
        return out

    elif field == 'M':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = dts.month
        return out

    elif field == 'D':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = dts.day
        return out

    elif field == 'h':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = dts.hour
        return out

    elif field == 'm':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = dts.min
        return out

    elif field == 's':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = dts.sec
        return out

    elif field == 'us':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = dts.us
        return out

    elif field == 'ns':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = dts.ps // 1000
        return out
    elif field == 'doy':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = get_day_of_year(dts.year, dts.month, dts.day)
        return out

    elif field == 'dow':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = dayofweek(dts.year, dts.month, dts.day)
        return out

    elif field == 'woy':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = get_week_of_year(dts.year, dts.month, dts.day)
        return out

    elif field == 'q':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = dts.month
                out[i] = ((out[i] - 1) // 3) + 1
        return out

    elif field == 'dim':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                out[i] = get_days_in_month(dts.year, dts.month)
        return out
    elif field == 'is_leap_year':
        return isleapyear_arr(get_date_field(dtindex, 'Y'))

    raise ValueError(f"Field {field} not supported")


@cython.wraparound(False)
@cython.boundscheck(False)
def get_timedelta_field(const int64_t[:] tdindex, str field):
    """
    Given a int64-based timedelta index, extract the days, hrs, sec.,
    field and return an array of these values.
    """
    cdef:
        Py_ssize_t i, count = len(tdindex)
        ndarray[int32_t] out
        pandas_timedeltastruct tds

    out = np.empty(count, dtype='i4')

    if field == 'days':
        with nogil:
            for i in range(count):
                if tdindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                td64_to_tdstruct(tdindex[i], &tds)
                out[i] = tds.days
        return out

    elif field == 'h':
        with nogil:
            for i in range(count):
                if tdindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                td64_to_tdstruct(tdindex[i], &tds)
                out[i] = tds.hrs
        return out

    elif field == 's':
        with nogil:
            for i in range(count):
                if tdindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                td64_to_tdstruct(tdindex[i], &tds)
                out[i] = tds.sec
        return out

    elif field == 'seconds':
        with nogil:
            for i in range(count):
                if tdindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                td64_to_tdstruct(tdindex[i], &tds)
                out[i] = tds.seconds
        return out

    elif field == 'ms':
        with nogil:
            for i in range(count):
                if tdindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                td64_to_tdstruct(tdindex[i], &tds)
                out[i] = tds.ms
        return out

    elif field == 'microseconds':
        with nogil:
            for i in range(count):
                if tdindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                td64_to_tdstruct(tdindex[i], &tds)
                out[i] = tds.microseconds
        return out

    elif field == 'us':
        with nogil:
            for i in range(count):
                if tdindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                td64_to_tdstruct(tdindex[i], &tds)
                out[i] = tds.us
        return out

    elif field == 'ns':
        with nogil:
            for i in range(count):
                if tdindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                td64_to_tdstruct(tdindex[i], &tds)
                out[i] = tds.ns
        return out

    elif field == 'nanoseconds':
        with nogil:
            for i in range(count):
                if tdindex[i] == NPY_NAT:
                    out[i] = -1
                    continue

                td64_to_tdstruct(tdindex[i], &tds)
                out[i] = tds.nanoseconds
        return out

    raise ValueError(f"Field {field} not supported")


cpdef isleapyear_arr(ndarray years):
    """vectorized version of isleapyear; NaT evaluates as False"""
    cdef:
        ndarray[int8_t] out

    out = np.zeros(len(years), dtype='int8')
    out[np.logical_or(years % 400 == 0,
                      np.logical_and(years % 4 == 0,
                                     years % 100 > 0))] = 1
    return out.view(bool)


@cython.wraparound(False)
@cython.boundscheck(False)
def build_isocalendar_sarray(const int64_t[:] dtindex):
    """
    Given a int64-based datetime array, return the ISO 8601 year, week, and day
    as a structured array.
    """
    cdef:
        Py_ssize_t i, count = len(dtindex)
        npy_datetimestruct dts
        ndarray[uint32_t] iso_years, iso_weeks, days
        iso_calendar_t ret_val

    sa_dtype = [
        ("year", "u4"),
        ("week", "u4"),
        ("day", "u4"),
    ]

    out = np.empty(count, dtype=sa_dtype)

    iso_years = out["year"]
    iso_weeks = out["week"]
    days = out["day"]

    with nogil:
        for i in range(count):
            if dtindex[i] == NPY_NAT:
                ret_val = 0, 0, 0
            else:
                dt64_to_dtstruct(dtindex[i], &dts)
                ret_val = get_iso_calendar(dts.year, dts.month, dts.day)

            iso_years[i] = ret_val[0]
            iso_weeks[i] = ret_val[1]
            days[i] = ret_val[2]
    return out


def get_locale_names(name_type: str, locale: object = None):
    """
    Returns an array of localized day or month names.

    Parameters
    ----------
    name_type : string, attribute of LocaleTime() in which to return localized
        names
    locale : string

    Returns
    -------
    list of locale names
    """
    with set_locale(locale, LC_TIME):
        return getattr(LocaleTime(), name_type)
