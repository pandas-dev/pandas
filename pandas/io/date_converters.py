"""This module is designed for community supported date conversion functions"""
from datetime import datetime, timedelta, time

from pandas.compat import range
import numpy as np
import pandas.lib as lib


def parse_date_time(date_col, time_col):
    date_col = _maybe_cast(date_col)
    time_col = _maybe_cast(time_col)
    return lib.try_parse_date_and_time(date_col, time_col)


def parse_date_fields(year_col, month_col, day_col):
    year_col = _maybe_cast(year_col)
    month_col = _maybe_cast(month_col)
    day_col = _maybe_cast(day_col)
    return lib.try_parse_year_month_day(year_col, month_col, day_col)


def parse_all_fields(year_col, month_col, day_col, hour_col, minute_col,
                     second_col):
    year_col = _maybe_cast(year_col)
    month_col = _maybe_cast(month_col)
    day_col = _maybe_cast(day_col)
    hour_col = _maybe_cast(hour_col)
    minute_col = _maybe_cast(minute_col)
    second_col = _maybe_cast(second_col)
    return lib.try_parse_datetime_components(year_col, month_col, day_col,
                    hour_col, minute_col, second_col)


def generic_parser(parse_func, *cols):
    N = _check_columns(cols)
    results = np.empty(N, dtype=object)

    for i in range(N):
        args = [c[i] for c in cols]
        results[i] = parse_func(*args)

    return results


def _maybe_cast(arr):
    if not arr.dtype.type == np.object_:
        arr = np.array(arr, dtype=object)
    return arr


def _check_columns(cols):
    if not ((len(cols) > 0)):
        raise AssertionError()

    N = len(cols[0])
    for c in cols[1:]:
        if not ((len(c) == N)):
            raise AssertionError()

    return N


## Datetime Conversion for date_parsers
## see also: create a community supported set of typical converters
##           https://github.com/pydata/pandas/issues/1180

def offset_datetime(dt_in, days=0, hours=0, minutes=0,
                    seconds=0, microseconds=0):
    '''appply corrective time offset using datetime.timedelta

    input
    -----
    dt_in : datetime.time or datetime.datetime object
    days : integer value (positive or negative) for days component of offset
    hours : integer value (positive or negative) for hours component of offset
    minutes : integer value (positive or negative) for
              minutes component of offset
    seconds : integer value (positive or negative) for
              seconds component of offset
    microseconds : integer value (positive or negative) for
                   microseconds component of offset

    output
    ------
    ti_corr : datetime.time or datetime.datetime object


    '''
    # if a excel time like '23.07.2013 24:00' they actually mean
    # in Python '23.07.2013 23:59', must be converted
#            offset = -10 # minutes
    delta = timedelta(days=days, hours=hours, minutes=minutes,
                      seconds=seconds, microseconds=microseconds)

    #check if offset it to me applied on datetime or time
    if type(dt_in) is time:
        #create psydo datetime
        dt_now = datetime.now()
        dt_base = datetime.combine(dt_now, dt_in)
    else:
        dt_base = dt_in

    dt_corr = (dt_base) + delta

    #if input is time, we return it.
    if type(dt_in) is time:
        dt_corr = dt_corr.time()

    return dt_corr


def dt2ti(dt_in):
    '''converts wrong datetime.datetime to datetime.time

    input
    -----
    dt_in : dt_in : datetime.time or datetime.datetime object

    output
    -------
    ti_corr : datetime.time object
    '''
    # so we correct those which are not of type :mod:datetime.time
    # impdt2tiortant hint:
    # http://stackoverflow.com/a/12906456
    if type(dt_in) is not time:
        dt_in = dt_in.time()
    elif type(dt_in) is datetime:
        dt_in = dt_in.time()
    else:
        pass

    return dt_in
