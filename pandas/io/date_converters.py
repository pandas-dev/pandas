"""This module is designed for community supported date conversion functions"""
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

    for i in xrange(N):
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
