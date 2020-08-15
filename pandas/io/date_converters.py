"""This module is designed for community supported date conversion functions"""
import warnings

import numpy as np

from pandas._libs.tslibs import parsing

from pandas import to_datetime


def parse_date_time(date_col, time_col):
    """Convert separate columns with dates and times into a single datetime column

    .. deprecated:: 1.1.0
       Use pd.to_datetime(date_col + " " + time_col) instead to get a Pandas Series.
       Use pd.to_datetime(date_col + " " + time_col).to_pydatetime() instead to get
        a Numpy array.
    """
    warnings.warn(
        """
        Use pd.to_datetime(date_col + " " + time_col) instead to get a Pandas Series.
        Use pd.to_datetime(date_col + " " + time_col).to_pydatetime() instead to get a Numpy array.
""",  # noqa: E501
        FutureWarning,
        stacklevel=2,
    )
    return to_datetime(date_col + " " + time_col).to_pydatetime()


def parse_date_fields(year_col, month_col, day_col):
    """Convert separate columns with years, months and days into a single date column

    .. deprecated:: 1.1.0
        Use pd.to_datetime({"year": year_col, "month": month_col, "day": day_col})
        instead to get a Pandas Series.
        Use ser = pd.to_datetime({"year": year_col, "month": month_col, "day": day_col})
        and np.array([s.to_pydatetime() for s in ser]) instead to get a Numpy array.
    """
    warnings.warn(
        """
        Use pd.to_datetime({"year": year_col, "month": month_col, "day": day_col}) instead to get a Pandas Series.
        Use ser = pd.to_datetime({"year": year_col, "month": month_col, "day": day_col}) and
        np.array([s.to_pydatetime() for s in ser]) instead to get a Numpy array.
""",  # noqa: E501
        FutureWarning,
        stacklevel=2,
    )

    ser = to_datetime({"year": year_col, "month": month_col, "day": day_col})
    return np.array([s.to_pydatetime() for s in ser])


def parse_all_fields(year_col, month_col, day_col, hour_col, minute_col, second_col):
    year_col = _maybe_cast(year_col)
    month_col = _maybe_cast(month_col)
    day_col = _maybe_cast(day_col)
    hour_col = _maybe_cast(hour_col)
    minute_col = _maybe_cast(minute_col)
    second_col = _maybe_cast(second_col)
    return parsing.try_parse_datetime_components(
        year_col, month_col, day_col, hour_col, minute_col, second_col
    )


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
    if not len(cols):
        raise AssertionError("There must be at least 1 column")

    head, tail = cols[0], cols[1:]

    N = len(head)

    for i, n in enumerate(map(len, tail)):
        if n != N:
            raise AssertionError(
                f"All columns must have the same length: {N}; column {i} has length {n}"
            )

    return N
