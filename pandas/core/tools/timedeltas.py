"""
timedelta support tools
"""

import numpy as np
import pandas as pd
from pandas._libs import tslibs
from pandas._libs.tslibs.timedeltas import (convert_to_timedelta64,
                                            parse_timedelta_unit)

from pandas.core.dtypes.common import is_list_like
from pandas.core.dtypes.generic import ABCSeries, ABCIndexClass

from pandas.core.arrays.timedeltas import sequence_to_td64ns


def to_timedelta(arg, unit='ns', box=True, errors='raise'):
    """
    Convert argument to timedelta

    Parameters
    ----------
    arg : string, timedelta, list, tuple, 1-d array, or Series
    unit : str, optional
        Denote the unit of the input, if input is an integer. Default 'ns'.
        Possible values:
        {'Y', 'M', 'W', 'D', 'days', 'day', 'hours', hour', 'hr', 'h',
        'm', 'minute', 'min', 'minutes', 'T', 'S', 'seconds', 'sec', 'second',
        'ms', 'milliseconds', 'millisecond', 'milli', 'millis', 'L',
        'us', 'microseconds', 'microsecond', 'micro', 'micros', 'U',
        'ns', 'nanoseconds', 'nano', 'nanos', 'nanosecond', 'N'}
    box : boolean, default True
        - If True returns a Timedelta/TimedeltaIndex of the results
        - if False returns a np.timedelta64 or ndarray of values of dtype
          timedelta64[ns]
    errors : {'ignore', 'raise', 'coerce'}, default 'raise'
        - If 'raise', then invalid parsing will raise an exception
        - If 'coerce', then invalid parsing will be set as NaT
        - If 'ignore', then invalid parsing will return the input

    Returns
    -------
    ret : timedelta64/arrays of timedelta64 if parsing succeeded

    Examples
    --------

    Parsing a single string to a Timedelta:

    >>> pd.to_timedelta('1 days 06:05:01.00003')
    Timedelta('1 days 06:05:01.000030')
    >>> pd.to_timedelta('15.5us')
    Timedelta('0 days 00:00:00.000015')

    Parsing a list or array of strings:

    >>> pd.to_timedelta(['1 days 06:05:01.00003', '15.5us', 'nan'])
    TimedeltaIndex(['1 days 06:05:01.000030', '0 days 00:00:00.000015', NaT],
                   dtype='timedelta64[ns]', freq=None)

    Converting numbers by specifying the `unit` keyword argument:

    >>> pd.to_timedelta(np.arange(5), unit='s')
    TimedeltaIndex(['00:00:00', '00:00:01', '00:00:02',
                    '00:00:03', '00:00:04'],
                   dtype='timedelta64[ns]', freq=None)
    >>> pd.to_timedelta(np.arange(5), unit='d')
    TimedeltaIndex(['0 days', '1 days', '2 days', '3 days', '4 days'],
                   dtype='timedelta64[ns]', freq=None)

    See Also
    --------
    pandas.DataFrame.astype : Cast argument to a specified dtype.
    pandas.to_datetime : Convert argument to datetime.
    """
    unit = parse_timedelta_unit(unit)

    if errors not in ('ignore', 'raise', 'coerce'):
        raise ValueError("errors must be one of 'ignore', "
                         "'raise', or 'coerce'}")

    if arg is None:
        return arg
    elif isinstance(arg, ABCSeries):
        from pandas import Series
        values = _convert_listlike(arg._values, unit=unit,
                                   box=False, errors=errors)
        return Series(values, index=arg.index, name=arg.name)
    elif isinstance(arg, ABCIndexClass):
        return _convert_listlike(arg, unit=unit, box=box,
                                 errors=errors, name=arg.name)
    elif isinstance(arg, np.ndarray) and arg.ndim == 0:
        # extract array scalar and process below
        arg = arg.item()
    elif is_list_like(arg) and getattr(arg, 'ndim', 1) == 1:
        return _convert_listlike(arg, unit=unit, box=box, errors=errors)
    elif getattr(arg, 'ndim', 1) > 1:
        raise TypeError('arg must be a string, timedelta, list, tuple, '
                        '1-d array, or Series')

    # ...so it must be a scalar value. Return scalar.
    return _coerce_scalar_to_timedelta_type(arg, unit=unit,
                                            box=box, errors=errors)


def _coerce_scalar_to_timedelta_type(r, unit='ns', box=True, errors='raise'):
    """Convert string 'r' to a timedelta object."""

    try:
        result = convert_to_timedelta64(r, unit)
    except ValueError:
        if errors == 'raise':
            raise
        elif errors == 'ignore':
            return r

        # coerce
        result = pd.NaT

    if box:
        result = tslibs.Timedelta(result)
    return result


def _convert_listlike(arg, unit='ns', box=True, errors='raise', name=None):
    """Convert a list of objects to a timedelta index object."""

    if isinstance(arg, (list, tuple)) or not hasattr(arg, 'dtype'):
        # This is needed only to ensure that in the case where we end up
        #  returning arg (errors == "ignore"), and where the input is a
        #  generator, we return a useful list-like instead of a
        #  used-up generator
        arg = np.array(list(arg), dtype=object)

    try:
        value = sequence_to_td64ns(arg, unit=unit,
                                   errors=errors, copy=False)[0]
    except ValueError:
        if errors == 'ignore':
            return arg
        else:
            # This else-block accounts for the cases when errors='raise'
            # and errors='coerce'. If errors == 'raise', these errors
            # should be raised. If errors == 'coerce', we shouldn't
            # expect any errors to be raised, since all parsing errors
            # cause coercion to pd.NaT. However, if an error / bug is
            # introduced that causes an Exception to be raised, we would
            # like to surface it.
            raise

    if box:
        from pandas import TimedeltaIndex
        value = TimedeltaIndex(value, unit='ns', name=name)
    return value
