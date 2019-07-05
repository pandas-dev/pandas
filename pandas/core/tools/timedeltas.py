"""
timedelta support tools
"""

import warnings

import numpy as np

from pandas._libs.tslibs import NaT
from pandas._libs.tslibs.timedeltas import Timedelta, parse_timedelta_unit
from pandas.util._decorators import deprecate_kwarg

from pandas.core.dtypes.common import is_list_like
from pandas.core.dtypes.generic import ABCIndexClass, ABCSeries

from pandas.core.arrays.timedeltas import sequence_to_td64ns


@deprecate_kwarg(old_arg_name="box", new_arg_name=None)
def to_timedelta(arg, unit="ns", box=True, errors="raise"):
    """
    Convert argument to timedelta.

    Timedeltas are absolute differences in times, expressed in difference
    units (e.g. days, hours, minutes, seconds). This method converts
    an argument from a recognized timedelta format / value into
    a Timedelta type.

    Parameters
    ----------
    arg : str, timedelta, list-like or Series
        The data to be converted to timedelta.
    unit : str, default 'ns'
        Denotes the unit of the arg. Possible values:
        ('Y', 'M', 'W', 'D', 'days', 'day', 'hours', hour', 'hr',
        'h', 'm', 'minute', 'min', 'minutes', 'T', 'S', 'seconds',
        'sec', 'second', 'ms', 'milliseconds', 'millisecond',
        'milli', 'millis', 'L', 'us', 'microseconds', 'microsecond',
        'micro', 'micros', 'U', 'ns', 'nanoseconds', 'nano', 'nanos',
        'nanosecond', 'N').
    box : bool, default True
        - If True returns a Timedelta/TimedeltaIndex of the results.
        - If False returns a numpy.timedelta64 or numpy.darray of
          values of dtype timedelta64[ns].

        .. deprecated:: 0.25.0
            Use :meth:`Series.to_numpy` or :meth:`Timedelta.to_timedelta64`
            instead to get an ndarray of values or numpy.timedelta64,
            respectively.

    errors : {'ignore', 'raise', 'coerce'}, default 'raise'
        - If 'raise', then invalid parsing will raise an exception.
        - If 'coerce', then invalid parsing will be set as NaT.
        - If 'ignore', then invalid parsing will return the input.

    Returns
    -------
    timedelta64 or numpy.array of timedelta64
        Output type returned if parsing succeeded.

    See Also
    --------
    DataFrame.astype : Cast argument to a specified dtype.
    to_datetime : Convert argument to datetime.

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

    Returning an ndarray by using the 'box' keyword argument:

    >>> pd.to_timedelta(np.arange(5), box=False)
    array([0, 1, 2, 3, 4], dtype='timedelta64[ns]')
    """
    unit = parse_timedelta_unit(unit)

    if errors not in ("ignore", "raise", "coerce"):
        raise ValueError("errors must be one of 'ignore', " "'raise', or 'coerce'}")

    if unit in {"Y", "y", "M"}:
        warnings.warn(
            "M and Y units are deprecated and " "will be removed in a future version.",
            FutureWarning,
            stacklevel=2,
        )

    if arg is None:
        return arg
    elif isinstance(arg, ABCSeries):
        values = _convert_listlike(arg._values, unit=unit, box=False, errors=errors)
        return arg._constructor(values, index=arg.index, name=arg.name)
    elif isinstance(arg, ABCIndexClass):
        return _convert_listlike(arg, unit=unit, box=box, errors=errors, name=arg.name)
    elif isinstance(arg, np.ndarray) and arg.ndim == 0:
        # extract array scalar and process below
        arg = arg.item()
    elif is_list_like(arg) and getattr(arg, "ndim", 1) == 1:
        return _convert_listlike(arg, unit=unit, box=box, errors=errors)
    elif getattr(arg, "ndim", 1) > 1:
        raise TypeError(
            "arg must be a string, timedelta, list, tuple, " "1-d array, or Series"
        )

    # ...so it must be a scalar value. Return scalar.
    return _coerce_scalar_to_timedelta_type(arg, unit=unit, box=box, errors=errors)


def _coerce_scalar_to_timedelta_type(r, unit="ns", box=True, errors="raise"):
    """Convert string 'r' to a timedelta object."""

    try:
        result = Timedelta(r, unit)
        if not box:
            # explicitly view as timedelta64 for case when result is pd.NaT
            result = result.asm8.view("timedelta64[ns]")
    except ValueError:
        if errors == "raise":
            raise
        elif errors == "ignore":
            return r

        # coerce
        result = NaT

    return result


def _convert_listlike(arg, unit="ns", box=True, errors="raise", name=None):
    """Convert a list of objects to a timedelta index object."""

    if isinstance(arg, (list, tuple)) or not hasattr(arg, "dtype"):
        # This is needed only to ensure that in the case where we end up
        #  returning arg (errors == "ignore"), and where the input is a
        #  generator, we return a useful list-like instead of a
        #  used-up generator
        arg = np.array(list(arg), dtype=object)

    try:
        value = sequence_to_td64ns(arg, unit=unit, errors=errors, copy=False)[0]
    except ValueError:
        if errors == "ignore":
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

        value = TimedeltaIndex(value, unit="ns", name=name)
    return value
