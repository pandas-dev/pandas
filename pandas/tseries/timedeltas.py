"""
timedelta support tools
"""

import numpy as np
import pandas.tslib as tslib
from pandas.core.common import (ABCSeries, is_integer_dtype,
                                is_timedelta64_dtype, is_list_like,
                                _ensure_object, ABCIndexClass)
from pandas.util.decorators import deprecate_kwarg


@deprecate_kwarg(old_arg_name='coerce', new_arg_name='errors',
                 mapping={True: 'coerce', False: 'raise'})
def to_timedelta(arg, unit='ns', box=True, errors='raise', coerce=None):
    """
    Convert argument to timedelta

    Parameters
    ----------
    arg : string, timedelta, list, tuple, 1-d array, or Series
    unit : unit of the arg (D,h,m,s,ms,us,ns) denote the unit, which is an
        integer/float number
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
    """
    unit = _validate_timedelta_unit(unit)

    def _convert_listlike(arg, box, unit, name=None):

        if isinstance(arg, (list, tuple)) or not hasattr(arg, 'dtype'):
            arg = np.array(list(arg), dtype='O')

        # these are shortcutable
        if is_timedelta64_dtype(arg):
            value = arg.astype('timedelta64[ns]')
        elif is_integer_dtype(arg):
            value = arg.astype('timedelta64[{0}]'.format(
                unit)).astype('timedelta64[ns]', copy=False)
        else:
            value = tslib.array_to_timedelta64(
                _ensure_object(arg), unit=unit, errors=errors)
            value = value.astype('timedelta64[ns]', copy=False)

        if box:
            from pandas import TimedeltaIndex
            value = TimedeltaIndex(value, unit='ns', name=name)
        return value

    if arg is None:
        return arg
    elif isinstance(arg, ABCSeries):
        from pandas import Series
        values = _convert_listlike(arg._values, box=False, unit=unit)
        return Series(values, index=arg.index, name=arg.name, dtype='m8[ns]')
    elif isinstance(arg, ABCIndexClass):
        return _convert_listlike(arg, box=box, unit=unit, name=arg.name)
    elif is_list_like(arg) and getattr(arg, 'ndim', 1) == 1:
        return _convert_listlike(arg, box=box, unit=unit)
    elif getattr(arg, 'ndim', 1) > 1:
        raise TypeError('arg must be a string, timedelta, list, tuple, '
                        '1-d array, or Series')

    # ...so it must be a scalar value. Return scalar.
    return _coerce_scalar_to_timedelta_type(arg, unit=unit,
                                            box=box, errors=errors)

_unit_map = {
    'Y': 'Y',
    'y': 'Y',
    'W': 'W',
    'w': 'W',
    'D': 'D',
    'd': 'D',
    'days': 'D',
    'Days': 'D',
    'day': 'D',
    'Day': 'D',
    'M': 'M',
    'H': 'h',
    'h': 'h',
    'm': 'm',
    'T': 'm',
    'S': 's',
    's': 's',
    'L': 'ms',
    'MS': 'ms',
    'ms': 'ms',
    'US': 'us',
    'us': 'us',
    'NS': 'ns',
    'ns': 'ns',
}


def _validate_timedelta_unit(arg):
    """ provide validation / translation for timedelta short units """
    try:
        return _unit_map[arg]
    except:
        if arg is None:
            return 'ns'
        raise ValueError("invalid timedelta unit {0} provided".format(arg))


def _coerce_scalar_to_timedelta_type(r, unit='ns', box=True, errors='raise'):
    """
    convert strings to timedelta; coerce to Timedelta (if box), else
    np.timedelta64
    """

    result = tslib.convert_to_timedelta(r, unit, errors)
    if box:
        result = tslib.Timedelta(result)

    return result
