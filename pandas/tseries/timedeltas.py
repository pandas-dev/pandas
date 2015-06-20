"""
timedelta support tools
"""

import re
import numpy as np
import pandas.tslib as tslib
from pandas import compat
from pandas.core.common import (ABCSeries, is_integer_dtype,
                                is_timedelta64_dtype, is_list_like,
                                isnull, _ensure_object)

def to_timedelta(arg, unit='ns', box=True, coerce=False):
    """
    Convert argument to timedelta

    Parameters
    ----------
    arg : string, timedelta, array of strings (with possible NAs)
    unit : unit of the arg (D,h,m,s,ms,us,ns) denote the unit, which is an integer/float number
    box : boolean, default True
        If True returns a Timedelta/TimedeltaIndex of the results
        if False returns a np.timedelta64 or ndarray of values of dtype timedelta64[ns]
    coerce : force errors to NaT (False by default)

    Returns
    -------
    ret : timedelta64/arrays of timedelta64 if parsing succeeded
    """
    unit = _validate_timedelta_unit(unit)

    def _convert_listlike(arg, box, unit):

        if isinstance(arg, (list,tuple)) or ((hasattr(arg,'__iter__') and not hasattr(arg,'dtype'))):
            arg = np.array(list(arg), dtype='O')

        # these are shortcutable
        if is_timedelta64_dtype(arg):
            value = arg.astype('timedelta64[ns]')
        elif is_integer_dtype(arg):
            value = arg.astype('timedelta64[{0}]'.format(unit)).astype('timedelta64[ns]', copy=False)
        else:
            value = tslib.array_to_timedelta64(_ensure_object(arg), unit=unit, coerce=coerce)
            value = value.astype('timedelta64[ns]', copy=False)

        if box:
            from pandas import TimedeltaIndex
            value = TimedeltaIndex(value,unit='ns')
        return value

    if arg is None:
        return arg
    elif isinstance(arg, ABCSeries):
        from pandas import Series
        values = _convert_listlike(arg.values, box=False, unit=unit)
        return Series(values, index=arg.index, name=arg.name, dtype='m8[ns]')
    elif is_list_like(arg):
        return _convert_listlike(arg, box=box, unit=unit)

    # ...so it must be a scalar value. Return scalar.
    return _coerce_scalar_to_timedelta_type(arg, unit=unit, box=box, coerce=coerce)

_unit_map = {
    'Y' : 'Y',
    'y' : 'Y',
    'W' : 'W',
    'w' : 'W',
    'D' : 'D',
    'd' : 'D',
    'days' : 'D',
    'Days' : 'D',
    'day'  : 'D',
    'Day'  : 'D',
    'M'    : 'M',
    'H'  : 'h',
    'h'  : 'h',
    'm'  : 'm',
    'T'  : 'm',
    'S'  : 's',
    's'  : 's',
    'L'  : 'ms',
    'MS' : 'ms',
    'ms' : 'ms',
    'US' : 'us',
    'us' : 'us',
    'NS' : 'ns',
    'ns' : 'ns',
    }

def _validate_timedelta_unit(arg):
    """ provide validation / translation for timedelta short units """
    try:
        return _unit_map[arg]
    except:
        if arg is None:
            return 'ns'
        raise ValueError("invalid timedelta unit {0} provided".format(arg))

def _coerce_scalar_to_timedelta_type(r, unit='ns', box=True, coerce=False):
    """ convert strings to timedelta; coerce to Timedelta (if box), else np.timedelta64"""

    result = tslib.convert_to_timedelta(r,unit,coerce)
    if box:
        result = tslib.Timedelta(result)

    return result
