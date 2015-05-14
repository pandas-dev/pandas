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

        if is_timedelta64_dtype(arg):
            value = arg.astype('timedelta64[ns]')
        elif is_integer_dtype(arg):

            # these are shortcutable
            value = arg.astype('timedelta64[{0}]'.format(unit)).astype('timedelta64[ns]')
        else:
            try:
                value = tslib.array_to_timedelta64(_ensure_object(arg), unit=unit, coerce=coerce)
            except:

                # try to process strings fast; may need to fallback
                try:
                    value = np.array([ _get_string_converter(r, unit=unit)() for r in arg ],dtype='m8[ns]')
                except:
                    value = np.array([ _coerce_scalar_to_timedelta_type(r, unit=unit, coerce=coerce) for r in arg ])
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
_unit_scale = {
    'd' : 86400*1e9,
    'h' : 3600*1e9,
    'm' : 60*1e9,
    's' : 1e9,
    'ms' : 1e6,
    'us' : 1e3,
    'ns' : 1,
    }

def _validate_timedelta_unit(arg):
    """ provide validation / translation for timedelta short units """
    try:
        return _unit_map[arg]
    except:
        if arg is None:
            return 'ns'
        raise ValueError("invalid timedelta unit {0} provided".format(arg))

_short_search = re.compile(
    "^\s*(?P<neg>-?)\s*(?P<value>\d*\.?\d*)\s*(?P<unit>d|s|ms|us|ns)?\s*$",re.IGNORECASE)
_full_search = re.compile(
    "^\s*(?P<neg>-?)\s*(?P<days>\d*?\.?\d*?)?\s*(days|d|day)?,?\s*\+?(?P<time>\d{1,2}:\d{2}:\d{2})?(?P<frac>\.\d+)?\s*$",re.IGNORECASE)
_nat_search = re.compile(
    "^\s*(nat|nan)\s*$",re.IGNORECASE)
_whitespace = re.compile('^\s*$')
_number_split = re.compile("^(\d+\.?\d*)")

# construct the full2_search
abbrevs = [('d' ,'days|d|day'),
           ('h' ,'hours|h|hour'),
           ('m' ,'minutes|min|minute|m'),
           ('s' ,'seconds|sec|second|s'),
           ('ms','milliseconds|milli|millis|millisecond|ms'),
           ('us','microseconds|micro|micros|microsecond|us'),
           ('ns','nanoseconds|nano|nanos|nanosecond|ns')]

_full_search2 = re.compile(''.join(
    ["^\s*(?P<neg>-?)\s*"] + [ "(?P<" + p + ">\\d+\.?\d*\s*(" + ss + "))?\\s*" for p, ss in abbrevs ] + ['$']))

def _coerce_scalar_to_timedelta_type(r, unit='ns', box=True, coerce=False):
    """ convert strings to timedelta; coerce to Timedelta (if box), else np.timedelta64"""

    if isinstance(r, compat.string_types):

        # we are already converting to nanoseconds
        converter = _get_string_converter(r, unit=unit)
        r = converter()
        unit='ns'

    result = tslib.convert_to_timedelta(r,unit,coerce)
    if box:
        result = tslib.Timedelta(result)

    return result

def _get_string_converter(r, unit='ns'):
    """ return a string converter for r to process the timedelta format """

    # treat as a nan
    if isnull(r):
        def convert(r=None, unit=None):
            return tslib.iNaT
        return convert

    if _whitespace.search(r):
        def convert(r=None, unit=None):
            return tslib.iNaT
        return convert

    m = _short_search.search(r)
    if m:
        def convert(r=None, unit=unit, m=m):
            if r is not None:
                m = _short_search.search(r)

            gd = m.groupdict()

            r = float(gd['value'])
            u = gd.get('unit')
            if u is not None:
                unit = u.lower()
            result = tslib.cast_from_unit(r, unit)
            if gd['neg']:
                result *= -1
            return result
        return convert

    m = _full_search.search(r)
    if m:
        def convert(r=None, unit=None, m=m):
            if r is not None:
                m = _full_search.search(r)

            gd = m.groupdict()

            # handle time
            value = 0
            time = gd['time']
            if time:
                (hh,mm,ss) = time.split(':')
                value += int((float(hh)*3600 + float(mm)*60 + float(ss))*1e9)

            # handle frac
            frac = gd['frac']
            if frac:
                value += round(float(frac)*1e9)

            # handle days (possibly negative)
            is_neg = gd['neg']
            if gd['days']:
                days = int((float(gd['days'] or 0) * 86400)*1e9)
                if is_neg:
                    days *= -1
                value += days
            else:
                if is_neg:
                    value *= -1
            return tslib.cast_from_unit(value, 'ns')
        return convert

    # look for combo strings
    m = _full_search2.search(r)
    if m:
        def convert(r=None, unit=None, m=m):
            if r is not None:
                m = _full_search2.search(r)

            gd = m.groupdict()

            # the parser
            def parse(k, v):
                if v is None:
                    return 0
                v = float(_number_split.search(v).group())
                return int(v*_unit_scale[k])

            # handle non-days
            days = gd.pop('days',None)
            neg = gd.pop('neg',None)
            value = 0
            for k, v in gd.items():
                value += parse(k,v)

            # parse days / neg
            if days:
                days = parse('days',days)
                if neg:
                    days *= -1
                value += days
            else:
                if neg:
                    value *= -1

            return tslib.cast_from_unit(value, 'ns')
        return convert

    m = _nat_search.search(r)
    if m:
        def convert(r=None, unit=None, m=m):
            return tslib.iNaT
        return convert

    # no converter
    raise ValueError("cannot create timedelta string converter for [{0}]".format(r))

