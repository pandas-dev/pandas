"""
timedelta support tools
"""

import re
from datetime import timedelta

import numpy as np
import pandas.tslib as tslib
from pandas import compat, _np_version_under1p7
from pandas.core.common import (ABCSeries, is_integer, is_timedelta64_dtype,
                                _values_from_object, is_list_like, isnull)

repr_timedelta = tslib.repr_timedelta64
repr_timedelta64 = tslib.repr_timedelta64

def to_timedelta(arg, box=True, unit='ns'):
    """
    Convert argument to timedelta

    Parameters
    ----------
    arg : string, timedelta, array of strings (with possible NAs)
    box : boolean, default True
        If True returns a Series of the results, if False returns ndarray of values
    unit : unit of the arg (D,s,ms,us,ns) denote the unit, which is an integer/float number

    Returns
    -------
    ret : timedelta64/arrays of timedelta64 if parsing succeeded
    """
    if _np_version_under1p7:
        raise ValueError("to_timedelta is not support for numpy < 1.7")

    def _convert_listlike(arg, box):

        if isinstance(arg, (list,tuple)):
            arg = np.array(arg, dtype='O')

        if is_timedelta64_dtype(arg):
            if box:
                from pandas import Series
                return Series(arg,dtype='m8[ns]')
            return arg

        value = np.array([ _coerce_scalar_to_timedelta_type(r, unit=unit) for r in arg ])
        if box:
            from pandas import Series
            value = Series(value,dtype='m8[ns]')
        return value

    if arg is None:
        return arg
    elif isinstance(arg, ABCSeries):
        from pandas import Series
        values = _convert_listlike(arg.values, box=False)
        return Series(values, index=arg.index, name=arg.name, dtype='m8[ns]')
    elif is_list_like(arg):
        return _convert_listlike(arg, box=box)

    # ...so it must be a scalar value. Return scalar.
    return _coerce_scalar_to_timedelta_type(arg, unit=unit)

_short_search = re.compile(
    "^\s*(?P<neg>-?)\s*(?P<value>\d*\.?\d*)\s*(?P<unit>d|s|ms|us|ns)?\s*$",re.IGNORECASE)
_full_search = re.compile(
    "^\s*(?P<neg>-?)\s*(?P<days>\d+)?\s*(days|d)?,?\s*(?P<time>\d{2}:\d{2}:\d{2})?(?P<frac>\.\d+)?\s*$",re.IGNORECASE)
_nat_search = re.compile(
    "^\s*(nat|nan)\s*$",re.IGNORECASE)
_whitespace = re.compile('^\s*$')

def _coerce_scalar_to_timedelta_type(r, unit='ns'):
    """ convert strings to timedelta; coerce to np.timedelta64"""

    if isinstance(r, compat.string_types):

        # we are already converting to nanoseconds
        converter = _get_string_converter(r, unit=unit)
        r = converter()
        unit='ns'

    return tslib.convert_to_timedelta(r,unit)

def _get_string_converter(r, unit='ns'):
    """ return a string converter for r to process the timedelta format """

    # treat as a nan
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
            if gd['neg']:
                r *= -1
            return tslib.cast_from_unit(r, unit)
        return convert

    m = _full_search.search(r)
    if m:
        def convert(r=None, unit=None, m=m):
            if r is not None:
                m = _full_search.search(r)

            gd = m.groupdict()

            # convert to seconds
            value = float(gd['days'] or 0) * 86400

            time = gd['time']
            if time:
                (hh,mm,ss) = time.split(':')
                value += float(hh)*3600 + float(mm)*60 + float(ss)

            frac = gd['frac']
            if frac:
                value += float(frac)

            if gd['neg']:
                value *= -1
            return tslib.cast_from_unit(value, 's')
        return convert

    m = _nat_search.search(r)
    if m:
        def convert(r=None, unit=None, m=m):
            return tslib.iNaT
        return convert

    # no converter
    raise ValueError("cannot create timedelta string converter")

def _possibly_cast_to_timedelta(value, coerce=True):
    """ try to cast to timedelta64, if already a timedeltalike, then make
        sure that we are [ns] (as numpy 1.6.2 is very buggy in this regards,
        don't force the conversion unless coerce is True

        if coerce='compat' force a compatibilty coercerion (to timedeltas) if needeed
        """

    # coercion compatability
    if coerce == 'compat' and _np_version_under1p7:

        def convert(td, dtype):

            # we have an array with a non-object dtype
            if hasattr(td,'item'):
                td = td.astype(np.int64).item()
                if td == tslib.iNaT:
                    return td
                if dtype == 'm8[us]':
                    td *= 1000
                return td

            if isnull(td) or td == tslib.compat_NaT or td == tslib.iNaT:
                return tslib.iNaT

            # convert td value to a nanosecond value
            d = td.days
            s = td.seconds
            us = td.microseconds

            if dtype == 'object' or dtype == 'm8[ns]':
                td = 1000*us + (s + d * 24 * 3600) * 10 ** 9
            else:
                raise ValueError("invalid conversion of dtype in np < 1.7 [%s]" % dtype)

            return td

        # < 1.7 coercion
        if not is_list_like(value):
            value = np.array([ value ])

        dtype = value.dtype
        return np.array([ convert(v,dtype) for v in value ], dtype='m8[ns]')

    # deal with numpy not being able to handle certain timedelta operations
    if isinstance(value, (ABCSeries, np.ndarray)) and value.dtype.kind == 'm':
        if value.dtype != 'timedelta64[ns]':
            value = value.astype('timedelta64[ns]')
        return value

    # we don't have a timedelta, but we want to try to convert to one (but
    # don't force it)
    if coerce:
        new_value = tslib.array_to_timedelta64(
            _values_from_object(value).astype(object), coerce=False)
        if new_value.dtype == 'i8':
            value = np.array(new_value, dtype='timedelta64[ns]')

    return value

