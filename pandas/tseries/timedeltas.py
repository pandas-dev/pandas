"""
timedelta support tools
"""

import re
from datetime import timedelta

import numpy as np
import pandas.tslib as tslib
from pandas import compat, _np_version_under1p7
from pandas.core.common import ABCSeries, is_integer, _values_from_object

timedelta_search = re.compile(
    "^(?P<value>-?\d*\.?\d*)(?P<unit>D|s|ms|us|ns)?$")

def _coerce_scalar_to_timedelta_type(r, unit='ns'):
    # kludgy here until we have a timedelta scalar
    # handle the numpy < 1.7 case

    if isinstance(r, compat.string_types):
        m = timedelta_search.search(r)
        if m:
            r = float(m.groupdict()['value'])
            u = m.groupdict().get('unit')
            if u is not None:
                unit = u
        else:
            raise ValueError("cannot convert timedelta scalar value!")

        r = tslib.cast_from_unit(unit, r)
        r = timedelta(microseconds=int(r)/1000)

    if is_integer(r):
        r = tslib.cast_from_unit(unit, r)
        r = timedelta(microseconds=int(r)/1000)

    if _np_version_under1p7:
        if not isinstance(r, timedelta):
            raise AssertionError("Invalid type for timedelta scalar: %s" % type(r))
        if compat.PY3:
            # convert to microseconds in timedelta64
            r = np.timedelta64(int(r.total_seconds()*1e9 + r.microseconds*1000))
        else:
            return r

    if isinstance(r, timedelta):
        r = np.timedelta64(r)
    elif not isinstance(r, np.timedelta64):
        raise AssertionError("Invalid type for timedelta scalar: %s" % type(r))
    return r.astype('timedelta64[ns]')

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

            if td == tslib.compat_NaT:
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

