# -*- coding: utf-8 -*-
from datetime import timedelta
from pandas.compat import zip
from pandas import compat
import re
import warnings

import numpy as np

from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.common import (
    is_period_arraylike,
    is_timedelta64_dtype,
    is_datetime64_dtype)

from pandas.tseries.offsets import DateOffset
from pandas.util._decorators import deprecate_kwarg
import pandas.tseries.offsets as offsets

from pandas._libs.tslib import Timedelta
from pandas._libs.tslibs.frequencies import (  # noqa
    get_freq_code, _base_and_stride, _period_str_to_code,
    _INVALID_FREQ_ERROR, opattern, _lite_rule_alias, _dont_uppercase,
    _period_code_map, _reverse_period_code_map)
from pandas._libs.tslibs.resolution import (Resolution,
                                            _FrequencyInferer,
                                            _TimedeltaFrequencyInferer)
from pandas._libs.tslibs.parsing import _get_rule_month
from pandas._libs.tslibs.ccalendar import MONTH_NUMBERS

from pytz import AmbiguousTimeError


class FreqGroup(object):
    FR_ANN = 1000
    FR_QTR = 2000
    FR_MTH = 3000
    FR_WK = 4000
    FR_BUS = 5000
    FR_DAY = 6000
    FR_HR = 7000
    FR_MIN = 8000
    FR_SEC = 9000
    FR_MS = 10000
    FR_US = 11000
    FR_NS = 12000


RESO_NS = 0
RESO_US = 1
RESO_MS = 2
RESO_SEC = 3
RESO_MIN = 4
RESO_HR = 5
RESO_DAY = 6


def get_to_timestamp_base(base):
    """
    Return frequency code group used for base of to_timestamp against
    frequency code.

    Example
    -------
    # Return day freq code against longer freq than day
    >>> get_to_timestamp_base(get_freq_code('D')[0])
    6000
    >>> get_to_timestamp_base(get_freq_code('W')[0])
    6000
    >>> get_to_timestamp_base(get_freq_code('M')[0])
    6000

    # Return second freq code against hour between second
    >>> get_to_timestamp_base(get_freq_code('H')[0])
    9000
    >>> get_to_timestamp_base(get_freq_code('S')[0])
    9000
    """
    if base < FreqGroup.FR_BUS:
        return FreqGroup.FR_DAY
    if FreqGroup.FR_HR <= base <= FreqGroup.FR_SEC:
        return FreqGroup.FR_SEC
    return base


def get_freq(freq):
    """
    Return frequency code of given frequency str.
    If input is not string, return input as it is.

    Example
    -------
    >>> get_freq('A')
    1000

    >>> get_freq('3A')
    1000
    """
    if isinstance(freq, compat.string_types):
        base, mult = get_freq_code(freq)
        freq = base
    return freq


def _get_freq_str(base, mult=1):
    code = _reverse_period_code_map.get(base)
    if mult == 1:
        return code
    return str(mult) + code


# ---------------------------------------------------------------------
# Offset names ("time rules") and related functions

from pandas._libs.tslibs.offsets import _offset_to_period_map  # noqa:E402
from pandas.tseries.offsets import (Nano, Micro, Milli, Second,  # noqa
                                    Minute, Hour,
                                    Day, BDay, CDay, Week, MonthBegin,
                                    MonthEnd, BMonthBegin, BMonthEnd,
                                    QuarterBegin, QuarterEnd, BQuarterBegin,
                                    BQuarterEnd, YearBegin, YearEnd,
                                    BYearBegin, BYearEnd, prefix_mapping)
try:
    cday = CDay()
except NotImplementedError:
    cday = None

#: cache of previously seen offsets
_offset_map = {}


def get_period_alias(offset_str):
    """ alias to closest period strings BQ->Q etc"""
    return _offset_to_period_map.get(offset_str, None)


_name_to_offset_map = {'days': Day(1),
                       'hours': Hour(1),
                       'minutes': Minute(1),
                       'seconds': Second(1),
                       'milliseconds': Milli(1),
                       'microseconds': Micro(1),
                       'nanoseconds': Nano(1)}


@deprecate_kwarg(old_arg_name='freqstr', new_arg_name='freq')
def to_offset(freq):
    """
    Return DateOffset object from string or tuple representation
    or datetime.timedelta object

    Parameters
    ----------
    freq : str, tuple, datetime.timedelta, DateOffset or None

    Returns
    -------
    delta : DateOffset
        None if freq is None

    Raises
    ------
    ValueError
        If freq is an invalid frequency

    See Also
    --------
    pandas.DateOffset

    Examples
    --------
    >>> to_offset('5min')
    <5 * Minutes>

    >>> to_offset('1D1H')
    <25 * Hours>

    >>> to_offset(('W', 2))
    <2 * Weeks: weekday=6>

    >>> to_offset((2, 'B'))
    <2 * BusinessDays>

    >>> to_offset(datetime.timedelta(days=1))
    <Day>

    >>> to_offset(Hour())
    <Hour>
    """
    if freq is None:
        return None

    if isinstance(freq, DateOffset):
        return freq

    if isinstance(freq, tuple):
        name = freq[0]
        stride = freq[1]
        if isinstance(stride, compat.string_types):
            name, stride = stride, name
        name, _ = _base_and_stride(name)
        delta = get_offset(name) * stride

    elif isinstance(freq, timedelta):
        delta = None
        freq = Timedelta(freq)
        try:
            for name in freq.components._fields:
                offset = _name_to_offset_map[name]
                stride = getattr(freq.components, name)
                if stride != 0:
                    offset = stride * offset
                    if delta is None:
                        delta = offset
                    else:
                        delta = delta + offset
        except Exception:
            raise ValueError(_INVALID_FREQ_ERROR.format(freq))

    else:
        delta = None
        stride_sign = None
        try:
            splitted = re.split(opattern, freq)
            if splitted[-1] != '' and not splitted[-1].isspace():
                # the last element must be blank
                raise ValueError('last element must be blank')
            for sep, stride, name in zip(splitted[0::4], splitted[1::4],
                                         splitted[2::4]):
                if sep != '' and not sep.isspace():
                    raise ValueError('separator must be spaces')
                prefix = _lite_rule_alias.get(name) or name
                if stride_sign is None:
                    stride_sign = -1 if stride.startswith('-') else 1
                if not stride:
                    stride = 1
                if prefix in Resolution._reso_str_bump_map.keys():
                    stride, name = Resolution.get_stride_from_decimal(
                        float(stride), prefix
                    )
                stride = int(stride)
                offset = get_offset(name)
                offset = offset * int(np.fabs(stride) * stride_sign)
                if delta is None:
                    delta = offset
                else:
                    delta = delta + offset
        except Exception:
            raise ValueError(_INVALID_FREQ_ERROR.format(freq))

    if delta is None:
        raise ValueError(_INVALID_FREQ_ERROR.format(freq))

    return delta


def get_base_alias(freqstr):
    """
    Returns the base frequency alias, e.g., '5D' -> 'D'
    """
    return _base_and_stride(freqstr)[0]


def get_offset(name):
    """
    Return DateOffset object associated with rule name

    Examples
    --------
    get_offset('EOM') --> BMonthEnd(1)
    """
    if name not in _dont_uppercase:
        name = name.upper()
        name = _lite_rule_alias.get(name, name)
        name = _lite_rule_alias.get(name.lower(), name)
    else:
        name = _lite_rule_alias.get(name, name)

    if name not in _offset_map:
        try:
            split = name.split('-')
            klass = prefix_mapping[split[0]]
            # handles case where there's no suffix (and will TypeError if too
            # many '-')
            offset = klass._from_name(*split[1:])
        except (ValueError, TypeError, KeyError):
            # bad prefix or suffix
            raise ValueError(_INVALID_FREQ_ERROR.format(name))
        # cache
        _offset_map[name] = offset
    # do not return cache because it's mutable
    return _offset_map[name].copy()


getOffset = get_offset


def get_standard_freq(freq):
    """
    Return the standardized frequency string
    """

    msg = ("get_standard_freq is deprecated. Use to_offset(freq).rule_code "
           "instead.")
    warnings.warn(msg, FutureWarning, stacklevel=2)
    return to_offset(freq).rule_code


# ---------------------------------------------------------------------
# Period codes


def infer_freq(index, warn=True):
    """
    Infer the most likely frequency given the input index. If the frequency is
    uncertain, a warning will be printed.

    Parameters
    ----------
    index : DatetimeIndex or TimedeltaIndex
      if passed a Series will use the values of the series (NOT THE INDEX)
    warn : boolean, default True

    Returns
    -------
    freq : string or None
        None if no discernible frequency
        TypeError if the index is not datetime-like
        ValueError if there are less than three values.
    """
    import pandas as pd

    if isinstance(index, ABCSeries):
        values = index._values
        if not (is_datetime64_dtype(values) or
                is_timedelta64_dtype(values) or
                values.dtype == object):
            raise TypeError("cannot infer freq from a non-convertible dtype "
                            "on a Series of {dtype}".format(dtype=index.dtype))
        index = values

    if is_period_arraylike(index):
        raise TypeError("PeriodIndex given. Check the `freq` attribute "
                        "instead of using infer_freq.")
    elif isinstance(index, pd.TimedeltaIndex):
        inferer = _TimedeltaFrequencyInferer(index, warn=warn)
        return inferer.get_freq()

    if isinstance(index, pd.Index) and not isinstance(index, pd.DatetimeIndex):
        if isinstance(index, (pd.Int64Index, pd.Float64Index)):
            raise TypeError("cannot infer freq from a non-convertible index "
                            "type {type}".format(type=type(index)))
        index = index.values

    if not isinstance(index, pd.DatetimeIndex):
        try:
            index = pd.DatetimeIndex(index)
        except AmbiguousTimeError:
            index = pd.DatetimeIndex(index.asi8)

    inferer = _FrequencyInferer(index, warn=warn)
    return inferer.get_freq()


def _maybe_coerce_freq(code):
    """ we might need to coerce a code to a rule_code
    and uppercase it

    Parameters
    ----------
    source : string
        Frequency converting from

    Returns
    -------
    string code
    """

    assert code is not None
    if isinstance(code, offsets.DateOffset):
        code = code.rule_code
    return code.upper()


def is_subperiod(source, target):
    """
    Returns True if downsampling is possible between source and target
    frequencies

    Parameters
    ----------
    source : string
        Frequency converting from
    target : string
        Frequency converting to

    Returns
    -------
    is_subperiod : boolean
    """

    if target is None or source is None:
        return False
    source = _maybe_coerce_freq(source)
    target = _maybe_coerce_freq(target)

    if _is_annual(target):
        if _is_quarterly(source):
            return _quarter_months_conform(_get_rule_month(source),
                                           _get_rule_month(target))
        return source in ['D', 'C', 'B', 'M', 'H', 'T', 'S', 'L', 'U', 'N']
    elif _is_quarterly(target):
        return source in ['D', 'C', 'B', 'M', 'H', 'T', 'S', 'L', 'U', 'N']
    elif _is_monthly(target):
        return source in ['D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N']
    elif _is_weekly(target):
        return source in [target, 'D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N']
    elif target == 'B':
        return source in ['B', 'H', 'T', 'S', 'L', 'U', 'N']
    elif target == 'C':
        return source in ['C', 'H', 'T', 'S', 'L', 'U', 'N']
    elif target == 'D':
        return source in ['D', 'H', 'T', 'S', 'L', 'U', 'N']
    elif target == 'H':
        return source in ['H', 'T', 'S', 'L', 'U', 'N']
    elif target == 'T':
        return source in ['T', 'S', 'L', 'U', 'N']
    elif target == 'S':
        return source in ['S', 'L', 'U', 'N']
    elif target == 'L':
        return source in ['L', 'U', 'N']
    elif target == 'U':
        return source in ['U', 'N']
    elif target == 'N':
        return source in ['N']


def is_superperiod(source, target):
    """
    Returns True if upsampling is possible between source and target
    frequencies

    Parameters
    ----------
    source : string
        Frequency converting from
    target : string
        Frequency converting to

    Returns
    -------
    is_superperiod : boolean
    """
    if target is None or source is None:
        return False
    source = _maybe_coerce_freq(source)
    target = _maybe_coerce_freq(target)

    if _is_annual(source):
        if _is_annual(target):
            return _get_rule_month(source) == _get_rule_month(target)

        if _is_quarterly(target):
            smonth = _get_rule_month(source)
            tmonth = _get_rule_month(target)
            return _quarter_months_conform(smonth, tmonth)
        return target in ['D', 'C', 'B', 'M', 'H', 'T', 'S', 'L', 'U', 'N']
    elif _is_quarterly(source):
        return target in ['D', 'C', 'B', 'M', 'H', 'T', 'S', 'L', 'U', 'N']
    elif _is_monthly(source):
        return target in ['D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N']
    elif _is_weekly(source):
        return target in [source, 'D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N']
    elif source == 'B':
        return target in ['D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N']
    elif source == 'C':
        return target in ['D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N']
    elif source == 'D':
        return target in ['D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N']
    elif source == 'H':
        return target in ['H', 'T', 'S', 'L', 'U', 'N']
    elif source == 'T':
        return target in ['T', 'S', 'L', 'U', 'N']
    elif source == 'S':
        return target in ['S', 'L', 'U', 'N']
    elif source == 'L':
        return target in ['L', 'U', 'N']
    elif source == 'U':
        return target in ['U', 'N']
    elif source == 'N':
        return target in ['N']


def _is_annual(rule):
    rule = rule.upper()
    return rule == 'A' or rule.startswith('A-')


def _quarter_months_conform(source, target):
    snum = MONTH_NUMBERS[source]
    tnum = MONTH_NUMBERS[target]
    return snum % 3 == tnum % 3


def _is_quarterly(rule):
    rule = rule.upper()
    return rule == 'Q' or rule.startswith('Q-') or rule.startswith('BQ')


def _is_monthly(rule):
    rule = rule.upper()
    return rule == 'M' or rule == 'BM'


def _is_weekly(rule):
    rule = rule.upper()
    return rule == 'W' or rule.startswith('W-')
