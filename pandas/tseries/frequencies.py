from datetime import timedelta
from pandas.compat import range, long, zip
from pandas import compat
import re
import warnings

import numpy as np

import pandas.core.algorithms as algos
from pandas.core.algorithms import unique
from pandas.tseries.offsets import DateOffset
from pandas.util.decorators import cache_readonly
import pandas.tseries.offsets as offsets
import pandas.core.common as com
import pandas.lib as lib
import pandas.tslib as tslib
from pandas.tslib import Timedelta
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


US_RESO = 0
MS_RESO = 1
S_RESO = 2
T_RESO = 3
H_RESO = 4
D_RESO = 5


class Resolution(object):

    # defined in period.pyx
    # note that these are different from freq codes
    RESO_US = US_RESO
    RESO_MS = MS_RESO
    RESO_SEC = S_RESO
    RESO_MIN = T_RESO
    RESO_HR = H_RESO
    RESO_DAY = D_RESO

    _reso_str_map = {
        RESO_US: 'microsecond',
        RESO_MS: 'millisecond',
        RESO_SEC: 'second',
        RESO_MIN: 'minute',
        RESO_HR: 'hour',
        RESO_DAY: 'day'}

    _str_reso_map = dict([(v, k) for k, v in compat.iteritems(_reso_str_map)])

    _reso_freq_map = {
        'year': 'A',
        'quarter': 'Q',
        'month': 'M',
        'day': 'D',
        'hour': 'H',
        'minute': 'T',
        'second': 'S',
        'millisecond': 'L',
        'microsecond': 'U',
        'nanosecond': 'N'}

    _freq_reso_map = dict([(v, k)
                           for k, v in compat.iteritems(_reso_freq_map)])

    @classmethod
    def get_str(cls, reso):
        """
        Return resolution str against resolution code.

        Example
        -------
        >>> Resolution.get_str(Resolution.RESO_SEC)
        'second'
        """
        return cls._reso_str_map.get(reso, 'day')

    @classmethod
    def get_reso(cls, resostr):
        """
        Return resolution str against resolution code.

        Example
        -------
        >>> Resolution.get_reso('second')
        2

        >>> Resolution.get_reso('second') == Resolution.RESO_SEC
        True
        """
        return cls._str_reso_map.get(resostr, cls.RESO_DAY)

    @classmethod
    def get_freq_group(cls, resostr):
        """
        Return frequency str against resolution str.

        Example
        -------
        >>> f.Resolution.get_freq_group('day')
        4000
        """
        return get_freq_group(cls.get_freq(resostr))

    @classmethod
    def get_freq(cls, resostr):
        """
        Return frequency str against resolution str.

        Example
        -------
        >>> f.Resolution.get_freq('day')
        'D'
        """
        return cls._reso_freq_map[resostr]

    @classmethod
    def get_str_from_freq(cls, freq):
        """
        Return resolution str against frequency str.

        Example
        -------
        >>> Resolution.get_str_from_freq('H')
        'hour'
        """
        return cls._freq_reso_map.get(freq, 'day')

    @classmethod
    def get_reso_from_freq(cls, freq):
        """
        Return resolution code against frequency str.

        Example
        -------
        >>> Resolution.get_reso_from_freq('H')
        4

        >>> Resolution.get_reso_from_freq('H') == Resolution.RESO_HR
        True
        """
        return cls.get_reso(cls.get_str_from_freq(freq))


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


def get_freq_group(freq):
    """
    Return frequency code group of given frequency str or offset.

    Example
    -------
    >>> get_freq_group('W-MON')
    4000

    >>> get_freq_group('W-FRI')
    4000
    """
    if isinstance(freq, offsets.DateOffset):
        freq = freq.rule_code

    if isinstance(freq, compat.string_types):
        base, mult = get_freq_code(freq)
        freq = base
    elif isinstance(freq, int):
        pass
    else:
        raise ValueError('input must be str, offset or int')
    return (freq // 1000) * 1000


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


def get_freq_code(freqstr):
    """
    Return freq str or tuple to freq code and stride (mult)

    Parameters
    ----------
    freqstr : str or tuple

    Returns
    -------
    return : tuple of base frequency code and stride (mult)

    Example
    -------
    >>> get_freq_code('3D')
    (6000, 3)

    >>> get_freq_code('D')
    (6000, 1)

    >>> get_freq_code(('D', 3))
    (6000, 3)
    """
    if isinstance(freqstr, DateOffset):
        freqstr = (freqstr.rule_code, freqstr.n)

    if isinstance(freqstr, tuple):
        if (com.is_integer(freqstr[0]) and
                com.is_integer(freqstr[1])):
            # e.g., freqstr = (2000, 1)
            return freqstr
        else:
            # e.g., freqstr = ('T', 5)
            try:
                code = _period_str_to_code(freqstr[0])
                stride = freqstr[1]
            except:
                if com.is_integer(freqstr[1]):
                    raise
                code = _period_str_to_code(freqstr[1])
                stride = freqstr[0]
            return code, stride

    if com.is_integer(freqstr):
        return (freqstr, 1)

    base, stride = _base_and_stride(freqstr)
    code = _period_str_to_code(base)

    return code, stride


def _get_freq_str(base, mult=1):
    code = _reverse_period_code_map.get(base)
    if mult == 1:
        return code
    return str(mult) + code


# ---------------------------------------------------------------------
# Offset names ("time rules") and related functions


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

_offset_to_period_map = {
    'WEEKDAY': 'D',
    'EOM': 'M',
    'BM': 'M',
    'BQS': 'Q',
    'QS': 'Q',
    'BQ': 'Q',
    'BA': 'A',
    'AS': 'A',
    'BAS': 'A',
    'MS': 'M',
    'D': 'D',
    'C': 'C',
    'B': 'B',
    'T': 'T',
    'S': 'S',
    'L': 'L',
    'U': 'U',
    'N': 'N',
    'H': 'H',
    'Q': 'Q',
    'A': 'A',
    'W': 'W',
    'M': 'M'
}

need_suffix = ['QS', 'BQ', 'BQS', 'AS', 'BA', 'BAS']
for __prefix in need_suffix:
    for _m in tslib._MONTHS:
        _offset_to_period_map['%s-%s' % (__prefix, _m)] = \
            _offset_to_period_map[__prefix]
for __prefix in ['A', 'Q']:
    for _m in tslib._MONTHS:
        _alias = '%s-%s' % (__prefix, _m)
        _offset_to_period_map[_alias] = _alias

_days = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']
for _d in _days:
    _offset_to_period_map['W-%s' % _d] = 'W-%s' % _d


def get_period_alias(offset_str):
    """ alias to closest period strings BQ->Q etc"""
    return _offset_to_period_map.get(offset_str, None)

_rule_aliases = {
    # Legacy rules that will continue to map to their original values
    # essentially for the rest of time
    'WEEKDAY': 'B',
    'EOM': 'BM',
    'W@MON': 'W-MON',
    'W@TUE': 'W-TUE',
    'W@WED': 'W-WED',
    'W@THU': 'W-THU',
    'W@FRI': 'W-FRI',
    'W@SAT': 'W-SAT',
    'W@SUN': 'W-SUN',
    'Q@JAN': 'BQ-JAN',
    'Q@FEB': 'BQ-FEB',
    'Q@MAR': 'BQ-MAR',
    'A@JAN': 'BA-JAN',
    'A@FEB': 'BA-FEB',
    'A@MAR': 'BA-MAR',
    'A@APR': 'BA-APR',
    'A@MAY': 'BA-MAY',
    'A@JUN': 'BA-JUN',
    'A@JUL': 'BA-JUL',
    'A@AUG': 'BA-AUG',
    'A@SEP': 'BA-SEP',
    'A@OCT': 'BA-OCT',
    'A@NOV': 'BA-NOV',
    'A@DEC': 'BA-DEC',
}

_lite_rule_alias = {
    'W': 'W-SUN',
    'Q': 'Q-DEC',

    'A': 'A-DEC',  # YearEnd(month=12),
    'AS': 'AS-JAN',  # YearBegin(month=1),
    'BA': 'BA-DEC',  # BYearEnd(month=12),
    'BAS': 'BAS-JAN',  # BYearBegin(month=1),

    'Min': 'T',
    'min': 'T',
    'ms': 'L',
    'us': 'U',
    'ns': 'N'
}

# TODO: Can this be killed?
for _i, _weekday in enumerate(['MON', 'TUE', 'WED', 'THU', 'FRI']):
    for _iweek in range(4):
        _name = 'WOM-%d%s' % (_iweek + 1, _weekday)
        _rule_aliases[_name.replace('-', '@')] = _name

# Note that _rule_aliases is not 1:1 (d[BA]==d[A@DEC]), and so traversal
# order matters when constructing an inverse. we pick one. #2331
# Used in get_legacy_offset_name
_legacy_reverse_map = dict((v, k) for k, v in
                           reversed(sorted(compat.iteritems(_rule_aliases))))

_name_to_offset_map = {'days': Day(1),
                       'hours': Hour(1),
                       'minutes': Minute(1),
                       'seconds': Second(1),
                       'milliseconds': Milli(1),
                       'microseconds': Micro(1),
                       'nanoseconds': Nano(1)}


def to_offset(freqstr):
    """
    Return DateOffset object from string representation or
    Timedelta object

    Examples
    --------
    >>> to_offset('5Min')
    Minute(5)
    """
    if freqstr is None:
        return None

    if isinstance(freqstr, DateOffset):
        return freqstr

    if isinstance(freqstr, tuple):
        name = freqstr[0]
        stride = freqstr[1]
        if isinstance(stride, compat.string_types):
            name, stride = stride, name
        name, _ = _base_and_stride(name)
        delta = get_offset(name) * stride

    elif isinstance(freqstr, timedelta):
        delta = None
        freqstr = Timedelta(freqstr)
        try:
            for name in freqstr.components._fields:
                offset = _name_to_offset_map[name]
                stride = getattr(freqstr.components, name)
                if stride != 0:
                    offset = stride * offset
                    if delta is None:
                        delta = offset
                    else:
                        delta = delta + offset
        except Exception:
            raise ValueError("Could not evaluate %s" % freqstr)

    else:
        delta = None
        stride_sign = None
        try:
            for stride, name, _ in opattern.findall(freqstr):
                offset = get_offset(name)
                if stride_sign is None:
                    stride_sign = -1 if stride.startswith('-') else 1
                if not stride:
                    stride = 1
                stride = int(stride)
                offset = offset * int(np.fabs(stride) * stride_sign)
                if delta is None:
                    delta = offset
                else:
                    delta = delta + offset
        except Exception:
            raise ValueError("Could not evaluate %s" % freqstr)

    if delta is None:
        raise ValueError('Unable to understand %s as a frequency' % freqstr)

    return delta


# hack to handle WOM-1MON
opattern = re.compile(r'([\-]?\d*)\s*([A-Za-z]+([\-@][\dA-Za-z\-]+)?)')


def _base_and_stride(freqstr):
    """
    Return base freq and stride info from string representation

    Examples
    --------
    _freq_and_stride('5Min') -> 'Min', 5
    """
    groups = opattern.match(freqstr)

    if not groups:
        raise ValueError("Could not evaluate %s" % freqstr)

    stride = groups.group(1)

    if len(stride):
        stride = int(stride)
    else:
        stride = 1

    base = groups.group(2)

    return (base, stride)


def get_base_alias(freqstr):
    """
    Returns the base frequency alias, e.g., '5D' -> 'D'
    """
    return _base_and_stride(freqstr)[0]


_dont_uppercase = set(('MS', 'ms'))


_LEGACY_FREQ_WARNING = 'Freq "{0}" is deprecated, use "{1}" as alternative.'


def get_offset(name):
    """
    Return DateOffset object associated with rule name

    Examples
    --------
    get_offset('EOM') --> BMonthEnd(1)
    """
    if name not in _dont_uppercase:
        name = name.upper()

        if name in _rule_aliases:
            new = _rule_aliases[name]
            warnings.warn(_LEGACY_FREQ_WARNING.format(name, new),
                          FutureWarning, stacklevel=2)
            name = new
        elif name.lower() in _rule_aliases:
            new = _rule_aliases[name.lower()]
            warnings.warn(_LEGACY_FREQ_WARNING.format(name, new),
                          FutureWarning, stacklevel=2)
            name = new

        name = _lite_rule_alias.get(name, name)
        name = _lite_rule_alias.get(name.lower(), name)

    else:
        if name in _rule_aliases:
            new = _rule_aliases[name]
            warnings.warn(_LEGACY_FREQ_WARNING.format(name, new),
                          FutureWarning, stacklevel=2)
            name = new
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
            raise ValueError('Bad rule name requested: %s.' % name)
        # cache
        _offset_map[name] = offset
    # do not return cache because it's mutable
    return _offset_map[name].copy()


getOffset = get_offset


def get_offset_name(offset):
    """
    Return rule name associated with a DateOffset object

    Examples
    --------
    get_offset_name(BMonthEnd(1)) --> 'EOM'
    """

    msg = "get_offset_name(offset) is deprecated. Use offset.freqstr instead"
    warnings.warn(msg, FutureWarning, stacklevel=2)
    return offset.freqstr


def get_legacy_offset_name(offset):
    """
    Return the pre pandas 0.8.0 name for the date offset
    """

    # This only used in test_timeseries_legacy.py

    name = offset.name
    return _legacy_reverse_map.get(name, name)


def get_standard_freq(freq):
    """
    Return the standardized frequency string
    """
    if freq is None:
        return None

    if isinstance(freq, DateOffset):
        return freq.rule_code

    code, stride = get_freq_code(freq)
    return _get_freq_str(code, stride)

# ---------------------------------------------------------------------
# Period codes

# period frequency constants corresponding to scikits timeseries
# originals
_period_code_map = {
    # Annual freqs with various fiscal year ends.
    # eg, 2005 for A-FEB runs Mar 1, 2004 to Feb 28, 2005
    "A-DEC": 1000,  # Annual - December year end
    "A-JAN": 1001,  # Annual - January year end
    "A-FEB": 1002,  # Annual - February year end
    "A-MAR": 1003,  # Annual - March year end
    "A-APR": 1004,  # Annual - April year end
    "A-MAY": 1005,  # Annual - May year end
    "A-JUN": 1006,  # Annual - June year end
    "A-JUL": 1007,  # Annual - July year end
    "A-AUG": 1008,  # Annual - August year end
    "A-SEP": 1009,  # Annual - September year end
    "A-OCT": 1010,  # Annual - October year end
    "A-NOV": 1011,  # Annual - November year end

    # Quarterly frequencies with various fiscal year ends.
    # eg, Q42005 for Q-OCT runs Aug 1, 2005 to Oct 31, 2005
    "Q-DEC": 2000,    # Quarterly - December year end
    "Q-JAN": 2001,    # Quarterly - January year end
    "Q-FEB": 2002,    # Quarterly - February year end
    "Q-MAR": 2003,    # Quarterly - March year end
    "Q-APR": 2004,    # Quarterly - April year end
    "Q-MAY": 2005,    # Quarterly - May year end
    "Q-JUN": 2006,    # Quarterly - June year end
    "Q-JUL": 2007,    # Quarterly - July year end
    "Q-AUG": 2008,    # Quarterly - August year end
    "Q-SEP": 2009,    # Quarterly - September year end
    "Q-OCT": 2010,    # Quarterly - October year end
    "Q-NOV": 2011,    # Quarterly - November year end

    "M": 3000,        # Monthly

    "W-SUN": 4000,    # Weekly - Sunday end of week
    "W-MON": 4001,    # Weekly - Monday end of week
    "W-TUE": 4002,    # Weekly - Tuesday end of week
    "W-WED": 4003,    # Weekly - Wednesday end of week
    "W-THU": 4004,    # Weekly - Thursday end of week
    "W-FRI": 4005,    # Weekly - Friday end of week
    "W-SAT": 4006,    # Weekly - Saturday end of week

    "B": 5000,        # Business days
    "D": 6000,        # Daily
    "H": 7000,        # Hourly
    "T": 8000,        # Minutely
    "S": 9000,        # Secondly
    "L": 10000,       # Millisecondly
    "U": 11000,       # Microsecondly
    "N": 12000,       # Nanosecondly
}

_reverse_period_code_map = {}
for _k, _v in compat.iteritems(_period_code_map):
    _reverse_period_code_map[_v] = _k

# Additional aliases
_period_code_map.update({
    "Q": 2000,  # Quarterly - December year end (default quarterly)
    "A": 1000,  # Annual
    "W": 4000,  # Weekly
    "C": 5000,  # Custom Business Day
})


def _period_alias_dictionary():
    """
    Build freq alias dictionary to support freqs from original c_dates.c file
    of the scikits.timeseries library.
    """
    alias_dict = {}

    M_aliases = ["M", "MTH", "MONTH", "MONTHLY"]
    B_aliases = ["B", "BUS", "BUSINESS", "BUSINESSLY", "WEEKDAY"]
    D_aliases = ["D", "DAY", "DLY", "DAILY"]
    H_aliases = ["H", "HR", "HOUR", "HRLY", "HOURLY"]
    T_aliases = ["T", "MIN", "MINUTE", "MINUTELY"]
    S_aliases = ["S", "SEC", "SECOND", "SECONDLY"]
    L_aliases = ["L", "ms", "MILLISECOND", "MILLISECONDLY"]
    U_aliases = ["U", "US", "MICROSECOND", "MICROSECONDLY"]
    N_aliases = ["N", "NS", "NANOSECOND", "NANOSECONDLY"]

    for k in M_aliases:
        alias_dict[k] = 'M'

    for k in B_aliases:
        alias_dict[k] = 'B'

    for k in D_aliases:
        alias_dict[k] = 'D'

    for k in H_aliases:
        alias_dict[k] = 'H'

    for k in T_aliases:
        alias_dict[k] = 'T'

    for k in S_aliases:
        alias_dict[k] = 'S'

    for k in L_aliases:
        alias_dict[k] = 'L'

    for k in U_aliases:
        alias_dict[k] = 'U'

    for k in N_aliases:
        alias_dict[k] = 'N'

    A_prefixes = ["A", "Y", "ANN", "ANNUAL", "ANNUALLY", "YR", "YEAR",
                  "YEARLY"]

    Q_prefixes = ["Q", "QTR", "QUARTER", "QUARTERLY", "Q-E",
                  "QTR-E", "QUARTER-E", "QUARTERLY-E"]

    month_names = [
        ["DEC", "DECEMBER"],
        ["JAN", "JANUARY"],
        ["FEB", "FEBRUARY"],
        ["MAR", "MARCH"],
        ["APR", "APRIL"],
        ["MAY", "MAY"],
        ["JUN", "JUNE"],
        ["JUL", "JULY"],
        ["AUG", "AUGUST"],
        ["SEP", "SEPTEMBER"],
        ["OCT", "OCTOBER"],
        ["NOV", "NOVEMBER"]]

    seps = ["@", "-"]

    for k in A_prefixes:
        alias_dict[k] = 'A'
        for m_tup in month_names:
            for sep in seps:
                m1, m2 = m_tup
                alias_dict[k + sep + m1] = 'A-' + m1
                alias_dict[k + sep + m2] = 'A-' + m1

    for k in Q_prefixes:
        alias_dict[k] = 'Q'
        for m_tup in month_names:
            for sep in seps:
                m1, m2 = m_tup
                alias_dict[k + sep + m1] = 'Q-' + m1
                alias_dict[k + sep + m2] = 'Q-' + m1

    W_prefixes = ["W", "WK", "WEEK", "WEEKLY"]

    day_names = [
        ["SUN", "SUNDAY"],
        ["MON", "MONDAY"],
        ["TUE", "TUESDAY"],
        ["WED", "WEDNESDAY"],
        ["THU", "THURSDAY"],
        ["FRI", "FRIDAY"],
        ["SAT", "SATURDAY"]]

    for k in W_prefixes:
        alias_dict[k] = 'W'
        for d_tup in day_names:
            for sep in ["@", "-"]:
                d1, d2 = d_tup
                alias_dict[k + sep + d1] = 'W-' + d1
                alias_dict[k + sep + d2] = 'W-' + d1

    return alias_dict


_period_alias_dict = _period_alias_dictionary()


def _period_str_to_code(freqstr):
    # hack
    if freqstr in _rule_aliases:
        new = _rule_aliases[freqstr]
        warnings.warn(_LEGACY_FREQ_WARNING.format(freqstr, new),
                      FutureWarning, stacklevel=3)
        freqstr = new
    freqstr = _lite_rule_alias.get(freqstr, freqstr)

    if freqstr not in _dont_uppercase:
        lower = freqstr.lower()
        if lower in _rule_aliases:
            new = _rule_aliases[lower]
            warnings.warn(_LEGACY_FREQ_WARNING.format(lower, new),
                          FutureWarning, stacklevel=3)
            freqstr = new
        freqstr = _lite_rule_alias.get(lower, freqstr)

    try:
        if freqstr not in _dont_uppercase:
            freqstr = freqstr.upper()
        return _period_code_map[freqstr]
    except KeyError:
        try:
            alias = _period_alias_dict[freqstr]
            warnings.warn(_LEGACY_FREQ_WARNING.format(freqstr, alias),
                          FutureWarning, stacklevel=3)
        except KeyError:
            raise ValueError("Unknown freqstr: %s" % freqstr)

        return _period_code_map[alias]


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

    if isinstance(index, com.ABCSeries):
        values = index._values
        if not (com.is_datetime64_dtype(values) or
                com.is_timedelta64_dtype(values) or
                values.dtype == object):
            raise TypeError("cannot infer freq from a non-convertible "
                            "dtype on a Series of {0}".format(index.dtype))
        index = values

    if com.is_period_arraylike(index):
        raise TypeError("PeriodIndex given. Check the `freq` attribute "
                        "instead of using infer_freq.")
    elif isinstance(index, pd.TimedeltaIndex):
        inferer = _TimedeltaFrequencyInferer(index, warn=warn)
        return inferer.get_freq()

    if isinstance(index, pd.Index) and not isinstance(index, pd.DatetimeIndex):
        if isinstance(index, (pd.Int64Index, pd.Float64Index)):
            raise TypeError("cannot infer freq from a non-convertible index "
                            "type {0}".format(type(index)))
        index = index.values

    if not isinstance(index, pd.DatetimeIndex):
        try:
            index = pd.DatetimeIndex(index)
        except AmbiguousTimeError:
            index = pd.DatetimeIndex(index.asi8)

    inferer = _FrequencyInferer(index, warn=warn)
    return inferer.get_freq()

_ONE_MICRO = long(1000)
_ONE_MILLI = _ONE_MICRO * 1000
_ONE_SECOND = _ONE_MILLI * 1000
_ONE_MINUTE = 60 * _ONE_SECOND
_ONE_HOUR = 60 * _ONE_MINUTE
_ONE_DAY = 24 * _ONE_HOUR


class _FrequencyInferer(object):
    """
    Not sure if I can avoid the state machine here
    """

    def __init__(self, index, warn=True):
        self.index = index
        self.values = np.asarray(index).view('i8')

        # This moves the values, which are implicitly in UTC, to the
        # the timezone so they are in local time
        if hasattr(index, 'tz'):
            if index.tz is not None:
                self.values = tslib.tz_convert(self.values, 'UTC', index.tz)

        self.warn = warn

        if len(index) < 3:
            raise ValueError('Need at least 3 dates to infer frequency')

        self.is_monotonic = (self.index.is_monotonic_increasing or
                             self.index.is_monotonic_decreasing)

    @cache_readonly
    def deltas(self):
        return tslib.unique_deltas(self.values)

    @cache_readonly
    def deltas_asi8(self):
        return tslib.unique_deltas(self.index.asi8)

    @cache_readonly
    def is_unique(self):
        return len(self.deltas) == 1

    @cache_readonly
    def is_unique_asi8(self):
        return len(self.deltas_asi8) == 1

    def get_freq(self):
        if not self.is_monotonic or not self.index.is_unique:
            return None

        delta = self.deltas[0]
        if _is_multiple(delta, _ONE_DAY):
            return self._infer_daily_rule()
        else:
            # Business hourly, maybe. 17: one day / 65: one weekend
            if self.hour_deltas in ([1, 17], [1, 65], [1, 17, 65]):
                return 'BH'
            # Possibly intraday frequency.  Here we use the
            # original .asi8 values as the modified values
            # will not work around DST transitions.  See #8772
            elif not self.is_unique_asi8:
                return None
            delta = self.deltas_asi8[0]
            if _is_multiple(delta, _ONE_HOUR):
                # Hours
                return _maybe_add_count('H', delta / _ONE_HOUR)
            elif _is_multiple(delta, _ONE_MINUTE):
                # Minutes
                return _maybe_add_count('T', delta / _ONE_MINUTE)
            elif _is_multiple(delta, _ONE_SECOND):
                # Seconds
                return _maybe_add_count('S', delta / _ONE_SECOND)
            elif _is_multiple(delta, _ONE_MILLI):
                # Milliseconds
                return _maybe_add_count('L', delta / _ONE_MILLI)
            elif _is_multiple(delta, _ONE_MICRO):
                # Microseconds
                return _maybe_add_count('U', delta / _ONE_MICRO)
            else:
                # Nanoseconds
                return _maybe_add_count('N', delta)

    @cache_readonly
    def day_deltas(self):
        return [x / _ONE_DAY for x in self.deltas]

    @cache_readonly
    def hour_deltas(self):
        return [x / _ONE_HOUR for x in self.deltas]

    @cache_readonly
    def fields(self):
        return tslib.build_field_sarray(self.values)

    @cache_readonly
    def rep_stamp(self):
        return lib.Timestamp(self.values[0])

    def month_position_check(self):
        # TODO: cythonize this, very slow
        calendar_end = True
        business_end = True
        calendar_start = True
        business_start = True

        years = self.fields['Y']
        months = self.fields['M']
        days = self.fields['D']
        weekdays = self.index.dayofweek

        from calendar import monthrange
        for y, m, d, wd in zip(years, months, days, weekdays):

            if calendar_start:
                calendar_start &= d == 1
            if business_start:
                business_start &= d == 1 or (d <= 3 and wd == 0)

            if calendar_end or business_end:
                _, daysinmonth = monthrange(y, m)
                cal = d == daysinmonth
                if calendar_end:
                    calendar_end &= cal
                if business_end:
                    business_end &= cal or (daysinmonth - d < 3 and wd == 4)
            elif not calendar_start and not business_start:
                break

        if calendar_end:
            return 'ce'
        elif business_end:
            return 'be'
        elif calendar_start:
            return 'cs'
        elif business_start:
            return 'bs'
        else:
            return None

    @cache_readonly
    def mdiffs(self):
        nmonths = self.fields['Y'] * 12 + self.fields['M']
        return tslib.unique_deltas(nmonths.astype('i8'))

    @cache_readonly
    def ydiffs(self):
        return tslib.unique_deltas(self.fields['Y'].astype('i8'))

    def _infer_daily_rule(self):
        annual_rule = self._get_annual_rule()
        if annual_rule:
            nyears = self.ydiffs[0]
            month = _month_aliases[self.rep_stamp.month]
            return _maybe_add_count('%s-%s' % (annual_rule, month), nyears)

        quarterly_rule = self._get_quarterly_rule()
        if quarterly_rule:
            nquarters = self.mdiffs[0] / 3
            mod_dict = {0: 12, 2: 11, 1: 10}
            month = _month_aliases[mod_dict[self.rep_stamp.month % 3]]
            return _maybe_add_count('%s-%s' % (quarterly_rule, month),
                                    nquarters)

        monthly_rule = self._get_monthly_rule()
        if monthly_rule:
            return _maybe_add_count(monthly_rule, self.mdiffs[0])

        if self.is_unique:
            days = self.deltas[0] / _ONE_DAY
            if days % 7 == 0:
                # Weekly
                alias = _weekday_rule_aliases[self.rep_stamp.weekday()]
                return _maybe_add_count('W-%s' % alias, days / 7)
            else:
                return _maybe_add_count('D', days)

        # Business daily. Maybe
        if self.day_deltas == [1, 3]:
            return 'B'

        wom_rule = self._get_wom_rule()
        if wom_rule:
            return wom_rule

    def _get_annual_rule(self):
        if len(self.ydiffs) > 1:
            return None

        if len(algos.unique(self.fields['M'])) > 1:
            return None

        pos_check = self.month_position_check()
        return {'cs': 'AS', 'bs': 'BAS',
                'ce': 'A', 'be': 'BA'}.get(pos_check)

    def _get_quarterly_rule(self):
        if len(self.mdiffs) > 1:
            return None

        if not self.mdiffs[0] % 3 == 0:
            return None

        pos_check = self.month_position_check()
        return {'cs': 'QS', 'bs': 'BQS',
                'ce': 'Q', 'be': 'BQ'}.get(pos_check)

    def _get_monthly_rule(self):
        if len(self.mdiffs) > 1:
            return None
        pos_check = self.month_position_check()
        return {'cs': 'MS', 'bs': 'BMS',
                'ce': 'M', 'be': 'BM'}.get(pos_check)

    def _get_wom_rule(self):
        #         wdiffs = unique(np.diff(self.index.week))
        # We also need -47, -49, -48 to catch index spanning year boundary
        #     if not lib.ismember(wdiffs, set([4, 5, -47, -49, -48])).all():
        #         return None

        weekdays = unique(self.index.weekday)
        if len(weekdays) > 1:
            return None

        week_of_months = unique((self.index.day - 1) // 7)
        # Only attempt to infer up to WOM-4. See #9425
        week_of_months = week_of_months[week_of_months < 4]
        if len(week_of_months) == 0 or len(week_of_months) > 1:
            return None

        # get which week
        week = week_of_months[0] + 1
        wd = _weekday_rule_aliases[weekdays[0]]

        return 'WOM-%d%s' % (week, wd)


class _TimedeltaFrequencyInferer(_FrequencyInferer):

    def _infer_daily_rule(self):
        if self.is_unique:
            days = self.deltas[0] / _ONE_DAY
            if days % 7 == 0:
                # Weekly
                alias = _weekday_rule_aliases[self.rep_stamp.weekday()]
                return _maybe_add_count('W-%s' % alias, days / 7)
            else:
                return _maybe_add_count('D', days)


def _maybe_add_count(base, count):
    if count != 1:
        return '%d%s' % (count, base)
    else:
        return base


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


_get_rule_month = tslib._get_rule_month


def _is_annual(rule):
    rule = rule.upper()
    return rule == 'A' or rule.startswith('A-')


def _quarter_months_conform(source, target):
    snum = _month_numbers[source]
    tnum = _month_numbers[target]
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


DAYS = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']

MONTHS = tslib._MONTHS
_month_numbers = tslib._MONTH_NUMBERS
_month_aliases = tslib._MONTH_ALIASES
_weekday_rule_aliases = dict((k, v) for k, v in enumerate(DAYS))


def _is_multiple(us, mult):
    return us % mult == 0
