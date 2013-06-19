from datetime import datetime
import re

import numpy as np

from pandas.core.algorithms import unique
from pandas.tseries.offsets import DateOffset
from pandas.util.decorators import cache_readonly
import pandas.tseries.offsets as offsets
import pandas.core.common as com
import pandas.lib as lib
import pandas.tslib as tslib


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


class Resolution(object):

    RESO_US = 0
    RESO_SEC = 1
    RESO_MIN = 2
    RESO_HR = 3
    RESO_DAY = 4

    @classmethod
    def get_str(cls, reso):
        return {cls.RESO_US: 'microsecond',
                cls.RESO_SEC: 'second',
                cls.RESO_MIN: 'minute',
                cls.RESO_HR: 'hour',
                cls.RESO_DAY: 'day'}.get(reso, 'day')


def get_reso_string(reso):
    return Resolution.get_str(reso)


def get_to_timestamp_base(base):
    if base < FreqGroup.FR_BUS:
        return FreqGroup.FR_DAY
    if FreqGroup.FR_HR <= base <= FreqGroup.FR_SEC:
        return FreqGroup.FR_SEC
    return base


def get_freq_group(freq):
    if isinstance(freq, basestring):
        base, mult = get_freq_code(freq)
        freq = base
    return (freq // 1000) * 1000


def get_freq(freq):
    if isinstance(freq, basestring):
        base, mult = get_freq_code(freq)
        freq = base
    return freq


def get_freq_code(freqstr):
    """

    Parameters
    ----------

    Returns
    -------
    """
    if isinstance(freqstr, DateOffset):
        freqstr = (get_offset_name(freqstr), freqstr.n)

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


#----------------------------------------------------------------------
# Offset names ("time rules") and related functions


from pandas.tseries.offsets import (Micro, Milli, Second, Minute, Hour,
                                    Day, BDay, CDay, Week, MonthBegin,
                                    MonthEnd, BMonthBegin, BMonthEnd,
                                    QuarterBegin, QuarterEnd, BQuarterBegin,
                                    BQuarterEnd, YearBegin, YearEnd,
                                    BYearBegin, BYearEnd,
                                    )
try:
    cday = CDay()
except NotImplementedError:
    cday = None

_offset_map = {
    'D': Day(),
    'C': cday,
    'B': BDay(),
    'H': Hour(),
    'T': Minute(),
    'S': Second(),
    'L': Milli(),
    'U': Micro(),
    None: None,

    # Monthly - Calendar
    'M': MonthEnd(),
    'MS': MonthBegin(),

    # Monthly - Business
    'BM': BMonthEnd(),
    'BMS': BMonthBegin(),

    # Annual - Calendar
    'A-JAN': YearEnd(month=1),
    'A-FEB': YearEnd(month=2),
    'A-MAR': YearEnd(month=3),
    'A-APR': YearEnd(month=4),
    'A-MAY': YearEnd(month=5),
    'A-JUN': YearEnd(month=6),
    'A-JUL': YearEnd(month=7),
    'A-AUG': YearEnd(month=8),
    'A-SEP': YearEnd(month=9),
    'A-OCT': YearEnd(month=10),
    'A-NOV': YearEnd(month=11),
    'A-DEC': YearEnd(month=12),

    # Annual - Calendar (start)
    'AS-JAN': YearBegin(month=1),
    'AS-FEB': YearBegin(month=2),
    'AS-MAR': YearBegin(month=3),
    'AS-APR': YearBegin(month=4),
    'AS-MAY': YearBegin(month=5),
    'AS-JUN': YearBegin(month=6),
    'AS-JUL': YearBegin(month=7),
    'AS-AUG': YearBegin(month=8),
    'AS-SEP': YearBegin(month=9),
    'AS-OCT': YearBegin(month=10),
    'AS-NOV': YearBegin(month=11),
    'AS-DEC': YearBegin(month=12),

    # Annual - Business
    'BA-JAN': BYearEnd(month=1),
    'BA-FEB': BYearEnd(month=2),
    'BA-MAR': BYearEnd(month=3),
    'BA-APR': BYearEnd(month=4),
    'BA-MAY': BYearEnd(month=5),
    'BA-JUN': BYearEnd(month=6),
    'BA-JUL': BYearEnd(month=7),
    'BA-AUG': BYearEnd(month=8),
    'BA-SEP': BYearEnd(month=9),
    'BA-OCT': BYearEnd(month=10),
    'BA-NOV': BYearEnd(month=11),
    'BA-DEC': BYearEnd(month=12),

    # Annual - Business (Start)
    'BAS-JAN': BYearBegin(month=1),
    'BAS-FEB': BYearBegin(month=2),
    'BAS-MAR': BYearBegin(month=3),
    'BAS-APR': BYearBegin(month=4),
    'BAS-MAY': BYearBegin(month=5),
    'BAS-JUN': BYearBegin(month=6),
    'BAS-JUL': BYearBegin(month=7),
    'BAS-AUG': BYearBegin(month=8),
    'BAS-SEP': BYearBegin(month=9),
    'BAS-OCT': BYearBegin(month=10),
    'BAS-NOV': BYearBegin(month=11),
    'BAS-DEC': BYearBegin(month=12),

    # Quarterly - Calendar
    # 'Q'     : QuarterEnd(startingMonth=3),
    'Q-JAN': QuarterEnd(startingMonth=1),
    'Q-FEB': QuarterEnd(startingMonth=2),
    'Q-MAR': QuarterEnd(startingMonth=3),
    'Q-APR': QuarterEnd(startingMonth=4),
    'Q-MAY': QuarterEnd(startingMonth=5),
    'Q-JUN': QuarterEnd(startingMonth=6),
    'Q-JUL': QuarterEnd(startingMonth=7),
    'Q-AUG': QuarterEnd(startingMonth=8),
    'Q-SEP': QuarterEnd(startingMonth=9),
    'Q-OCT': QuarterEnd(startingMonth=10),
    'Q-NOV': QuarterEnd(startingMonth=11),
    'Q-DEC': QuarterEnd(startingMonth=12),

    # Quarterly - Calendar (Start)
    'QS': QuarterBegin(startingMonth=1),
    'QS-JAN': QuarterBegin(startingMonth=1),
    'QS-FEB': QuarterBegin(startingMonth=2),
    'QS-MAR': QuarterBegin(startingMonth=3),
    'QS-APR': QuarterBegin(startingMonth=4),
    'QS-MAY': QuarterBegin(startingMonth=5),
    'QS-JUN': QuarterBegin(startingMonth=6),
    'QS-JUL': QuarterBegin(startingMonth=7),
    'QS-AUG': QuarterBegin(startingMonth=8),
    'QS-SEP': QuarterBegin(startingMonth=9),
    'QS-OCT': QuarterBegin(startingMonth=10),
    'QS-NOV': QuarterBegin(startingMonth=11),
    'QS-DEC': QuarterBegin(startingMonth=12),

    # Quarterly - Business
    'BQ-JAN': BQuarterEnd(startingMonth=1),
    'BQ-FEB': BQuarterEnd(startingMonth=2),
    'BQ-MAR': BQuarterEnd(startingMonth=3),

    'BQ': BQuarterEnd(startingMonth=12),
    'BQ-APR': BQuarterEnd(startingMonth=4),
    'BQ-MAY': BQuarterEnd(startingMonth=5),
    'BQ-JUN': BQuarterEnd(startingMonth=6),
    'BQ-JUL': BQuarterEnd(startingMonth=7),
    'BQ-AUG': BQuarterEnd(startingMonth=8),
    'BQ-SEP': BQuarterEnd(startingMonth=9),
    'BQ-OCT': BQuarterEnd(startingMonth=10),
    'BQ-NOV': BQuarterEnd(startingMonth=11),
    'BQ-DEC': BQuarterEnd(startingMonth=12),

    # Quarterly - Business (Start)
    'BQS-JAN': BQuarterBegin(startingMonth=1),
    'BQS': BQuarterBegin(startingMonth=1),
    'BQS-FEB': BQuarterBegin(startingMonth=2),
    'BQS-MAR': BQuarterBegin(startingMonth=3),
    'BQS-APR': BQuarterBegin(startingMonth=4),
    'BQS-MAY': BQuarterBegin(startingMonth=5),
    'BQS-JUN': BQuarterBegin(startingMonth=6),
    'BQS-JUL': BQuarterBegin(startingMonth=7),
    'BQS-AUG': BQuarterBegin(startingMonth=8),
    'BQS-SEP': BQuarterBegin(startingMonth=9),
    'BQS-OCT': BQuarterBegin(startingMonth=10),
    'BQS-NOV': BQuarterBegin(startingMonth=11),
    'BQS-DEC': BQuarterBegin(startingMonth=12),

    # Weekly
    'W-MON': Week(weekday=0),
    'W-TUE': Week(weekday=1),
    'W-WED': Week(weekday=2),
    'W-THU': Week(weekday=3),
    'W-FRI': Week(weekday=4),
    'W-SAT': Week(weekday=5),
    'W-SUN': Week(weekday=6),

}

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
    'H': 'H',
    'Q': 'Q',
    'A': 'A',
    'W': 'W',
    'M': 'M'
}

need_suffix = ['QS', 'BQ', 'BQS', 'AS', 'BA', 'BAS']
_months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP',
           'OCT', 'NOV', 'DEC']
for __prefix in need_suffix:
    for _m in _months:
        _offset_to_period_map['%s-%s' % (__prefix, _m)] = \
            _offset_to_period_map[__prefix]
for __prefix in ['A', 'Q']:
    for _m in _months:
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
    'W': 'W-SUN',

    'Q@JAN': 'BQ-JAN',
    'Q@FEB': 'BQ-FEB',
    'Q@MAR': 'BQ-MAR',
    'Q': 'Q-DEC',

    'A': 'A-DEC',  # YearEnd(month=12),
    'AS': 'AS-JAN',  # YearBegin(month=1),
    'BA': 'BA-DEC',  # BYearEnd(month=12),
    'BAS': 'BAS-JAN',  # BYearBegin(month=1),

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

    # lite aliases
    'Min': 'T',
    'min': 'T',
    'ms': 'L',
    'us': 'U'
}

for _i, _weekday in enumerate(['MON', 'TUE', 'WED', 'THU', 'FRI']):
    for _iweek in xrange(4):
        _name = 'WOM-%d%s' % (_iweek + 1, _weekday)
        _offset_map[_name] = offsets.WeekOfMonth(week=_iweek, weekday=_i)
        _rule_aliases[_name.replace('-', '@')] = _name

# Note that _rule_aliases is not 1:1 (d[BA]==d[A@DEC]), and so traversal
# order matters when constructing an inverse. we pick one. #2331
_legacy_reverse_map = dict((v, k) for k, v in
                           reversed(sorted(_rule_aliases.iteritems())))

# for helping out with pretty-printing and name-lookups

_offset_names = {}
for name, offset in _offset_map.iteritems():
    if offset is None:
        continue
    offset.name = name
    _offset_names[offset] = name


def inferTimeRule(index):
    from pandas.tseries.index import DatetimeIndex
    import warnings
    warnings.warn("This method is deprecated, use infer_freq or inferred_freq"
                  " attribute of DatetimeIndex", FutureWarning)

    freq = DatetimeIndex(index).inferred_freq
    if freq is None:
        raise Exception('Unable to infer time rule')

    offset = to_offset(freq)
    return get_legacy_offset_name(offset)


def to_offset(freqstr):
    """
    Return DateOffset object from string representation

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
        if isinstance(stride, basestring):
            name, stride = stride, name
        name, _ = _base_and_stride(name)
        delta = get_offset(name) * stride
    else:
        delta = None
        stride_sign = None
        try:
            for stride, name, _ in opattern.findall(freqstr):
                offset = get_offset(name)
                if not stride:
                    stride = 1
                stride = int(stride)
                if stride_sign is None:
                    stride_sign = np.sign(stride)
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
opattern = re.compile(r'([\-]?\d*)\s*([A-Za-z]+([\-@]\d*[A-Za-z]+)?)')


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

_dont_uppercase = ['MS', 'ms']


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
            name = _rule_aliases[name]
        elif name.lower() in _rule_aliases:
            name = _rule_aliases[name.lower()]
    else:
        if name in _rule_aliases:
            name = _rule_aliases[name]

    offset = _offset_map.get(name)

    if offset is not None:
        return offset
    else:
        raise ValueError('Bad rule name requested: %s.' % name)


getOffset = get_offset


def hasOffsetName(offset):
    return offset in _offset_names


def get_offset_name(offset):
    """
    Return rule name associated with a DateOffset object

    Examples
    --------
    get_offset_name(BMonthEnd(1)) --> 'EOM'
    """
    name = _offset_names.get(offset)

    if name is not None:
        return name
    else:
        raise ValueError('Bad rule given: %s.' % offset)


def get_legacy_offset_name(offset):
    """
    Return the pre pandas 0.8.0 name for the date offset
    """
    name = _offset_names.get(offset)
    return _legacy_reverse_map.get(name, name)

get_offset_name = get_offset_name


def get_standard_freq(freq):
    """
    Return the standardized frequency string
    """
    if freq is None:
        return None

    if isinstance(freq, DateOffset):
        return get_offset_name(freq)

    code, stride = get_freq_code(freq)
    return _get_freq_str(code, stride)

#----------------------------------------------------------------------
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
}

_reverse_period_code_map = {}
for _k, _v in _period_code_map.iteritems():
    _reverse_period_code_map[_v] = _k

# Additional aliases
_period_code_map.update({
    "Q": 2000,  # Quarterly - December year end (default quarterly)
    "A": 1000,  # Annual
    "W": 4000,  # Weekly
})


def _period_alias_dictionary():
    """
    Build freq alias dictionary to support freqs from original c_dates.c file
    of the scikits.timeseries library.
    """
    alias_dict = {}

    M_aliases = ["M", "MTH", "MONTH", "MONTHLY"]
    B_aliases = ["B", "BUS", "BUSINESS", "BUSINESSLY", 'WEEKDAY']
    D_aliases = ["D", "DAY", "DLY", "DAILY"]
    H_aliases = ["H", "HR", "HOUR", "HRLY", "HOURLY"]
    T_aliases = ["T", "MIN", "MINUTE", "MINUTELY"]
    S_aliases = ["S", "SEC", "SECOND", "SECONDLY"]

    for k in M_aliases:
        alias_dict[k] = 'M'

    for k in B_aliases:
        alias_dict[k] = 'B'

    for k in D_aliases:
        alias_dict[k] = 'D'

    for k in H_aliases:
        alias_dict[k] = 'H'

    for k in T_aliases:
        alias_dict[k] = 'Min'

    for k in S_aliases:
        alias_dict[k] = 'S'

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

_reso_period_map = {
    "year": "A",
    "quarter": "Q",
    "month": "M",
    "day": "D",
    "hour": "H",
    "minute": "T",
    "second": "S",
}


def _infer_period_group(freqstr):
    return _period_group(_reso_period_map[freqstr])


def _period_group(freqstr):
    base, mult = get_freq_code(freqstr)
    return base // 1000 * 1000

_period_alias_dict = _period_alias_dictionary()


def _period_str_to_code(freqstr):
    # hack
    freqstr = _rule_aliases.get(freqstr, freqstr)
    freqstr = _rule_aliases.get(freqstr.lower(), freqstr)

    try:
        freqstr = freqstr.upper()
        return _period_code_map[freqstr]
    except:
        alias = _period_alias_dict[freqstr]
        return _period_code_map[alias]


def infer_freq(index, warn=True):
    """
    Infer the most likely frequency given the input index. If the frequency is
    uncertain, a warning will be printed

    Parameters
    ----------
    index : DatetimeIndex
    warn : boolean, default True

    Returns
    -------
    freq : string or None
        None if no discernible frequency
    """
    from pandas.tseries.index import DatetimeIndex

    if not isinstance(index, DatetimeIndex):
        index = DatetimeIndex(index)

    inferer = _FrequencyInferer(index, warn=warn)
    return inferer.get_freq()

_ONE_MICRO = 1000L
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
        self.warn = warn

        if len(index) < 3:
            raise ValueError('Need at least 3 dates to infer frequency')

        self.is_monotonic = self.index.is_monotonic

    @cache_readonly
    def deltas(self):
        return tslib.unique_deltas(self.values)

    @cache_readonly
    def is_unique(self):
        return len(self.deltas) == 1

    def get_freq(self):
        if not self.is_monotonic or not self.index.is_unique:
            return None

        delta = self.deltas[0]
        if _is_multiple(delta, _ONE_DAY):
            return self._infer_daily_rule()
        else:
            # Possibly intraday frequency
            if not self.is_unique:
                return None
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
            wd = datetime(y, m, d).weekday()

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
            return monthly_rule

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
        wdiffs = unique(np.diff(self.index.week))
        if not lib.ismember(wdiffs, set([4, 5])).all():
            return None

        weekdays = unique(self.index.weekday)
        if len(weekdays) > 1:
            return None

        # get which week
        week = (self.index[0].day - 1) // 7 + 1
        wd = _weekday_rule_aliases[weekdays[0]]

        return 'WOM-%d%s' % (week, wd)

import pandas.core.algorithms as algos


def _maybe_add_count(base, count):
    if count > 1:
        return '%d%s' % (count, base)
    else:
        return base


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
    if isinstance(source, offsets.DateOffset):
        source = source.rule_code

    if isinstance(target, offsets.DateOffset):
        target = target.rule_code

    target = target.upper()
    source = source.upper()
    if _is_annual(target):
        if _is_quarterly(source):
            return _quarter_months_conform(_get_rule_month(source),
                                           _get_rule_month(target))
        return source in ['D', 'C', 'B', 'M', 'H', 'T', 'S']
    elif _is_quarterly(target):
        return source in ['D', 'C', 'B', 'M', 'H', 'T', 'S']
    elif target == 'M':
        return source in ['D', 'C', 'B', 'H', 'T', 'S']
    elif _is_weekly(target):
        return source in [target, 'D', 'C', 'B', 'H', 'T', 'S']
    elif target == 'B':
        return source in ['B', 'H', 'T', 'S']
    elif target == 'C':
        return source in ['C', 'H', 'T', 'S']
    elif target == 'D':
        return source in ['D', 'H', 'T', 'S']
    elif target == 'H':
        return source in ['H', 'T', 'S']
    elif target == 'T':
        return source in ['T', 'S']
    elif target == 'S':
        return source in ['S']


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
    if isinstance(source, offsets.DateOffset):
        source = source.rule_code

    if isinstance(target, offsets.DateOffset):
        target = target.rule_code

    target = target.upper()
    source = source.upper()
    if _is_annual(source):
        if _is_annual(target):
            return _get_rule_month(source) == _get_rule_month(target)

        if _is_quarterly(target):
            smonth = _get_rule_month(source)
            tmonth = _get_rule_month(target)
            return _quarter_months_conform(smonth, tmonth)
        return target in ['D', 'C', 'B', 'M', 'H', 'T', 'S']
    elif _is_quarterly(source):
        return target in ['D', 'C', 'B', 'M', 'H', 'T', 'S']
    elif source == 'M':
        return target in ['D', 'C', 'B', 'H', 'T', 'S']
    elif _is_weekly(source):
        return target in [source, 'D', 'C', 'B', 'H', 'T', 'S']
    elif source == 'B':
        return target in ['D', 'C', 'B', 'H', 'T', 'S']
    elif source == 'C':
        return target in ['D', 'C', 'B', 'H', 'T', 'S']
    elif source == 'D':
        return target in ['D', 'C', 'B', 'H', 'T', 'S']
    elif source == 'H':
        return target in ['H', 'T', 'S']
    elif source == 'T':
        return target in ['T', 'S']
    elif source == 'S':
        return target in ['S']


def _get_rule_month(source, default='DEC'):
    source = source.upper()
    if '-' not in source:
        return default
    else:
        return source.split('-')[1]


def _is_annual(rule):
    rule = rule.upper()
    return rule == 'A' or rule.startswith('A-')


def _quarter_months_conform(source, target):
    snum = _month_numbers[source]
    tnum = _month_numbers[target]
    return snum % 3 == tnum % 3


def _is_quarterly(rule):
    rule = rule.upper()
    return rule == 'Q' or rule.startswith('Q-')


def _is_weekly(rule):
    rule = rule.upper()
    return rule == 'W' or rule.startswith('W-')


DAYS = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']

MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
          'AUG', 'SEP', 'OCT', 'NOV', 'DEC']

_month_numbers = dict((k, i) for i, k in enumerate(MONTHS))


_weekday_rule_aliases = dict((k, v) for k, v in enumerate(DAYS))
_month_aliases = dict((k + 1, v) for k, v in enumerate(MONTHS))


def _is_multiple(us, mult):
    return us % mult == 0
