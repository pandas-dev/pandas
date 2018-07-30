# -*- coding: utf-8 -*-
# cython: profile=False
import re

cimport numpy as cnp
cnp.import_array()

from util cimport is_integer_object, is_string_object

from ccalendar import MONTH_NUMBERS

# ----------------------------------------------------------------------
# Constants

# hack to handle WOM-1MON
opattern = re.compile(
    r'([+\-]?\d*|[+\-]?\d*\.\d*)\s*([A-Za-z]+([\-][\dA-Za-z\-]+)?)'
)

INVALID_FREQ_ERR_MSG = "Invalid frequency: {0}"

# ---------------------------------------------------------------------
# Period codes


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
    "N": 12000}       # Nanosecondly


_reverse_period_code_map = {
    _period_code_map[key]: key for key in _period_code_map}

# Yearly aliases; careful not to put these in _reverse_period_code_map
_period_code_map.update({'Y' + key[1:]: _period_code_map[key]
                         for key in _period_code_map
                         if key.startswith('A-')})

_period_code_map.update({
    "Q": 2000,   # Quarterly - December year end (default quarterly)
    "A": 1000,   # Annual
    "W": 4000,   # Weekly
    "C": 5000})  # Custom Business Day

_lite_rule_alias = {
    'W': 'W-SUN',
    'Q': 'Q-DEC',

    'A': 'A-DEC',      # YearEnd(month=12),
    'Y': 'A-DEC',
    'AS': 'AS-JAN',    # YearBegin(month=1),
    'YS': 'AS-JAN',
    'BA': 'BA-DEC',    # BYearEnd(month=12),
    'BY': 'BA-DEC',
    'BAS': 'BAS-JAN',  # BYearBegin(month=1),
    'BYS': 'BAS-JAN',

    'Min': 'T',
    'min': 'T',
    'ms': 'L',
    'us': 'U',
    'ns': 'N'}

_dont_uppercase = set(('MS', 'ms'))

# ----------------------------------------------------------------------

cpdef get_freq_code(freqstr):
    """
    Return freq str or tuple to freq code and stride (mult)

    Parameters
    ----------
    freqstr : str or tuple

    Returns
    -------
    return : tuple of base frequency code and stride (mult)

    Examples
    --------
    >>> get_freq_code('3D')
    (6000, 3)

    >>> get_freq_code('D')
    (6000, 1)

    >>> get_freq_code(('D', 3))
    (6000, 3)
    """
    if getattr(freqstr, '_typ', None) == 'dateoffset':
        freqstr = (freqstr.rule_code, freqstr.n)

    if isinstance(freqstr, tuple):
        if (is_integer_object(freqstr[0]) and
                is_integer_object(freqstr[1])):
            # e.g., freqstr = (2000, 1)
            return freqstr
        else:
            # e.g., freqstr = ('T', 5)
            try:
                code = _period_str_to_code(freqstr[0])
                stride = freqstr[1]
            except:
                if is_integer_object(freqstr[1]):
                    raise
                code = _period_str_to_code(freqstr[1])
                stride = freqstr[0]
            return code, stride

    if is_integer_object(freqstr):
        return (freqstr, 1)

    base, stride = _base_and_stride(freqstr)
    code = _period_str_to_code(base)

    return code, stride


cpdef _base_and_stride(freqstr):
    """
    Return base freq and stride info from string representation

    Examples
    --------
    _freq_and_stride('5Min') -> 'Min', 5
    """
    groups = opattern.match(freqstr)

    if not groups:
        raise ValueError("Could not evaluate {freq}".format(freq=freqstr))

    stride = groups.group(1)

    if len(stride):
        stride = int(stride)
    else:
        stride = 1

    base = groups.group(2)

    return (base, stride)


cpdef _period_str_to_code(freqstr):
    freqstr = _lite_rule_alias.get(freqstr, freqstr)

    if freqstr not in _dont_uppercase:
        lower = freqstr.lower()
        freqstr = _lite_rule_alias.get(lower, freqstr)

    if freqstr not in _dont_uppercase:
        freqstr = freqstr.upper()
    try:
        return _period_code_map[freqstr]
    except KeyError:
        raise ValueError(INVALID_FREQ_ERR_MSG.format(freqstr))


cpdef str get_freq_str(base, mult=1):
    """
    Return the summary string associated with this offset code, possibly
    adjusted by a multiplier.

    Parameters
    ----------
    base : int (member of FreqGroup)

    Returns
    -------
    freq_str : str

    Examples
    --------
    >>> get_freq_str(1000)
    'A-DEC'

    >>> get_freq_str(2000, 2)
    '2Q-DEC'

    >>> get_freq_str("foo")
    """
    code = _reverse_period_code_map.get(base)
    if mult == 1:
        return code
    return str(mult) + code


cpdef str get_base_alias(freqstr):
    """
    Returns the base frequency alias, e.g., '5D' -> 'D'

    Parameters
    ----------
    freqstr : str

    Returns
    -------
    base_alias : str
    """
    return _base_and_stride(freqstr)[0]


cpdef int get_to_timestamp_base(int base):
    """
    Return frequency code group used for base of to_timestamp against
    frequency code.

    Parameters
    ----------
    base : int (member of FreqGroup)

    Returns
    -------
    base : int

    Examples
    --------
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
    elif FreqGroup.FR_HR <= base <= FreqGroup.FR_SEC:
        return FreqGroup.FR_SEC
    return base


cpdef object get_freq(object freq):
    """
    Return frequency code of given frequency str.
    If input is not string, return input as it is.

    Examples
    --------
    >>> get_freq('A')
    1000

    >>> get_freq('3A')
    1000
    """
    if is_string_object(freq):
        base, mult = get_freq_code(freq)
        freq = base
    return freq


# ----------------------------------------------------------------------
# Frequency comparison

cpdef bint is_subperiod(source, target):
    """
    Returns True if downsampling is possible between source and target
    frequencies

    Parameters
    ----------
    source : string or DateOffset
        Frequency converting from
    target : string or DateOffset
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
            return _quarter_months_conform(get_rule_month(source),
                                           get_rule_month(target))
        return source in {'D', 'C', 'B', 'M', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif _is_quarterly(target):
        return source in {'D', 'C', 'B', 'M', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif _is_monthly(target):
        return source in {'D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif _is_weekly(target):
        return source in {target, 'D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif target == 'B':
        return source in {'B', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif target == 'C':
        return source in {'C', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif target == 'D':
        return source in {'D', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif target == 'H':
        return source in {'H', 'T', 'S', 'L', 'U', 'N'}
    elif target == 'T':
        return source in {'T', 'S', 'L', 'U', 'N'}
    elif target == 'S':
        return source in {'S', 'L', 'U', 'N'}
    elif target == 'L':
        return source in {'L', 'U', 'N'}
    elif target == 'U':
        return source in {'U', 'N'}
    elif target == 'N':
        return source in {'N'}


cpdef bint is_superperiod(source, target):
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
            return get_rule_month(source) == get_rule_month(target)

        if _is_quarterly(target):
            smonth = get_rule_month(source)
            tmonth = get_rule_month(target)
            return _quarter_months_conform(smonth, tmonth)
        return target in {'D', 'C', 'B', 'M', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif _is_quarterly(source):
        return target in {'D', 'C', 'B', 'M', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif _is_monthly(source):
        return target in {'D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif _is_weekly(source):
        return target in {source, 'D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif source == 'B':
        return target in {'D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif source == 'C':
        return target in {'D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif source == 'D':
        return target in {'D', 'C', 'B', 'H', 'T', 'S', 'L', 'U', 'N'}
    elif source == 'H':
        return target in {'H', 'T', 'S', 'L', 'U', 'N'}
    elif source == 'T':
        return target in {'T', 'S', 'L', 'U', 'N'}
    elif source == 'S':
        return target in {'S', 'L', 'U', 'N'}
    elif source == 'L':
        return target in {'L', 'U', 'N'}
    elif source == 'U':
        return target in {'U', 'N'}
    elif source == 'N':
        return target in {'N'}


cdef str _maybe_coerce_freq(code):
    """ we might need to coerce a code to a rule_code
    and uppercase it

    Parameters
    ----------
    source : string or DateOffset
        Frequency converting from

    Returns
    -------
    code : string
    """
    assert code is not None
    if getattr(code, '_typ', None) == 'dateoffset':
        # i.e. isinstance(code, ABCDateOffset):
        code = code.rule_code
    return code.upper()


cdef bint _quarter_months_conform(str source, str target):
    snum = MONTH_NUMBERS[source]
    tnum = MONTH_NUMBERS[target]
    return snum % 3 == tnum % 3


cdef bint _is_annual(str rule):
    rule = rule.upper()
    return rule == 'A' or rule.startswith('A-')


cdef bint _is_quarterly(str rule):
    rule = rule.upper()
    return rule == 'Q' or rule.startswith('Q-') or rule.startswith('BQ')


cdef bint _is_monthly(str rule):
    rule = rule.upper()
    return rule == 'M' or rule == 'BM'


cdef bint _is_weekly(str rule):
    rule = rule.upper()
    return rule == 'W' or rule.startswith('W-')


# ----------------------------------------------------------------------

cpdef object get_rule_month(object source, object default='DEC'):
    """
    Return starting month of given freq, default is December.

    Parameters
    ----------
    source : object
    default : object (default "DEC")

    Returns
    -------
    rule_month: object (usually string)

    Examples
    --------
    >>> get_rule_month('D')
    'DEC'

    >>> get_rule_month('A-JAN')
    'JAN'
    """
    if hasattr(source, 'freqstr'):
        source = source.freqstr
    source = source.upper()
    if '-' not in source:
        return default
    else:
        return source.split('-')[1]
