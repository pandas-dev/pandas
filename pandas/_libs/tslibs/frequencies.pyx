# -*- coding: utf-8 -*-
# cython: profile=False
import re

cimport cython

import numpy as np
cimport numpy as np
np.import_array()

from util cimport is_integer_object

# ----------------------------------------------------------------------
# Constants

# hack to handle WOM-1MON
opattern = re.compile(
    r'([+\-]?\d*|[+\-]?\d*\.\d*)\s*([A-Za-z]+([\-][\dA-Za-z\-]+)?)'
)

_INVALID_FREQ_ERROR = "Invalid frequency: {0}"

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

    Example
    -------
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
        raise ValueError(_INVALID_FREQ_ERROR.format(freqstr))
