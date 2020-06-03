cimport numpy as cnp
cnp.import_array()

from pandas._libs.tslibs.util cimport is_integer_object

from pandas._libs.tslibs.offsets cimport is_offset_object
from pandas._libs.tslibs.offsets import (
    INVALID_FREQ_ERR_MSG,
    _dont_uppercase,
    _lite_rule_alias,
    base_and_stride,
    opattern,
)

# ---------------------------------------------------------------------
# Period codes


class FreqGroup:
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

# Map attribute-name resolutions to resolution abbreviations
_attrname_to_abbrevs = {
    "year": "A",
    "quarter": "Q",
    "month": "M",
    "day": "D",
    "hour": "H",
    "minute": "T",
    "second": "S",
    "millisecond": "L",
    "microsecond": "U",
    "nanosecond": "N",
}
cdef dict attrname_to_abbrevs = _attrname_to_abbrevs


# ----------------------------------------------------------------------

def get_freq_group(freq) -> int:
    """
    Return frequency code group of given frequency str or offset.

    Examples
    --------
    >>> get_freq_group('W-MON')
    4000

    >>> get_freq_group('W-FRI')
    4000
    """
    if is_offset_object(freq):
        freq = freq.rule_code

    if isinstance(freq, str):
        freq = attrname_to_abbrevs.get(freq, freq)
        base, mult = get_freq_code(freq)
        freq = base
    elif isinstance(freq, int):
        pass
    else:
        raise ValueError('input must be str, offset or int')
    return (freq // 1000) * 1000


cpdef get_freq_code(freqstr):
    """
    Return freq str or tuple to freq code and stride (mult)

    Parameters
    ----------
    freqstr : str or tuple

    Returns
    -------
    return : tuple of base frequency code and stride (mult)

    Raises
    ------
    TypeError : if passed a tuple witth incorrect types

    Examples
    --------
    >>> get_freq_code('3D')
    (6000, 3)

    >>> get_freq_code('D')
    (6000, 1)

    >>> get_freq_code(('D', 3))
    (6000, 3)
    """
    if is_offset_object(freqstr):
        freqstr = (freqstr.rule_code, freqstr.n)

    if isinstance(freqstr, tuple):
        if is_integer_object(freqstr[0]) and is_integer_object(freqstr[1]):
            # e.g., freqstr = (2000, 1)
            return freqstr
        elif is_integer_object(freqstr[0]):
            # Note: passing freqstr[1] below will raise TypeError if that
            #  is not a str
            code = _period_str_to_code(freqstr[1])
            stride = freqstr[0]
            return code, stride
        else:
            # e.g., freqstr = ('T', 5)
            code = _period_str_to_code(freqstr[0])
            stride = freqstr[1]
            return code, stride

    if is_integer_object(freqstr):
        return freqstr, 1

    base, stride = base_and_stride(freqstr)
    code = _period_str_to_code(base)

    return code, stride


cpdef _period_str_to_code(str freqstr):
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
