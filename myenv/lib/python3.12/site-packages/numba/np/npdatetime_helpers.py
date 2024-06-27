"""
Helper functions for np.timedelta64 and np.datetime64.
For now, multiples-of-units (for example timedeltas expressed in tens
of seconds) are not supported.
"""


import numpy as np


DATETIME_UNITS = {
    'Y': 0,   # Years
    'M': 1,   # Months
    'W': 2,   # Weeks
    # Yes, there's a gap here
    'D': 4,   # Days
    'h': 5,   # Hours
    'm': 6,   # Minutes
    's': 7,   # Seconds
    'ms': 8,  # Milliseconds
    'us': 9,  # Microseconds
    'ns': 10,  # Nanoseconds
    'ps': 11,  # Picoseconds
    'fs': 12,  # Femtoseconds
    'as': 13,  # Attoseconds
    '': 14,   # "generic", i.e. unit-less
}

NAT = np.timedelta64('nat').astype(np.int64)

# NOTE: numpy has several inconsistent functions for timedelta casting:
# - can_cast_timedelta64_{metadata,units}() disallows "safe" casting
#   to and from generic units
# - cast_timedelta_to_timedelta() allows casting from (but not to)
#   generic units
# - compute_datetime_metadata_greatest_common_divisor() allows casting from
#   generic units (used for promotion)


def same_kind(src, dest):
    """
    Whether the *src* and *dest* units are of the same kind.
    """
    return (DATETIME_UNITS[src] < 5) == (DATETIME_UNITS[dest] < 5)


def can_cast_timedelta_units(src, dest):
    # Mimic NumPy's "safe" casting and promotion
    # `dest` must be more precise than `src` and they must be compatible
    # for conversion.
    # XXX should we switch to enforcing "same-kind" for Numpy 1.10+ ?
    src = DATETIME_UNITS[src]
    dest = DATETIME_UNITS[dest]
    if src == dest:
        return True
    if src == 14:
        return True
    if src > dest:
        return False
    if dest == 14:
        # unit-less timedelta64 is not compatible with anything else
        return False
    if src <= 1 and dest > 1:
        # Cannot convert between months or years and other units
        return False
    return True


# Exact conversion factors from one unit to the immediately more precise one
_factors = {
    0: (1, 12),   # Years -> Months
    2: (4, 7),    # Weeks -> Days
    4: (5, 24),   # Days -> Hours
    5: (6, 60),   # Hours -> Minutes
    6: (7, 60),   # Minutes -> Seconds
    7: (8, 1000),
    8: (9, 1000),
    9: (10, 1000),
    10: (11, 1000),
    11: (12, 1000),
    12: (13, 1000),
}


def _get_conversion_multiplier(big_unit_code, small_unit_code):
    """
    Return an integer multiplier allowing to convert from *big_unit_code*
    to *small_unit_code*.
    None is returned if the conversion is not possible through a
    simple integer multiplication.
    """
    # Mimics get_datetime_units_factor() in NumPy's datetime.c,
    # with a twist to allow no-op conversion from generic units.
    if big_unit_code == 14:
        return 1
    c = big_unit_code
    factor = 1
    while c < small_unit_code:
        try:
            c, mult = _factors[c]
        except KeyError:
            # No possible conversion
            return None
        factor *= mult
    if c == small_unit_code:
        return factor
    else:
        return None


def get_timedelta_conversion_factor(src_unit, dest_unit):
    """
    Return an integer multiplier allowing to convert from timedeltas
    of *src_unit* to *dest_unit*.
    """
    return _get_conversion_multiplier(DATETIME_UNITS[src_unit],
                                      DATETIME_UNITS[dest_unit])


def get_datetime_timedelta_conversion(datetime_unit, timedelta_unit):
    """
    Compute a possible conversion for combining *datetime_unit* and
    *timedelta_unit* (presumably for adding or subtracting).
    Return (result unit, integer datetime multiplier, integer timedelta
    multiplier). RuntimeError is raised if the combination is impossible.
    """
    # XXX now unused (I don't know where / how Numpy uses this)
    dt_unit_code = DATETIME_UNITS[datetime_unit]
    td_unit_code = DATETIME_UNITS[timedelta_unit]
    if td_unit_code == 14 or dt_unit_code == 14:
        return datetime_unit, 1, 1
    if td_unit_code < 2 and dt_unit_code >= 2:
        # Cannot combine Y or M timedelta64 with a finer-grained datetime64
        raise RuntimeError("cannot combine datetime64(%r) and timedelta64(%r)"
                           % (datetime_unit, timedelta_unit))
    dt_factor, td_factor = 1, 1

    # If years or months, the datetime unit is first scaled to weeks or days,
    # then conversion continues below.  This is the same algorithm as used
    # in Numpy's get_datetime_conversion_factor() (src/multiarray/datetime.c):
    # """Conversions between years/months and other units use
    # the factor averaged over the 400 year leap year cycle."""
    if dt_unit_code == 0:
        if td_unit_code >= 4:
            dt_factor = 97 + 400 * 365
            td_factor = 400
            dt_unit_code = 4
        elif td_unit_code == 2:
            dt_factor = 97 + 400 * 365
            td_factor = 400 * 7
            dt_unit_code = 2
    elif dt_unit_code == 1:
        if td_unit_code >= 4:
            dt_factor = 97 + 400 * 365
            td_factor = 400 * 12
            dt_unit_code = 4
        elif td_unit_code == 2:
            dt_factor = 97 + 400 * 365
            td_factor = 400 * 12 * 7
            dt_unit_code = 2

    if td_unit_code >= dt_unit_code:
        factor = _get_conversion_multiplier(dt_unit_code, td_unit_code)
        assert factor is not None, (dt_unit_code, td_unit_code)
        return timedelta_unit, dt_factor * factor, td_factor
    else:
        factor = _get_conversion_multiplier(td_unit_code, dt_unit_code)
        assert factor is not None, (dt_unit_code, td_unit_code)
        return datetime_unit, dt_factor, td_factor * factor


def combine_datetime_timedelta_units(datetime_unit, timedelta_unit):
    """
    Return the unit result of combining *datetime_unit* with *timedelta_unit*
    (e.g. by adding or subtracting).  None is returned if combining
    those units is forbidden.
    """
    dt_unit_code = DATETIME_UNITS[datetime_unit]
    td_unit_code = DATETIME_UNITS[timedelta_unit]
    if dt_unit_code == 14:
        return timedelta_unit
    elif td_unit_code == 14:
        return datetime_unit
    if td_unit_code < 2 and dt_unit_code >= 2:
        return None
    if dt_unit_code > td_unit_code:
        return datetime_unit
    else:
        return timedelta_unit


def get_best_unit(unit_a, unit_b):
    """
    Get the best (i.e. finer-grained) of two units.
    """
    a = DATETIME_UNITS[unit_a]
    b = DATETIME_UNITS[unit_b]
    if a == 14:
        return unit_b
    if b == 14:
        return unit_a
    if b > a:
        return unit_b
    return unit_a


def datetime_minimum(a, b):
    pass


def datetime_maximum(a, b):
    pass
