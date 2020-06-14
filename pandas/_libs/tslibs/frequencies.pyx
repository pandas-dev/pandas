
from .dtypes import FreqGroup

# ----------------------------------------------------------------------


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
