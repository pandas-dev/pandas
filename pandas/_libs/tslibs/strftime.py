"""Strftime-related classes and functions.
"""
class UnsupportedDatetimeDirective(ValueError):
    """The strftime format contains a directive that is not supported in this context."""


def convert_dtformat(
    strftime_fmt: str,
    new_style_fmt: bool = False,
) -> str:
    """Convert a strftime formatting string into a string formatting template string.

    Parameters
    ----------
    strftime_fmt : str
        The strftime format string specification, e.g. `"%Y-%m-%d %H:%M:%S"`. Note
        that not all directives are eligible to successful usage of string formatting.
        Unsupported directives will lead to an `UnsupportedDatetimeDirective` being raised.
    new_style_fmt : bool, default: False
        Whether the output string should be new-style
        e.g. "{year}-{month:02d}-{day:02d} {hour:02d}:{min:02d}:{sec:02d}"
        or old-style e.g. "%(year)s-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d"

    Returns
    -------
    fmt_out : str
        A string that may be used to format a `datetime` variable. The style of this string
        is either old-style or new-style depending on `new_style_formatting`.
        For old-style, it may be used as `fmt_out % get_datetime_fmt_dct(dt)`.
        For new-style, it may be used as `fmt_out.format(**get_datetime_fmt_dct(dt))`

    Raises
    ------
    UnsupportedDatetimeDirective
        Raised when the received `strftime_fmt` format contains a directive for which the output
        can not currently be created using string formatting.

    See Also
    --------
    - `strftime format codes reference <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes>`_
    - `Stackoverflow post <https://stackoverflow.com/a/43495629/7262247>`_ explaining
       how old-style formatting is faster than new-style formatting, itself faster than
       `datetime.strftime`.
    """
    non_supported = (
        # All of these below are names and therefore are not in the numpy or datetime attr representation
        "%a",  # Weekday as locale’s abbreviated name.
        "%A",  # Weekday as locale’s full name.
        "%w",  # Weekday as a decimal number, where 0 is Sunday and 6 is Saturday.
        "%b",  # Month as locale’s abbreviated name.
        "%B",  # Month as locale’s full name.
        # TODO All of those below can probably be derived easily from the numbers on the instance
        "%y",  # Year without century as a zero-padded decimal number.  >> year % 100
        "%I",  # Hour (12-hour clock) as a zero-padded decimal number. >> hour % 12
        "%p",  # Locale’s equivalent of either AM or PM.  >> "pm" if (hour // 12) else "am"
        # TODO Below Time offset and timezone information ... but may be hard
        "%z",  # UTC offset in the form ±HHMM[SS[.ffffff]] (empty string if the object is naive).
        "%Z",  # Time zone name (empty string if the object is naive).
        # We do not want to enter into these below, we do not want to re-create the datetime implementation
        "%j",  # Day of the year as a zero-padded decimal number.
        "%U",  # Week number of the year (Sunday as the first day of the week) as a zero-padded decimal number. All days in a new year preceding the first Sunday are considered to be in week 0.
        "%W",  # Week number of the year (Monday as the first day of the week) as a zero-padded decimal number. All days in a new year preceding the first Monday are considered to be in week 0.
        "%c",  # Locale’s appropriate date and time representation.
        "%x",  # Locale’s appropriate date representation.
        "%X",  # Locale’s appropriate time representation.
    )

    # Raise if unsupported directive found in `strftime_fmt`
    for key in non_supported:
        if key in strftime_fmt:
            raise UnsupportedDatetimeDirective(f"Unsupported datetime formatting directive: {key!r}")

    # Mapping between strftime and string formatting, according to both styles
    if new_style_fmt:
        supported = {
            "%d": "{day:02d}",    # Day of the month as a zero-padded decimal number.
            "%m": "{month:02d}",  # Month as a zero-padded decimal number.
            "%Y": "{year}",       # Year with century as a decimal number.
            "%H": "{hour:02d}",   # Hour (24-hour clock) as a zero-padded decimal number.
            "%M": "{min:02d}",    # Minute as a zero-padded decimal number.
            "%S": "{sec:02d}",    # Second as a zero-padded decimal number.
            "%f": "{us:06d}",     # Microsecond as a decimal number, zero-padded to 6 digits.
        }
    else:
        supported = {
            "%d": "%(day)02d",    # Day of the month as a zero-padded decimal number.
            "%m": "%(month)02d",  # Month as a zero-padded decimal number.
            "%Y": "%(year)s",     # Year with century as a decimal number.
            "%H": "%(hour)02d",   # Hour (24-hour clock) as a zero-padded decimal number.
            "%M": "%(min)02d",    # Minute as a zero-padded decimal number.
            "%S": "%(sec)02d",    # Second as a zero-padded decimal number.
            "%f": "%(us)06d",     # Microsecond as a decimal number, zero-padded to 6 digits.
        }

    # Create the output by replacing all directives
    for key, replacement in supported.items():
        strftime_fmt = strftime_fmt.replace(key, replacement)

    return strftime_fmt
