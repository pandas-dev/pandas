"""Strftime-related classes and functions.
"""


class UnsupportedStrFmtDirective(ValueError):
    """The format contains a directive that is not supported in this context.
    """


_COMMON_UNSUPPORTED = (
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


_COMMON_MAP = {
    "%d": ("day", "02d"),    # Day of the month as a zero-padded decimal number.
    "%m": ("month", "02d"),  # Month as a zero-padded decimal number.
    "%Y": ("year", "d"),     # Year with century as a decimal number.
    "%H": ("hour", "02d"),   # Hour (24-hour clock) as a zero-padded decimal number.
    "%M": ("min", "02d"),    # Minute as a zero-padded decimal number.
    "%S": ("sec", "02d"),    # Second as a zero-padded decimal number.
}

_DATETIME_MAP = {
    "%f": ("us", "06d"),     # Microsecond as a decimal number, zero-padded to 6 digits.
}

_PERIOD_MAP = {
    "%f": ("fyear", "02d"),  # 'Fiscal' year without century as zero-padded decimal number [00,99]
    "%F": ("Fyear", "d"),    # 'Fiscal' year with century as a decimal number
    "%q": ("q", "d"),        # Quarter as a decimal number [1,4]
    "%l": ("ms", "03d"),     # Microsecond as a decimal number, zero-padded to 3 digits.
    "%u": ("us", "06d"),     # Microsecond as a decimal number, zero-padded to 6 digits.
    "%n": ("ns", "09d"),     # Microsecond as a decimal number, zero-padded to 9 digits.
}


def convert_strftime_format(
    strftime_fmt: str,
    target: str = "datetime",
    new_style_fmt: bool = False,
) -> str:
    """Convert a strftime formatting string into a formatting template string.

    The set of supported directives varies according to the `target`.

    This method can be tested on a single instance of

     - `datetime` or `Timestamp`, through
       `pandas.core.tools.datetimes.fast_strftime`. The
        result may be compared with `datetime.strftime` or `Timestamp.strftime`

     - `Period` through `Period.fast_strftime`. The result may be compared
        with `Period.strftime`.

    On array-like objects, this method is used in several places:

     - Subclasses of `DatelikeOps` now rely on this method in their
       `self.strftime(fmt, fast_strftime=True)` default implementation, which
       delegates to `_format_native_types`.

        - `DatetimeArray._format_native_types` relies on
          `tslib.format_array_from_datetime` which relies on this function
        - `PeriodArray._format_native_types` directly relies on this function.
        - `TimedeltaArray._format_native_types` does not currently support
          custom formats.

    In addition, `Datetime64Formatter` and `Datetime64TZFormatter` also
    rely on this when their attribute `fast_strftime` is `True` (default).

    Parameters
    ----------
    strftime_fmt : str
        The strftime format string specification, e.g. `"%Y-%m-%d %H:%M:%S"`. Note
        that not all directives are eligible to successful usage of string formatting.
        Unsupported directives will lead to an `UnsupportedStrFmtDirective` being raised.
    target : { "datetime", "date", "time", "period" }, default: "datetime"
        The kind of data that will be formatted using this template.
    new_style_fmt : bool, default: False
        Whether the output string should be new-style
        e.g. "{year}-{month:02d}-{day:02d} {hour:02d}:{min:02d}:{sec:02d}"
        or old-style e.g. "%(year)s-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d"

    Returns
    -------
    fmt_out : str
        A string that may be used to format a `datetime` variable. The style of this string
        is either old-style or new-style depending on `new_style_formatting`.
        For old-style, it may be used as `fmt_out % fmt_dct`.
        For new-style, it may be used as `fmt_out.format(**fmt_dct)`

    Raises
    ------
    UnsupportedStrFmtDirective
        Raised when the received `strftime_fmt` format contains a directive for which the output
        can not currently be created using string formatting.

    See Also
    --------
    `strftime format codes reference <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes>`_

    `Stackoverflow post <https://stackoverflow.com/a/43495629/7262247>`_
    explaining how old-style formatting is faster than new-style formatting,
    itself faster than datetime.strftime`.

    See `Period.strftime` doc for all the supported period directives (same
    directives as the :func:`time.strftime` function of the standard Python
    distribution, as well as specific additional directives ``%f``, ``%F``,
    ``%q``, ``%l``, ``%u``, ``%n``).
    """
    if target in ("datetime", "date", "time"):
        directive_maps = (_COMMON_MAP, _DATETIME_MAP)
    elif target == "period":
        directive_maps = (_COMMON_MAP, _PERIOD_MAP)
    else:
        raise ValueError(f"Invalid target: {target!r}")

    # Raise if unsupported directive found in `strftime_fmt`
    for key in _COMMON_UNSUPPORTED:
        if key in strftime_fmt:
            raise UnsupportedStrFmtDirective(f"Unsupported directive: {key!r}")

    # Mapping between strftime and string formatting, according to both styles
    if new_style_fmt:
        esc = "/_+\\"

        # Escape the %% before searching for directives, same as strftime
        strftime_fmt = strftime_fmt.replace("%%", esc)

        esc_l = "+^_\\"
        esc_r = "/_^+"

        # Create the output by replacing all directives
        for _map in directive_maps:
            for key, (_name, _fmt) in _map.items():
                # for example replace "%d" by "{day:02d}" but with escaped { and }
                strftime_fmt = strftime_fmt.replace(key, f"{esc_l}{_name}:{_fmt}{esc_r}")

        # Restore the %% into %
        strftime_fmt = strftime_fmt.replace(esc, "%")

        # Escape remaining curly braces
        strftime_fmt = strftime_fmt.replace("{", "{{").replace("}", "}}")

        # Finally replace our placeholders
        strftime_fmt = strftime_fmt.replace(esc_l, "{").replace(esc_r, "}")

    else:
        esc = "/_^+"

        # Escape the %% before searching for directives, same as strftime
        strftime_fmt = strftime_fmt.replace("%%", esc * 2)

        # Create the output by replacing all directives
        for _map in directive_maps:
            for key, (_name, _fmt) in _map.items():
                # for example replace "%d" by "%(day)02d" but with escaped %
                strftime_fmt = strftime_fmt.replace(key, f"{esc}({_name}){_fmt}")

        # Escape remaining percent signs
        strftime_fmt = strftime_fmt.replace("%", "%%")

        # Finally replace our placeholder
        strftime_fmt = strftime_fmt.replace(esc, "%")

    return strftime_fmt
