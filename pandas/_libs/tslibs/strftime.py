"""Strftime-related classes and functions."""

from __future__ import annotations

from datetime import time
import locale


class UnsupportedStrFmtDirective(ValueError):
    """The format contains a directive that is not supported in this context."""


_COMMON_UNSUPPORTED = (
    # 1- Names not in the numpy or datetime attr representation
    "%a",  # Weekday as locale's abbreviated name.
    "%A",  # Weekday as locale's full name.
    "%w",  # Weekday as a decimal number, where 0 is Sunday and 6 is Saturday.
    "%b",  # Month as locale's abbreviated name.
    "%B",  # Month as locale's full name.
    # 2- TODO Below Time offset and timezone information ... but may be hard
    "%z",  # UTC offset in the form ±HHMM[SS[.ffffff]] ("" if tz naive).
    "%Z",  # Time zone name ("" if tz naive).
    # 3- Probably too complex ones for now
    "%j",  # Day of the year as a zero-padded decimal number.
    "%U",  # Week number of the year (Sunday as the first day of the week) as
    # a zero-padded decimal number. All days in a new year preceding the first
    # Sunday are considered to be in week 0.
    "%W",  # Week number of the year (Monday as the first day of the week) as
    # a zero-padded decimal number. All days in a new year preceding the first
    # Monday are considered to be in week 0.
    "%c",  # Locale's appropriate date and time representation.
    "%x",  # Locale's appropriate date representation.
    "%X",  # Locale's appropriate time representation.
)


_COMMON_MAP = {
    "%d": ("day", "02d"),  # Day of the month as a zero-padded decimal number.
    "%m": ("month", "02d"),  # Month as a zero-padded decimal number.
    "%Y": ("year", "04d"),  # Year with century as a 0-padded decimal number.
    "%y": ("shortyear", "02d"),  # Year without century as 0-padded decimal nb.
    "%H": ("hour", "02d"),  # Hour (24-hour clock) as 0-padded decimal number.
    "%I": ("hour12", "02d"),  # Hour (12-hour clock) as a 0-padded decimal nb.
    "%p": ("ampm", "s"),  # Locale's equivalent of either AM or PM.
    "%M": ("min", "02d"),  # Minute as a zero-padded decimal number.
    "%S": ("sec", "02d"),  # Second as a zero-padded decimal number.
}

_DATETIME_MAP = {
    "%f": ("us", "06d"),  # Microsecond as decimal number, 0-padded to 6 digits
}
_DATETIME_UNSUPPORTED = (
    "%F",
    "%q",
    "%l",
    "%u",
    "%n",
)

_PERIOD_MAP = {
    "%f": (
        "fyear",
        "02d",
    ),  # 'Fiscal' year without century as zero-padded decimal number [00,99]
    "%F": ("Fyear", "d"),  # 'Fiscal' year with century as a decimal number
    "%q": ("q", "d"),  # Quarter as a decimal number [1,4]
    "%l": ("ms", "03d"),  # Millisecond as decimal number, 0-padded 3 digits
    "%u": ("us", "06d"),  # Microsecond as decimal number, 0-padded 6 digits
    "%n": ("ns", "09d"),  # Nanosecond as decimal number, 0-padded 9 digits
}
_PERIOD_UNSUPPORTED = ()


class LocaleSpecificDtStrings:
    """A container for date/time strings used in a specific locale.

    We will use these when formatting datetime as string using string templates, which
    is faster than strftime when executed on arrays.

    `get_current_locale_specific_string()` is the recommended way to get an instance,
    as it provides caching.

    Attributes
    ----------
    am : str
        Used in the %p strftime directive. Locale's equivalent of AM.
    pm : str
        Used in the %p strftime directive. Locale's equivalent of PM.
    """

    __slots__ = ("am", "pm")

    def __init__(self, am: str, pm: str) -> None:
        self.am = am
        self.pm = pm

    def __repr__(self) -> str:
        attrs = ", ".join([f"{k}={getattr(self, k)!r}" for k in type(self).__slots__])
        return f"{type(self).__name__}({attrs})"

    @classmethod
    def get_current(cls):
        return LocaleSpecificDtStrings(
            am=time(1).strftime("%p"),
            pm=time(13).strftime("%p"),
        )


_locale_specifics: dict[str, LocaleSpecificDtStrings] = {}


def get_current_locale_specific_string() -> LocaleSpecificDtStrings:
    """Return a `LocaleSpecificDtStrings` for the current locale.

    This function caches results in the `_locale_specifics` dict.
    """

    # Get current locale
    current_locale = locale.setlocale(locale.LC_ALL)

    try:
        # Any entry in cache for current locale ?
        return _locale_specifics[current_locale]
    except KeyError:
        # Create it using current locale, and cache it
        o = LocaleSpecificDtStrings.get_current()
        _locale_specifics[current_locale] = o
        return o


def convert_strftime_format(
    strftime_fmt: str,
    target: str,
    new_style_fmt: bool = False,
) -> tuple[str, LocaleSpecificDtStrings]:
    """Convert a strftime formatting string into a formatting template string.

    The set of supported directives varies according to the `target`.

    This method can be tested on a single instance of

     - `Timestamp`, through `Timestamp._strftime_pystr`. The result may be compared
       with `Timestamp.strftime`

     - `Period` through `Period._strftime_pystr`. The result may be compared
        with `Period.strftime`.

    On array-like objects, this method is used in several places:

     - Subclasses of `DatelikeOps` now rely on this method in their
       `self.strftime(fmt)` default implementation, which delegates to
       `_format_native_types`.

        - `DatetimeArray._format_native_types` relies on
          `tslib.format_array_from_datetime` which relies on this function
        - `PeriodArray._format_native_types` directly relies on this function.
        - `TimedeltaArray._format_native_types` does not currently support
          custom formats.

    In addition, `Datetime64Formatter` and `Datetime64TZFormatter` rely on this
    too.

    Parameters
    ----------
    strftime_fmt : str
        The strftime format string specification, e.g. `"%Y-%m-%d %H:%M:%S"`.
        Note that not all directives are eligible to successful usage of string
        formatting. Unsupported directives will lead to an
        `UnsupportedStrFmtDirective` being raised.
    target : { "datetime", "date", "time", "period" }, default: "datetime"
        The kind of data that will be formatted using this template.
    new_style_fmt : bool, default: False
        Whether the output string should be new-style
        e.g. "{year}-{month:02d}-{day:02d} {hour:02d}:{min:02d}:{sec:02d}"
        or old-style
        e.g. "%(year)s-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d"

    Returns
    -------
    fmt_out : str
        A string that may be used to format a `datetime` variable. The style of
        this string is either old-style or new-style depending on
        `new_style_formatting`.
        For old-style, it may be used as `fmt_out % fmt_dct`.
        For new-style, it may be used as `fmt_out.format(**fmt_dct)`
    locale_dt_strings : LocaleSpecificDtStrings
        An object containing the locale-specific strings needed for some of the
        directives. For example locale_dt_strings.am and locale_dt_strings.pm should be
        used to fill the "ampm" part of the template, induced by directive %p.

    Raises
    ------
    UnsupportedStrFmtDirective
        Raised when the received `strftime_fmt` format contains a directive for
        which the output can not currently be created using string formatting.

    See Also
    --------
    `strftime format codes reference
    <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format\
-codes>`_

    `Stackoverflow post <https://stackoverflow.com/a/43495629/7262247>`_
    explaining how old-style formatting is faster than new-style formatting,
    itself faster than datetime.strftime`.

    See `Period.strftime` doc for all the supported period directives (same
    directives as the :func:`time.strftime` function of the standard Python
    distribution, as well as specific additional directives ``%f``, ``%F``,
    ``%q``, ``%l``, ``%u``, ``%n``).
    """
    unsupported: tuple[tuple[str, ...], ...]
    if target in ("datetime", "date", "time"):
        directive_maps = (_COMMON_MAP, _DATETIME_MAP)
        unsupported = (_COMMON_UNSUPPORTED, _DATETIME_UNSUPPORTED)
    elif target == "period":
        directive_maps = (_COMMON_MAP, _PERIOD_MAP)
        unsupported = (_COMMON_UNSUPPORTED, _PERIOD_UNSUPPORTED)
    else:
        raise ValueError(f"Invalid target: {target!r}")

    # Raise if unsupported directive found in `strftime_fmt`
    for _u in unsupported:
        for key in _u:
            if key in strftime_fmt:
                raise UnsupportedStrFmtDirective(f"Unsupported directive: '{key}'")

    # Find an escape sequence, that we will use to replace all '%' signs
    esc = _create_escape_sequence(strftime_fmt, init_esc="-+", prefix="-")

    # Escape the %% before searching for directives (we will put them back at the end)
    strftime_fmt = strftime_fmt.replace("%%", esc * 2)

    # Mapping between strftime and string formatting, according to both styles
    if new_style_fmt:
        # Escape single curly braces
        strftime_fmt = strftime_fmt.replace("{", "{{").replace("}", "}}")

        # Create the output by replacing all directives
        for _map in directive_maps:
            for key, (_name, _fmt) in _map.items():
                # for example replace "%d" by "{day:02d}"
                strftime_fmt = strftime_fmt.replace(key, f"{{{_name}:{_fmt}}}")

        # If there are remaining percent signs, be conservative and fallback
        if "%" in strftime_fmt:
            raise UnsupportedStrFmtDirective("Unsupported directive found")

    else:
        # Create the output by replacing all directives
        for _map in directive_maps:
            for key, (_name, _fmt) in _map.items():
                # for example replace "%d" by "%(day)02d" but with escaped %
                strftime_fmt = strftime_fmt.replace(key, f"{esc}({_name}){_fmt}")

        # If there are remaining percent signs, raise 'unsupported directive' so that
        # the caller can fallback to OS C strftime engine.
        if "%" in strftime_fmt:
            raise UnsupportedStrFmtDirective("Unsupported directive found")

    # Restore the escaped %%
    strftime_fmt = strftime_fmt.replace(esc, "%")

    return strftime_fmt, get_current_locale_specific_string()


def _create_escape_sequence(txt: str, init_esc: str = "+", prefix: str = "-") -> str:
    """Return a unique string that does not exist in txt, by prepending as many
    `prefix` as necessary to the initial proposal `init_esc`."""

    if init_esc in prefix:
        raise ValueError("`ini_esc` must not be a subset of `prefix`")

    # Prepend `ini_esc` with `str_to_add` as many times as necessary
    while init_esc in txt:
        init_esc = prefix + init_esc
    return init_esc
