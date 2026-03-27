from __future__ import annotations

from datetime import (
    datetime,
    time,
)
from typing import TYPE_CHECKING

import numpy as np

from pandas._libs.lib import is_list_like

from pandas.core.dtypes.generic import (
    ABCIndex,
    ABCSeries,
)
from pandas.core.dtypes.missing import notna

if TYPE_CHECKING:
    from pandas._typing import DateTimeErrorChoices


def to_time(
    arg,
    format: str | None = None,
    infer_time_format: bool = False,
    errors: DateTimeErrorChoices = "raise",
):
    """
    Parse time strings to time objects using fixed strptime formats ("%H:%M",
    "%H%M", "%I:%M%p", "%I%M%p", "%H:%M:%S", "%H%M%S", "%I:%M:%S%p",
    "%I%M%S%p")

    Use infer_time_format if all the strings are in the same format to speed
    up conversion.

    Parameters
    ----------
    arg : string in time format, datetime.time, list, tuple, 1-d array,  Series
    format : str, default None
        Format used to convert arg into a time object.  If None, fixed formats
        are used.
    infer_time_format: bool, default False
        Infer the time format based on the first non-NaN element.  If all
        strings are in the same format, this will speed up conversion.
    errors : {'raise', 'coerce'}, default 'raise'
        - If 'raise', then invalid parsing will raise an exception
        - If 'coerce', then invalid parsing will be set as None

    Returns
    -------
    datetime.time
    """
    if errors not in ("raise", "coerce"):
        raise ValueError("errors must be one of 'raise', or 'coerce'.")

    def _convert_listlike(arg, format):
        if isinstance(arg, (list, tuple)):
            arg = np.array(arg, dtype="O")

        elif getattr(arg, "ndim", 1) > 1:
            raise TypeError(
                "arg must be a string, datetime, list, tuple, 1-d array, or Series"
            )

        arg = np.asarray(arg, dtype="O")

        if infer_time_format and format is None:
            format = _guess_time_format_for_array(arg)

        times: list[time | None] = []
        # On non-English locales, strptime doesn't recognize English AM/PM.
        # Check once and use a fallback parser for %p formats if needed.
        need_ampm_fallback = not _LOCALE_SUPPORTS_ENGLISH_AMPM

        if format is not None:
            if need_ampm_fallback and "%p" in format:
                _parse = _strptime_ampm_fallback
            else:
                _parse = lambda elem, fmt: datetime.strptime(elem, fmt).time()
            for element in arg:
                try:
                    times.append(_parse(element, format))
                except (ValueError, TypeError) as err:
                    if errors == "raise":
                        msg = (
                            f"Cannot convert {element} to a time with given "
                            f"format {format}"
                        )
                        raise ValueError(msg) from err
                    times.append(None)
        else:
            formats = _time_formats[:]
            format_found = False
            for element in arg:
                time_object = None
                try:
                    time_object = time.fromisoformat(element)
                except (ValueError, TypeError):
                    for time_format in formats:
                        try:
                            if need_ampm_fallback and "%p" in time_format:
                                time_object = _strptime_ampm_fallback(
                                    element, time_format
                                )
                            else:
                                time_object = datetime.strptime(
                                    element, time_format
                                ).time()
                            if not format_found:
                                # Put the found format in front
                                fmt = formats.pop(formats.index(time_format))
                                formats.insert(0, fmt)
                                format_found = True
                            break
                        except (ValueError, TypeError):
                            continue

                if time_object is not None:
                    times.append(time_object)
                elif errors == "raise":
                    raise ValueError(f"Cannot convert arg {arg} to a time")
                else:
                    times.append(None)

        return times

    if arg is None:
        return arg
    elif isinstance(arg, time):
        return arg
    elif isinstance(arg, ABCSeries):
        values = _convert_listlike(arg._values, format)
        return arg._constructor(values, index=arg.index, name=arg.name)
    elif isinstance(arg, ABCIndex):
        return _convert_listlike(arg, format)
    elif is_list_like(arg):
        return _convert_listlike(arg, format)

    return _convert_listlike(np.array([arg]), format)[0]


# Fixed time formats for time parsing
_time_formats = [
    "%H:%M",
    "%H%M",
    "%I:%M%p",
    "%I%M%p",
    "%H:%M:%S",
    "%H%M%S",
    "%I:%M:%S%p",
    "%I%M%S%p",
]


def _locale_supports_english_ampm() -> bool:
    """Check whether the current locale's strptime recognizes English AM/PM."""
    try:
        datetime.strptime("01PM", "%I%p")
        return True
    except ValueError:
        return False


_LOCALE_SUPPORTS_ENGLISH_AMPM: bool = _locale_supports_english_ampm()


def _strptime_ampm_fallback(element: str, time_format: str) -> time:
    """
    Parse a time string with %p when the locale doesn't recognize English AM/PM.

    Strips the AM/PM suffix, converts %I→%H / %p→empty, and adjusts
    the hour manually.
    """
    element_lower = element.lower()
    if element_lower.endswith("am"):
        is_pm = False
        base = element[:-2]
    elif element_lower.endswith("pm"):
        is_pm = True
        base = element[:-2]
    else:
        raise ValueError(
            f"Cannot convert {element} to a time with given format {time_format}"
        )

    alt_format = time_format.replace("%I", "%H").replace("%p", "")
    parsed = datetime.strptime(base.rstrip(), alt_format.rstrip())
    hour = parsed.hour
    if is_pm and hour != 12:
        hour += 12
    elif not is_pm and hour == 12:
        hour = 0
    return time(hour, parsed.minute, parsed.second)


def _guess_time_format_for_array(arr):
    # Try to guess the format based on the first non-NaN element
    need_ampm_fallback = not _LOCALE_SUPPORTS_ENGLISH_AMPM
    non_nan_elements = notna(arr).nonzero()[0]
    if len(non_nan_elements):
        element = arr[non_nan_elements[0]]
        for time_format in _time_formats:
            try:
                if need_ampm_fallback and "%p" in time_format:
                    _strptime_ampm_fallback(element, time_format)
                else:
                    datetime.strptime(element, time_format)
                return time_format
            except ValueError:
                pass

    return None
