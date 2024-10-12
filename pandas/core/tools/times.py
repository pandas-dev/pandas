from __future__ import annotations
from datetime import datetime, time
from typing import TYPE_CHECKING, List, Union
import numpy as np
from pandas._libs.lib import is_list_like
from pandas.core.dtypes.generic import ABCIndex, ABCSeries
from pandas.core.dtypes.missing import notna

if TYPE_CHECKING:
    from pandas._typing import DateTimeErrorChoices


def to_time(
        arg,
        format: Union[str, List[str], None] = None,
        infer_time_format: bool = False,
        errors: DateTimeErrorChoices = "raise",
        custom_formats: List[str] = None,
):
    """
Parse time strings to time objects using fixed strptime formats ("%H:%M",
"%H%M", "%I:%M%p", "%I%M%p", "%H:%M:%S", "%H%M%S", "%I:%M:%S%p",
"%I%M%S%p") and additional custom formats.

Use infer_time_format if all the strings are in the same format to speed
up conversion.

Parameters
----------
arg : string in time format, datetime.time, list, tuple, 1-d array, Series
format : str or list of str, default None
    Format(s) used to convert arg into a time object. If None, fixed
    formats are used.
infer_time_format: bool, default False
    Infer the time format based on the first non-NaN element. If all
    strings are in the same format, this will speed up conversion.
errors : {'raise', 'coerce'}, default 'raise'
    - If 'raise', then invalid parsing will raise an exception
    - If 'coerce', then invalid parsing will be set as None
custom_formats : list of str, default None
    Additional custom time formats to use.
Returns
-------
datetime.time or list of datetime.time
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
        times = []
        if format is not None:
            if isinstance(format, list):
                for element in arg:
                    for fmt in format:
                        try:
                            times.append(datetime.strptime(element, fmt).time())
                            break
                        except (ValueError, TypeError):
                            continue
                    else:
                        if errors == "raise":
                            msg = (
                                f"Cannot convert {element} to a time with given "f"formats {format}")
                            raise ValueError(msg)
                        times.append(None)
            else:
                for element in arg:
                    try:
                        times.append(datetime.strptime(element, format).time())
                    except (ValueError, TypeError) as err:
                        if errors == "raise":
                            msg = (f"Cannot convert {element} to a time withgiven "f"format {format}")
                            raise ValueError(msg) from err
                        times.append(None)
        else:
            formats = _time_formats + (custom_formats or [])
            for element in arg:
                time_object = None
                try:
                    time_object = time.fromisoformat(element)
                except (ValueError, TypeError):
                    for time_format in formats:
                        try:
                            time_object = datetime.strptime(element, time_format).time()
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


def _guess_time_format_for_array(arr):
    # Try to guess the format based on the first non-NaN element
    non_nan_elements = notna(arr).nonzero()[0]
    if len(non_nan_elements):
        element = arr[non_nan_elements[0]]
        for time_format in _time_formats:
            try:
                datetime.strptime(element, time_format)
                return time_format
            except ValueError:
                pass
    return None
