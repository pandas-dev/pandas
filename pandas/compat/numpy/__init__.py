""" support numpy compatibility across versions """

from distutils.version import LooseVersion
import re

import numpy as np

# numpy versioning
_np_version = np.__version__
_nlv = LooseVersion(_np_version)
np_version_under1p18 = _nlv < LooseVersion("1.18")
np_version_under1p19 = _nlv < LooseVersion("1.19")
np_version_under1p20 = _nlv < LooseVersion("1.20")
is_numpy_dev = ".dev" in str(_nlv)
_min_numpy_ver = "1.17.3"


if _nlv < _min_numpy_ver:
    raise ImportError(
        f"this version of pandas is incompatible with numpy < {_min_numpy_ver}\n"
        f"your numpy version is {_np_version}.\n"
        f"Please upgrade numpy to >= {_min_numpy_ver} to use this pandas version"
    )


_tz_regex = re.compile("[+-]0000$")


def _tz_replacer(tstring):
    if isinstance(tstring, str):
        if tstring.endswith("Z"):
            tstring = tstring[:-1]
        elif _tz_regex.search(tstring):
            tstring = tstring[:-5]
    return tstring


def np_datetime64_compat(tstring: str, unit: str = "ns"):
    """
    provide compat for construction of strings to numpy datetime64's with
    tz-changes in 1.11 that make '2015-01-01 09:00:00Z' show a deprecation
    warning, when need to pass '2015-01-01 09:00:00'
    """
    tstring = _tz_replacer(tstring)
    return np.datetime64(tstring, unit)


def np_array_datetime64_compat(arr, dtype="M8[ns]"):
    """
    provide compat for construction of an array of strings to a
    np.array(..., dtype=np.datetime64(..))
    tz-changes in 1.11 that make '2015-01-01 09:00:00Z' show a deprecation
    warning, when need to pass '2015-01-01 09:00:00'
    """
    # is_list_like; can't import as it would be circular
    if hasattr(arr, "__iter__") and not isinstance(arr, (str, bytes)):
        arr = [_tz_replacer(s) for s in arr]
    else:
        arr = _tz_replacer(arr)

    return np.array(arr, dtype=dtype)


__all__ = [
    "np",
    "_np_version",
    "is_numpy_dev",
]
