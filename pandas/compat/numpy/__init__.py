""" support numpy compatibility across versions """

import re

import numpy as np

from pandas.util.version import Version

# numpy versioning
_np_version = np.__version__
_nlv = Version(_np_version)
np_version_under1p19 = _nlv < Version("1.19")
np_version_under1p20 = _nlv < Version("1.20")
np_version_under1p22 = _nlv < Version("1.22")
is_numpy_dev = _nlv.dev is not None
_min_numpy_ver = "1.18.5"

if is_numpy_dev or not np_version_under1p22:
    np_percentile_argname = "method"
else:
    np_percentile_argname = "interpolation"


if _nlv < Version(_min_numpy_ver):
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


__all__ = [
    "np",
    "_np_version",
    "is_numpy_dev",
]
