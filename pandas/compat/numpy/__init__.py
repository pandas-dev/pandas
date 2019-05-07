""" support numpy compatiblitiy across versions """

import re
import numpy as np
from distutils.version import LooseVersion


# numpy versioning
_np_version = np.__version__
_nlv = LooseVersion(_np_version)
_np_version_under1p14 = _nlv < LooseVersion('1.14')
_np_version_under1p15 = _nlv < LooseVersion('1.15')
_np_version_under1p16 = _nlv < LooseVersion('1.16')
_np_version_under1p17 = _nlv < LooseVersion('1.17')


if _nlv < '1.13.3':
    raise ImportError('this version of pandas is incompatible with '
                      'numpy < 1.13.3\n'
                      'your numpy version is {0}.\n'
                      'Please upgrade numpy to >= 1.13.3 to use '
                      'this pandas version'.format(_np_version))


_tz_regex = re.compile('[+-]0000$')


def tz_replacer(s):
    if isinstance(s, str):
        if s.endswith('Z'):
            s = s[:-1]
        elif _tz_regex.search(s):
            s = s[:-5]
    return s


def np_datetime64_compat(s, *args, **kwargs):
    """
    provide compat for construction of strings to numpy datetime64's with
    tz-changes in 1.11 that make '2015-01-01 09:00:00Z' show a deprecation
    warning, when need to pass '2015-01-01 09:00:00'
    """
    s = tz_replacer(s)
    return np.datetime64(s, *args, **kwargs)


def np_array_datetime64_compat(arr, *args, **kwargs):
    """
    provide compat for construction of an array of strings to a
    np.array(..., dtype=np.datetime64(..))
    tz-changes in 1.11 that make '2015-01-01 09:00:00Z' show a deprecation
    warning, when need to pass '2015-01-01 09:00:00'
    """
    # is_list_like
    if (hasattr(arr, '__iter__') and not isinstance(arr, (str, bytes))):
        arr = [tz_replacer(s) for s in arr]
    else:
        arr = tz_replacer(arr)

    return np.array(arr, *args, **kwargs)


__all__ = ['np',
           '_np_version_under1p14',
           '_np_version_under1p15',
           '_np_version_under1p16',
           '_np_version_under1p17'
           ]
