"""support numpy compatibility across versions"""

import warnings

import numpy as np

from pandas.util.version import Version

# numpy versioning
_np_version = np.__version__
_nlv = Version(_np_version)
np_version_gt2 = _nlv >= Version("2.0.0")
is_numpy_dev = _nlv.dev is not None
_min_numpy_ver = "1.26.0"


if _nlv < Version(_min_numpy_ver):
    raise ImportError(
        f"Please upgrade numpy to >= {_min_numpy_ver} to use this pandas version.\n"
        f"Your numpy version is {_np_version}."
    )


np_long: type
np_ulong: type

if np_version_gt2:
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                r".*In the future `np\.long` will be defined as.*",
                FutureWarning,
            )
            np_long = np.long
            np_ulong = np.ulong
    except AttributeError:
        np_long = np.int_
        np_ulong = np.uint
else:
    np_long = np.int_
    np_ulong = np.uint


__all__ = [
    "_np_version",
    "is_numpy_dev",
]
