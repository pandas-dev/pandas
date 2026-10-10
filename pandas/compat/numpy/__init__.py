"""support numpy compatibility across versions"""

import numpy as np

from pandas.util.version import Version

# numpy versioning
_np_version = np.__version__
_nlv = Version(_np_version)
np_version_gt2_2 = _nlv >= Version("2.2.0")
np_version_gt2_3 = _nlv >= Version("2.3.0")
np_version_gt2_5 = _nlv >= Version("2.5.0")
np_version_gt2_6 = _nlv >= Version("2.6.0.dev0")
is_numpy_dev = _nlv.dev is not None
_min_numpy_ver = "2.0.2"


if _nlv < Version(_min_numpy_ver):
    raise ImportError(
        f"Please upgrade numpy to >= {_min_numpy_ver} to use this pandas version.\n"
        f"Your numpy version is {_np_version}."
    )


__all__ = [
    "_np_version",
    "is_numpy_dev",
]
