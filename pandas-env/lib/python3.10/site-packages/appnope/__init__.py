__version__ = "0.1.4"

import re
import sys
import platform


def _v(version_s):
    return tuple(int(s) for s in re.findall(r"\d+", version_s))


if sys.platform != "darwin" or _v(platform.mac_ver()[0]) < _v("10.9"):
    from ._dummy import *  # noqa
else:
    from ._nope import *  # noqa
