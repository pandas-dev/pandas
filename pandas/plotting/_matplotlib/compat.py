# being a bit too dynamic
from __future__ import annotations

from pandas.util.version import Version


def _mpl_version(version, op):
    def inner():
        try:
            import matplotlib as mpl
        except ImportError:
            return False
        return op(Version(mpl.__version__), Version(version))

    return inner
