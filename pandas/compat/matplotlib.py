from __future__ import annotations

from pandas.util.version import Version

try:
    import matplotlib

    _matplotlib_version = matplotlib.__version__
    _matplotliblv = Version(_matplotlib_version)
    is_at_least_matplotlib_36 = _matplotliblv >= Version("3.6.0")

except ImportError:
    is_at_least_matplotlib_36 = True
