from __future__ import annotations

import warnings

from pandas.core.ndframe import *  # noqa F401

warnings.warn(
    "pandas.core.generic is deprecated "
    "and will be removed from pandas in a future version. "
    "Use pandas.core.ndframe with instead.",
    FutureWarning,
    stacklevel=2,
)
