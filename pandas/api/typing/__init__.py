"""
Public API classes that store intermediate results useful for type-hinting.
"""

from pandas.core.groupby import (
    DataFrameGroupBy,
    SeriesGroupBy,
)
from pandas.core.resample import (
    Resampler,
    TimeGrouper,
)
from pandas.core.window import (
    Expanding,
    ExponentialMovingWindow,
    Rolling,
    Window,
)

# TODO: Can't import Styler without importing jinja2
# from pandas.io.formats.style import Styler
from pandas.io.json._json import JsonReader
from pandas.io.stata import StataReader

__all__ = [
    "DataFrameGroupBy",
    "Expanding",
    "ExponentialMovingWindow",
    "JsonReader",
    "Resampler",
    "Rolling",
    "SeriesGroupBy",
    "StataReader",
    # See TODO above
    # "Styler",
    "TimeGrouper",
    "Window",
]
