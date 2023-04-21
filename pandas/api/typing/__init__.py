# pylint: disable=undefined-all-variable
"""
Public API classes that store intermediate results useful for type-hinting.
"""

from pandas.core.groupby import (
    DataFrameGroupBy,
    SeriesGroupBy,
)
from pandas.core.resample import (
    DatetimeIndexResamplerGroupby,
    PeriodIndexResamplerGroupby,
    Resampler,
    TimedeltaIndexResamplerGroupby,
    TimeGrouper,
)
from pandas.core.window import (
    Expanding,
    ExpandingGroupby,
    ExponentialMovingWindow,
    ExponentialMovingWindowGroupby,
    Rolling,
    RollingGroupby,
    Window,
)

# TODO: Can't import Styler without importing jinja2
# from pandas.io.formats.style import Styler


def __getattr__(key: str):
    if key == "JsonReader":
        from pandas.io.json._json import JsonReader

        return JsonReader
    elif key == "StataReader":
        from pandas.io.stata import StataReader

        return StataReader
    else:
        raise AttributeError(f"module 'pandas.api.typing' has no attribute '{key}'")


def __dir__() -> list[str]:
    # include lazy imports defined in __getattr__ in dir()
    base = list(globals().keys())
    result = base + ["JsonReader", "StataReader"]
    return sorted(result)


__all__ = [
    "DataFrameGroupBy",
    "DatetimeIndexResamplerGroupby",
    "Expanding",
    "ExpandingGroupby",
    "ExponentialMovingWindow",
    "ExponentialMovingWindowGroupby",
    "JsonReader",  # pyright: ignore[reportUnsupportedDunderAll]
    "PeriodIndexResamplerGroupby",
    "Resampler",
    "Rolling",
    "RollingGroupby",
    "SeriesGroupBy",
    "StataReader",  # pyright: ignore[reportUnsupportedDunderAll]
    # See TODO above
    # "Styler",
    "TimedeltaIndexResamplerGroupby",
    "TimeGrouper",
    "Window",
]
