"""
Public API classes that store intermediate results useful for type-hinting.
"""

from typing import TYPE_CHECKING

from pandas._libs import NaTType
from pandas._libs.lib import NoDefault
from pandas._libs.missing import NAType

from pandas.core.col import Expression
from pandas.core.groupby import (
    DataFrameGroupBy,
    SeriesGroupBy,
)
from pandas.core.indexes.frozen import FrozenList
from pandas.core.resample import (
    DatetimeIndexResamplerGroupBy,
    PeriodIndexResamplerGroupBy,
    Resampler,
    TimedeltaIndexResamplerGroupBy,
    TimeGrouper,
)
from pandas.core.window import (
    Expanding,
    ExpandingGroupBy,
    ExponentialMovingWindow,
    ExponentialMovingWindowGroupBy,
    Rolling,
    RollingGroupBy,
    Window,
)

# TODO: Can't import Styler without importing jinja2
# from pandas.io.formats.style import Styler
from pandas.io.json._json import JsonReader
from pandas.io.parquet import ParquetFileReader
from pandas.io.sas.sasreader import SASReader
from pandas.io.stata import StataReader

__all__ = [
    "DataFrameGroupBy",
    "DatetimeIndexResamplerGroupBy",
    "Expanding",
    "ExpandingGroupBy",
    "ExponentialMovingWindow",
    "ExponentialMovingWindowGroupBy",
    "Expression",
    "FrozenList",
    "JsonReader",
    "NAType",
    "NaTType",
    "NoDefault",
    "ParquetFileReader",
    "PeriodIndexResamplerGroupBy",
    "Resampler",
    "Rolling",
    "RollingGroupBy",
    "SASReader",
    "SeriesGroupBy",
    "StataReader",
    "TimeGrouper",
    "TimedeltaIndexResamplerGroupBy",
    "Window",
]


if TYPE_CHECKING:
    # GH#49578 this module exists for type-hinting, so annotating with a
    # deprecated spelling has to keep working until it is removed. Binding these
    # unconditionally would instead put them back in the public namespace.
    DatetimeIndexResamplerGroupby = DatetimeIndexResamplerGroupBy
    ExpandingGroupby = ExpandingGroupBy
    ExponentialMovingWindowGroupby = ExponentialMovingWindowGroupBy
    PeriodIndexResamplerGroupby = PeriodIndexResamplerGroupBy
    RollingGroupby = RollingGroupBy
    TimedeltaIndexResamplerGroupby = TimedeltaIndexResamplerGroupBy
else:
    # would otherwise land in the public namespace; see test_api_typing
    del TYPE_CHECKING

    def __getattr__(name: str) -> object:
        deprecated = {
            "DatetimeIndexResamplerGroupby": DatetimeIndexResamplerGroupBy,
            "ExpandingGroupby": ExpandingGroupBy,
            "ExponentialMovingWindowGroupby": ExponentialMovingWindowGroupBy,
            "PeriodIndexResamplerGroupby": PeriodIndexResamplerGroupBy,
            "RollingGroupby": RollingGroupBy,
            "TimedeltaIndexResamplerGroupby": TimedeltaIndexResamplerGroupBy,
        }
        if name in deprecated:
            # imported here rather than at module level, which would expose
            # them as pandas.api.typing attributes
            import warnings

            from pandas.errors import Pandas4Warning

            new = deprecated[name]
            warnings.warn(
                f"pandas.api.typing.{name} is deprecated and will be removed in "
                f"a future version. Use pandas.api.typing.{new.__name__} instead.",
                Pandas4Warning,
                stacklevel=2,
            )
            return new
        raise AttributeError(f"module 'pandas.api.typing' has no attribute '{name}'")
