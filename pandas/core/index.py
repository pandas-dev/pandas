import warnings

from pandas.core.indexes.api import (
    CategoricalIndex,
    DatetimeIndex,
    Float64Index,
    Index,
    Int64Index,
    IntervalIndex,
    InvalidIndexError,
    MultiIndex,
    NaT,
    NumericIndex,
    PeriodIndex,
    RangeIndex,
    TimedeltaIndex,
    UInt64Index,
    _new_Index,
    ensure_index,
    ensure_index_from_sequences,
    get_objs_combined_axis,
)

__all__ = [
    "CategoricalIndex",
    "DatetimeIndex",
    "Float64Index",
    "Index",
    "Int64Index",
    "IntervalIndex",
    "InvalidIndexError",
    "MultiIndex",
    "NaT",
    "NumericIndex",
    "PeriodIndex",
    "RangeIndex",
    "TimedeltaIndex",
    "UInt64Index",
    "_new_Index",
    "ensure_index",
    "ensure_index_from_sequences",
    "get_objs_combined_axis",
]


# GH#30193
warnings.warn(
    "pandas.core.index is deprecated and will be removed in a future version.  "
    "The public classes are available in the top-level namespace.",
    FutureWarning,
    stacklevel=2,
)
