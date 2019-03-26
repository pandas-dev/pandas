
# pylint: disable=W0614,W0401,W0611
# flake8: noqa

import numpy as np

from pandas.core.arrays.integer import (
    Int8Dtype,
    Int16Dtype,
    Int32Dtype,
    Int64Dtype,
    UInt8Dtype,
    UInt16Dtype,
    UInt32Dtype,
    UInt64Dtype,
)
from pandas.core.algorithms import factorize, unique, value_counts
from pandas.core.dtypes.missing import isna, isnull, notna, notnull
from pandas.core.dtypes.dtypes import (
    CategoricalDtype,
    PeriodDtype,
    IntervalDtype,
    DatetimeTZDtype,
)
from pandas.core.arrays import Categorical, array
from pandas.core.groupby import Grouper
from pandas.io.formats.format import set_eng_float_format
from pandas.core.index import (Index, CategoricalIndex, Int64Index,
                               UInt64Index, RangeIndex, Float64Index,
                               MultiIndex, IntervalIndex,
                               TimedeltaIndex, DatetimeIndex,
                               PeriodIndex, NaT)
from pandas.core.indexes.period import Period, period_range
from pandas.core.indexes.timedeltas import Timedelta, timedelta_range
from pandas.core.indexes.datetimes import Timestamp, date_range, bdate_range
from pandas.core.indexes.interval import Interval, interval_range

from pandas.core.series import Series
from pandas.core.frame import DataFrame
from pandas.core.panel import Panel

# TODO: Remove import when statsmodels updates #18264
from pandas.core.reshape.reshape import get_dummies

from pandas.core.indexing import IndexSlice
from pandas.core.tools.numeric import to_numeric
from pandas.tseries.offsets import DateOffset
from pandas.core.tools.datetimes import to_datetime
from pandas.core.tools.timedeltas import to_timedelta


# Deprecation: xref gh-16747
class TimeGrouper(object):

    def __new__(cls, *args, **kwargs):
        from pandas.core.resample import TimeGrouper
        import warnings
        warnings.warn("pd.TimeGrouper is deprecated and will be removed; "
                      "Please use pd.Grouper(freq=...)",
                      FutureWarning, stacklevel=2)
        return TimeGrouper(*args, **kwargs)
