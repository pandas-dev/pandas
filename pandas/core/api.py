
# pylint: disable=W0614,W0401,W0611
# flake8: noqa

import numpy as np

from pandas.core.algorithms import factorize, unique, value_counts
from pandas.core.typed.missing import isnull, notnull
from pandas.core.categorical import Categorical
from pandas.core.groupby import Grouper
from pandas.formats.format import set_eng_float_format
from pandas.core.index import (Index, CategoricalIndex, Int64Index,
                               UInt64Index, RangeIndex, Float64Index,
                               MultiIndex, IntervalIndex)
from pandas.indexes.interval import Interval, interval_range

from pandas.core.series import Series
from pandas.core.frame import DataFrame
from pandas.core.panel import Panel, WidePanel
from pandas.core.panel4d import Panel4D
from pandas.core.reshape import (pivot_simple as pivot, get_dummies,
                                 lreshape, wide_to_long)

from pandas.core.indexing import IndexSlice
from pandas.tseries.offsets import DateOffset
from pandas.tseries.tools import to_datetime
from pandas.tseries.index import (DatetimeIndex, Timestamp,
                                  date_range, bdate_range)
from pandas.tseries.tdi import TimedeltaIndex, Timedelta
from pandas.tseries.period import Period, PeriodIndex

# see gh-14094.
from pandas.util.depr_module import _DeprecatedModule

_removals = ['day', 'bday', 'businessDay', 'cday', 'customBusinessDay',
             'customBusinessMonthEnd', 'customBusinessMonthBegin',
             'monthEnd', 'yearEnd', 'yearBegin', 'bmonthEnd', 'bmonthBegin',
             'cbmonthEnd', 'cbmonthBegin', 'bquarterEnd', 'quarterEnd',
             'byearEnd', 'week']
datetools = _DeprecatedModule(deprmod='pandas.core.datetools',
                              removals=_removals)

from pandas.core.config import (get_option, set_option, reset_option,
                                describe_option, option_context, options)


# deprecation, xref #13790
def match(*args, **kwargs):
    import warnings

    warnings.warn("pd.match() is deprecated and will be removed "
                  "in a future version",
                  FutureWarning, stacklevel=2)
    from pandas.core.algorithms import match
    return match(*args, **kwargs)


def groupby(*args, **kwargs):
    import warnings

    warnings.warn("pd.groupby() is deprecated and will be removed "
                  "Please use the Series.groupby() or "
                  "DataFrame.groupby() methods",
                  FutureWarning, stacklevel=2)
    return args[0].groupby(*args[1:], **kwargs)
