# pylint: disable=E1101

import numpy as np

from pandas import Series
from pandas.core.indexes.datetimes import date_range
from pandas.core.indexes.period import period_range
from pandas.tseries.offsets import BDay

bday = BDay()

# The various methods we support
downsample_methods = ['min', 'max', 'first', 'last', 'sum', 'mean', 'sem',
                      'median', 'prod', 'var', 'ohlc']
upsample_methods = ['count', 'size']
series_methods = ['nunique']
resample_methods = downsample_methods + upsample_methods + series_methods


def _simple_ts(start, end, freq='D'):
    rng = date_range(start, end, freq=freq)
    return Series(np.random.randn(len(rng)), index=rng)


def _simple_pts(start, end, freq='D'):
    rng = period_range(start, end, freq=freq)
    return Series(np.random.randn(len(rng)), index=rng)
