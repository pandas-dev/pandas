import pandas as pd
import pandas.util.testing as tm
from pandas.util._decorators import cache_readonly

_ts = tm.makeTimeSeries()


class TestData(object):

    @cache_readonly
    def ts(self):
        ts = _ts.copy()
        ts.name = 'ts'
        return ts

    @cache_readonly
    def series(self):
        series = tm.makeStringSeries()
        series.name = 'series'
        return series

    @cache_readonly
    def objSeries(self):
        objSeries = tm.makeObjectSeries()
        objSeries.name = 'objects'
        return objSeries

    @cache_readonly
    def empty(self):
        return pd.Series([], index=[])
