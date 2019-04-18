# -*- coding: utf-8 -*-

import pandas as pd

# Segregated collection of methods that require the BlockManager internal data
# structure


class TestSeriesBlockInternals(object):

    def test_setitem_invalidates_datetime_index_freq(self):
        # GH#24096 altering a datetime64tz Series inplace invalidates the
        #  `freq` attribute on the underlying DatetimeIndex

        dti = pd.date_range('20130101', periods=3, tz='US/Eastern')
        ts = dti[1]
        ser = pd.Series(dti)
        assert ser._values is not dti
        assert ser._values._data.base is not dti._data._data.base
        assert dti.freq == 'D'
        ser.iloc[1] = pd.NaT
        assert ser._values.freq is None

        # check that the DatetimeIndex was not altered in place
        assert ser._values is not dti
        assert ser._values._data.base is not dti._data._data.base
        assert dti[1] == ts
        assert dti.freq == 'D'

    def test_dt64tz_setitem_does_not_mutate_dti(self):
        # GH#21907, GH#24096
        dti = pd.date_range('2016-01-01', periods=10, tz='US/Pacific')
        ts = dti[0]
        ser = pd.Series(dti)
        assert ser._values is not dti
        assert ser._values._data.base is not dti._data._data.base
        assert ser._data.blocks[0].values is not dti
        assert (ser._data.blocks[0].values._data.base
                is not dti._data._data.base)

        ser[::3] = pd.NaT
        assert ser[0] is pd.NaT
        assert dti[0] == ts
