from datetime import datetime

import numpy as np

from pandas import Series, DataFrame

from pandas.tseries.index import date_range
from pandas.tseries.offsets import Minute
from pandas.tseries.period import period_range
from pandas.tseries.resample import DatetimeIndex, TimeGrouper

import unittest
import nose

from pandas.util.testing import assert_series_equal, assert_almost_equal

class TestResample(unittest.TestCase):

    def setUp(self):
        dti = DatetimeIndex(start=datetime(2005,1,1),
                            end=datetime(2005,1,10), freq='Min')

        self.series = Series(np.random.rand(len(dti)), dti)

    def test_custom_grouper(self):

        dti = DatetimeIndex(freq='Min', start=datetime(2005,1,1),
                            end=datetime(2005,1,10))

        data = np.array([1]*len(dti))
        s = Series(data, index=dti)

        b = TimeGrouper(Minute(5))
        g = s.groupby(b)

        self.assertEquals(g.ngroups, 2593)

        # construct expected val
        arr = [5] * 2592
        arr.append(1)
        idx = dti[0:-1:5]
        idx = idx.append(DatetimeIndex([np.datetime64(dti[-1])]))
        expect = Series(arr, index=idx)

        # cython returns float for now
        result = g.agg(np.sum)
        assert_series_equal(result, expect.astype(float))

        data = np.random.rand(len(dti), 10)
        df = DataFrame(data, index=dti)
        r = df.groupby(b).agg(np.sum)

        self.assertEquals(len(r.columns), 10)
        self.assertEquals(len(r.index), 2593)

    def test_convert_basic(self):
        s = self.series

        result = s.convert('5Min', how='last')

        grouper = TimeGrouper(Minute(5), closed='right', label='right')
        expect = s.groupby(grouper).agg(lambda x: x[-1])

        assert_series_equal(result, expect)

        # from daily
        dti = DatetimeIndex(start=datetime(2005,1,1), end=datetime(2005,1,10),
                            freq='D')

        s = Series(np.random.rand(len(dti)), dti)

        # to weekly
        result = s.convert('w-sun', how='last')

        self.assertEquals(len(result), 3)
        self.assert_((result.index.dayofweek == [6,6,6]).all())
        self.assertEquals(result.irow(0), s['1/2/2005'])
        self.assertEquals(result.irow(1), s['1/9/2005'])
        self.assertEquals(result.irow(2), s.irow(-1))

        result = s.convert('W-MON', how='last')
        self.assertEquals(len(result), 2)
        self.assert_((result.index.dayofweek == [0,0]).all())
        self.assertEquals(result.irow(0), s['1/3/2005'])
        self.assertEquals(result.irow(1), s['1/10/2005'])

        result = s.convert('W-TUE', how='last')
        self.assertEquals(len(result), 2)
        self.assert_((result.index.dayofweek == [1,1]).all())
        self.assertEquals(result.irow(0), s['1/4/2005'])
        self.assertEquals(result.irow(1), s['1/10/2005'])

        result = s.convert('W-WED', how='last')
        self.assertEquals(len(result), 2)
        self.assert_((result.index.dayofweek == [2,2]).all())
        self.assertEquals(result.irow(0), s['1/5/2005'])
        self.assertEquals(result.irow(1), s['1/10/2005'])

        result = s.convert('W-THU', how='last')
        self.assertEquals(len(result), 2)
        self.assert_((result.index.dayofweek == [3,3]).all())
        self.assertEquals(result.irow(0), s['1/6/2005'])
        self.assertEquals(result.irow(1), s['1/10/2005'])

        result = s.convert('W-FRI', how='last')
        self.assertEquals(len(result), 2)
        self.assert_((result.index.dayofweek == [4,4]).all())
        self.assertEquals(result.irow(0), s['1/7/2005'])
        self.assertEquals(result.irow(1), s['1/10/2005'])

        # to biz day
        result = s.convert('B', how='last')
        self.assertEquals(len(result), 6)
        self.assert_((result.index.dayofweek == [0,1,2,3,4,0]).all())
        self.assertEquals(result.irow(0), s['1/3/2005'])
        self.assertEquals(result.irow(1), s['1/4/2005'])
        self.assertEquals(result.irow(5), s['1/10/2005'])

    def test_convert_upsample(self):
        # from daily
        dti = DatetimeIndex(start=datetime(2005,1,1), end=datetime(2005,1,10),
                            freq='D')

        s = Series(np.random.rand(len(dti)), dti)

        # to minutely, by padding
        result = s.convert('Min', method='pad')
        self.assertEquals(len(result), 12961)
        self.assertEquals(result[0], s[0])
        self.assertEquals(result[-1], s[-1])

    def test_convert_ohlc(self):
        s = self.series

        grouper = TimeGrouper(Minute(5), closed='right', label='right')
        expect = s.groupby(grouper).agg(lambda x: x[-1])
        result = s.convert('5Min', how='ohlc')

        self.assertEquals(len(result), len(expect))
        self.assertEquals(len(result.columns), 4)

        xs = result.irow(-1)
        self.assertEquals(xs['open'], s[-5])
        self.assertEquals(xs['high'], s[-5:].max())
        self.assertEquals(xs['low'], s[-5:].min())
        self.assertEquals(xs['close'], s[-1])

        xs = result.irow(1)
        self.assertEquals(xs['open'], s[1])
        self.assertEquals(xs['high'], s[1:6].max())
        self.assertEquals(xs['low'], s[1:6].min())
        self.assertEquals(xs['close'], s[5])

    def test_convert_reconvert(self):
        dti = DatetimeIndex(start=datetime(2005,1,1), end=datetime(2005,1,10),
                            freq='D')
        s = Series(np.random.rand(len(dti)), dti)
        result = s.convert('B').convert('8H')
        self.assertEquals(len(result), 22)

    def test_resample_timestamp_to_period(self):
        ts = _simple_ts('1/1/1990', '1/1/2000')

        result = ts.convert('A-DEC', kind='period')
        expected = ts.convert('A-DEC')
        expected.index = period_range('1990', '2000', freq='a-dec')
        assert_series_equal(result, expected)

        result = ts.convert('M', kind='period')
        expected = ts.convert('M')
        expected.index = period_range('1990-01', '2000-01', freq='M')
        assert_series_equal(result, expected)

        result = ts.convert('BM', kind='period')
        expected = ts.convert('BM')
        expected.index = period_range('1990-01', '2000-01', freq='M')
        assert_series_equal(result, expected)


def _simple_ts(start, end, freq='D'):
    rng = date_range(start, end, freq=freq)
    return Series(np.random.randn(len(rng)), index=rng)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

