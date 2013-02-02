import unittest
import nose
from datetime import datetime

from pandas.util.py3compat import StringIO, BytesIO

import pandas as pd
import pandas.io.data as web
from pandas.util.testing import (network, assert_frame_equal,
                                 assert_series_equal,
                                 assert_almost_equal)
from numpy.testing.decorators import slow

import urllib2


class TestYahoo(unittest.TestCase):

    @slow
    @network
    def test_yahoo(self):
        # asserts that yahoo is minimally working and that it throws
        # an excecption when DataReader can't get a 200 response from
        # yahoo
        start = datetime(2010, 1, 1)
        end = datetime(2013, 01, 27)

        try:
            self.assertEquals(
                web.DataReader("F", 'yahoo', start, end)['Close'][-1],
                13.68)

            self.assertRaises(
                Exception,
                lambda: web.DataReader("NON EXISTENT TICKER", 'yahoo',
                                      start, end))
        except urllib2.URLError:
            try:
                urllib2.urlopen('http://www.google.com')
            except urllib2.URLError:
                raise nose.SkipTest
            else:
                raise

    @slow
    @network
    def test_get_quote(self):
        df = web.get_quote_yahoo(pd.Series(['GOOG', 'AAPL', 'GOOG']))
        assert_series_equal(df.ix[0], df.ix[2])

    @slow
    @network
    def test_get_components(self):

        df = web.get_components_yahoo() #Dow Jones (default)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 30

        df = web.get_components_yahoo('^GDAXI') #DAX
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 30
        assert df[df.name.str.contains('adidas', case=False)].index == 'ADS.DE'

        df = web.get_components_yahoo('^NDX') #NASDAQ-100
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 100
        #Usual culprits, should be around for a while
        assert 'AAPL' in df.index
        assert 'GOOG' in df.index
        assert 'AMZN' in df.index

    @slow
    @network
    def test_get_data(self):
        #single symbol
        #http://finance.yahoo.com/q/hp?s=GOOG&a=09&b=08&c=2010&d=09&e=10&f=2010&g=d
        df = web.get_data_yahoo('GOOG')
        assert df.Volume.ix['OCT-08-2010'] == 2859200

        sl = ['AAPL', 'AMZN', 'GOOG']
        pan = web.get_data_yahoo(sl, '2012')
        ts = pan.Close.GOOG.index[pan.Close.AAPL > pan.Close.GOOG]
        assert ts[0].dayofyear == 96

        dfi = web.get_components_yahoo('^DJI')
        pan = web.get_data_yahoo(dfi, 'JAN-01-12', 'JAN-31-13')
        expected = [19.02, 28.23, 25.39]
        result = pan.Close.ix['01-18-12'][['GE', 'MSFT', 'INTC']].tolist()
        assert result == expected

        pan = web.get_data_yahoo(dfi, 'JAN-01-12', 'JAN-31-13',
                                 adjust_price=True)
        expected = [18.38, 27.45, 24.54]
        result = pan.Close.ix['01-18-12'][['GE', 'MSFT', 'INTC']].tolist()
        assert result == expected

        pan = web.get_data_yahoo(dfi, '2011', ret_index=True)
        d = [[ 1.31810193,  1.08170606,  1.05281026],
             [ 1.31810193,  1.09352518,  1.05658242],
             [ 1.30228471,  1.09815005,  1.05054696],
             [ 1.30521383,  1.08119219,  1.03545832]]

        expected = pd.DataFrame(d)
        result = pan.Ret_Index[['GE', 'INTC', 'MSFT']].ix[-5:-1]
        assert_almost_equal(result.values, expected.values)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
