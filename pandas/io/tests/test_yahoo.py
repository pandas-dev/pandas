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

        df = web.get_components_yahoo('^DJI') #Dow Jones
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 30

        df = web.get_components_yahoo('^GDAXI') #DAX
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 30
        assert df[df.name.str.contains('adidas', case=False)].index == 'ADS.DE'

        df = web.get_components_yahoo('^NDX') #NASDAQ-100
        assert isinstance(df, pd.DataFrame)
        #assert len(df) == 100
        #Usual culprits, should be around for a while
        assert 'AAPL' in df.index
        assert 'GOOG' in df.index
        assert 'AMZN' in df.index

    @slow
    @network
    def test_get_data(self):
        import numpy as np
        #single symbol
        #http://finance.yahoo.com/q/hp?s=GOOG&a=09&b=08&c=2010&d=09&e=10&f=2010&g=d
        df = web.get_data_yahoo('GOOG')
        assert df.Volume.ix['OCT-08-2010'] == 2859200

        sl = ['AAPL', 'AMZN', 'GOOG']
        pan = web.get_data_yahoo(sl, '2012')
        ts = pan.Close.GOOG.index[pan.Close.AAPL > pan.Close.GOOG]
        assert ts[0].dayofyear == 96

        #dfi = web.get_components_yahoo('^DJI')
        #pan = web.get_data_yahoo(dfi, 'JAN-01-12', 'JAN-31-12')
        pan = web.get_data_yahoo(['GE', 'MSFT', 'INTC'], 'JAN-01-12', 'JAN-31-12')
        expected = [19.02, 28.23, 25.39]
        result = pan.Close.ix['01-18-12'][['GE', 'MSFT', 'INTC']].tolist()
        assert result == expected

        # sanity checking
        t= np.array(result)
        assert     np.issubdtype(t.dtype, np.floating)
        assert     t.shape == (3,)

        expected = [[ 18.99,  28.4 ,  25.18],
                    [ 18.58,  28.31,  25.13],
                    [ 19.03,  28.16,  25.52],
                    [ 18.81,  28.82,  25.87]]
        result = pan.Open.ix['Jan-15-12':'Jan-20-12'][['GE', 'MSFT', 'INTC']].values
        assert (result == expected).all()

        #Check ret_index
        pan = web.get_data_yahoo(['GE', 'INTC', 'IBM'], '1977', '1987',
                                 ret_index=True)
        tstamp = pan.Ret_Index.INTC.first_valid_index()
        result = pan.Ret_Index.ix[tstamp]['INTC']
        expected = 1.0
        assert result == expected

        # sanity checking
        t= np.array(pan)
        assert     np.issubdtype(t.dtype, np.floating)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
