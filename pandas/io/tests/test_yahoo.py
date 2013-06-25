import unittest
import nose
from datetime import datetime
import warnings

import pandas as pd
import pandas.io.data as web
from pandas.util.testing import network, assert_series_equal, with_connectivity_check


class TestYahoo(unittest.TestCase):

    @with_connectivity_check("http://www.google.com")
    def test_yahoo(self):
        # asserts that yahoo is minimally working and that it throws
        # an exception when DataReader can't get a 200 response from
        # yahoo
        start = datetime(2010, 1, 1)
        end = datetime(2013, 01, 27)

        self.assertEquals(
            web.DataReader("F", 'yahoo', start, end)['Close'][-1],
            13.68)

        self.assertRaises(
            Exception,
            lambda: web.DataReader("NON EXISTENT TICKER", 'yahoo',
                                      start, end))

    @network
    def test_get_quote(self):
        df = web.get_quote_yahoo(pd.Series(['GOOG', 'AAPL', 'GOOG']))
        assert_series_equal(df.ix[0], df.ix[2])


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

    @network
    def test_options(self):
        try:
            import lxml
        except ImportError:
            raise nose.SkipTest

        ##### FAILING #####
        raise nose.SkipTest('this test is currently failing')

        # aapl has monthlies
        aapl = web.Options('aapl', 'yahoo')
        today = datetime.today()
        year = today.year
        month = today.month+1
        if (month>12):
            year = year +1
            month = 1
        expiry=datetime(year, month, 1)
        (calls, puts) = aapl.get_options_data(expiry=expiry)
        assert len(calls)>1
        assert len(puts)>1
        (calls, puts) = aapl.get_near_stock_price(call=True, put=True, expiry=expiry)
        assert len(calls)==5
        assert len(puts)==5
        calls = aapl.get_call_data(expiry=expiry)
        assert len(calls)>1
        puts = aapl.get_put_data(expiry=expiry)
        assert len(puts)>1

    @network
    def test_options_warnings(self):
        try:
            import lxml
        except ImportError:
            raise nose.SkipTest
        with warnings.catch_warnings(record=True) as w:
            warnings.resetwarnings()
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # aapl has monthlies
            aapl = web.Options('aapl')
            today = datetime.today()
            year = today.year
            month = today.month+1
            if (month>12):
                year = year +1
                month = 1
            (calls, puts) = aapl.get_options_data(month=month, year=year)
            (calls, puts) = aapl.get_near_stock_price(call=True, put=True, month=month, year=year)
            calls = aapl.get_call_data(month=month, year=year)
            puts = aapl.get_put_data(month=month, year=year)
            print(w)
            assert len(w) == 5
            assert "deprecated" in str(w[0].message)
            assert "deprecated" in str(w[1].message)
            assert "deprecated" in str(w[2].message)
            assert "deprecated" in str(w[3].message)
            assert "deprecated" in str(w[4].message)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
