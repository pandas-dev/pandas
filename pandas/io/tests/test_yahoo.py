import unittest
import nose
from datetime import datetime

import pandas as pd
import numpy as np
import pandas.io.data as web
from pandas.util.testing import (network, assert_series_equal,
                                 assert_produces_warning)
from numpy.testing import assert_array_equal


class TestYahoo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import lxml
        except ImportError:
            raise nose.SkipTest

    @network
    def test_yahoo(self):
        # asserts that yahoo is minimally working and that it throws
        # an exception when DataReader can't get a 200 response from
        # yahoo
        start = datetime(2010, 1, 1)
        end = datetime(2013, 01, 27)

        self.assertEquals( web.DataReader("F", 'yahoo', start,
                                          end)['Close'][-1], 13.68)

    @network
    def test_yahoo_fails(self):
        start = datetime(2010, 1, 1)
        end = datetime(2013, 01, 27)
        self.assertRaises(Exception, web.DataReader, "NON EXISTENT TICKER",
                          'yahoo', start, end)

    @network
    def test_get_quote(self):
        df = web.get_quote_yahoo(pd.Series(['GOOG', 'AAPL', 'GOOG']))
        assert_series_equal(df.ix[0], df.ix[2])

    @network
    def test_get_components_dow_jones(self):
        df = web.get_components_yahoo('^DJI') #Dow Jones
        assert isinstance(df, pd.DataFrame)
        self.assertEqual(len(df), 30)

    @network
    def test_get_components_dax(self):
        df = web.get_components_yahoo('^GDAXI') #DAX
        assert isinstance(df, pd.DataFrame)
        self.assertEqual(len(df), 30)
        self.assertEqual(df[df.name.str.contains('adidas', case=False)].index,
                         'ADS.DE')

    @network
    def test_get_components_nasdaq_100(self):
        df = web.get_components_yahoo('^NDX') #NASDAQ-100
        assert isinstance(df, pd.DataFrame)
        # Usual culprits, should be around for a while
        assert 'AAPL' in df.index
        assert 'GOOG' in df.index
        assert 'AMZN' in df.index

    @network
    def test_get_data_single_symbol(self):
        #single symbol
        #http://finance.yahoo.com/q/hp?s=GOOG&a=09&b=08&c=2010&d=09&e=10&f=2010&g=d
        df = web.get_data_yahoo('GOOG')
        self.assertEqual(df.Volume.ix['OCT-08-2010'], 2859200)

    @network
    def test_get_data_multiple_symbols(self):
        sl = ['AAPL', 'AMZN', 'GOOG']
        pan = web.get_data_yahoo(sl, '2012')

        def testit():
            ts = pan.Close.GOOG.index[pan.Close.AAPL > pan.Close.GOOG]
            self.assertEquals(ts[0].dayofyear, 96)

        if hasattr(pan.Close, 'GOOG') and hasattr(pan.Close, 'AAPL'):
            testit()
        else:
            self.assertRaises(AttributeError, testit)

    @network
    def test_get_data_multiple_symbols_two_dates(self):
        pan = web.get_data_yahoo(['GE', 'MSFT', 'INTC'], 'JAN-01-12', 'JAN-31-12')
        result = pan.Close.ix['01-18-12']
        self.assertEqual(len(result), 3)

        # sanity checking
        assert np.issubdtype(result.dtype, np.floating)

        expected = np.array([[ 18.99,  28.4 ,  25.18],
                             [ 18.58,  28.31,  25.13],
                             [ 19.03,  28.16,  25.52],
                             [ 18.81,  28.82,  25.87]])
        result = pan.Open.ix['Jan-15-12':'Jan-20-12']
        assert_array_equal(np.array(expected).shape, result.shape)

    @network
    def test_get_date_ret_index(self):
        pan = web.get_data_yahoo(['GE', 'INTC', 'IBM'], '1977', '1987',
                                 ret_index=True)
        self.assert_(hasattr(pan, 'Ret_Index'))
        if hasattr(pan, 'Ret_Index') and hasattr(pan.Ret_Index, 'INTC'):
            tstamp = pan.Ret_Index.INTC.first_valid_index()
            result = pan.Ret_Index.ix[tstamp]['INTC']
            self.assertEqual(result, 1.0)

        # sanity checking
        assert np.issubdtype(pan.values.dtype, np.floating)


class TestYahooOptions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import lxml
        except ImportError:
            raise nose.SkipTest

        # aapl has monthlies
        cls.aapl = web.Options('aapl', 'yahoo')
        today = datetime.today()
        year = today.year
        month = today.month + 1
        if month > 12:
            year = year + 1
            month = 1
        cls.expiry = datetime(year, month, 1)

    @classmethod
    def tearDownClass(cls):
        del cls.aapl, cls.expiry

    @network
    def test_get_options_data(self):
        calls, puts = self.aapl.get_options_data(expiry=self.expiry)
        assert len(calls)>1
        assert len(puts)>1

    @network
    def test_get_near_stock_price(self):
        calls, puts = self.aapl.get_near_stock_price(call=True, put=True,
                                                     expiry=self.expiry)
        self.assertEqual(len(calls), 5)
        self.assertEqual(len(puts), 5)

    @network
    def test_get_call_data(self):
        calls = self.aapl.get_call_data(expiry=self.expiry)
        assert len(calls)>1

    @network
    def test_get_put_data(self):
        puts = self.aapl.get_put_data(expiry=self.expiry)
        assert len(puts)>1


class TestOptionsWarnings(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import lxml
        except ImportError:
            raise nose.SkipTest

        with assert_produces_warning(FutureWarning):
            cls.aapl = web.Options('aapl')

        today = datetime.today()
        cls.year = today.year
        cls.month = today.month + 1
        if cls.month > 12:
            cls.year += 1
            cls.month = 1

    @classmethod
    def tearDownClass(cls):
        del cls.aapl, cls.year, cls.month

    @network
    def test_get_options_data_warning(self):
        with assert_produces_warning(FutureWarning):
            print('month: {0}, year: {1}'.format(self.month, self.year))
            self.aapl.get_options_data(month=self.month, year=self.year)

    @network
    def test_get_near_stock_price_warning(self):
        with assert_produces_warning(FutureWarning):
            print('month: {0}, year: {1}'.format(self.month, self.year))
            calls_near, puts_near = self.aapl.get_near_stock_price(call=True,
                                                                   put=True,
                                                                   month=self.month,
                                                                   year=self.year)

    @network
    def test_get_call_data_warning(self):
        with assert_produces_warning(FutureWarning):
            print('month: {0}, year: {1}'.format(self.month, self.year))
            self.aapl.get_call_data(month=self.month, year=self.year)

    @network
    def test_get_put_data_warning(self):
        with assert_produces_warning(FutureWarning):
            print('month: {0}, year: {1}'.format(self.month, self.year))
            self.aapl.get_put_data(month=self.month, year=self.year)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
