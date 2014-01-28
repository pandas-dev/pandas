from __future__ import print_function
from pandas import compat
import warnings
import nose
from nose.tools import assert_equal
from datetime import datetime

import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas.io import data as web
from pandas.io.data import DataReader, SymbolWarning
from pandas.util.testing import (assert_series_equal, assert_produces_warning,
                                 network, assert_frame_equal)
import pandas.util.testing as tm
from numpy.testing import assert_array_equal

if compat.PY3:
    from urllib.error import HTTPError
else:
    from urllib2 import HTTPError

def _skip_if_no_lxml():
    try:
        import lxml
    except ImportError:
        raise nose.SkipTest("no lxml")


def assert_n_failed_equals_n_null_columns(wngs, obj, cls=SymbolWarning):
    all_nan_cols = pd.Series(dict((k, pd.isnull(v).all()) for k, v in
                                  compat.iteritems(obj)))
    n_all_nan_cols = all_nan_cols.sum()
    valid_warnings = pd.Series([wng for wng in wngs if isinstance(wng, cls)])
    assert_equal(len(valid_warnings), n_all_nan_cols)
    failed_symbols = all_nan_cols[all_nan_cols].index
    msgs = valid_warnings.map(lambda x: x.message)
    assert msgs.str.contains('|'.join(failed_symbols)).all()


class TestGoogle(tm.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestGoogle, cls).setUpClass()
        cls.locales = tm.get_locales(prefix='en_US')
        if not cls.locales:
            raise nose.SkipTest("US English locale not available for testing")

    @classmethod
    def tearDownClass(cls):
        super(TestGoogle, cls).tearDownClass()
        del cls.locales

    @network
    def test_google(self):
        # asserts that google is minimally working and that it throws
        # an exception when DataReader can't get a 200 response from
        # google
        start = datetime(2010, 1, 1)
        end = datetime(2013, 1, 27)

        for locale in self.locales:
            with tm.set_locale(locale):
                panel = web.DataReader("F", 'google', start, end)
            self.assertEquals(panel.Close[-1], 13.68)

        self.assertRaises(Exception, web.DataReader, "NON EXISTENT TICKER",
                          'google', start, end)

    @network
    def test_get_quote_fails(self):
        self.assertRaises(NotImplementedError, web.get_quote_google,
                          pd.Series(['GOOG', 'AAPL', 'GOOG']))

    @network
    def test_get_goog_volume(self):
        for locale in self.locales:
            with tm.set_locale(locale):
                df = web.get_data_google('GOOG').sort_index()
            self.assertEqual(df.Volume.ix['OCT-08-2010'], 2863473)

    @network
    def test_get_multi1(self):
        for locale in self.locales:
            sl = ['AAPL', 'AMZN', 'GOOG']
            with tm.set_locale(locale):
                pan = web.get_data_google(sl, '2012')
            ts = pan.Close.GOOG.index[pan.Close.AAPL > pan.Close.GOOG]
            if (hasattr(pan, 'Close') and hasattr(pan.Close, 'GOOG') and
                hasattr(pan.Close, 'AAPL')):
                self.assertEquals(ts[0].dayofyear, 96)
            else:
                self.assertRaises(AttributeError, lambda: pan.Close)

    @network
    def test_get_multi2(self):
        with warnings.catch_warnings(record=True) as w:
            for locale in self.locales:
                with tm.set_locale(locale):
                    pan = web.get_data_google(['GE', 'MSFT', 'INTC'],
                                              'JAN-01-12', 'JAN-31-12')
                result = pan.Close.ix['01-18-12']
                assert_n_failed_equals_n_null_columns(w, result)

                # sanity checking

                assert np.issubdtype(result.dtype, np.floating)
                result = pan.Open.ix['Jan-15-12':'Jan-20-12']
                self.assertEqual((4, 3), result.shape)
                assert_n_failed_equals_n_null_columns(w, result)


class TestYahoo(tm.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestYahoo, cls).setUpClass()
        _skip_if_no_lxml()

    @network
    def test_yahoo(self):
        # asserts that yahoo is minimally working and that it throws
        # an exception when DataReader can't get a 200 response from
        # yahoo
        start = datetime(2010, 1, 1)
        end = datetime(2013, 1, 27)

        self.assertEquals( web.DataReader("F", 'yahoo', start,
                                          end)['Close'][-1], 13.68)

    @network
    def test_yahoo_fails(self):
        start = datetime(2010, 1, 1)
        end = datetime(2013, 1, 27)
        self.assertRaises(Exception, web.DataReader, "NON EXISTENT TICKER",
                          'yahoo', start, end)

    @network
    def test_get_quote_series(self):
        df = web.get_quote_yahoo(pd.Series(['GOOG', 'AAPL', 'GOOG']))
        assert_series_equal(df.ix[0], df.ix[2])

    @network
    def test_get_quote_string(self):
        df = web.get_quote_yahoo('GOOG')

    @network
    def test_get_quote_stringlist(self):
        df = web.get_quote_yahoo(['GOOG', 'AAPL', 'GOOG'])
        assert_series_equal(df.ix[0], df.ix[2])

    @network
    def test_get_components_dow_jones(self):
        raise nose.SkipTest('unreliable test, receive partial components back for dow_jones')

        df = web.get_components_yahoo('^DJI') #Dow Jones
        assert isinstance(df, pd.DataFrame)
        self.assertEqual(len(df), 30)

    @network
    def test_get_components_dax(self):
        raise nose.SkipTest('unreliable test, receive partial components back for dax')

        df = web.get_components_yahoo('^GDAXI') #DAX
        assert isinstance(df, pd.DataFrame)
        self.assertEqual(len(df), 30)
        self.assertEqual(df[df.name.str.contains('adidas', case=False)].index,
                         'ADS.DE')

    @network
    def test_get_components_nasdaq_100(self):
        # as of 7/12/13 the conditional will test false because the link is invalid
        raise nose.SkipTest('unreliable test, receive partial components back for nasdaq_100')

        df = web.get_components_yahoo('^NDX') #NASDAQ-100
        assert isinstance(df, pd.DataFrame)

        if len(df) > 1:
            # Usual culprits, should be around for a while
            assert 'AAPL' in df.index
            assert 'GOOG' in df.index
            assert 'AMZN' in df.index
        else:
            expected = DataFrame({'exchange': 'N/A', 'name': '@^NDX'},
                                 index=['@^NDX'])
            assert_frame_equal(df, expected)

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
        self.assertEqual(expected.shape, result.shape)

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


class TestYahooOptions(tm.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestYahooOptions, cls).setUpClass()
        _skip_if_no_lxml()

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
        super(TestYahooOptions, cls).tearDownClass()
        del cls.aapl, cls.expiry

    @network
    def test_get_options_data(self):
        try:
            calls, puts = self.aapl.get_options_data(expiry=self.expiry)
        except IndexError:
            warnings.warn("IndexError thrown no tables found")
        else:
            assert len(calls)>1
            assert len(puts)>1

    def test_get_options_data(self):

        # regression test GH6105
        self.assertRaises(ValueError,self.aapl.get_options_data,month=3)
        self.assertRaises(ValueError,self.aapl.get_options_data,year=1992)

    @network
    def test_get_near_stock_price(self):
        try:
            calls, puts = self.aapl.get_near_stock_price(call=True, put=True,
                                                        expiry=self.expiry)
        except IndexError:
            warnings.warn("IndexError thrown no tables found")
        else:
            self.assertEqual(len(calls), 5)
            self.assertEqual(len(puts), 5)

    @network
    def test_get_call_data(self):
        try:
            calls = self.aapl.get_call_data(expiry=self.expiry)
        except IndexError:
            warnings.warn("IndexError thrown no tables found")
        else:
            assert len(calls)>1

    @network
    def test_get_put_data(self):
        try:
            puts = self.aapl.get_put_data(expiry=self.expiry)
        except IndexError:
            warnings.warn("IndexError thrown no tables found")
        else:
            assert len(puts)>1


class TestOptionsWarnings(tm.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestOptionsWarnings, cls).setUpClass()
        _skip_if_no_lxml()

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
        super(TestOptionsWarnings, cls).tearDownClass()
        del cls.aapl, cls.year, cls.month

    @network
    def test_get_options_data_warning(self):
        with assert_produces_warning():
            print('month: {0}, year: {1}'.format(self.month, self.year))
            try:
                self.aapl.get_options_data(month=self.month, year=self.year)
            except IndexError:
                warnings.warn("IndexError thrown no tables found")

    @network
    def test_get_near_stock_price_warning(self):
        with assert_produces_warning():
            print('month: {0}, year: {1}'.format(self.month, self.year))
            try:
                calls_near, puts_near = self.aapl.get_near_stock_price(call=True,
                                                                    put=True,
                                                                    month=self.month,
                                                                    year=self.year)
            except IndexError:
                warnings.warn("IndexError thrown no tables found")

    @network
    def test_get_call_data_warning(self):
        with assert_produces_warning():
            print('month: {0}, year: {1}'.format(self.month, self.year))
            try:
                self.aapl.get_call_data(month=self.month, year=self.year)
            except IndexError:
                warnings.warn("IndexError thrown no tables found")

    @network
    def test_get_put_data_warning(self):
        with assert_produces_warning():
            print('month: {0}, year: {1}'.format(self.month, self.year))
            try:
                self.aapl.get_put_data(month=self.month, year=self.year)
            except IndexError:
                warnings.warn("IndexError thrown no tables found")


class TestDataReader(tm.TestCase):
    def test_is_s3_url(self):
        from pandas.io.common import _is_s3_url
        self.assert_(_is_s3_url("s3://pandas/somethingelse.com"))

    @network
    def test_read_yahoo(self):
        gs = DataReader("GS", "yahoo")
        assert isinstance(gs, DataFrame)

    @network
    def test_read_google(self):
        gs = DataReader("GS", "google")
        assert isinstance(gs, DataFrame)

    @network
    def test_read_fred(self):
        vix = DataReader("VIXCLS", "fred")
        assert isinstance(vix, DataFrame)

    @network
    def test_read_famafrench(self):
        for name in ("F-F_Research_Data_Factors",
                     "F-F_Research_Data_Factors_weekly", "6_Portfolios_2x3",
                     "F-F_ST_Reversal_Factor"):
            ff = DataReader(name, "famafrench")
            assert ff
            assert isinstance(ff, dict)


class TestFred(tm.TestCase):
    @network
    def test_fred(self):

        # Throws an exception when DataReader can't get a 200 response from
        # FRED.

        start = datetime(2010, 1, 1)
        end = datetime(2013, 1, 27)

        received = web.DataReader("GDP", "fred", start, end)['GDP'].tail(1)[0]
        self.assertEquals(int(received), 16535)

        self.assertRaises(Exception, web.DataReader, "NON EXISTENT SERIES",
                          'fred', start, end)

    @network
    def test_fred_nan(self):
        start = datetime(2010, 1, 1)
        end = datetime(2013, 1, 27)
        df = web.DataReader("DFII5", "fred", start, end)
        assert pd.isnull(df.ix['2010-01-01'][0])

    @network
    def test_fred_parts(self):
        start = datetime(2010, 1, 1)
        end = datetime(2013, 1, 27)
        df = web.get_data_fred("CPIAUCSL", start, end)
        self.assertEqual(df.ix['2010-05-01'][0], 217.23)

        t = df.CPIAUCSL.values
        assert np.issubdtype(t.dtype, np.floating)
        self.assertEqual(t.shape, (37,))

    @network
    def test_fred_part2(self):
        expected = [[576.7],
                    [962.9],
                    [684.7],
                    [848.3],
                    [933.3]]
        result = web.get_data_fred("A09024USA144NNBR", start="1915").ix[:5]
        assert_array_equal(result.values, np.array(expected))

    @network
    def test_invalid_series(self):
        name = "NOT A REAL SERIES"
        self.assertRaises(Exception, web.get_data_fred, name)

    @network
    def test_fred_multi(self):
        names = ['CPIAUCSL', 'CPALTT01USQ661S', 'CPILFESL']
        start = datetime(2010, 1, 1)
        end = datetime(2013, 1, 27)

        received = web.DataReader(names, "fred", start, end).head(1)
        expected = DataFrame([[217.478, 0.99701529, 220.544]], columns=names,
                             index=[pd.tslib.Timestamp('2010-01-01 00:00:00')])
        expected.index.rename('DATE', inplace=True)
        assert_frame_equal(received, expected, check_less_precise=True)

    @network
    def test_fred_multi_bad_series(self):

        names = ['NOTAREALSERIES', 'CPIAUCSL', "ALSO FAKE"]
        with tm.assertRaises(HTTPError):
            DataReader(names, data_source="fred")

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
