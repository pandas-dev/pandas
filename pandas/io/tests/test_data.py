from __future__ import print_function
from pandas import compat
import warnings
import nose
from nose.tools import assert_equal
from datetime import datetime
import os

import numpy as np
import pandas as pd
from pandas import DataFrame, Timestamp
from pandas.util.testing import (assert_series_equal, assert_produces_warning,
                                 network, assert_frame_equal)
import pandas.util.testing as tm

with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
    from pandas.io import data as web

from pandas.io.data import DataReader, SymbolWarning, RemoteDataError, _yahoo_codes

if compat.PY3:
    from urllib.error import HTTPError
else:
    from urllib2 import HTTPError


def _skip_if_no_lxml():
    try:
        import lxml
    except ImportError:
        raise nose.SkipTest("no lxml")

def _skip_if_no_bs():
    try:
        import bs4
        import html5lib
    except ImportError:
        raise nose.SkipTest("no html5lib/bs4")


def assert_n_failed_equals_n_null_columns(wngs, obj, cls=SymbolWarning):
    all_nan_cols = pd.Series(dict((k, pd.isnull(v).all()) for k, v in
                                  compat.iteritems(obj)))
    n_all_nan_cols = all_nan_cols.sum()
    valid_warnings = pd.Series([wng for wng in wngs if wng.category == cls])
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
            self.assertEqual(panel.Close[-1], 13.68)

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
            self.assertEqual(df.Volume.ix['JAN-02-2015'], 1446662)

    @network
    def test_get_multi1(self):
        for locale in self.locales:
            sl = ['AAPL', 'AMZN', 'GOOG']
            with tm.set_locale(locale):
                pan = web.get_data_google(sl, '2012', '2013')
            ts = pan.Close.GOOG.index[pan.Close.AAPL < pan.Close.GOOG]
            if (hasattr(pan, 'Close') and hasattr(pan.Close, 'GOOG') and
                hasattr(pan.Close, 'AAPL')):
                self.assertEqual(ts[0].dayofyear, 3)
            else:
                self.assertRaises(AttributeError, lambda: pan.Close)

    @network
    def test_get_multi_invalid(self):
        sl = ['AAPL', 'AMZN', 'INVALID']
        with tm.assert_produces_warning(SymbolWarning):
            pan = web.get_data_google(sl, '2012')
            self.assertIn('INVALID', pan.minor_axis)

    @network
    def test_get_multi_all_invalid(self):
        sl = ['INVALID', 'INVALID2', 'INVALID3']
        with tm.assert_produces_warning(SymbolWarning):
            self.assertRaises(RemoteDataError, web.get_data_google, sl, '2012')

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

                self.assertTrue(np.issubdtype(result.dtype, np.floating))
                result = pan.Open.ix['Jan-15-12':'Jan-20-12']
                self.assertEqual((4, 3), result.shape)
                assert_n_failed_equals_n_null_columns(w, result)

    @network
    def test_dtypes(self):
        #GH3995, #GH8980
        data = web.get_data_google('F', start='JAN-01-10', end='JAN-27-13')
        self.assertTrue(np.issubdtype(data.Open.dtype, np.number))
        self.assertTrue(np.issubdtype(data.Close.dtype, np.number))
        self.assertTrue(np.issubdtype(data.Low.dtype, np.number))
        self.assertTrue(np.issubdtype(data.High.dtype, np.number))
        self.assertTrue(np.issubdtype(data.Volume.dtype, np.number))

    @network
    def test_unicode_date(self):
        #GH8967
        data = web.get_data_google('F', start='JAN-01-10', end='JAN-27-13')
        self.assertEqual(data.index.name, 'Date')


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

        self.assertEqual(web.DataReader("F", 'yahoo', start, end)['Close'][-1],
                         13.68)

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
    def test_get_quote_string(self):
        _yahoo_codes.update({'MarketCap': 'j1'})
        df = web.get_quote_yahoo('GOOG')
        self.assertFalse(pd.isnull(df['MarketCap'][0]))

    @network
    def test_get_quote_stringlist(self):
        df = web.get_quote_yahoo(['GOOG', 'AAPL', 'GOOG'])
        assert_series_equal(df.ix[0], df.ix[2])

    @network
    def test_get_components_dow_jones(self):
        raise nose.SkipTest('unreliable test, receive partial components back for dow_jones')

        df = web.get_components_yahoo('^DJI') #Dow Jones
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 30)

    @network
    def test_get_components_dax(self):
        raise nose.SkipTest('unreliable test, receive partial components back for dax')

        df = web.get_components_yahoo('^GDAXI') #DAX
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 30)
        self.assertEqual(df[df.name.str.contains('adidas', case=False)].index,
                         'ADS.DE')

    @network
    def test_get_components_nasdaq_100(self):
        # as of 7/12/13 the conditional will test false because the link is invalid
        raise nose.SkipTest('unreliable test, receive partial components back for nasdaq_100')

        df = web.get_components_yahoo('^NDX') #NASDAQ-100
        self.assertIsInstance(df, pd.DataFrame)

        if len(df) > 1:
            # Usual culprits, should be around for a while
            self.assertTrue('AAPL' in df.index)
            self.assertTrue('GOOG' in df.index)
            self.assertTrue('AMZN' in df.index)
        else:
            expected = DataFrame({'exchange': 'N/A', 'name': '@^NDX'},
                                 index=['@^NDX'])
            assert_frame_equal(df, expected)

    @network
    def test_get_data_single_symbol(self):
        #single symbol
        #http://finance.yahoo.com/q/hp?s=GOOG&a=09&b=08&c=2010&d=09&e=10&f=2010&g=d
        # just test that we succeed
        web.get_data_yahoo('GOOG')

    @network
    def test_get_data_interval(self):
        # daily interval data
        pan = web.get_data_yahoo('XOM', '2013-01-01', '2013-12-31', interval='d')
        self.assertEqual(len(pan), 252)

        # weekly interval data
        pan = web.get_data_yahoo('XOM', '2013-01-01', '2013-12-31', interval='w')
        self.assertEqual(len(pan), 53)

        # montly interval data
        pan = web.get_data_yahoo('XOM', '2013-01-01', '2013-12-31', interval='m')
        self.assertEqual(len(pan), 12)

        # dividend data
        pan = web.get_data_yahoo('XOM', '2013-01-01', '2013-12-31', interval='v')
        self.assertEqual(len(pan), 4)

        # test fail on invalid interval
        self.assertRaises(ValueError, web.get_data_yahoo, 'XOM', interval='NOT VALID')

    @network
    def test_get_data_multiple_symbols(self):
        # just test that we succeed
        sl = ['AAPL', 'AMZN', 'GOOG']
        web.get_data_yahoo(sl, '2012')

    @network
    def test_get_data_multiple_symbols_two_dates(self):
        pan = web.get_data_yahoo(['GE', 'MSFT', 'INTC'], 'JAN-01-12',
                                 'JAN-31-12')
        result = pan.Close.ix['01-18-12']
        self.assertEqual(len(result), 3)

        # sanity checking
        self.assertTrue(np.issubdtype(result.dtype, np.floating))

        expected = np.array([[18.99,  28.4, 25.18],
                             [18.58, 28.31, 25.13],
                             [19.03, 28.16, 25.52],
                             [18.81, 28.82, 25.87]])
        result = pan.Open.ix['Jan-15-12':'Jan-20-12']
        self.assertEqual(expected.shape, result.shape)

    @network
    def test_get_date_ret_index(self):
        pan = web.get_data_yahoo(['GE', 'INTC', 'IBM'], '1977', '1987',
                                 ret_index=True)
        self.assertTrue(hasattr(pan, 'Ret_Index'))
        if hasattr(pan, 'Ret_Index') and hasattr(pan.Ret_Index, 'INTC'):
            tstamp = pan.Ret_Index.INTC.first_valid_index()
            result = pan.Ret_Index.ix[tstamp]['INTC']
            self.assertEqual(result, 1.0)

        # sanity checking
        self.assertTrue(np.issubdtype(pan.values.dtype, np.floating))


class TestYahooOptions(tm.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestYahooOptions, cls).setUpClass()
        _skip_if_no_lxml()
        _skip_if_no_bs()

        # aapl has monthlies
        cls.aapl = web.Options('aapl', 'yahoo')
        d = (Timestamp.today() + pd.offsets.MonthBegin(1)).normalize()
        cls.year = d.year
        cls.month = d.month
        cls.expiry = d
        cls.expiry2 = d + pd.offsets.MonthBegin(1)
        cls.dirpath = tm.get_data_path()
        cls.html1 = os.path.join(cls.dirpath, 'yahoo_options1.html')
        cls.html2 = os.path.join(cls.dirpath, 'yahoo_options2.html')
        cls.html3 = os.path.join(cls.dirpath, 'yahoo_options3.html') #Empty table GH#22
        cls.data1 = cls.aapl._option_frames_from_url(cls.html1)['puts']

    @classmethod
    def tearDownClass(cls):
        super(TestYahooOptions, cls).tearDownClass()
        del cls.aapl, cls.expiry

    @network
    def test_get_options_data(self):
        # regression test GH6105
        self.assertRaises(ValueError, self.aapl.get_options_data, month=3)
        self.assertRaises(ValueError, self.aapl.get_options_data, year=1992)

        try:
            options = self.aapl.get_options_data(expiry=self.expiry)
        except RemoteDataError as e:
            raise nose.SkipTest(e)
        self.assertTrue(len(options) > 1)

    @network
    def test_get_near_stock_price(self):
        try:
            options = self.aapl.get_near_stock_price(call=True, put=True,
                                                     expiry=[self.expiry,self.expiry2])
        except RemoteDataError as e:
            raise nose.SkipTest(e)
        self.assertTrue(len(options) > 1)

    @network
    def test_get_call_data(self):
        try:
            calls = self.aapl.get_call_data(expiry=self.expiry)
        except RemoteDataError as e:
            raise nose.SkipTest(e)
        self.assertTrue(len(calls) > 1)

    @network
    def test_get_put_data(self):
        try:
            puts = self.aapl.get_put_data(expiry=self.expiry)
        except RemoteDataError as e:
            raise nose.SkipTest(e)
        self.assertTrue(len(puts) > 1)

    @network
    def test_get_expiry_dates(self):
        try:
            dates, _ = self.aapl._get_expiry_dates_and_links()
        except RemoteDataError as e:
            raise nose.SkipTest(e)
        self.assertTrue(len(dates) > 1)

    @network
    def test_get_all_data(self):
        try:
            data = self.aapl.get_all_data(put=True)
        except RemoteDataError as e:
            raise nose.SkipTest(e)
        self.assertTrue(len(data) > 1)

    @network
    def test_get_data_with_list(self):
        try:
            data = self.aapl.get_call_data(expiry=self.aapl.expiry_dates)
        except RemoteDataError as e:
            raise nose.SkipTest(e)
        self.assertTrue(len(data) > 1)

    @network
    def test_get_all_data_calls_only(self):
        try:
            data = self.aapl.get_all_data(call=True, put=False)
        except RemoteDataError as e:
            raise nose.SkipTest(e)
        self.assertTrue(len(data) > 1)

    @network
    def test_get_underlying_price(self):
        #GH7
        try:
            options_object = web.Options('^spxpm', 'yahoo')
            url = options_object._yahoo_url_from_expiry(options_object.expiry_dates[0])
            root = options_object._parse_url(url)
            quote_price = options_object._underlying_price_from_root(root)
        except RemoteDataError as e:
            raise nose.SkipTest(e)
        self.assertIsInstance(quote_price, float)

    def test_sample_page_price_quote_time1(self):
        #Tests the weekend quote time format
        price, quote_time = self.aapl._underlying_price_and_time_from_url(self.html1)
        self.assertIsInstance(price, (int, float, complex))
        self.assertIsInstance(quote_time, (datetime, Timestamp))

    def test_chop(self):
        #regression test for #7625
        self.aapl.chop_data(self.data1, above_below=2, underlying_price=np.nan)
        chopped = self.aapl.chop_data(self.data1, above_below=2, underlying_price=100)
        self.assertIsInstance(chopped, DataFrame)
        self.assertTrue(len(chopped) > 1)

    def test_chop_out_of_strike_range(self):
        #regression test for #7625
        self.aapl.chop_data(self.data1, above_below=2, underlying_price=np.nan)
        chopped = self.aapl.chop_data(self.data1, above_below=2, underlying_price=100000)
        self.assertIsInstance(chopped, DataFrame)
        self.assertTrue(len(chopped) > 1)


    @network
    def test_sample_page_price_quote_time2(self):
        #Tests the EDT page format
        #regression test for #8741
        price, quote_time = self.aapl._underlying_price_and_time_from_url(self.html2)
        self.assertIsInstance(price, (int, float, complex))
        self.assertIsInstance(quote_time, (datetime, Timestamp))

    @network
    def test_sample_page_chg_float(self):
        #Tests that numeric columns with comma's are appropriately dealt with
        self.assertEqual(self.data1['Chg'].dtype, 'float64')

    @network
    def test_month_year(self):
        try:
            data = self.aapl.get_call_data(month=self.month, year=self.year)
        except RemoteDataError as e:
            raise nose.SkipTest(e)

        self.assertTrue(len(data) > 1)

    @network
    def test_empty_table(self):
        #GH22
        empty = self.aapl._option_frames_from_url(self.html3)['puts']
        self.assertTrue(len(empty) == 0)


class TestOptionsWarnings(tm.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestOptionsWarnings, cls).setUpClass()

    @classmethod
    def tearDownClass(cls):
        super(TestOptionsWarnings, cls).tearDownClass()

    @network
    def test_options_source_warning(self):
        with assert_produces_warning():
            aapl = web.Options('aapl')


class TestDataReader(tm.TestCase):
    def test_is_s3_url(self):
        from pandas.io.common import _is_s3_url
        self.assertTrue(_is_s3_url("s3://pandas/somethingelse.com"))

    @network
    def test_read_yahoo(self):
        gs = DataReader("GS", "yahoo")
        self.assertIsInstance(gs, DataFrame)

    @network
    def test_read_google(self):
        gs = DataReader("GS", "google")
        self.assertIsInstance(gs, DataFrame)

    @network
    def test_read_fred(self):
        vix = DataReader("VIXCLS", "fred")
        self.assertIsInstance(vix, DataFrame)

    @network
    def test_read_famafrench(self):
        for name in ("F-F_Research_Data_Factors",
                     "F-F_Research_Data_Factors_weekly", "6_Portfolios_2x3",
                     "F-F_ST_Reversal_Factor", "F-F_Momentum_Factor"):
            ff = DataReader(name, "famafrench")
            self.assertTrue(ff is not None)
            self.assertIsInstance(ff, dict)


class TestFred(tm.TestCase):
    @network
    def test_fred(self):

        # Throws an exception when DataReader can't get a 200 response from
        # FRED.

        start = datetime(2010, 1, 1)
        end = datetime(2013, 1, 27)

        received = web.DataReader("GDP", "fred", start, end)['GDP'].tail(1)[0]
        self.assertTrue(int(received) > 10000)

        self.assertRaises(Exception, web.DataReader, "NON EXISTENT SERIES",
                          'fred', start, end)

    @network
    def test_fred_nan(self):
        start = datetime(2010, 1, 1)
        end = datetime(2013, 1, 27)
        df = web.DataReader("DFII5", "fred", start, end)
        self.assertTrue(pd.isnull(df.ix['2010-01-01'][0]))

    @network
    def test_fred_parts(self):
        raise nose.SkipTest('buggy as of 2/18/14; maybe a data revision?')

        start = datetime(2010, 1, 1)
        end = datetime(2013, 1, 27)
        df = web.get_data_fred("CPIAUCSL", start, end)
        self.assertEqual(df.ix['2010-05-01'][0], 217.23)

        t = df.CPIAUCSL.values
        self.assertTrue(np.issubdtype(t.dtype, np.floating))
        self.assertEqual(t.shape, (37,))

    @network
    def test_fred_part2(self):
        expected = [[576.7],
                    [962.9],
                    [684.7],
                    [848.3],
                    [933.3]]
        result = web.get_data_fred("A09024USA144NNBR", start="1915").ix[:5]
        tm.assert_numpy_array_equal(result.values, np.array(expected))

    @network
    def test_invalid_series(self):
        name = "NOT A REAL SERIES"
        self.assertRaises(Exception, web.get_data_fred, name)

    @network
    def test_fred_multi(self):
        raise nose.SkipTest('buggy as of 2/18/14; maybe a data revision?')

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
