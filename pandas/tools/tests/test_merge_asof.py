import nose
import os

import numpy as np
import pandas as pd
from pandas import (merge_asof, read_csv,
                    to_datetime, Timedelta)
from pandas.tools.merge import MergeError
from pandas.util import testing as tm
from pandas.util.testing import assert_frame_equal


class TestAsOfMerge(tm.TestCase):
    _multiprocess_can_split_ = True

    def read_data(self, name, dedupe=False):
        path = os.path.join(tm.get_data_path(), name)
        x = read_csv(path)
        if dedupe:
            x = (x.drop_duplicates(['time', 'ticker'], keep='last')
                  .reset_index(drop=True)
                 )
        x.time = to_datetime(x.time)
        return x

    def setUp(self):

        self.trades = self.read_data('trades.csv')
        self.quotes = self.read_data('quotes.csv', dedupe=True)
        self.asof = self.read_data('asof.csv')
        self.tolerance = self.read_data('tolerance.csv')
        self.allow_exact_matches = self.read_data('allow_exact_matches.csv')
        self.allow_exact_matches_and_tolerance = self.read_data(
            'allow_exact_matches_and_tolerance.csv')

    def test_examples1(self):
        """ doc-string examples """

        left = pd.DataFrame({'a': [1, 5, 10],
                             'left_val': ['a', 'b', 'c']})
        right = pd.DataFrame({'a': [1, 2, 3, 6, 7],
                              'right_val': [1, 2, 3, 6, 7]})

        pd.merge_asof(left, right, on='a')

    def test_examples2(self):
        """ doc-string examples """

        trades = pd.DataFrame({
            'time': pd.to_datetime(['20160525 13:30:00.023',
                                    '20160525 13:30:00.038',
                                    '20160525 13:30:00.048',
                                    '20160525 13:30:00.048',
                                    '20160525 13:30:00.048']),
            'ticker': ['MSFT', 'MSFT',
                       'GOOG', 'GOOG', 'AAPL'],
            'price': [51.95, 51.95,
                      720.77, 720.92, 98.00],
            'quantity': [75, 155,
                         100, 100, 100]},
            columns=['time', 'ticker', 'price', 'quantity'])

        quotes = pd.DataFrame({
            'time': pd.to_datetime(['20160525 13:30:00.023',
                                    '20160525 13:30:00.023',
                                    '20160525 13:30:00.030',
                                    '20160525 13:30:00.041',
                                    '20160525 13:30:00.048',
                                    '20160525 13:30:00.049',
                                    '20160525 13:30:00.072',
                                    '20160525 13:30:00.075']),
            'ticker': ['GOOG', 'MSFT', 'MSFT',
                       'MSFT', 'GOOG', 'AAPL', 'GOOG',
                       'MSFT'],
            'bid': [720.50, 51.95, 51.97, 51.99,
                    720.50, 97.99, 720.50, 52.01],
            'ask': [720.93, 51.96, 51.98, 52.00,
                    720.93, 98.01, 720.88, 52.03]},
            columns=['time', 'ticker', 'bid', 'ask'])

        pd.merge_asof(trades, quotes,
                      on='time',
                      by='ticker')

        pd.merge_asof(trades, quotes,
                      on='time',
                      by='ticker',
                      tolerance=pd.Timedelta('2ms'))

        pd.merge_asof(trades, quotes,
                      on='time',
                      by='ticker',
                      tolerance=pd.Timedelta('10ms'),
                      allow_exact_matches=False)

    def test_basic(self):

        expected = self.asof
        trades = self.trades
        quotes = self.quotes

        result = merge_asof(trades, quotes,
                            on='time',
                            by='ticker')
        assert_frame_equal(result, expected)

    def test_basic_categorical(self):

        expected = self.asof
        trades = self.trades.copy()
        trades.ticker = trades.ticker.astype('category')
        quotes = self.quotes.copy()
        quotes.ticker = quotes.ticker.astype('category')

        result = merge_asof(trades, quotes,
                            on='time',
                            by='ticker')
        assert_frame_equal(result, expected)

    def test_missing_right_by(self):

        expected = self.asof
        trades = self.trades
        quotes = self.quotes

        q = quotes[quotes.ticker != 'MSFT']
        result = merge_asof(trades, q,
                            on='time',
                            by='ticker')
        expected.loc[expected.ticker == 'MSFT', ['bid', 'ask']] = np.nan
        assert_frame_equal(result, expected)

    def test_basic2(self):

        expected = self.read_data('asof2.csv')
        trades = self.read_data('trades2.csv')
        quotes = self.read_data('quotes2.csv', dedupe=True)

        result = merge_asof(trades, quotes,
                            on='time',
                            by='ticker')
        assert_frame_equal(result, expected)

    def test_basic_no_by(self):
        f = lambda x: x[x.ticker == 'MSFT'].drop('ticker', axis=1) \
            .reset_index(drop=True)

        # just use a single ticker
        expected = f(self.asof)
        trades = f(self.trades)
        quotes = f(self.quotes)

        result = merge_asof(trades, quotes,
                            on='time')
        assert_frame_equal(result, expected)

    def test_valid_join_keys(self):

        trades = self.trades
        quotes = self.quotes

        with self.assertRaises(MergeError):
            merge_asof(trades, quotes,
                       left_on='time',
                       right_on='bid',
                       by='ticker')

        with self.assertRaises(MergeError):
            merge_asof(trades, quotes,
                       on=['time', 'ticker'],
                       by='ticker')

        with self.assertRaises(MergeError):
            merge_asof(trades, quotes,
                       by='ticker')

    def test_with_duplicates(self):

        q = pd.concat([self.quotes, self.quotes]).sort_values(
            ['time', 'ticker']).reset_index(drop=True)
        result = merge_asof(self.trades, q,
                            on='time',
                            by='ticker')
        expected = self.read_data('asof.csv')
        assert_frame_equal(result, expected)

    def test_with_duplicates_no_on(self):

        df1 = pd.DataFrame({'key': [1, 1, 3],
                            'left_val': [1, 2, 3]})
        df2 = pd.DataFrame({'key': [1, 2, 2],
                            'right_val': [1, 2, 3]})
        result = merge_asof(df1, df2, on='key')
        expected = pd.DataFrame({'key': [1, 1, 3],
                                 'left_val': [1, 2, 3],
                                 'right_val': [1, 1, 3]})
        assert_frame_equal(result, expected)

    def test_valid_allow_exact_matches(self):

        trades = self.trades
        quotes = self.quotes

        with self.assertRaises(MergeError):
            merge_asof(trades, quotes,
                       on='time',
                       by='ticker',
                       allow_exact_matches='foo')

    def test_valid_tolerance(self):

        trades = self.trades
        quotes = self.quotes

        # dti
        merge_asof(trades, quotes,
                   on='time',
                   by='ticker',
                   tolerance=Timedelta('1s'))

        # integer
        merge_asof(trades.reset_index(), quotes.reset_index(),
                   on='index',
                   by='ticker',
                   tolerance=1)

        # incompat
        with self.assertRaises(MergeError):
            merge_asof(trades, quotes,
                       on='time',
                       by='ticker',
                       tolerance=1)

        # invalid
        with self.assertRaises(MergeError):
            merge_asof(trades.reset_index(), quotes.reset_index(),
                       on='index',
                       by='ticker',
                       tolerance=1.0)

        # invalid negative
        with self.assertRaises(MergeError):
            merge_asof(trades, quotes,
                       on='time',
                       by='ticker',
                       tolerance=-Timedelta('1s'))

        with self.assertRaises(MergeError):
            merge_asof(trades.reset_index(), quotes.reset_index(),
                       on='index',
                       by='ticker',
                       tolerance=-1)

    def test_non_sorted(self):

        trades = self.trades.sort_values('time', ascending=False)
        quotes = self.quotes.sort_values('time', ascending=False)

        # we require that we are already sorted on time & quotes
        self.assertFalse(trades.time.is_monotonic)
        self.assertFalse(quotes.time.is_monotonic)
        with self.assertRaises(ValueError):
            merge_asof(trades, quotes,
                       on='time',
                       by='ticker')

        trades = self.trades.sort_values('time')
        self.assertTrue(trades.time.is_monotonic)
        self.assertFalse(quotes.time.is_monotonic)
        with self.assertRaises(ValueError):
            merge_asof(trades, quotes,
                       on='time',
                       by='ticker')

        quotes = self.quotes.sort_values('time')
        self.assertTrue(trades.time.is_monotonic)
        self.assertTrue(quotes.time.is_monotonic)

        # ok, though has dupes
        merge_asof(trades, self.quotes,
                   on='time',
                   by='ticker')

    def test_tolerance(self):

        trades = self.trades
        quotes = self.quotes

        result = merge_asof(trades, quotes,
                            on='time',
                            by='ticker',
                            tolerance=Timedelta('1day'))
        expected = self.tolerance
        assert_frame_equal(result, expected)

    def test_allow_exact_matches(self):

        result = merge_asof(self.trades, self.quotes,
                            on='time',
                            by='ticker',
                            allow_exact_matches=False)
        expected = self.allow_exact_matches
        assert_frame_equal(result, expected)

    def test_allow_exact_matches_and_tolerance(self):

        result = merge_asof(self.trades, self.quotes,
                            on='time',
                            by='ticker',
                            tolerance=Timedelta('100ms'),
                            allow_exact_matches=False)
        expected = self.allow_exact_matches_and_tolerance
        assert_frame_equal(result, expected)

    def test_allow_exact_matches_and_tolerance2(self):
        # GH 13695
        df1 = pd.DataFrame({
            'time': pd.to_datetime(['2016-07-15 13:30:00.030']),
            'username': ['bob']})
        df2 = pd.DataFrame({
            'time': pd.to_datetime(['2016-07-15 13:30:00.000',
                                    '2016-07-15 13:30:00.030']),
            'version': [1, 2]})

        result = pd.merge_asof(df1, df2, on='time')
        expected = pd.DataFrame({
            'time': pd.to_datetime(['2016-07-15 13:30:00.030']),
            'username': ['bob'],
            'version': [2]})
        assert_frame_equal(result, expected)

        result = pd.merge_asof(df1, df2, on='time', allow_exact_matches=False)
        expected = pd.DataFrame({
            'time': pd.to_datetime(['2016-07-15 13:30:00.030']),
            'username': ['bob'],
            'version': [1]})
        assert_frame_equal(result, expected)

        result = pd.merge_asof(df1, df2, on='time', allow_exact_matches=False,
                               tolerance=pd.Timedelta('10ms'))
        expected = pd.DataFrame({
            'time': pd.to_datetime(['2016-07-15 13:30:00.030']),
            'username': ['bob'],
            'version': [np.nan]})
        assert_frame_equal(result, expected)

    def test_allow_exact_matches_and_tolerance3(self):
        # GH 13709
        df1 = pd.DataFrame({
            'time': pd.to_datetime(['2016-07-15 13:30:00.030',
                                   '2016-07-15 13:30:00.030']),
            'username': ['bob', 'charlie']})
        df2 = pd.DataFrame({
            'time': pd.to_datetime(['2016-07-15 13:30:00.000',
                                    '2016-07-15 13:30:00.030']),
            'version': [1, 2]})

        result = pd.merge_asof(df1, df2, on='time', allow_exact_matches=False,
                               tolerance=pd.Timedelta('10ms'))
        expected = pd.DataFrame({
            'time': pd.to_datetime(['2016-07-15 13:30:00.030',
                                    '2016-07-15 13:30:00.030']),
            'username': ['bob', 'charlie'],
            'version': [np.nan, np.nan]})
        assert_frame_equal(result, expected)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
