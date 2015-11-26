import os
from datetime import datetime

import warnings
import nose
import pandas as pd
from pandas import compat
from pandas.util.testing import network, assert_frame_equal, with_connectivity_check
from numpy.testing.decorators import slow
import pandas.util.testing as tm

if compat.PY3:
    raise nose.SkipTest("python-gflags does not support Python 3 yet")

try:
    import httplib2
    import apiclient

    # deprecated
    with warnings.catch_warnings(record=True):
        import pandas.io.ga as ga

    from pandas.io.ga import GAnalytics, read_ga
    from pandas.io.auth import AuthenticationConfigError, reset_default_token_store
    from pandas.io import auth
except ImportError:
    raise nose.SkipTest("need httplib2 and auth libs")


class TestGoogle(tm.TestCase):

    _multiprocess_can_split_ = True

    def test_remove_token_store(self):
        auth.DEFAULT_TOKEN_FILE = 'test.dat'
        with open(auth.DEFAULT_TOKEN_FILE, 'w') as fh:
            fh.write('test')

        reset_default_token_store()
        self.assertFalse(os.path.exists(auth.DEFAULT_TOKEN_FILE))

    @with_connectivity_check("http://www.google.com")
    def test_getdata(self):
        try:
            end_date = datetime.now()
            start_date = end_date - pd.offsets.Day() * 5
            end_date = end_date.strftime('%Y-%m-%d')
            start_date = start_date.strftime('%Y-%m-%d')

            reader = GAnalytics()
            df = reader.get_data(
                metrics=['avgTimeOnSite', 'visitors', 'newVisits',
                         'pageviewsPerVisit'],
                start_date=start_date,
                end_date=end_date,
                dimensions=['date', 'hour'],
                parse_dates={'ts': ['date', 'hour']},
                index_col=0)

            self.assertIsInstance(df, pd.DataFrame)
            self.assertIsInstance(df.index, pd.DatetimeIndex)
            self.assertGreater(len(df), 1)
            self.assertTrue('date' not in df)
            self.assertTrue('hour' not in df)
            self.assertEqual(df.index.name, 'ts')
            self.assertTrue('avgTimeOnSite' in df)
            self.assertTrue('visitors' in df)
            self.assertTrue('newVisits' in df)
            self.assertTrue('pageviewsPerVisit' in df)

            df2 = read_ga(
                metrics=['avgTimeOnSite', 'visitors', 'newVisits',
                         'pageviewsPerVisit'],
                start_date=start_date,
                end_date=end_date,
                dimensions=['date', 'hour'],
                parse_dates={'ts': ['date', 'hour']},
                index_col=0)

            assert_frame_equal(df, df2)

        except AuthenticationConfigError:
            raise nose.SkipTest("authentication error")

    @with_connectivity_check("http://www.google.com")
    def test_iterator(self):
        try:
            reader = GAnalytics()

            it = reader.get_data(
                metrics='visitors',
                start_date='2005-1-1',
                dimensions='date',
                max_results=10, chunksize=5,
                index_col=0)

            df1 = next(it)
            df2 = next(it)

            for df in [df1, df2]:
                self.assertIsInstance(df, pd.DataFrame)
                self.assertIsInstance(df.index, pd.DatetimeIndex)
                self.assertEqual(len(df), 5)
                self.assertTrue('date' not in df)
                self.assertEqual(df.index.name, 'date')
                self.assertTrue('visitors' in df)

            self.assertTrue((df2.index > df1.index).all())

        except AuthenticationConfigError:
            raise nose.SkipTest("authentication error")

    def test_v2_advanced_segment_format(self):
        advanced_segment_id = 1234567
        query = ga.format_query('google_profile_id', ['visits'], '2013-09-01', segment=advanced_segment_id)
        self.assertEqual(query['segment'], 'gaid::' + str(advanced_segment_id), "An integer value should be formatted as an advanced segment.")

    def test_v2_dynamic_segment_format(self):
        dynamic_segment_id = 'medium==referral'
        query = ga.format_query('google_profile_id', ['visits'], '2013-09-01', segment=dynamic_segment_id)
        self.assertEqual(query['segment'], 'dynamic::ga:' + str(dynamic_segment_id), "A string value with more than just letters and numbers should be formatted as a dynamic segment.")

    def test_v3_advanced_segment_common_format(self):
        advanced_segment_id = 'aZwqR234'
        query = ga.format_query('google_profile_id', ['visits'], '2013-09-01', segment=advanced_segment_id)
        self.assertEqual(query['segment'], 'gaid::' + str(advanced_segment_id), "A string value with just letters and numbers should be formatted as an advanced segment.")

    def test_v3_advanced_segment_weird_format(self):
        advanced_segment_id = '_aZwqR234-s1'
        query = ga.format_query('google_profile_id', ['visits'], '2013-09-01', segment=advanced_segment_id)
        self.assertEqual(query['segment'], 'gaid::' + str(advanced_segment_id), "A string value with just letters, numbers, and hyphens should be formatted as an advanced segment.")

    def test_v3_advanced_segment_with_underscore_format(self):
        advanced_segment_id = 'aZwqR234_s1'
        query = ga.format_query('google_profile_id', ['visits'], '2013-09-01', segment=advanced_segment_id)
        self.assertEqual(query['segment'], 'gaid::' + str(advanced_segment_id), "A string value with just letters, numbers, and underscores should be formatted as an advanced segment.")

    @with_connectivity_check("http://www.google.com")
    def test_segment(self):
        try:
            end_date = datetime.now()
            start_date = end_date - pd.offsets.Day() * 5
            end_date = end_date.strftime('%Y-%m-%d')
            start_date = start_date.strftime('%Y-%m-%d')

            reader = GAnalytics()
            df = reader.get_data(
                metrics=['avgTimeOnSite', 'visitors', 'newVisits',
                         'pageviewsPerVisit'],
                start_date=start_date,
                end_date=end_date,
                segment=-2,
                dimensions=['date', 'hour'],
                parse_dates={'ts': ['date', 'hour']},
                index_col=0)

            self.assertIsInstance(df, pd.DataFrame)
            self.assertIsInstance(df.index, pd.DatetimeIndex)
            self.assertGreater(len(df), 1)
            self.assertTrue('date' not in df)
            self.assertTrue('hour' not in df)
            self.assertEqual(df.index.name, 'ts')
            self.assertTrue('avgTimeOnSite' in df)
            self.assertTrue('visitors' in df)
            self.assertTrue('newVisits' in df)
            self.assertTrue('pageviewsPerVisit' in df)

            # dynamic
            df = read_ga(
                metrics=['avgTimeOnSite', 'visitors', 'newVisits',
                         'pageviewsPerVisit'],
                start_date=start_date,
                end_date=end_date,
                segment="source=~twitter",
                dimensions=['date', 'hour'],
                parse_dates={'ts': ['date', 'hour']},
                index_col=0)

            assert isinstance(df, pd.DataFrame)
            assert isinstance(df.index, pd.DatetimeIndex)
            self.assertGreater(len(df), 1)
            self.assertTrue('date' not in df)
            self.assertTrue('hour' not in df)
            self.assertEqual(df.index.name, 'ts')
            self.assertTrue('avgTimeOnSite' in df)
            self.assertTrue('visitors' in df)
            self.assertTrue('newVisits' in df)
            self.assertTrue('pageviewsPerVisit' in df)

        except AuthenticationConfigError:
            raise nose.SkipTest("authentication error")


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
