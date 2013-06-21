import os
import unittest
from datetime import datetime

import nose
import pandas as pd
from pandas import DataFrame
from pandas.util.testing import network, assert_frame_equal, with_connectivity_check
from numpy.testing.decorators import slow

try:
    import httplib2
    from pandas.io.ga import GAnalytics, read_ga
    from pandas.io.auth import AuthenticationConfigError, reset_token_store
    from pandas.io import auth
except ImportError:
    raise nose.SkipTest

class TestGoogle(unittest.TestCase):

    _multiprocess_can_split_ = True

    def test_remove_token_store(self):
        auth.DEFAULT_TOKEN_FILE = 'test.dat'
        with open(auth.DEFAULT_TOKEN_FILE, 'w') as fh:
            fh.write('test')

        reset_token_store()
        self.assert_(not os.path.exists(auth.DEFAULT_TOKEN_FILE))

    @slow
    @network
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
                parse_dates={'ts': ['date', 'hour']})

            assert isinstance(df, DataFrame)
            assert isinstance(df.index, pd.DatetimeIndex)
            assert len(df) > 1
            assert 'date' not in df
            assert 'hour' not in df
            assert df.index.name == 'ts'
            assert 'avgTimeOnSite' in df
            assert 'visitors' in df
            assert 'newVisits' in df
            assert 'pageviewsPerVisit' in df

            df2 = read_ga(
                metrics=['avgTimeOnSite', 'visitors', 'newVisits',
                         'pageviewsPerVisit'],
                start_date=start_date,
                end_date=end_date,
                dimensions=['date', 'hour'],
                parse_dates={'ts': ['date', 'hour']})

            assert_frame_equal(df, df2)

        except AuthenticationConfigError:
            raise nose.SkipTest

    @slow
    @with_connectivity_check("http://www.google.com")
    def test_iterator(self):
        try:
            reader = GAnalytics()

            it = reader.get_data(
                metrics='visitors',
                start_date='2005-1-1',
                dimensions='date',
                max_results=10, chunksize=5)

            df1 = it.next()
            df2 = it.next()

            for df in [df1, df2]:
                assert isinstance(df, DataFrame)
                assert isinstance(df.index, pd.DatetimeIndex)
                assert len(df) == 5
                assert 'date' not in df
                assert df.index.name == 'date'
                assert 'visitors' in df

            assert (df2.index > df1.index).all()

        except AuthenticationConfigError:
            raise nose.SkipTest

    @slow
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
                parse_dates={'ts': ['date', 'hour']})

            assert isinstance(df, DataFrame)
            assert isinstance(df.index, pd.DatetimeIndex)
            assert len(df) > 1
            assert 'date' not in df
            assert 'hour' not in df
            assert df.index.name == 'ts'
            assert 'avgTimeOnSite' in df
            assert 'visitors' in df
            assert 'newVisits' in df
            assert 'pageviewsPerVisit' in df

            #dynamic
            df = read_ga(
                metrics=['avgTimeOnSite', 'visitors', 'newVisits',
                         'pageviewsPerVisit'],
                start_date=start_date,
                end_date=end_date,
                segment="source=~twitter",
                dimensions=['date', 'hour'],
                parse_dates={'ts': ['date', 'hour']})

            assert isinstance(df, DataFrame)
            assert isinstance(df.index, pd.DatetimeIndex)
            assert len(df) > 1
            assert 'date' not in df
            assert 'hour' not in df
            assert df.index.name == 'ts'
            assert 'avgTimeOnSite' in df
            assert 'visitors' in df
            assert 'newVisits' in df
            assert 'pageviewsPerVisit' in df

        except AuthenticationConfigError:
            raise nose.SkipTest

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
