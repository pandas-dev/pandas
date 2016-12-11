import numpy as np

import pandas.util.testing as tm
from pandas import DataFrame, DatetimeIndex

from pandas.util.testing import assert_frame_equal, assert_series_equal


class TestTimeSeriesPartialSlices(tm.TestCase):
    _multiprocess_can_split_ = True

    def assert_exact(self, df, ts, value):
        element = df['a'][ts]

        # Series should return scalar
        self.assertIsInstance(element, np.int64)
        self.assertEqual(element, value)

        # Frame should raise (exact match)
        self.assertRaises(KeyError, df.__getitem__, ts)

        # TODO: test falling to column selection

    def assert_slice(self, df, ts, the_slice):
        # Series should return slice
        expected = df['a'][the_slice]
        assert_series_equal(df['a'][ts], expected)

        # Frame should return slice as well
        expected = df[the_slice]
        assert_frame_equal(df[ts], expected)

    def assert_key_error(self, df, ts):
        self.assertRaises(KeyError, df['a'].__getitem__, ts)
        self.assertRaises(KeyError, df.__getitem__, ts)

    def test_partial_slices_day(self):
        df = DataFrame({'a': [1, 2, 3]}, DatetimeIndex(['2011-12-31',
                                                        '2012-01-01',
                                                        '2012-01-02']),
                       dtype=np.int64)

        self.assertEqual(df.index.resolution, 'day')

        # Timestamp with resolution 'day'
        self.assert_exact(df, '2011-12-31', 1)
        self.assert_exact(df, '2012-01-01', 2)
        self.assert_exact(df, '2012-01-02', 3)

        # Timestamp with resolution less precise than 'day'
        for ts in ['2011', '2011-12']:
            self.assert_slice(df, ts, slice(None, 1))

        # The same as previous but several elements in the slice
        for ts in ['2012', '2012-01']:
            self.assert_slice(df, ts, slice(1, None))

        # Timestamp with resolution more precise than 'day'
        # Compatible with existing key
        for ts in ['2012-01-01 00', '2012-01-01 00:00',
                   '2012-01-01 00:00:00']:
            self.assert_exact(df, ts, 2)

        # Timestamp with resolution more precise than 'day'
        # Not compatible with existing key
        for ts in ['2012-01-01 01', '2012-01-01 00:01',
                   '2012-01-01 00:00:01']:
            self.assert_key_error(df, ts)

    def test_partial_slice_hour(self):
        df = DataFrame({'a': [1, 2, 3]}, DatetimeIndex(['2011-12-31 23',
                                                        '2012-01-01 00',
                                                        '2012-01-01 01']),
                       dtype=np.int64)

        self.assertEqual(df.index.resolution, 'hour')

        # Timestamp with resolution 'hour'
        self.assert_exact(df, '2011-12-31 23', 1)
        self.assert_exact(df, '2012-01-01 00', 2)
        self.assert_exact(df, '2012-01-01 01', 3)

        # Timestamp with resolution less precise than 'hour'
        for ts in ['2011', '2011-12', '2011-12-31']:
            self.assert_slice(df, ts, slice(None, 1))

        # The same as previous but several elements in the slice
        for ts in ['2012', '2012-01', '2012-01-01']:
            self.assert_slice(df, ts, slice(1, None))

        # Timestamp with resolution more precise than 'hour'
        # Compatible with existing key
        for ts in ['2012-01-01 00:00',
                   '2012-01-01 00:00:00']:
            self.assert_exact(df, ts, 2)

        # Timestamp with resolution more precise than 'hour'
        # Not compatible with existing key
        for ts in ['2012-01-01 00:01',
                   '2012-01-01 00:00:01']:
            self.assert_key_error(df, ts)

    def test_partial_slice_minute(self):
        df = DataFrame({'a': [1, 2, 3]},
                       DatetimeIndex(['2011-12-31 23:59',
                                      '2012-01-01 00:00',
                                      '2012-01-01 00:01']),
                       dtype=np.int64)

        self.assertEqual(df.index.resolution, 'minute')

        # Timestamp with resolution 'minute'
        self.assert_exact(df, '2011-12-31 23:59', 1)
        self.assert_exact(df, '2012-01-01 00:00', 2)
        self.assert_exact(df, '2012-01-01 00:01', 3)

        # Timestamp with resolution less precise than 'minute'
        for ts in ['2011', '2011-12', '2011-12-31',
                   '2011-12-31 23']:
            self.assert_slice(df, ts, slice(None, 1))

        # The same as previous but several elements in the slice
        for ts in ['2012', '2012-01', '2012-01-01',
                   '2012-01-01 00']:
            self.assert_slice(df, ts, slice(1, None))

        # Timestamp with resolution more precise than 'minute'
        # Compatible with existing key
        for ts in ['2012-01-01 00:00:00']:
            self.assert_exact(df, ts, 2)

        # Timestamp with resolution more precise than 'minute'
        # Not compatible with existing key
        for ts in ['2012-01-01 00:00:01']:
            self.assert_key_error(df, ts)

    def test_partial_slice_second(self):
        # See GH14826
        df = DataFrame({'a': [1, 2, 3]},
                       DatetimeIndex(['2011-12-31 23:59:59',
                                      '2012-01-01 00:00:00',
                                      '2012-01-01 00:00:01']),
                       dtype=np.int64)

        self.assertEqual(df.index.resolution, 'second')

        # Timestamp with resolution 'second'
        self.assert_exact(df, '2011-12-31 23:59:59', 1)
        self.assert_exact(df, '2012-01-01 00:00:00', 2)
        self.assert_exact(df, '2012-01-01 00:00:01', 3)

        # Timestamp with resolution less precise than 'second'
        for ts in ['2011', '2011-12', '2011-12-31',
                   '2011-12-31 23', '2011-12-31 23:59']:
            self.assert_slice(df, ts, slice(None, 1))

        # The same as previous but several elements in the slice
        for ts in ['2012', '2012-01', '2012-01-01',
                   '2012-01-01 00', '2012-01-01 00:00']:
            self.assert_slice(df, ts, slice(1, None))

        # Not possible to create a string that represents timestamp
        # that is more exact then 'second'
