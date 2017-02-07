import pandas as pd
import pandas.util.testing as tm


class TestDatetimeIndex(tm.TestCase):

    def test_fillna_datetime64(self):
        # GH 11343
        for tz in ['US/Eastern', 'Asia/Tokyo']:
            idx = pd.DatetimeIndex(['2011-01-01 09:00', pd.NaT,
                                    '2011-01-01 11:00'])

            exp = pd.DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00',
                                    '2011-01-01 11:00'])
            self.assert_index_equal(
                idx.fillna(pd.Timestamp('2011-01-01 10:00')), exp)

            # tz mismatch
            exp = pd.Index([pd.Timestamp('2011-01-01 09:00'),
                            pd.Timestamp('2011-01-01 10:00', tz=tz),
                            pd.Timestamp('2011-01-01 11:00')], dtype=object)
            self.assert_index_equal(
                idx.fillna(pd.Timestamp('2011-01-01 10:00', tz=tz)), exp)

            # object
            exp = pd.Index([pd.Timestamp('2011-01-01 09:00'), 'x',
                            pd.Timestamp('2011-01-01 11:00')], dtype=object)
            self.assert_index_equal(idx.fillna('x'), exp)

            idx = pd.DatetimeIndex(['2011-01-01 09:00', pd.NaT,
                                    '2011-01-01 11:00'], tz=tz)

            exp = pd.DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00',
                                    '2011-01-01 11:00'], tz=tz)
            self.assert_index_equal(
                idx.fillna(pd.Timestamp('2011-01-01 10:00', tz=tz)), exp)

            exp = pd.Index([pd.Timestamp('2011-01-01 09:00', tz=tz),
                            pd.Timestamp('2011-01-01 10:00'),
                            pd.Timestamp('2011-01-01 11:00', tz=tz)],
                           dtype=object)
            self.assert_index_equal(
                idx.fillna(pd.Timestamp('2011-01-01 10:00')), exp)

            # object
            exp = pd.Index([pd.Timestamp('2011-01-01 09:00', tz=tz),
                            'x',
                            pd.Timestamp('2011-01-01 11:00', tz=tz)],
                           dtype=object)
            self.assert_index_equal(idx.fillna('x'), exp)
