import nose
import unittest
import pandas as pd
from datetime import datetime


class TestTsUtil(unittest.TestCase):
    def test_min_valid(self):
        # Ensure that Timestamp.min is a valid Timestamp
        pd.Timestamp(pd.Timestamp.min)

    def test_max_valid(self):
        # Ensure that Timestamp.max is a valid Timestamp
        pd.Timestamp(pd.Timestamp.max)

    def test_datetime_convert_roundtrip(self):
        # Ensure that min/max can be converted to python datetime and back without loss
        self.assertEqual(pd.Timestamp.max, pd.Timestamp(pd.Timestamp.max.to_pydatetime()))
        self.assertEqual(pd.Timestamp.min, pd.Timestamp(pd.Timestamp.min.to_pydatetime()))

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)