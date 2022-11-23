import pandas as pd
import unittest

class MyTestCase(unittest.TestCase):

    def test_1_infer(self):
        with self.assertRaises(ValueError):
            ts = pd.Timestamp("2020-10-25 00:00", tz='Europe/Berlin')
            # This works:
            ts.floor('H')
            # But this gives an error:
            ts.floor('H', ambiguous='infer')

    def test_2_True(self):
        ts = pd.Timestamp("2020-10-25 00:00", tz='Europe/Berlin')
        actual = ts.floor('H', ambiguous=True)
        expected = ts.floor('H')
        self.assertEqual(actual, expected)

    def test_3_False(self):
        ts = pd.Timestamp("2020-10-25 00:00", tz='Europe/Berlin')
        actual = ts.floor('H', ambiguous=False)
        expected = ts.floor('H')
        self.assertEqual(actual, expected)

    def test_4_NaT(self):
        ts = pd.Timestamp("2020-10-25 00:00", tz='Europe/Berlin')
        actual = ts.floor('H', ambiguous='NaT')
        expected = ts.floor('H')
        self.assertEqual(actual, expected)

    def test_5_None(self):
        ts = pd.Timestamp("2020-10-25 00:00", tz='Europe/Berlin')
        actual = str(ts.floor('H'))
        expected = '2020-10-25 00:00:00+02:00'
        self.assertEqual(actual, expected)

    def test_6_other(self):
        with self.assertRaises(ValueError):
            ts = pd.Timestamp("2020-10-25 00:00", tz='Europe/Berlin')
            ts.floor('H', ambiguous='other')


if __name__ == '__main__':
    unittest.main()
