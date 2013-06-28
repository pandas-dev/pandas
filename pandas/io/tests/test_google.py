import unittest
import nose
from datetime import datetime

import numpy as np
import pandas as pd
import pandas.io.data as web
from pandas.util.testing import network, with_connectivity_check


class TestGoogle(unittest.TestCase):

    @network
    def test_google(self):
        # asserts that google is minimally working and that it throws
        # an exception when DataReader can't get a 200 response from
        # google
        start = datetime(2010, 1, 1)
        end = datetime(2013, 01, 27)

        self.assertEquals(
            web.DataReader("F", 'google', start, end)['Close'][-1],
            13.68)

        self.assertRaises(Exception, web.DataReader, "NON EXISTENT TICKER",
                          'google', start, end)

    @network
    def test_get_quote_fails(self):
        self.assertRaises(NotImplementedError, web.get_quote_google,
                          pd.Series(['GOOG', 'AAPL', 'GOOG']))

    @network
    def test_get_goog_volume(self):
        df = web.get_data_google('GOOG')
        self.assertEqual(df.Volume.ix['OCT-08-2010'], 2863473)

    @network
    def test_get_multi1(self):
        sl = ['AAPL', 'AMZN', 'GOOG']
        pan = web.get_data_google(sl, '2012')

        def testit():
            ts = pan.Close.GOOG.index[pan.Close.AAPL > pan.Close.GOOG]
            self.assertEquals(ts[0].dayofyear, 96)

        if (hasattr(pan, 'Close') and hasattr(pan.Close, 'GOOG') and
            hasattr(pan.Close, 'AAPL')):
            testit()
        else:
            self.assertRaises(AttributeError, testit)

    @network
    def test_get_multi2(self):
        pan = web.get_data_google(['GE', 'MSFT', 'INTC'], 'JAN-01-12',
                                  'JAN-31-12')
        result = pan.Close.ix['01-18-12']
        self.assertEqual(len(result), 3)

        # sanity checking
        assert np.issubdtype(result.dtype, np.floating)

        expected = np.array([[ 18.99,  28.4 ,  25.18],
                             [ 18.58,  28.31,  25.13],
                             [ 19.03,  28.16,  25.52],
                             [ 18.81,  28.82,  25.87]])
        result = pan.Open.ix['Jan-15-12':'Jan-20-12']
        self.assertEqual(np.array(expected).shape, result.shape)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
