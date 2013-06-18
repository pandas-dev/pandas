import unittest
import nose
from datetime import datetime

from pandas.util.py3compat import StringIO, BytesIO

import pandas as pd
import pandas.io.data as web
from pandas.util.testing import (network, assert_frame_equal,
                                 assert_series_equal,
                                 assert_almost_equal, with_connectivity_check)
from numpy.testing.decorators import slow

import urllib2


class TestFred(unittest.TestCase):

    @slow
    @with_connectivity_check("http://www.google.com")
    def test_fred(self):
        """
        Throws an exception when DataReader can't get a 200 response from
        FRED.
        """
        start = datetime(2010, 1, 1)
        end = datetime(2013, 01, 27)

        self.assertEquals(
            web.DataReader("GDP", "fred", start, end)['GDP'].tail(1),
            16004.5)

        self.assertRaises(
            Exception,
            lambda: web.DataReader("NON EXISTENT SERIES", 'fred',
                                   start, end))

    @slow
    @network
    def test_fred_nan(self):
        start = datetime(2010, 1, 1)
        end = datetime(2013, 01, 27)
        df = web.DataReader("DFII5", "fred", start, end)
        assert pd.isnull(df.ix['2010-01-01'])

    @slow
    @network
    def test_fred_parts(self):
        import numpy as np

        start = datetime(2010, 1, 1)
        end = datetime(2013, 01, 27)
        df = web.get_data_fred("CPIAUCSL", start, end)
        assert df.ix['2010-05-01'] == 217.23

        t = np.array(df.CPIAUCSL.tolist())
        assert np.issubdtype(t.dtype, np.floating)
        assert t.shape == (37,)

        # Test some older ones:
        expected = [[576.7],
                    [962.9],
                    [684.7],
                    [848.3],
                    [933.3]]
        result = web.get_data_fred("A09024USA144NNBR", start="1915").ix[:5]
        assert (result.values == expected).all()

    @slow
    @network
    def test_invalid_series(self):
        name = "NOT A REAL SERIES"
        self.assertRaises(Exception, web.get_data_fred, name)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
