import unittest
import nose
from datetime import datetime

import pandas as pd
import numpy as np
import pandas.io.data as web
from pandas.util.testing import network
from numpy.testing import assert_array_equal


class TestFred(unittest.TestCase):
    @network
    def test_fred(self):
        """
        Throws an exception when DataReader can't get a 200 response from
        FRED.
        """
        start = datetime(2010, 1, 1)
        end = datetime(2013, 01, 27)

        self.assertEquals(
            web.DataReader("GDP", "fred", start, end)['GDP'].tail(1),
            15984.1)

        self.assertRaises(Exception, web.DataReader, "NON EXISTENT SERIES",
                          'fred', start, end)

    @network
    def test_fred_nan(self):
        start = datetime(2010, 1, 1)
        end = datetime(2013, 01, 27)
        df = web.DataReader("DFII5", "fred", start, end)
        assert pd.isnull(df.ix['2010-01-01'])

    @network
    def test_fred_parts(self):
        start = datetime(2010, 1, 1)
        end = datetime(2013, 01, 27)
        df = web.get_data_fred("CPIAUCSL", start, end)
        self.assertEqual(df.ix['2010-05-01'], 217.23)

        t = df.CPIAUCSL.values
        assert np.issubdtype(t.dtype, np.floating)
        self.assertEqual(t.shape, (37,))

    @network
    def test_fred_part2(self):
        expected = [[576.7],
                    [962.9],
                    [684.7],
                    [848.3],
                    [933.3]]
        result = web.get_data_fred("A09024USA144NNBR", start="1915").ix[:5]
        assert_array_equal(result.values, np.array(expected))

    @network
    def test_invalid_series(self):
        name = "NOT A REAL SERIES"
        self.assertRaises(Exception, web.get_data_fred, name)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
