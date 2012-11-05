from pandas.util.py3compat import StringIO, BytesIO
from datetime import datetime
import csv
import os
import sys
import re
import unittest
import pandas.io.data as pd
import nose
from pandas.util.testing import network
import urllib2

class TestYahoo(unittest.TestCase):

    @network
    def test_yahoo(self):
        # asserts that yahoo is minimally working and that it throws
        # an excecption when DataReader can't get a 200 response from
        # yahoo
        start = datetime(2010,1,1)
        end = datetime(2012,1,24)

        try:
            self.assertEquals(
                pd.DataReader("F", 'yahoo', start, end)['Close'][-1],
                12.82)

            self.assertRaises(
                Exception,
                lambda: pd.DataReader("NON EXISTENT TICKER", 'yahoo',
                                      start, end))
        except urllib2.URLError:
            try:
                urllib2.urlopen('http://www.google.com')
            except urllib2.URLError:
                raise nose.SkipTest
            else:
                raise

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
