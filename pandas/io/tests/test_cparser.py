"""
C/Cython ascii file parser tests
"""

from pandas.util.py3compat import StringIO, BytesIO
from datetime import datetime
import csv
import os
import sys
import re
import unittest

import nose

from numpy import nan
import numpy as np

from pandas import DataFrame, Series, Index, isnull, MultiIndex
import pandas.io.parsers as parsers
from pandas.io.parsers import (read_csv, read_table, read_fwf,
                               ExcelFile, TextParser)
from pandas.util.testing import (assert_almost_equal, assert_frame_equal,
                                 assert_series_equal, network)
import pandas.lib as lib
from pandas.util import py3compat
from pandas.lib import Timestamp
from pandas.tseries.index import date_range


import pandas._parser as parser


def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth

class TestCParser(unittest.TestCase):

    def setUp(self):
        self.dirpath = curpath()
        self.csv1 = os.path.join(self.dirpath, 'test1.csv')
        self.csv2 = os.path.join(self.dirpath, 'test2.csv')
        self.xls1 = os.path.join(self.dirpath, 'test.xls')

    def test_string_filename(self):
        reader = parser.TextReader(self.csv1)
        result = reader.read()

    def test_file_handle(self):
        try:
            f = open(self.csv1, 'rb')
            reader = parser.TextReader(f)
            result = reader.read()
        finally:
            f.close()

    # def test_StringIO(self):
    #     text = open(self.csv1, 'rb').read()

    #     reader = parser.TextReader(BytesIO(text))
    #     result = reader.read()

if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

