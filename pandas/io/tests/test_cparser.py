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

import pandas.util.testing as tm

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

    def test_file_handle_mmap(self):
        try:
            f = open(self.csv1, 'rb')
            reader = parser.TextReader(f, memory_map=True)
            result = reader.read()
        finally:
            f.close()

    def test_StringIO(self):
        text = open(self.csv1, 'rb').read()
        reader = parser.TextReader(BytesIO(text))
        result = reader.read()

    def test_string_factorize(self):
        # should this be optional?
        data = 'a\nb\na\nb\na'
        reader = parser.TextReader(StringIO(data))
        result = reader.read()
        self.assert_(len(set(map(id, result[0]))) == 2)

    def test_skipinitialspace(self):
        data = ('a,   b\n'
                'a,   b\n'
                'a,   b\n'
                'a,   b')

        reader = parser.TextReader(StringIO(data), skipinitialspace=True)
        result = reader.read()

        self.assert_(np.array_equal(result[0], ['a', 'a', 'a', 'a']))
        self.assert_(np.array_equal(result[1], ['b', 'b', 'b', 'b']))

    def test_parse_booleans(self):
        data = 'True\nFalse\nTrue\nTrue'

        reader = parser.TextReader(StringIO(data))
        result = reader.read()

        self.assert_(result[0].dtype == np.bool_)

    def test_delimit_whitespace(self):
        data = 'a  b\na\t\t "b"\n"a"\t \t b'

        reader = parser.TextReader(StringIO(data), delim_whitespace=True)
        result = reader.read()

        self.assert_(np.array_equal(result[0], ['a', 'a', 'a']))
        self.assert_(np.array_equal(result[1], ['b', 'b', 'b']))

    def test_embedded_newline(self):
        data = 'a\n"hello\nthere"\nthis'

        reader = parser.TextReader(StringIO(data))
        result = reader.read()

        expected = ['a', 'hello\nthere', 'this']
        self.assert_(np.array_equal(result[0], expected))

    def test_euro_decimal(self):
        data = '12345,67\n345,678'

        reader = parser.TextReader(StringIO(data), delimiter=':',
                                   decimal=',')
        result = reader.read()

        expected = [12345.67, 345.678]
        tm.assert_almost_equal(result[0], expected)

    def test_integer_thousands(self):
        data = '123,456\n12,500'

        reader = parser.TextReader(StringIO(data), delimiter=':',
                                   thousands=',')
        result = reader.read()

        expected = [123456, 12500]
        tm.assert_almost_equal(result[0], expected)

    def test_na_substitution(self):
        pass

if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

