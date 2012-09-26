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

from pandas._parser import TextReader
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

    def test_file_handle(self):
        try:
            f = open(self.csv1, 'rb')
            reader = TextReader(f)
            result = reader.read()
        finally:
            f.close()

    def test_string_filename(self):
        reader = TextReader(self.csv1, header=None)
        result = reader.read()

    def test_file_handle_mmap(self):
        try:
            f = open(self.csv1, 'rb')
            reader = TextReader(f, memory_map=True, header=None)
            result = reader.read()
        finally:
            f.close()

    def test_StringIO(self):
        text = open(self.csv1, 'rb').read()
        reader = TextReader(BytesIO(text), header=None)
        result = reader.read()

    def test_string_factorize(self):
        # should this be optional?
        data = 'a\nb\na\nb\na'
        reader = TextReader(StringIO(data), header=None)
        result = reader.read()
        self.assert_(len(set(map(id, result[0]))) == 2)

    def test_skipinitialspace(self):
        data = ('a,   b\n'
                'a,   b\n'
                'a,   b\n'
                'a,   b')

        reader = TextReader(StringIO(data), skipinitialspace=True,
                            header=None)
        result = reader.read()

        self.assert_(np.array_equal(result[0], ['a', 'a', 'a', 'a']))
        self.assert_(np.array_equal(result[1], ['b', 'b', 'b', 'b']))

    def test_parse_booleans(self):
        data = 'True\nFalse\nTrue\nTrue'

        reader = TextReader(StringIO(data), header=None)
        result = reader.read()

        self.assert_(result[0].dtype == np.bool_)

    def test_delimit_whitespace(self):
        data = 'a  b\na\t\t "b"\n"a"\t \t b'

        reader = TextReader(StringIO(data), delim_whitespace=True,
                            header=None)
        result = reader.read()

        self.assert_(np.array_equal(result[0], ['a', 'a', 'a']))
        self.assert_(np.array_equal(result[1], ['b', 'b', 'b']))

    def test_embedded_newline(self):
        data = 'a\n"hello\nthere"\nthis'

        reader = TextReader(StringIO(data), header=None)
        result = reader.read()

        expected = ['a', 'hello\nthere', 'this']
        self.assert_(np.array_equal(result[0], expected))

    def test_euro_decimal(self):
        data = '12345,67\n345,678'

        reader = TextReader(StringIO(data), delimiter=':',
                            decimal=',', header=None)
        result = reader.read()

        expected = [12345.67, 345.678]
        tm.assert_almost_equal(result[0], expected)

    def test_integer_thousands(self):
        data = '123,456\n12,500'

        reader = TextReader(StringIO(data), delimiter=':',
                            thousands=',', header=None)
        result = reader.read()

        expected = [123456, 12500]
        tm.assert_almost_equal(result[0], expected)

    def test_skip_bad_lines(self):
        data = ('a:b:c\n'
                'd:e:f\n'
                'g:h:i\n'
                'j:k\n'
                'l:m:n')

        reader = TextReader(StringIO(data), delimiter=':',
                            header=None)
        self.assertRaises(parser.CParserError, reader.read)

        reader = TextReader(StringIO(data), delimiter=':',
                            header=None,
                            error_bad_lines=False,
                            warn_bad_lines=False)
        result = reader.read()
        expected = {0: ['a', 'd', 'g', 'l'],
                    1: ['b', 'e', 'h', 'm'],
                    2: ['c', 'f', 'i', 'n']}
        assert_array_dicts_equal(result, expected)

    def test_header_not_enough_lines(self):
        data = ('skip this\n'
                'skip this\n'
                'a,b,c\n'
                '1,2,3\n'
                '4,5,6')

        reader = TextReader(StringIO(data), delimiter=',', header=2,
                            as_recarray=True)
        header = reader.header
        expected = ['a', 'b', 'c']
        self.assertEquals(header, expected)

        recs = reader.read()
        expected = {'a': [1, 4], 'b': [2, 5], 'c': [3, 6]}
        assert_array_dicts_equal(expected, recs)

        # not enough rows
        self.assertRaises(parser.CParserError, TextReader, StringIO(data),
                          delimiter=',', header=5, as_recarray=True)

    def test_escapechar(self):
        data = ('\\"hello world\"\n'
                '\\"hello world\"\n'
                '\\"hello world\"')

        reader = TextReader(StringIO(data), delimiter=',', header=None,
                            escapechar='\\')
        result = reader.read()
        expected = {0: ['"hello world"'] * 3}
        assert_array_dicts_equal(result, expected)

    def test_eof_has_eol(self):
        # handling of new line at EOF
        pass

    def test_na_substitution(self):
        pass

def assert_array_dicts_equal(left, right):
    for k, v in left.iteritems():
        assert(np.array_equal(v, right[k]))

if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

