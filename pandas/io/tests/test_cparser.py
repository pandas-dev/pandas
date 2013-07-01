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
                               TextParser)
from pandas.util.testing import (assert_almost_equal, assert_frame_equal,
                                 assert_series_equal, network)
import pandas.lib as lib
from pandas.util import py3compat
from pandas.lib import Timestamp

import pandas.util.testing as tm

from pandas.parser import TextReader
import pandas.parser as parser


class TestCParser(unittest.TestCase):

    def setUp(self):
        self.dirpath = tm.get_data_path()
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
        src = BytesIO(text)
        reader = TextReader(src, header=None)
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
        # too many lines, see #2430 for why
        data = ('a:b:c\n'
                'd:e:f\n'
                'g:h:i\n'
                'j:k:l:m\n'
                'l:m:n\n'
                'o:p:q:r')

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

        stderr = sys.stderr
        sys.stderr = StringIO()
        try:
            reader = TextReader(StringIO(data), delimiter=':',
                                header=None,
                                error_bad_lines=False,
                                warn_bad_lines=True)
            reader.read()
            val = sys.stderr.getvalue()
            self.assertTrue('Skipping line 4' in val)
            self.assertTrue('Skipping line 6' in val)
        finally:
            sys.stderr = stderr

    def test_header_not_enough_lines(self):
        data = ('skip this\n'
                'skip this\n'
                'a,b,c\n'
                '1,2,3\n'
                '4,5,6')

        reader = TextReader(StringIO(data), delimiter=',', header=2,
                            as_recarray=True)
        header = reader.header
        expected = [['a', 'b', 'c']]
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

    def test_numpy_string_dtype(self):
        data = """\
a,1
aa,2
aaa,3
aaaa,4
aaaaa,5"""

        def _make_reader(**kwds):
            return TextReader(StringIO(data), delimiter=',', header=None,
                              **kwds)

        reader = _make_reader(dtype='S5,i4')
        result = reader.read()

        self.assert_(result[0].dtype == 'S5')

        ex_values = np.array(['a', 'aa', 'aaa', 'aaaa', 'aaaaa'], dtype='S5')
        self.assert_((result[0] == ex_values).all())
        self.assert_(result[1].dtype == 'i4')

        reader = _make_reader(dtype='S4')
        result = reader.read()
        self.assert_(result[0].dtype == 'S4')
        ex_values = np.array(['a', 'aa', 'aaa', 'aaaa', 'aaaa'], dtype='S4')
        self.assert_((result[0] == ex_values).all())
        self.assert_(result[1].dtype == 'S4')

        reader = _make_reader(dtype='S4', as_recarray=True)
        result = reader.read()
        self.assert_(result['0'].dtype == 'S4')
        ex_values = np.array(['a', 'aa', 'aaa', 'aaaa', 'aaaa'], dtype='S4')
        self.assert_((result['0'] == ex_values).all())
        self.assert_(result['1'].dtype == 'S4')

    def test_pass_dtype(self):
        data = """\
one,two
1,a
2,b
3,c
4,d"""

        def _make_reader(**kwds):
            return TextReader(StringIO(data), delimiter=',', **kwds)

        reader = _make_reader(dtype={'one': 'u1', 1: 'S1'})
        result = reader.read()
        self.assert_(result[0].dtype == 'u1')
        self.assert_(result[1].dtype == 'S1')

        reader = _make_reader(dtype={'one': np.uint8, 1: object})
        result = reader.read()
        self.assert_(result[0].dtype == 'u1')
        self.assert_(result[1].dtype == 'O')

        reader = _make_reader(dtype={'one': np.dtype('u1'),
                                     1: np.dtype('O')})
        result = reader.read()
        self.assert_(result[0].dtype == 'u1')
        self.assert_(result[1].dtype == 'O')

    def test_usecols(self):
        data = """\
a,b,c
1,2,3
4,5,6
7,8,9
10,11,12"""

        def _make_reader(**kwds):
            return TextReader(StringIO(data), delimiter=',', **kwds)

        reader = _make_reader(usecols=(1, 2))
        result = reader.read()

        exp = _make_reader().read()
        self.assertEquals(len(result), 2)
        self.assertTrue((result[1] == exp[1]).all())
        self.assertTrue((result[2] == exp[2]).all())

    def test_cr_delimited(self):
        def _test(text, **kwargs):
            nice_text = text.replace('\r', '\r\n')
            result = TextReader(StringIO(text), **kwargs).read()
            expected = TextReader(StringIO(nice_text), **kwargs).read()
            assert_array_dicts_equal(result, expected)

        data = 'a,b,c\r1,2,3\r4,5,6\r7,8,9\r10,11,12'
        _test(data, delimiter=',')

        data = 'a  b  c\r1  2  3\r4  5  6\r7  8  9\r10  11  12'
        _test(data, delim_whitespace=True)

        data = 'a,b,c\r1,2,3\r4,5,6\r,88,9\r10,11,12'
        _test(data, delimiter=',')

        sample = ('A,B,C,D,E,F,G,H,I,J,K,L,M,N,O\r'
                  'AAAAA,BBBBB,0,0,0,0,0,0,0,0,0,0,0,0,0\r'
                  ',BBBBB,0,0,0,0,0,0,0,0,0,0,0,0,0')
        _test(sample, delimiter=',')

        data = 'A  B  C\r  2  3\r4  5  6'
        _test(data, delim_whitespace=True)

    def test_empty_field_eof(self):
        data = 'a,b,c\n1,2,3\n4,,'

        result = TextReader(StringIO(data), delimiter=',').read()

        expected = {0: np.array([1, 4]),
                    1: np.array(['2', ''], dtype=object),
                    2: np.array(['3', ''], dtype=object)}
        assert_array_dicts_equal(result, expected)


def assert_array_dicts_equal(left, right):
    for k, v in left.iteritems():
        assert(np.array_equal(v, right[k]))

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
