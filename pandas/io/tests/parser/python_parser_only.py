# -*- coding: utf-8 -*-

"""
Tests that apply specifically to the Python parser. Unless specifically
stated as a Python-specific issue, the goal is to eventually move as many of
these tests out of this module as soon as the C parser can accept further
arguments when parsing.
"""

import sys
import nose

import pandas.util.testing as tm
from pandas import DataFrame, Index
from pandas import compat
from pandas.compat import StringIO, BytesIO, u


class PythonParserTests(object):
    def test_negative_skipfooter_raises(self):
        text = """#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
1/1/2000,1.,2.,3.
1/2/2000,4,5,6
1/3/2000,7,8,9
"""

        with tm.assertRaisesRegexp(
                ValueError, 'skip footer cannot be negative'):
            self.read_csv(StringIO(text), skipfooter=-1)

    def test_sniff_delimiter(self):
        text = """index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
"""
        data = self.read_csv(StringIO(text), index_col=0, sep=None)
        self.assertTrue(data.index.equals(Index(['foo', 'bar', 'baz'])))

        data2 = self.read_csv(StringIO(text), index_col=0, delimiter='|')
        tm.assert_frame_equal(data, data2)

        text = """ignore this
ignore this too
index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
"""
        data3 = self.read_csv(StringIO(text), index_col=0,
                              sep=None, skiprows=2)
        tm.assert_frame_equal(data, data3)

        text = u("""ignore this
ignore this too
index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
""").encode('utf-8')

        s = BytesIO(text)
        if compat.PY3:
            # somewhat False since the code never sees bytes
            from io import TextIOWrapper
            s = TextIOWrapper(s, encoding='utf-8')

        data4 = self.read_csv(s, index_col=0, sep=None, skiprows=2,
                              encoding='utf-8')
        tm.assert_frame_equal(data, data4)

    def test_BytesIO_input(self):
        if not compat.PY3:
            raise nose.SkipTest(
                "Bytes-related test - only needs to work on Python 3")

        data = BytesIO("שלום::1234\n562::123".encode('cp1255'))
        result = self.read_table(data, sep="::", encoding='cp1255')
        expected = DataFrame([[562, 123]], columns=["שלום", "1234"])
        tm.assert_frame_equal(result, expected)

    def test_single_line(self):
        # see gh-6607: sniff separator

        buf = StringIO()
        sys.stdout = buf

        try:
            df = self.read_csv(StringIO('1,2'), names=['a', 'b'],
                               header=None, sep=None)
            tm.assert_frame_equal(DataFrame({'a': [1], 'b': [2]}), df)
        finally:
            sys.stdout = sys.__stdout__

    def test_skip_footer(self):
        # see gh-6607
        data = """A,B,C
1,2,3
4,5,6
7,8,9
want to skip this
also also skip this
"""
        result = self.read_csv(StringIO(data), skip_footer=2)
        no_footer = '\n'.join(data.split('\n')[:-3])
        expected = self.read_csv(StringIO(no_footer))
        tm.assert_frame_equal(result, expected)

        result = self.read_csv(StringIO(data), nrows=3)
        tm.assert_frame_equal(result, expected)

        # skipfooter alias
        result = self.read_csv(StringIO(data), skipfooter=2)
        no_footer = '\n'.join(data.split('\n')[:-3])
        expected = self.read_csv(StringIO(no_footer))
        tm.assert_frame_equal(result, expected)

    def test_decompression_regex_sep(self):
        # see gh-6607

        try:
            import gzip
            import bz2
        except ImportError:
            raise nose.SkipTest('need gzip and bz2 to run')

        data = open(self.csv1, 'rb').read()
        data = data.replace(b',', b'::')
        expected = self.read_csv(self.csv1)

        with tm.ensure_clean() as path:
            tmp = gzip.GzipFile(path, mode='wb')
            tmp.write(data)
            tmp.close()

            result = self.read_csv(path, sep='::', compression='gzip')
            tm.assert_frame_equal(result, expected)

        with tm.ensure_clean() as path:
            tmp = bz2.BZ2File(path, mode='wb')
            tmp.write(data)
            tmp.close()

            result = self.read_csv(path, sep='::', compression='bz2')
            tm.assert_frame_equal(result, expected)

            self.assertRaises(ValueError, self.read_csv,
                              path, compression='bz3')

    def test_read_table_buglet_4x_multiindex(self):
        # see gh-6607
        text = """                      A       B       C       D        E
one two three   four
a   b   10.0032 5    -0.5109 -2.3358 -0.4645  0.05076  0.3640
a   q   20      4     0.4473  1.4152  0.2834  1.00661  0.1744
x   q   30      3    -0.6662 -0.5243 -0.3580  0.89145  2.5838"""

        df = self.read_table(StringIO(text), sep='\s+')
        self.assertEqual(df.index.names, ('one', 'two', 'three', 'four'))

        # see gh-6893
        data = '      A B C\na b c\n1 3 7 0 3 6\n3 1 4 1 5 9'
        expected = DataFrame.from_records(
            [(1, 3, 7, 0, 3, 6), (3, 1, 4, 1, 5, 9)],
            columns=list('abcABC'), index=list('abc'))
        actual = self.read_table(StringIO(data), sep='\s+')
        tm.assert_frame_equal(actual, expected)
