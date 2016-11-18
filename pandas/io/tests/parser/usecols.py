# -*- coding: utf-8 -*-

"""
Tests the usecols functionality during parsing
for all of the parsers defined in parsers.py
"""

import nose

import numpy as np
import pandas.util.testing as tm

from pandas import DataFrame, Index
from pandas.lib import Timestamp
from pandas.compat import StringIO


class UsecolsTests(object):

    def test_raise_on_mixed_dtype_usecols(self):
        # See gh-12678
        data = """a,b,c
        1000,2000,3000
        4000,5000,6000
        """
        msg = ("The elements of 'usecols' must "
               "either be all strings, all unicode, or all integers")
        usecols = [0, 'b', 2]

        with tm.assertRaisesRegexp(ValueError, msg):
            self.read_csv(StringIO(data), usecols=usecols)

    def test_usecols(self):
        data = """\
a,b,c
1,2,3
4,5,6
7,8,9
10,11,12"""

        result = self.read_csv(StringIO(data), usecols=(1, 2))
        result2 = self.read_csv(StringIO(data), usecols=('b', 'c'))
        exp = self.read_csv(StringIO(data))

        self.assertEqual(len(result.columns), 2)
        self.assertTrue((result['b'] == exp['b']).all())
        self.assertTrue((result['c'] == exp['c']).all())

        tm.assert_frame_equal(result, result2)

        result = self.read_csv(StringIO(data), usecols=[1, 2], header=0,
                               names=['foo', 'bar'])
        expected = self.read_csv(StringIO(data), usecols=[1, 2])
        expected.columns = ['foo', 'bar']
        tm.assert_frame_equal(result, expected)

        data = """\
1,2,3
4,5,6
7,8,9
10,11,12"""
        result = self.read_csv(StringIO(data), names=['b', 'c'],
                               header=None, usecols=[1, 2])

        expected = self.read_csv(StringIO(data), names=['a', 'b', 'c'],
                                 header=None)
        expected = expected[['b', 'c']]
        tm.assert_frame_equal(result, expected)

        result2 = self.read_csv(StringIO(data), names=['a', 'b', 'c'],
                                header=None, usecols=['b', 'c'])
        tm.assert_frame_equal(result2, result)

        # see gh-5766
        result = self.read_csv(StringIO(data), names=['a', 'b'],
                               header=None, usecols=[0, 1])

        expected = self.read_csv(StringIO(data), names=['a', 'b', 'c'],
                                 header=None)
        expected = expected[['a', 'b']]
        tm.assert_frame_equal(result, expected)

        # length conflict, passed names and usecols disagree
        self.assertRaises(ValueError, self.read_csv, StringIO(data),
                          names=['a', 'b'], usecols=[1], header=None)

    def test_usecols_index_col_False(self):
        # see gh-9082
        s = "a,b,c,d\n1,2,3,4\n5,6,7,8"
        s_malformed = "a,b,c,d\n1,2,3,4,\n5,6,7,8,"
        cols = ['a', 'c', 'd']
        expected = DataFrame({'a': [1, 5], 'c': [3, 7], 'd': [4, 8]})
        df = self.read_csv(StringIO(s), usecols=cols, index_col=False)
        tm.assert_frame_equal(expected, df)
        df = self.read_csv(StringIO(s_malformed),
                           usecols=cols, index_col=False)
        tm.assert_frame_equal(expected, df)

    def test_usecols_index_col_conflict(self):
        # see gh-4201: test that index_col as integer reflects usecols
        data = 'a,b,c,d\nA,a,1,one\nB,b,2,two'
        expected = DataFrame({'c': [1, 2]}, index=Index(
            ['a', 'b'], name='b'))

        df = self.read_csv(StringIO(data), usecols=['b', 'c'],
                           index_col=0)
        tm.assert_frame_equal(expected, df)

        df = self.read_csv(StringIO(data), usecols=['b', 'c'],
                           index_col='b')
        tm.assert_frame_equal(expected, df)

        df = self.read_csv(StringIO(data), usecols=[1, 2],
                           index_col='b')
        tm.assert_frame_equal(expected, df)

        df = self.read_csv(StringIO(data), usecols=[1, 2],
                           index_col=0)
        tm.assert_frame_equal(expected, df)

        expected = DataFrame(
            {'b': ['a', 'b'], 'c': [1, 2], 'd': ('one', 'two')})
        expected = expected.set_index(['b', 'c'])
        df = self.read_csv(StringIO(data), usecols=['b', 'c', 'd'],
                           index_col=['b', 'c'])
        tm.assert_frame_equal(expected, df)

    def test_usecols_implicit_index_col(self):
        # see gh-2654
        data = 'a,b,c\n4,apple,bat,5.7\n8,orange,cow,10'

        result = self.read_csv(StringIO(data), usecols=['a', 'b'])
        expected = DataFrame({'a': ['apple', 'orange'],
                              'b': ['bat', 'cow']}, index=[4, 8])

        tm.assert_frame_equal(result, expected)

    def test_usecols_regex_sep(self):
        # see gh-2733
        data = 'a  b  c\n4  apple  bat  5.7\n8  orange  cow  10'

        df = self.read_csv(StringIO(data), sep=r'\s+', usecols=('a', 'b'))

        expected = DataFrame({'a': ['apple', 'orange'],
                              'b': ['bat', 'cow']}, index=[4, 8])
        tm.assert_frame_equal(df, expected)

    def test_usecols_with_whitespace(self):
        data = 'a  b  c\n4  apple  bat  5.7\n8  orange  cow  10'

        result = self.read_csv(StringIO(data), delim_whitespace=True,
                               usecols=('a', 'b'))
        expected = DataFrame({'a': ['apple', 'orange'],
                              'b': ['bat', 'cow']}, index=[4, 8])

        tm.assert_frame_equal(result, expected)

    def test_usecols_with_integer_like_header(self):
        data = """2,0,1
        1000,2000,3000
        4000,5000,6000
        """

        usecols = [0, 1]  # column selection by index
        expected = DataFrame(data=[[1000, 2000],
                                   [4000, 5000]],
                             columns=['2', '0'])
        df = self.read_csv(StringIO(data), usecols=usecols)
        tm.assert_frame_equal(df, expected)

        usecols = ['0', '1']  # column selection by name
        expected = DataFrame(data=[[2000, 3000],
                                   [5000, 6000]],
                             columns=['0', '1'])
        df = self.read_csv(StringIO(data), usecols=usecols)
        tm.assert_frame_equal(df, expected)

    def test_usecols_with_parse_dates(self):
        # See gh-9755
        s = """a,b,c,d,e
        0,1,20140101,0900,4
        0,1,20140102,1000,4"""
        parse_dates = [[1, 2]]

        cols = {
            'a': [0, 0],
            'c_d': [
                Timestamp('2014-01-01 09:00:00'),
                Timestamp('2014-01-02 10:00:00')
            ]
        }
        expected = DataFrame(cols, columns=['c_d', 'a'])

        df = self.read_csv(StringIO(s), usecols=[0, 2, 3],
                           parse_dates=parse_dates)
        tm.assert_frame_equal(df, expected)

        df = self.read_csv(StringIO(s), usecols=[3, 0, 2],
                           parse_dates=parse_dates)
        tm.assert_frame_equal(df, expected)

    def test_usecols_with_parse_dates_and_full_names(self):
        # See gh-9755
        s = """0,1,20140101,0900,4
        0,1,20140102,1000,4"""
        parse_dates = [[1, 2]]
        names = list('abcde')

        cols = {
            'a': [0, 0],
            'c_d': [
                Timestamp('2014-01-01 09:00:00'),
                Timestamp('2014-01-02 10:00:00')
            ]
        }
        expected = DataFrame(cols, columns=['c_d', 'a'])

        df = self.read_csv(StringIO(s), names=names,
                           usecols=[0, 2, 3],
                           parse_dates=parse_dates)
        tm.assert_frame_equal(df, expected)

        df = self.read_csv(StringIO(s), names=names,
                           usecols=[3, 0, 2],
                           parse_dates=parse_dates)
        tm.assert_frame_equal(df, expected)

    def test_usecols_with_parse_dates_and_usecol_names(self):
        # See gh-9755
        s = """0,1,20140101,0900,4
        0,1,20140102,1000,4"""
        parse_dates = [[1, 2]]
        names = list('acd')

        cols = {
            'a': [0, 0],
            'c_d': [
                Timestamp('2014-01-01 09:00:00'),
                Timestamp('2014-01-02 10:00:00')
            ]
        }
        expected = DataFrame(cols, columns=['c_d', 'a'])

        df = self.read_csv(StringIO(s), names=names,
                           usecols=[0, 2, 3],
                           parse_dates=parse_dates)
        tm.assert_frame_equal(df, expected)

        df = self.read_csv(StringIO(s), names=names,
                           usecols=[3, 0, 2],
                           parse_dates=parse_dates)
        tm.assert_frame_equal(df, expected)

    def test_usecols_with_unicode_strings(self):
        # see gh-13219

        s = '''AAA,BBB,CCC,DDD
        0.056674973,8,True,a
        2.613230982,2,False,b
        3.568935038,7,False,a
        '''

        data = {
            'AAA': {
                0: 0.056674972999999997,
                1: 2.6132309819999997,
                2: 3.5689350380000002
            },
            'BBB': {0: 8, 1: 2, 2: 7}
        }
        expected = DataFrame(data)

        df = self.read_csv(StringIO(s), usecols=[u'AAA', u'BBB'])
        tm.assert_frame_equal(df, expected)

    def test_usecols_with_single_byte_unicode_strings(self):
        # see gh-13219

        s = '''A,B,C,D
        0.056674973,8,True,a
        2.613230982,2,False,b
        3.568935038,7,False,a
        '''

        data = {
            'A': {
                0: 0.056674972999999997,
                1: 2.6132309819999997,
                2: 3.5689350380000002
            },
            'B': {0: 8, 1: 2, 2: 7}
        }
        expected = DataFrame(data)

        df = self.read_csv(StringIO(s), usecols=[u'A', u'B'])
        tm.assert_frame_equal(df, expected)

    def test_usecols_with_mixed_encoding_strings(self):
        s = '''AAA,BBB,CCC,DDD
        0.056674973,8,True,a
        2.613230982,2,False,b
        3.568935038,7,False,a
        '''

        msg = ("The elements of 'usecols' must "
               "either be all strings, all unicode, or all integers")

        with tm.assertRaisesRegexp(ValueError, msg):
            self.read_csv(StringIO(s), usecols=[u'AAA', b'BBB'])

        with tm.assertRaisesRegexp(ValueError, msg):
            self.read_csv(StringIO(s), usecols=[b'AAA', u'BBB'])

    def test_usecols_with_multibyte_characters(self):
        s = '''あああ,いい,ううう,ええええ
        0.056674973,8,True,a
        2.613230982,2,False,b
        3.568935038,7,False,a
        '''
        data = {
            'あああ': {
                0: 0.056674972999999997,
                1: 2.6132309819999997,
                2: 3.5689350380000002
            },
            'いい': {0: 8, 1: 2, 2: 7}
        }
        expected = DataFrame(data)

        df = self.read_csv(StringIO(s), usecols=['あああ', 'いい'])
        tm.assert_frame_equal(df, expected)

    def test_usecols_with_multibyte_unicode_characters(self):
        raise nose.SkipTest('TODO: see gh-13253')

        s = '''あああ,いい,ううう,ええええ
        0.056674973,8,True,a
        2.613230982,2,False,b
        3.568935038,7,False,a
        '''
        data = {
            'あああ': {
                0: 0.056674972999999997,
                1: 2.6132309819999997,
                2: 3.5689350380000002
            },
            'いい': {0: 8, 1: 2, 2: 7}
        }
        expected = DataFrame(data)

        df = self.read_csv(StringIO(s), usecols=[u'あああ', u'いい'])
        tm.assert_frame_equal(df, expected)

    def test_empty_usecols(self):
        # should not raise
        data = 'a,b,c\n1,2,3\n4,5,6'
        expected = DataFrame()
        result = self.read_csv(StringIO(data), usecols=set([]))
        tm.assert_frame_equal(result, expected)

    def test_np_array_usecols(self):
        # See gh-12546
        data = 'a,b,c\n1,2,3'
        usecols = np.array(['a', 'b'])

        expected = DataFrame([[1, 2]], columns=usecols)
        result = self.read_csv(StringIO(data), usecols=usecols)
        tm.assert_frame_equal(result, expected)
