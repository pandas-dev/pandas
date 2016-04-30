# -*- coding: utf-8 -*-

"""
Tests that apply specifically to the CParser. Unless specifically stated
as a CParser-specific issue, the goal is to eventually move as many of
these tests out of this module as soon as the Python parser can accept
further arguments when parsing.
"""

import nose
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import DataFrame, Series, Index, MultiIndex
from pandas import compat
from pandas.compat import StringIO, range, lrange


class CParserTests(object):
    def test_buffer_overflow(self):
        # see gh-9205: test certain malformed input files that cause
        # buffer overflows in tokenizer.c

        malfw = "1\r1\r1\r 1\r 1\r"         # buffer overflow in words pointer
        malfs = "1\r1\r1\r 1\r 1\r11\r"     # buffer overflow in stream pointer
        malfl = "1\r1\r1\r 1\r 1\r11\r1\r"  # buffer overflow in lines pointer

        cperr = 'Buffer overflow caught - possible malformed input file.'

        for malf in (malfw, malfs, malfl):
            try:
                self.read_table(StringIO(malf))
            except Exception as err:
                self.assertIn(cperr, str(err))

    def test_buffer_rd_bytes(self):
        # see gh-12098: src->buffer in the C parser can be freed twice leading
        # to a segfault if a corrupt gzip file is read with 'read_csv' and the
        # buffer is filled more than once before gzip throws an exception

        data = '\x1F\x8B\x08\x00\x00\x00\x00\x00\x00\x03\xED\xC3\x41\x09' \
               '\x00\x00\x08\x00\xB1\xB7\xB6\xBA\xFE\xA5\xCC\x21\x6C\xB0' \
               '\xA6\x4D' + '\x55' * 267 + \
               '\x7D\xF7\x00\x91\xE0\x47\x97\x14\x38\x04\x00' \
               '\x1f\x8b\x08\x00VT\x97V\x00\x03\xed]\xefO'
        for i in range(100):
            try:
                self.read_csv(StringIO(data),
                              compression='gzip',
                              delim_whitespace=True)
            except Exception:
                pass

    def test_delim_whitespace_custom_terminator(self):
        # See gh-12912
        data = """a b c~1 2 3~4 5 6~7 8 9"""
        df = self.read_csv(StringIO(data), lineterminator='~',
                           delim_whitespace=True)
        expected = DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                             columns=['a', 'b', 'c'])
        tm.assert_frame_equal(df, expected)

    def test_parse_dates_empty_string(self):
        # see gh-2263
        s = StringIO("Date, test\n2012-01-01, 1\n,2")
        result = self.read_csv(s, parse_dates=["Date"], na_filter=False)
        self.assertTrue(result['Date'].isnull()[1])

    def test_dtype_and_names_error(self):
        # see gh-8833: passing both dtype and names
        # resulting in an error reporting issue
        data = """
1.0 1
2.0 2
3.0 3
"""
        # base cases
        result = self.read_csv(StringIO(data), sep='\s+', header=None)
        expected = DataFrame([[1.0, 1], [2.0, 2], [3.0, 3]])
        tm.assert_frame_equal(result, expected)

        result = self.read_csv(StringIO(data), sep='\s+',
                               header=None, names=['a', 'b'])
        expected = DataFrame(
            [[1.0, 1], [2.0, 2], [3.0, 3]], columns=['a', 'b'])
        tm.assert_frame_equal(result, expected)

        # fallback casting
        result = self.read_csv(StringIO(
            data), sep='\s+', header=None,
            names=['a', 'b'], dtype={'a': np.int32})
        expected = DataFrame([[1, 1], [2, 2], [3, 3]],
                             columns=['a', 'b'])
        expected['a'] = expected['a'].astype(np.int32)
        tm.assert_frame_equal(result, expected)

        data = """
1.0 1
nan 2
3.0 3
"""
        # fallback casting, but not castable
        with tm.assertRaisesRegexp(ValueError, 'cannot safely convert'):
            self.read_csv(StringIO(data), sep='\s+', header=None,
                          names=['a', 'b'], dtype={'a': np.int32})

    def test_passing_dtype(self):
        # see gh-6607
        df = DataFrame(np.random.rand(5, 2), columns=list(
            'AB'), index=['1A', '1B', '1C', '1D', '1E'])

        with tm.ensure_clean('__passing_str_as_dtype__.csv') as path:
            df.to_csv(path)

            # see gh-3795: passing 'str' as the dtype
            result = self.read_csv(path, dtype=str, index_col=0)
            tm.assert_series_equal(result.dtypes, Series(
                {'A': 'object', 'B': 'object'}))

            # we expect all object columns, so need to
            # convert to test for equivalence
            result = result.astype(float)
            tm.assert_frame_equal(result, df)

            # invalid dtype
            self.assertRaises(TypeError, self.read_csv, path,
                              dtype={'A': 'foo', 'B': 'float64'},
                              index_col=0)

            # valid but we don't support it (date)
            self.assertRaises(TypeError, self.read_csv, path,
                              dtype={'A': 'datetime64', 'B': 'float64'},
                              index_col=0)
            self.assertRaises(TypeError, self.read_csv, path,
                              dtype={'A': 'datetime64', 'B': 'float64'},
                              index_col=0, parse_dates=['B'])

            # valid but we don't support it
            self.assertRaises(TypeError, self.read_csv, path,
                              dtype={'A': 'timedelta64', 'B': 'float64'},
                              index_col=0)

        # see gh-12048: empty frame
        actual = self.read_csv(StringIO('A,B'), dtype=str)
        expected = DataFrame({'A': [], 'B': []}, index=[], dtype=str)
        tm.assert_frame_equal(actual, expected)

    def test_precise_conversion(self):
        # see gh-8002
        tm._skip_if_32bit()
        from decimal import Decimal

        normal_errors = []
        precise_errors = []

        # test numbers between 1 and 2
        for num in np.linspace(1., 2., num=500):
            # 25 decimal digits of precision
            text = 'a\n{0:.25}'.format(num)

            normal_val = float(self.read_csv(StringIO(text))['a'][0])
            precise_val = float(self.read_csv(
                StringIO(text), float_precision='high')['a'][0])
            roundtrip_val = float(self.read_csv(
                StringIO(text), float_precision='round_trip')['a'][0])
            actual_val = Decimal(text[2:])

            def error(val):
                return abs(Decimal('{0:.100}'.format(val)) - actual_val)

            normal_errors.append(error(normal_val))
            precise_errors.append(error(precise_val))

            # round-trip should match float()
            self.assertEqual(roundtrip_val, float(text[2:]))

        self.assertTrue(sum(precise_errors) <= sum(normal_errors))
        self.assertTrue(max(precise_errors) <= max(normal_errors))

    def test_compact_ints(self):
        if compat.is_platform_windows() and not self.low_memory:
            raise nose.SkipTest(
                "segfaults on win-64, only when all tests are run")

        data = ('0,1,0,0\n'
                '1,1,0,0\n'
                '0,1,0,1')

        result = self.read_csv(StringIO(data), delimiter=',', header=None,
                               compact_ints=True, as_recarray=True)
        ex_dtype = np.dtype([(str(i), 'i1') for i in range(4)])
        self.assertEqual(result.dtype, ex_dtype)

        result = self.read_csv(StringIO(data), delimiter=',', header=None,
                               as_recarray=True, compact_ints=True,
                               use_unsigned=True)
        ex_dtype = np.dtype([(str(i), 'u1') for i in range(4)])
        self.assertEqual(result.dtype, ex_dtype)

    def test_compact_ints_as_recarray(self):
        if compat.is_platform_windows() and self.low_memory:
            raise nose.SkipTest(
                "segfaults on win-64, only when all tests are run")

        data = ('0,1,0,0\n'
                '1,1,0,0\n'
                '0,1,0,1')

        result = self.read_csv(StringIO(data), delimiter=',', header=None,
                               compact_ints=True, as_recarray=True)
        ex_dtype = np.dtype([(str(i), 'i1') for i in range(4)])
        self.assertEqual(result.dtype, ex_dtype)

        result = self.read_csv(StringIO(data), delimiter=',', header=None,
                               as_recarray=True, compact_ints=True,
                               use_unsigned=True)
        ex_dtype = np.dtype([(str(i), 'u1') for i in range(4)])
        self.assertEqual(result.dtype, ex_dtype)

    def test_pass_dtype(self):
        data = """\
one,two
1,2.5
2,3.5
3,4.5
4,5.5"""

        result = self.read_csv(StringIO(data), dtype={'one': 'u1', 1: 'S1'})
        self.assertEqual(result['one'].dtype, 'u1')
        self.assertEqual(result['two'].dtype, 'object')

    def test_pass_dtype_as_recarray(self):
        if compat.is_platform_windows() and self.low_memory:
            raise nose.SkipTest(
                "segfaults on win-64, only when all tests are run")

        data = """\
one,two
1,2.5
2,3.5
3,4.5
4,5.5"""

        result = self.read_csv(StringIO(data), dtype={'one': 'u1', 1: 'S1'},
                               as_recarray=True)
        self.assertEqual(result['one'].dtype, 'u1')
        self.assertEqual(result['two'].dtype, 'S1')

    def test_empty_pass_dtype(self):
        data = 'one,two'
        result = self.read_csv(StringIO(data), dtype={'one': 'u1'})

        expected = DataFrame({'one': np.empty(0, dtype='u1'),
                              'two': np.empty(0, dtype=np.object)})
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_empty_with_index_pass_dtype(self):
        data = 'one,two'
        result = self.read_csv(StringIO(data), index_col=['one'],
                               dtype={'one': 'u1', 1: 'f'})

        expected = DataFrame({'two': np.empty(0, dtype='f')},
                             index=Index([], dtype='u1', name='one'))
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_empty_with_multiindex_pass_dtype(self):
        data = 'one,two,three'
        result = self.read_csv(StringIO(data), index_col=['one', 'two'],
                               dtype={'one': 'u1', 1: 'f8'})

        exp_idx = MultiIndex.from_arrays([np.empty(0, dtype='u1'),
                                          np.empty(0, dtype='O')],
                                         names=['one', 'two'])
        expected = DataFrame(
            {'three': np.empty(0, dtype=np.object)}, index=exp_idx)
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_empty_with_mangled_column_pass_dtype_by_names(self):
        data = 'one,one'
        result = self.read_csv(StringIO(data), dtype={
            'one': 'u1', 'one.1': 'f'})

        expected = DataFrame(
            {'one': np.empty(0, dtype='u1'), 'one.1': np.empty(0, dtype='f')})
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_empty_with_mangled_column_pass_dtype_by_indexes(self):
        data = 'one,one'
        result = self.read_csv(StringIO(data), dtype={0: 'u1', 1: 'f'})

        expected = DataFrame(
            {'one': np.empty(0, dtype='u1'), 'one.1': np.empty(0, dtype='f')})
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_empty_with_dup_column_pass_dtype_by_names(self):
        data = 'one,one'
        result = self.read_csv(
            StringIO(data), mangle_dupe_cols=False, dtype={'one': 'u1'})
        expected = pd.concat([Series([], name='one', dtype='u1')] * 2, axis=1)
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_empty_with_dup_column_pass_dtype_by_indexes(self):
        # FIXME in gh-9424
        raise nose.SkipTest(
            "gh-9424; known failure read_csv with duplicate columns")

        data = 'one,one'
        result = self.read_csv(
            StringIO(data), mangle_dupe_cols=False, dtype={0: 'u1', 1: 'f'})
        expected = pd.concat([Series([], name='one', dtype='u1'),
                              Series([], name='one', dtype='f')], axis=1)
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_usecols_dtypes(self):
        data = """\
1,2,3
4,5,6
7,8,9
10,11,12"""

        result = self.read_csv(StringIO(data), usecols=(0, 1, 2),
                               names=('a', 'b', 'c'),
                               header=None,
                               converters={'a': str},
                               dtype={'b': int, 'c': float},
                               )
        result2 = self.read_csv(StringIO(data), usecols=(0, 2),
                                names=('a', 'b', 'c'),
                                header=None,
                                converters={'a': str},
                                dtype={'b': int, 'c': float},
                                )
        self.assertTrue((result.dtypes == [object, np.int, np.float]).all())
        self.assertTrue((result2.dtypes == [object, np.float]).all())

    def test_memory_map(self):
        # it works!
        self.read_csv(self.csv1, memory_map=True)

    def test_disable_bool_parsing(self):
        # #2090

        data = """A,B,C
Yes,No,Yes
No,Yes,Yes
Yes,,Yes
No,No,No"""

        result = self.read_csv(StringIO(data), dtype=object)
        self.assertTrue((result.dtypes == object).all())

        result = self.read_csv(StringIO(data), dtype=object, na_filter=False)
        self.assertEqual(result['B'][2], '')

    def test_euro_decimal_format(self):
        data = """Id;Number1;Number2;Text1;Text2;Number3
1;1521,1541;187101,9543;ABC;poi;4,738797819
2;121,12;14897,76;DEF;uyt;0,377320872
3;878,158;108013,434;GHI;rez;2,735694704"""

        df2 = self.read_csv(StringIO(data), sep=';', decimal=',')
        self.assertEqual(df2['Number1'].dtype, float)
        self.assertEqual(df2['Number2'].dtype, float)
        self.assertEqual(df2['Number3'].dtype, float)

    def test_custom_lineterminator(self):
        data = 'a,b,c~1,2,3~4,5,6'

        result = self.read_csv(StringIO(data), lineterminator='~')
        expected = self.read_csv(StringIO(data.replace('~', '\n')))

        tm.assert_frame_equal(result, expected)

    def test_raise_on_passed_int_dtype_with_nas(self):
        # see gh-2631
        data = """YEAR, DOY, a
2001,106380451,10
2001,,11
2001,106380451,67"""
        self.assertRaises(ValueError, self.read_csv, StringIO(data),
                          sep=",", skipinitialspace=True,
                          dtype={'DOY': np.int64})

    def test_na_trailing_columns(self):
        data = """Date,Currenncy,Symbol,Type,Units,UnitPrice,Cost,Tax
2012-03-14,USD,AAPL,BUY,1000
2012-05-12,USD,SBUX,SELL,500"""

        result = self.read_csv(StringIO(data))
        self.assertEqual(result['Date'][1], '2012-05-12')
        self.assertTrue(result['UnitPrice'].isnull().all())

    def test_parse_ragged_csv(self):
        data = """1,2,3
1,2,3,4
1,2,3,4,5
1,2
1,2,3,4"""

        nice_data = """1,2,3,,
1,2,3,4,
1,2,3,4,5
1,2,,,
1,2,3,4,"""
        result = self.read_csv(StringIO(data), header=None,
                               names=['a', 'b', 'c', 'd', 'e'])

        expected = self.read_csv(StringIO(nice_data), header=None,
                                 names=['a', 'b', 'c', 'd', 'e'])

        tm.assert_frame_equal(result, expected)

        # too many columns, cause segfault if not careful
        data = "1,2\n3,4,5"

        result = self.read_csv(StringIO(data), header=None,
                               names=lrange(50))
        expected = self.read_csv(StringIO(data), header=None,
                                 names=lrange(3)).reindex(columns=lrange(50))

        tm.assert_frame_equal(result, expected)

    def test_tokenize_CR_with_quoting(self):
        # see gh-3453

        data = ' a,b,c\r"a,b","e,d","f,f"'

        result = self.read_csv(StringIO(data), header=None)
        expected = self.read_csv(StringIO(data.replace('\r', '\n')),
                                 header=None)
        tm.assert_frame_equal(result, expected)

        result = self.read_csv(StringIO(data))
        expected = self.read_csv(StringIO(data.replace('\r', '\n')))
        tm.assert_frame_equal(result, expected)

    def test_raise_on_no_columns(self):
        # single newline
        data = "\n"
        self.assertRaises(ValueError, self.read_csv, StringIO(data))

        # test with more than a single newline
        data = "\n\n\n"
        self.assertRaises(ValueError, self.read_csv, StringIO(data))

    def test_1000_sep_with_decimal(self):
        data = """A|B|C
1|2,334.01|5
10|13|10.
"""
        expected = DataFrame({
            'A': [1, 10],
            'B': [2334.01, 13],
            'C': [5, 10.]
        })

        tm.assert_equal(expected.A.dtype, 'int64')
        tm.assert_equal(expected.B.dtype, 'float')
        tm.assert_equal(expected.C.dtype, 'float')

        df = self.read_csv(StringIO(data), sep='|', thousands=',', decimal='.')
        tm.assert_frame_equal(df, expected)

        df = self.read_table(StringIO(data), sep='|',
                             thousands=',', decimal='.')
        tm.assert_frame_equal(df, expected)

        data_with_odd_sep = """A|B|C
1|2.334,01|5
10|13|10,
"""
        df = self.read_csv(StringIO(data_with_odd_sep),
                           sep='|', thousands='.', decimal=',')
        tm.assert_frame_equal(df, expected)

        df = self.read_table(StringIO(data_with_odd_sep),
                             sep='|', thousands='.', decimal=',')
        tm.assert_frame_equal(df, expected)

    def test_grow_boundary_at_cap(self):
        # See gh-12494
        #
        # Cause of error was that the C parser
        # was not increasing the buffer size when
        # the desired space would fill the buffer
        # to capacity, which would later cause a
        # buffer overflow error when checking the
        # EOF terminator of the CSV stream
        def test_empty_header_read(count):
            s = StringIO(',' * count)
            expected = DataFrame(columns=[
                'Unnamed: {i}'.format(i=i)
                for i in range(count + 1)])
            df = self.read_csv(s)
            tm.assert_frame_equal(df, expected)

        for count in range(1, 101):
            test_empty_header_read(count)

    def test_inf_parsing(self):
        data = """\
,A
a,inf
b,-inf
c,Inf
d,-Inf
e,INF
f,-INF
g,INf
h,-INf
i,inF
j,-inF"""
        inf = float('inf')
        expected = Series([inf, -inf] * 5)

        df = self.read_csv(StringIO(data), index_col=0)
        tm.assert_almost_equal(df['A'].values, expected.values)

        df = self.read_csv(StringIO(data), index_col=0, na_filter=False)
        tm.assert_almost_equal(df['A'].values, expected.values)
