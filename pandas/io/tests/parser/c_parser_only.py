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
from pandas import DataFrame, Series, Index, MultiIndex, Categorical
from pandas import compat
from pandas.compat import StringIO, range, lrange
from pandas.types.dtypes import CategoricalDtype


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

    def test_dtype_and_names_error(self):
        # see gh-8833: passing both dtype and names
        # resulting in an error reporting issue
        data = """
1.0 1
2.0 2
3.0 3
"""
        # base cases
        result = self.read_csv(StringIO(data), sep=r'\s+', header=None)
        expected = DataFrame([[1.0, 1], [2.0, 2], [3.0, 3]])
        tm.assert_frame_equal(result, expected)

        result = self.read_csv(StringIO(data), sep=r'\s+',
                               header=None, names=['a', 'b'])
        expected = DataFrame(
            [[1.0, 1], [2.0, 2], [3.0, 3]], columns=['a', 'b'])
        tm.assert_frame_equal(result, expected)

        # fallback casting
        result = self.read_csv(StringIO(
            data), sep=r'\s+', header=None,
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
            self.read_csv(StringIO(data), sep=r'\s+', header=None,
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

            # valid but unsupported - fixed width unicode string
            self.assertRaises(TypeError, self.read_csv, path,
                              dtype={'A': 'U8'},
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

    def test_categorical_dtype(self):
        # GH 10153
        data = """a,b,c
1,a,3.4
1,a,3.4
2,b,4.5"""
        expected = pd.DataFrame({'a': Categorical(['1', '1', '2']),
                                 'b': Categorical(['a', 'a', 'b']),
                                 'c': Categorical(['3.4', '3.4', '4.5'])})
        actual = self.read_csv(StringIO(data), dtype='category')
        tm.assert_frame_equal(actual, expected)

        actual = self.read_csv(StringIO(data), dtype=CategoricalDtype())
        tm.assert_frame_equal(actual, expected)

        actual = self.read_csv(StringIO(data), dtype={'a': 'category',
                                                      'b': 'category',
                                                      'c': CategoricalDtype()})
        tm.assert_frame_equal(actual, expected)

        actual = self.read_csv(StringIO(data), dtype={'b': 'category'})
        expected = pd.DataFrame({'a': [1, 1, 2],
                                 'b': Categorical(['a', 'a', 'b']),
                                 'c': [3.4, 3.4, 4.5]})
        tm.assert_frame_equal(actual, expected)

        actual = self.read_csv(StringIO(data), dtype={1: 'category'})
        tm.assert_frame_equal(actual, expected)

        # unsorted
        data = """a,b,c
1,b,3.4
1,b,3.4
2,a,4.5"""
        expected = pd.DataFrame({'a': Categorical(['1', '1', '2']),
                                 'b': Categorical(['b', 'b', 'a']),
                                 'c': Categorical(['3.4', '3.4', '4.5'])})
        actual = self.read_csv(StringIO(data), dtype='category')
        tm.assert_frame_equal(actual, expected)

        # missing
        data = """a,b,c
1,b,3.4
1,nan,3.4
2,a,4.5"""
        expected = pd.DataFrame({'a': Categorical(['1', '1', '2']),
                                 'b': Categorical(['b', np.nan, 'a']),
                                 'c': Categorical(['3.4', '3.4', '4.5'])})
        actual = self.read_csv(StringIO(data), dtype='category')
        tm.assert_frame_equal(actual, expected)

    def test_categorical_dtype_encoding(self):
        # GH 10153
        pth = tm.get_data_path('unicode_series.csv')
        encoding = 'latin-1'
        expected = self.read_csv(pth, header=None, encoding=encoding)
        expected[1] = Categorical(expected[1])
        actual = self.read_csv(pth, header=None, encoding=encoding,
                               dtype={1: 'category'})
        tm.assert_frame_equal(actual, expected)

        pth = tm.get_data_path('utf16_ex.txt')
        encoding = 'utf-16'
        expected = self.read_table(pth, encoding=encoding)
        expected = expected.apply(Categorical)
        actual = self.read_table(pth, encoding=encoding, dtype='category')
        tm.assert_frame_equal(actual, expected)

    def test_categorical_dtype_chunksize(self):
        # GH 10153
        data = """a,b
1,a
1,b
1,b
2,c"""
        expecteds = [pd.DataFrame({'a': [1, 1],
                                   'b': Categorical(['a', 'b'])}),
                     pd.DataFrame({'a': [1, 2],
                                   'b': Categorical(['b', 'c'])},
                                  index=[2, 3])]
        actuals = self.read_csv(StringIO(data), dtype={'b': 'category'},
                                chunksize=2)

        for actual, expected in zip(actuals, expecteds):
            tm.assert_frame_equal(actual, expected)

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

        with tm.assert_produces_warning(
                FutureWarning, check_stacklevel=False):
            result = self.read_csv(StringIO(data), dtype={
                'one': 'u1', 1: 'S1'}, as_recarray=True)
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

    def test_empty_with_dup_column_pass_dtype_by_indexes(self):
        # see gh-9424
        expected = pd.concat([Series([], name='one', dtype='u1'),
                              Series([], name='one.1', dtype='f')], axis=1)

        data = 'one,one'
        result = self.read_csv(StringIO(data), dtype={0: 'u1', 1: 'f'})
        tm.assert_frame_equal(result, expected, check_index_type=False)

        data = ''
        result = self.read_csv(StringIO(data), names=['one', 'one'],
                               dtype={0: 'u1', 1: 'f'})
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

    def test_parse_trim_buffers(self):
        # This test is part of a bugfix for issue #13703. It attmepts to
        # to stress the system memory allocator, to cause it to move the
        # stream buffer and either let the OS reclaim the region, or let
        # other memory requests of parser otherwise modify the contents
        # of memory space, where it was formely located.
        # This test is designed to cause a `segfault` with unpatched
        # `tokenizer.c`. Sometimes the test fails on `segfault`, other
        # times it fails due to memory corruption, which causes the
        # loaded DataFrame to differ from the expected one.

        # Generate a large mixed-type CSV file on-the-fly (one record is
        # approx 1.5KiB).
        record_ = \
            """9999-9,99:99,,,,ZZ,ZZ,,,ZZZ-ZZZZ,.Z-ZZZZ,-9.99,,,9.99,Z""" \
            """ZZZZ,,-99,9,ZZZ-ZZZZ,ZZ-ZZZZ,,9.99,ZZZ-ZZZZZ,ZZZ-ZZZZZ,""" \
            """ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,9""" \
            """99,ZZZ-ZZZZ,,ZZ-ZZZZ,,,,,ZZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZ,,,9,9,""" \
            """9,9,99,99,999,999,ZZZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZ,9,ZZ-ZZZZ,9.""" \
            """99,ZZ-ZZZZ,ZZ-ZZZZ,,,,ZZZZ,,,ZZ,ZZ,,,,,,,,,,,,,9,,,999.""" \
            """99,999.99,,,ZZZZZ,,,Z9,,,,,,,ZZZ,ZZZ,,,,,,,,,,,ZZZZZ,ZZ""" \
            """ZZZ,ZZZ-ZZZZZZ,ZZZ-ZZZZZZ,ZZ-ZZZZ,ZZ-ZZZZ,ZZ-ZZZZ,ZZ-ZZ""" \
            """ZZ,,,999999,999999,ZZZ,ZZZ,,,ZZZ,ZZZ,999.99,999.99,,,,Z""" \
            """ZZ-ZZZ,ZZZ-ZZZ,-9.99,-9.99,9,9,,99,,9.99,9.99,9,9,9.99,""" \
            """9.99,,,,9.99,9.99,,99,,99,9.99,9.99,,,ZZZ,ZZZ,,999.99,,""" \
            """999.99,ZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,,,ZZZZZ,ZZZZZ,ZZZ,ZZZ,9,9,""" \
            """,,,,,ZZZ-ZZZZ,ZZZ999Z,,,999.99,,999.99,ZZZ-ZZZZ,,,9.999""" \
            """,9.999,9.999,9.999,-9.999,-9.999,-9.999,-9.999,9.999,9.""" \
            """999,9.999,9.999,9.999,9.999,9.999,9.999,99999,ZZZ-ZZZZ,""" \
            """,9.99,ZZZ,,,,,,,,ZZZ,,,,,9,,,,9,,,,,,,,,,ZZZ-ZZZZ,ZZZ-Z""" \
            """ZZZ,,ZZZZZ,ZZZZZ,ZZZZZ,ZZZZZ,,,9.99,,ZZ-ZZZZ,ZZ-ZZZZ,ZZ""" \
            """,999,,,,ZZ-ZZZZ,ZZZ,ZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,,,99.99,99.99""" \
            """,,,9.99,9.99,9.99,9.99,ZZZ-ZZZZ,,,ZZZ-ZZZZZ,,,,,-9.99,-""" \
            """9.99,-9.99,-9.99,,,,,,,,,ZZZ-ZZZZ,,9,9.99,9.99,99ZZ,,-9""" \
            """.99,-9.99,ZZZ-ZZZZ,,,,,,,ZZZ-ZZZZ,9.99,9.99,9999,,,,,,,""" \
            """,,,-9.9,Z/Z-ZZZZ,999.99,9.99,,999.99,ZZ-ZZZZ,ZZ-ZZZZ,9.""" \
            """99,9.99,9.99,9.99,9.99,9.99,,ZZZ-ZZZZZ,ZZZ-ZZZZZ,ZZZ-ZZ""" \
            """ZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZZ,ZZZ,ZZZ,ZZZ,ZZZ,9.99,,,-9.99,ZZ""" \
            """-ZZZZ,-999.99,,-9999,,999.99,,,,999.99,99.99,,,ZZ-ZZZZZ""" \
            """ZZZ,ZZ-ZZZZ-ZZZZZZZ,,,,ZZ-ZZ-ZZZZZZZZ,ZZZZZZZZ,ZZZ-ZZZZ""" \
            """,9999,999.99,ZZZ-ZZZZ,-9.99,-9.99,ZZZ-ZZZZ,99:99:99,,99""" \
            """,99,,9.99,,-99.99,,,,,,9.99,ZZZ-ZZZZ,-9.99,-9.99,9.99,9""" \
            """.99,,ZZZ,,,,,,,ZZZ,ZZZ,,,,,"""

        # Set the number of lines so that a call to `parser_trim_buffers`
        # is triggered: after a couple of full chunks are consumed a
        # relatively small 'residual' chunk would cause reallocation
        # within the parser.
        chunksize, n_lines = 128, 2 * 128 + 15
        csv_data = "\n".join([record_] * n_lines) + "\n"

        # We will use StringIO to load the CSV from this text buffer.
        # pd.read_csv() will iterate over the file in chunks and will
        # finally read a residual chunk of really small size.

        # Generate the expected output: manually create the dataframe
        # by splitting by comma and repeating the `n_lines` times.
        row = tuple(val_ if val_ else float("nan")
                    for val_ in record_.split(","))
        expected = pd.DataFrame([row for _ in range(n_lines)],
                                dtype=object, columns=None, index=None)

        # Iterate over the CSV file in chunks of `chunksize` lines
        chunks_ = self.read_csv(StringIO(csv_data), header=None,
                                dtype=object, chunksize=chunksize)
        result = pd.concat(chunks_, axis=0, ignore_index=True)

        # Check for data corruption if there was no segfault
        tm.assert_frame_equal(result, expected)

    def test_internal_null_byte(self):
        # see gh-14012
        #
        # The null byte ('\x00') should not be used as a
        # true line terminator, escape character, or comment
        # character, only as a placeholder to indicate that
        # none was specified.
        #
        # This test should be moved to common.py ONLY when
        # Python's csv class supports parsing '\x00'.
        names = ['a', 'b', 'c']
        data = "1,2,3\n4,\x00,6\n7,8,9"
        expected = pd.DataFrame([[1, 2.0, 3], [4, np.nan, 6],
                                 [7, 8, 9]], columns=names)

        result = self.read_csv(StringIO(data), names=names)
        tm.assert_frame_equal(result, expected)

    def test_empty_dtype(self):
        # see gh-14712
        data = 'a,b'

        expected = pd.DataFrame(columns=['a', 'b'], dtype=np.float64)
        result = self.read_csv(StringIO(data), header=0, dtype=np.float64)
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame({'a': pd.Categorical([]),
                                 'b': pd.Categorical([])},
                                index=[])
        result = self.read_csv(StringIO(data), header=0,
                               dtype='category')
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(columns=['a', 'b'], dtype='datetime64[ns]')
        result = self.read_csv(StringIO(data), header=0,
                               dtype='datetime64[ns]')
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame({'a': pd.Series([], dtype='timedelta64[ns]'),
                                 'b': pd.Series([], dtype='timedelta64[ns]')},
                                index=[])
        result = self.read_csv(StringIO(data), header=0,
                               dtype='timedelta64[ns]')
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(columns=['a', 'b'])
        expected['a'] = expected['a'].astype(np.float64)
        result = self.read_csv(StringIO(data), header=0,
                               dtype={'a': np.float64})
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(columns=['a', 'b'])
        expected['a'] = expected['a'].astype(np.float64)
        result = self.read_csv(StringIO(data), header=0,
                               dtype={0: np.float64})
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(columns=['a', 'b'])
        expected['a'] = expected['a'].astype(np.int32)
        expected['b'] = expected['b'].astype(np.float64)
        result = self.read_csv(StringIO(data), header=0,
                               dtype={'a': np.int32, 1: np.float64})
        tm.assert_frame_equal(result, expected)

    def test_read_nrows_large(self):
        # gh-7626 - Read only nrows of data in for large inputs (>262144b)
        header_narrow = '\t'.join(['COL_HEADER_' + str(i)
                                  for i in range(10)]) + '\n'
        data_narrow = '\t'.join(['somedatasomedatasomedata1'
                                for i in range(10)]) + '\n'
        header_wide = '\t'.join(['COL_HEADER_' + str(i)
                                for i in range(15)]) + '\n'
        data_wide = '\t'.join(['somedatasomedatasomedata2'
                              for i in range(15)]) + '\n'
        test_input = (header_narrow + data_narrow * 1050 +
                      header_wide + data_wide * 2)

        df = self.read_csv(StringIO(test_input), sep='\t', nrows=1010)

        self.assertTrue(df.size == 1010 * 10)
