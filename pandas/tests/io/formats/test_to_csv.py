# -*- coding: utf-8 -*-

import os
import sys

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, compat
from pandas.util import testing as tm


class TestToCSV(object):

    @pytest.mark.xfail((3, 6, 5) > sys.version_info >= (3, 5),
                       reason=("Python csv library bug "
                               "(see https://bugs.python.org/issue32255)"))
    def test_to_csv_with_single_column(self):
        # see gh-18676, https://bugs.python.org/issue32255
        #
        # Python's CSV library adds an extraneous '""'
        # before the newline when the NaN-value is in
        # the first row. Otherwise, only the newline
        # character is added. This behavior is inconsistent
        # and was patched in https://bugs.python.org/pull_request4672.
        df1 = DataFrame([None, 1])
        expected1 = """\
""
1.0
"""
        with tm.ensure_clean('test.csv') as path:
            df1.to_csv(path, header=None, index=None)
            with open(path, 'r') as f:
                assert f.read() == expected1

        df2 = DataFrame([1, None])
        expected2 = """\
1.0
""
"""
        with tm.ensure_clean('test.csv') as path:
            df2.to_csv(path, header=None, index=None)
            with open(path, 'r') as f:
                assert f.read() == expected2

    def test_to_csv_defualt_encoding(self):
        # GH17097
        df = DataFrame({'col': ["AAAAA", "ÄÄÄÄÄ", "ßßßßß", "聞聞聞聞聞"]})

        with tm.ensure_clean('test.csv') as path:
            # the default to_csv encoding in Python 2 is ascii, and that in
            # Python 3 is uft-8.
            df.to_csv(path)
            tm.assert_frame_equal(pd.read_csv(path, index_col=0), df)

    def test_to_csv_quotechar(self):
        df = DataFrame({'col': [1, 2]})
        expected = """\
"","col"
"0","1"
"1","2"
"""

        with tm.ensure_clean('test.csv') as path:
            df.to_csv(path, quoting=1)  # 1=QUOTE_ALL
            with open(path, 'r') as f:
                assert f.read() == expected

        expected = """\
$$,$col$
$0$,$1$
$1$,$2$
"""

        with tm.ensure_clean('test.csv') as path:
            df.to_csv(path, quoting=1, quotechar="$")
            with open(path, 'r') as f:
                assert f.read() == expected

        with tm.ensure_clean('test.csv') as path:
            with pytest.raises(TypeError, match='quotechar'):
                df.to_csv(path, quoting=1, quotechar=None)

    def test_to_csv_doublequote(self):
        df = DataFrame({'col': ['a"a', '"bb"']})
        expected = '''\
"","col"
"0","a""a"
"1","""bb"""
'''

        with tm.ensure_clean('test.csv') as path:
            df.to_csv(path, quoting=1, doublequote=True)  # QUOTE_ALL
            with open(path, 'r') as f:
                assert f.read() == expected

        from _csv import Error
        with tm.ensure_clean('test.csv') as path:
            with pytest.raises(Error, match='escapechar'):
                df.to_csv(path, doublequote=False)  # no escapechar set

    def test_to_csv_escapechar(self):
        df = DataFrame({'col': ['a"a', '"bb"']})
        expected = '''\
"","col"
"0","a\\"a"
"1","\\"bb\\""
'''

        with tm.ensure_clean('test.csv') as path:  # QUOTE_ALL
            df.to_csv(path, quoting=1, doublequote=False, escapechar='\\')
            with open(path, 'r') as f:
                assert f.read() == expected

        df = DataFrame({'col': ['a,a', ',bb,']})
        expected = """\
,col
0,a\\,a
1,\\,bb\\,
"""

        with tm.ensure_clean('test.csv') as path:
            df.to_csv(path, quoting=3, escapechar='\\')  # QUOTE_NONE
            with open(path, 'r') as f:
                assert f.read() == expected

    def test_csv_to_string(self):
        df = DataFrame({'col': [1, 2]})
        expected_rows = [',col',
                         '0,1',
                         '1,2']
        expected = tm.convert_rows_list_to_csv_str(expected_rows)
        assert df.to_csv() == expected

    def test_to_csv_decimal(self):
        # see gh-781
        df = DataFrame({'col1': [1], 'col2': ['a'], 'col3': [10.1]})

        expected_rows = [',col1,col2,col3',
                         '0,1,a,10.1']
        expected_default = tm.convert_rows_list_to_csv_str(expected_rows)
        assert df.to_csv() == expected_default

        expected_rows = [';col1;col2;col3',
                         '0;1;a;10,1']
        expected_european_excel = tm.convert_rows_list_to_csv_str(
            expected_rows)
        assert df.to_csv(decimal=',', sep=';') == expected_european_excel

        expected_rows = [',col1,col2,col3',
                         '0,1,a,10.10']
        expected_float_format_default = tm.convert_rows_list_to_csv_str(
            expected_rows)
        assert df.to_csv(float_format='%.2f') == expected_float_format_default

        expected_rows = [';col1;col2;col3',
                         '0;1;a;10,10']
        expected_float_format = tm.convert_rows_list_to_csv_str(expected_rows)
        assert df.to_csv(decimal=',', sep=';',
                         float_format='%.2f') == expected_float_format

        # see gh-11553: testing if decimal is taken into account for '0.0'
        df = pd.DataFrame({'a': [0, 1.1], 'b': [2.2, 3.3], 'c': 1})

        expected_rows = ['a,b,c',
                         '0^0,2^2,1',
                         '1^1,3^3,1']
        expected = tm.convert_rows_list_to_csv_str(expected_rows)
        assert df.to_csv(index=False, decimal='^') == expected

        # same but for an index
        assert df.set_index('a').to_csv(decimal='^') == expected

        # same for a multi-index
        assert df.set_index(['a', 'b']).to_csv(decimal="^") == expected

    def test_to_csv_float_format(self):
        # testing if float_format is taken into account for the index
        # GH 11553
        df = pd.DataFrame({'a': [0, 1], 'b': [2.2, 3.3], 'c': 1})

        expected_rows = ['a,b,c',
                         '0,2.20,1',
                         '1,3.30,1']
        expected = tm.convert_rows_list_to_csv_str(expected_rows)
        assert df.set_index('a').to_csv(float_format='%.2f') == expected

        # same for a multi-index
        assert df.set_index(['a', 'b']).to_csv(
            float_format='%.2f') == expected

    def test_to_csv_na_rep(self):
        # see gh-11553
        #
        # Testing if NaN values are correctly represented in the index.
        df = DataFrame({'a': [0, np.NaN], 'b': [0, 1], 'c': [2, 3]})
        expected_rows = ['a,b,c',
                         '0.0,0,2',
                         '_,1,3']
        expected = tm.convert_rows_list_to_csv_str(expected_rows)

        assert df.set_index('a').to_csv(na_rep='_') == expected
        assert df.set_index(['a', 'b']).to_csv(na_rep='_') == expected

        # now with an index containing only NaNs
        df = DataFrame({'a': np.NaN, 'b': [0, 1], 'c': [2, 3]})
        expected_rows = ['a,b,c',
                         '_,0,2',
                         '_,1,3']
        expected = tm.convert_rows_list_to_csv_str(expected_rows)

        assert df.set_index('a').to_csv(na_rep='_') == expected
        assert df.set_index(['a', 'b']).to_csv(na_rep='_') == expected

        # check if na_rep parameter does not break anything when no NaN
        df = DataFrame({'a': 0, 'b': [0, 1], 'c': [2, 3]})
        expected_rows = ['a,b,c',
                         '0,0,2',
                         '0,1,3']
        expected = tm.convert_rows_list_to_csv_str(expected_rows)

        assert df.set_index('a').to_csv(na_rep='_') == expected
        assert df.set_index(['a', 'b']).to_csv(na_rep='_') == expected

    def test_to_csv_date_format(self):
        # GH 10209
        df_sec = DataFrame({'A': pd.date_range('20130101', periods=5, freq='s')
                            })
        df_day = DataFrame({'A': pd.date_range('20130101', periods=5, freq='d')
                            })

        expected_rows = [',A',
                         '0,2013-01-01 00:00:00',
                         '1,2013-01-01 00:00:01',
                         '2,2013-01-01 00:00:02',
                         '3,2013-01-01 00:00:03',
                         '4,2013-01-01 00:00:04']
        expected_default_sec = tm.convert_rows_list_to_csv_str(expected_rows)
        assert df_sec.to_csv() == expected_default_sec

        expected_rows = [',A',
                         '0,2013-01-01 00:00:00',
                         '1,2013-01-02 00:00:00',
                         '2,2013-01-03 00:00:00',
                         '3,2013-01-04 00:00:00',
                         '4,2013-01-05 00:00:00']
        expected_ymdhms_day = tm.convert_rows_list_to_csv_str(expected_rows)
        assert (df_day.to_csv(date_format='%Y-%m-%d %H:%M:%S') ==
                expected_ymdhms_day)

        expected_rows = [',A',
                         '0,2013-01-01',
                         '1,2013-01-01',
                         '2,2013-01-01',
                         '3,2013-01-01',
                         '4,2013-01-01']
        expected_ymd_sec = tm.convert_rows_list_to_csv_str(expected_rows)
        assert df_sec.to_csv(date_format='%Y-%m-%d') == expected_ymd_sec

        expected_rows = [',A',
                         '0,2013-01-01',
                         '1,2013-01-02',
                         '2,2013-01-03',
                         '3,2013-01-04',
                         '4,2013-01-05']
        expected_default_day = tm.convert_rows_list_to_csv_str(expected_rows)
        assert df_day.to_csv() == expected_default_day
        assert df_day.to_csv(date_format='%Y-%m-%d') == expected_default_day

        # see gh-7791
        #
        # Testing if date_format parameter is taken into account
        # for multi-indexed DataFrames.
        df_sec['B'] = 0
        df_sec['C'] = 1

        expected_rows = ['A,B,C',
                         '2013-01-01,0,1']
        expected_ymd_sec = tm.convert_rows_list_to_csv_str(expected_rows)

        df_sec_grouped = df_sec.groupby([pd.Grouper(key='A', freq='1h'), 'B'])
        assert (df_sec_grouped.mean().to_csv(date_format='%Y-%m-%d') ==
                expected_ymd_sec)

    def test_to_csv_multi_index(self):
        # see gh-6618
        df = DataFrame([1], columns=pd.MultiIndex.from_arrays([[1], [2]]))

        exp_rows = [',1',
                    ',2',
                    '0,1']
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv() == exp

        exp_rows = ['1', '2', '1']
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv(index=False) == exp

        df = DataFrame([1], columns=pd.MultiIndex.from_arrays([[1], [2]]),
                       index=pd.MultiIndex.from_arrays([[1], [2]]))

        exp_rows = [',,1', ',,2', '1,2,1']
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv() == exp

        exp_rows = ['1', '2', '1']
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv(index=False) == exp

        df = DataFrame(
            [1], columns=pd.MultiIndex.from_arrays([['foo'], ['bar']]))

        exp_rows = [',foo', ',bar', '0,1']
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv() == exp

        exp_rows = ['foo', 'bar', '1']
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv(index=False) == exp

    @pytest.mark.parametrize("ind,expected", [
        (pd.MultiIndex(levels=[[1.0]],
                       codes=[[0]],
                       names=["x"]),
         "x,data\n1.0,1\n"),
        (pd.MultiIndex(levels=[[1.], [2.]],
                       codes=[[0], [0]],
                       names=["x", "y"]),
         "x,y,data\n1.0,2.0,1\n")
    ])
    @pytest.mark.parametrize("klass", [
        pd.DataFrame, pd.Series
    ])
    def test_to_csv_single_level_multi_index(self, ind, expected, klass):
        # see gh-19589
        result = klass(pd.Series([1], ind, name="data")).to_csv(
            line_terminator="\n", header=True)
        assert result == expected

    def test_to_csv_string_array_ascii(self):
        # GH 10813
        str_array = [{'names': ['foo', 'bar']}, {'names': ['baz', 'qux']}]
        df = pd.DataFrame(str_array)
        expected_ascii = '''\
,names
0,"['foo', 'bar']"
1,"['baz', 'qux']"
'''
        with tm.ensure_clean('str_test.csv') as path:
            df.to_csv(path, encoding='ascii')
            with open(path, 'r') as f:
                assert f.read() == expected_ascii

    @pytest.mark.xfail(strict=False)
    def test_to_csv_string_array_utf8(self):
        # GH 10813
        str_array = [{'names': ['foo', 'bar']}, {'names': ['baz', 'qux']}]
        df = pd.DataFrame(str_array)
        expected_utf8 = '''\
,names
0,"['foo', 'bar']"
1,"['baz', 'qux']"
'''
        with tm.ensure_clean('unicode_test.csv') as path:
            df.to_csv(path, encoding='utf-8')
            with open(path, 'r') as f:
                assert f.read() == expected_utf8

    def test_to_csv_string_with_lf(self):
        # GH 20353
        data = {
            'int': [1, 2, 3],
            'str_lf': ['abc', 'd\nef', 'g\nh\n\ni']
        }
        df = pd.DataFrame(data)
        with tm.ensure_clean('lf_test.csv') as path:
            # case 1: The default line terminator(=os.linesep)(PR 21406)
            os_linesep = os.linesep.encode('utf-8')
            expected_noarg = (
                b'int,str_lf' + os_linesep +
                b'1,abc' + os_linesep +
                b'2,"d\nef"' + os_linesep +
                b'3,"g\nh\n\ni"' + os_linesep
            )
            df.to_csv(path, index=False)
            with open(path, 'rb') as f:
                assert f.read() == expected_noarg
        with tm.ensure_clean('lf_test.csv') as path:
            # case 2: LF as line terminator
            expected_lf = (
                b'int,str_lf\n'
                b'1,abc\n'
                b'2,"d\nef"\n'
                b'3,"g\nh\n\ni"\n'
            )
            df.to_csv(path, line_terminator='\n', index=False)
            with open(path, 'rb') as f:
                assert f.read() == expected_lf
        with tm.ensure_clean('lf_test.csv') as path:
            # case 3: CRLF as line terminator
            # 'line_terminator' should not change inner element
            expected_crlf = (
                b'int,str_lf\r\n'
                b'1,abc\r\n'
                b'2,"d\nef"\r\n'
                b'3,"g\nh\n\ni"\r\n'
            )
            df.to_csv(path, line_terminator='\r\n', index=False)
            with open(path, 'rb') as f:
                assert f.read() == expected_crlf

    def test_to_csv_string_with_crlf(self):
        # GH 20353
        data = {
            'int': [1, 2, 3],
            'str_crlf': ['abc', 'd\r\nef', 'g\r\nh\r\n\r\ni']
        }
        df = pd.DataFrame(data)
        with tm.ensure_clean('crlf_test.csv') as path:
            # case 1: The default line terminator(=os.linesep)(PR 21406)
            os_linesep = os.linesep.encode('utf-8')
            expected_noarg = (
                b'int,str_crlf' + os_linesep +
                b'1,abc' + os_linesep +
                b'2,"d\r\nef"' + os_linesep +
                b'3,"g\r\nh\r\n\r\ni"' + os_linesep
            )
            df.to_csv(path, index=False)
            with open(path, 'rb') as f:
                assert f.read() == expected_noarg
        with tm.ensure_clean('crlf_test.csv') as path:
            # case 2: LF as line terminator
            expected_lf = (
                b'int,str_crlf\n'
                b'1,abc\n'
                b'2,"d\r\nef"\n'
                b'3,"g\r\nh\r\n\r\ni"\n'
            )
            df.to_csv(path, line_terminator='\n', index=False)
            with open(path, 'rb') as f:
                assert f.read() == expected_lf
        with tm.ensure_clean('crlf_test.csv') as path:
            # case 3: CRLF as line terminator
            # 'line_terminator' should not change inner element
            expected_crlf = (
                b'int,str_crlf\r\n'
                b'1,abc\r\n'
                b'2,"d\r\nef"\r\n'
                b'3,"g\r\nh\r\n\r\ni"\r\n'
            )
            df.to_csv(path, line_terminator='\r\n', index=False)
            with open(path, 'rb') as f:
                assert f.read() == expected_crlf

    def test_to_csv_stdout_file(self, capsys):
        # GH 21561
        df = pd.DataFrame([['foo', 'bar'], ['baz', 'qux']],
                          columns=['name_1', 'name_2'])
        expected_rows = [',name_1,name_2',
                         '0,foo,bar',
                         '1,baz,qux']
        expected_ascii = tm.convert_rows_list_to_csv_str(expected_rows)

        df.to_csv(sys.stdout, encoding='ascii')
        captured = capsys.readouterr()

        assert captured.out == expected_ascii
        assert not sys.stdout.closed

    @pytest.mark.xfail(
        compat.is_platform_windows(),
        reason=("Especially in Windows, file stream should not be passed"
                "to csv writer without newline='' option."
                "(https://docs.python.org/3.6/library/csv.html#csv.writer)"))
    def test_to_csv_write_to_open_file(self):
        # GH 21696
        df = pd.DataFrame({'a': ['x', 'y', 'z']})
        expected = '''\
manual header
x
y
z
'''
        with tm.ensure_clean('test.txt') as path:
            with open(path, 'w') as f:
                f.write('manual header\n')
                df.to_csv(f, header=None, index=None)
            with open(path, 'r') as f:
                assert f.read() == expected

    def test_to_csv_write_to_open_file_with_newline_py3(self):
        # see gh-21696
        # see gh-20353
        df = pd.DataFrame({'a': ['x', 'y', 'z']})
        expected_rows = ["x",
                         "y",
                         "z"]
        expected = ("manual header\n" +
                    tm.convert_rows_list_to_csv_str(expected_rows))
        with tm.ensure_clean('test.txt') as path:
            with open(path, 'w', newline='') as f:
                f.write('manual header\n')
                df.to_csv(f, header=None, index=None)

            with open(path, 'rb') as f:
                assert f.read() == bytes(expected, 'utf-8')

    @pytest.mark.parametrize("to_infer", [True, False])
    @pytest.mark.parametrize("read_infer", [True, False])
    def test_to_csv_compression(self, compression_only,
                                read_infer, to_infer):
        # see gh-15008
        compression = compression_only

        if compression == "zip":
            pytest.skip("{compression} is not supported "
                        "for to_csv".format(compression=compression))

        # We'll complete file extension subsequently.
        filename = "test."

        if compression == "gzip":
            filename += "gz"
        else:
            # xz --> .xz
            # bz2 --> .bz2
            filename += compression

        df = DataFrame({"A": [1]})

        to_compression = "infer" if to_infer else compression
        read_compression = "infer" if read_infer else compression

        with tm.ensure_clean(filename) as path:
            df.to_csv(path, compression=to_compression)
            result = pd.read_csv(path, index_col=0,
                                 compression=read_compression)
            tm.assert_frame_equal(result, df)
