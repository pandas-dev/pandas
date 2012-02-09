from StringIO import StringIO
import sys
import unittest

from numpy import nan
from numpy.random import randn
import numpy as np

from pandas import DataFrame, Series
import pandas.core.format as fmt
import pandas.util.testing as tm

_frame = DataFrame(tm.getSeriesData())

class TestDataFrameFormatting(unittest.TestCase):

    def setUp(self):
        self.frame = _frame.copy()

    def test_repr_embedded_ndarray(self):
        arr = np.empty(10, dtype=[('err', object)])
        for i in range(len(arr)):
            arr['err'][i] = np.random.randn(i)

        df = DataFrame(arr)
        repr(df['err'])
        repr(df)
        df.to_string()

    def test_eng_float_formatter(self):
        self.frame.ix[5] = 0

        fmt.set_eng_float_format()
        result = repr(self.frame)

        fmt.set_eng_float_format(use_eng_prefix=True)
        repr(self.frame)

        fmt.set_eng_float_format(accuracy=0)
        repr(self.frame)

        fmt.reset_printoptions()

    def test_repr_tuples(self):
        buf = StringIO()

        df = DataFrame({'tups' : zip(range(10), range(10))})
        repr(df)
        df.to_string(col_space=10, buf=buf)

    def test_to_string_repr_unicode(self):
        buf = StringIO()

        unicode_values = [u'\u03c3'] * 10
        unicode_values = np.array(unicode_values, dtype=object)
        df = DataFrame({'unicode' : unicode_values})
        df.to_string(col_space=10, buf=buf)

        # it works!
        repr(df)

        # it works even if sys.stdin in None
        sys.stdin = None
        repr(df)
        sys.stdin = sys.__stdin__

    def test_to_string_unicode_columns(self):
        df = DataFrame({u'\u03c3' : np.arange(10.)})

        buf = StringIO()
        df.to_string(buf=buf)
        buf.getvalue()

        buf = StringIO()
        df.info(buf=buf)
        buf.getvalue()

        result = self.frame.to_string(force_unicode=True)
        self.assert_(isinstance(result, unicode))

    def test_to_string_unicode_two(self):
        dm = DataFrame({u'c/\u03c3': []})
        buf = StringIO()
        dm.to_string(buf)

    def test_to_string_with_formatters_unicode(self):
        df = DataFrame({u'c/\u03c3':[1,2,3]})
        result = df.to_string(formatters={u'c/\u03c3': lambda x: '%s' % x})
        self.assertEqual(result, (u'  c/\u03c3\n'
                                   '0   1\n'
                                   '1   2\n'
                                   '2   3'))

    def test_to_string_buffer_all_unicode(self):
        buf = StringIO()

        empty = DataFrame({u'c/\u03c3':Series()})
        nonempty = DataFrame({u'c/\u03c3':Series([1,2,3])})

        print >>buf, empty
        print >>buf, nonempty

        # this should work
        buf.getvalue()

    def test_unicode_problem_decoding_as_ascii(self):
        dm = DataFrame({u'c/\u03c3': Series({'test':np.NaN})})
        unicode(dm.to_string())

    def test_repr_corner(self):
        # representing infs poses no problems
        df = DataFrame({'foo' : np.inf * np.empty(10)})
        foo = repr(df)

    def test_to_string(self):
        from pandas import read_table
        import re

        # big mixed
        biggie = DataFrame({'A' : randn(200),
                            'B' : tm.makeStringIndex(200)},
                            index=range(200))

        biggie['A'][:20] = nan
        biggie['B'][:20] = nan
        s = biggie.to_string()

        buf = StringIO()
        retval = biggie.to_string(buf=buf)
        self.assert_(retval is None)
        self.assertEqual(buf.getvalue(), s)

        self.assert_(isinstance(s, basestring))

        # print in right order
        result = biggie.to_string(columns=['B', 'A'], col_space=17,
                                  float_format='%.5f'.__mod__)
        lines = result.split('\n')
        header = lines[0].strip().split()
        joined = '\n'.join([re.sub('\s+', ' ', x).strip() for x in lines[1:]])
        recons = read_table(StringIO(joined), names=header, sep=' ')
        tm.assert_series_equal(recons['B'], biggie['B'])
        self.assertEqual(recons['A'].count(), biggie['A'].count())
        self.assert_((np.abs(recons['A'].dropna() -
                             biggie['A'].dropna()) < 0.1).all())

        # expected = ['B', 'A']
        # self.assertEqual(header, expected)

        result = biggie.to_string(columns=['A'], col_space=17)
        header = result.split('\n')[0].strip().split()
        expected = ['A']
        self.assertEqual(header, expected)

        biggie.to_string(columns=['B', 'A'],
                         formatters={'A' : lambda x: '%.1f' % x})

        biggie.to_string(columns=['B', 'A'], float_format=str)
        biggie.to_string(columns=['B', 'A'], col_space=12,
                         float_format=str)

        frame = DataFrame(index=np.arange(200))
        frame.to_string()

    def test_to_string_no_header(self):
        df = DataFrame({'x' : [1, 2, 3],
                        'y' : [4, 5, 6]})

        df_s = df.to_string(header=False)
        expected = "0  1  4\n1  2  5\n2  3  6"

        assert(df_s == expected)

    def test_to_string_no_index(self):
        df = DataFrame({'x' : [1, 2, 3],
                        'y' : [4, 5, 6]})

        df_s = df.to_string(index=False)
        expected = " x  y\n 1  4\n 2  5\n 3  6"

        assert(df_s == expected)

    def test_to_string_float_formatting(self):
        fmt.reset_printoptions()
        fmt.set_printoptions(precision=6, column_space=12)

        df = DataFrame({'x' : [0, 0.25, 3456.000, 12e+45, 1.64e+6,
                               1.7e+8, 1.253456, np.pi, -1e6]})

        df_s = df.to_string()

        # Python 2.5 just wants me to be sad. And debian 32-bit
        #sys.version_info[0] == 2 and sys.version_info[1] < 6:
        if '%.4g' % 1.7e8 == '1.7e+008':
            expected = ('              x\n0  0.00000e+000\n1  2.50000e-001\n'
                        '2  3.45600e+003\n3  1.20000e+046\n4  1.64000e+006\n'
                        '5  1.70000e+008\n6  1.25346e+000\n7  3.14159e+000\n'
                        '8 -1.00000e+006')
        else:
            expected = ('             x\n0  0.00000e+00\n1  2.50000e-01\n'
                        '2  3.45600e+03\n3  1.20000e+46\n4  1.64000e+06\n'
                        '5  1.70000e+08\n6  1.25346e+00\n7  3.14159e+00\n'
                        '8 -1.00000e+06')
        assert(df_s == expected)

        df = DataFrame({'x' : [3234, 0.253]})
        df_s = df.to_string()

        expected = ('          x\n'
                    '0  3234.000\n'
                    '1     0.253')
        assert(df_s == expected)

        fmt.reset_printoptions()
        self.assertEqual(fmt.print_config.precision, 7)

        df = DataFrame({'x': [1e9, 0.2512]})
        df_s = df.to_string()
        # Python 2.5 just wants me to be sad. And debian 32-bit
        #sys.version_info[0] == 2 and sys.version_info[1] < 6:
        if '%.4g' % 1.7e8 == '1.7e+008':
            expected = ('               x\n'
                        '0  1.000000e+009\n'
                        '1  2.512000e-001')
        else:
            expected = ('              x\n'
                        '0  1.000000e+09\n'
                        '1  2.512000e-01')
        assert(df_s == expected)

    def test_to_string_left_justify_cols(self):
        fmt.reset_printoptions()
        df = DataFrame({'x' : [3234, 0.253]})
        df_s = df.to_string(justify='left')
        expected = ('   x       \n'
                    '0  3234.000\n'
                    '1     0.253')
        assert(df_s == expected)

    def test_to_string_format_na(self):
        fmt.reset_printoptions()
        df = DataFrame({'A' : [np.nan, -1, -2.1234, 3, 4],
                        'B' : [np.nan, 'foo', 'foooo', 'fooooo', 'bar']})
        result = df.to_string()

        expected = ('        A       B\n'
                    '0     NaN     NaN\n'
                    '1 -1.0000     foo\n'
                    '2 -2.1234   foooo\n'
                    '3  3.0000  fooooo\n'
                    '4  4.0000     bar')
        self.assertEqual(result, expected)

        df = DataFrame({'A' : [np.nan, -1., -2., 3., 4.],
                        'B' : [np.nan, 'foo', 'foooo', 'fooooo', 'bar']})
        result = df.to_string()

        expected = ('    A       B\n'
                    '0 NaN     NaN\n'
                    '1  -1     foo\n'
                    '2  -2   foooo\n'
                    '3   3  fooooo\n'
                    '4   4     bar')
        self.assertEqual(result, expected)

    def test_to_html(self):
        # big mixed
        biggie = DataFrame({'A' : randn(200),
                            'B' : tm.makeStringIndex(200)},
                            index=range(200))

        biggie['A'][:20] = nan
        biggie['B'][:20] = nan
        s = biggie.to_html()

        buf = StringIO()
        retval = biggie.to_html(buf=buf)
        self.assert_(retval is None)
        self.assertEqual(buf.getvalue(), s)

        self.assert_(isinstance(s, basestring))

        biggie.to_html(columns=['B', 'A'], col_space=17)
        biggie.to_html(columns=['B', 'A'],
                       formatters={'A' : lambda x: '%.1f' % x})

        biggie.to_html(columns=['B', 'A'], float_format=str)
        biggie.to_html(columns=['B', 'A'], col_space=12,
                       float_format=str)

        frame = DataFrame(index=np.arange(200))
        frame.to_html()

    def test_to_html_with_no_bold(self):
        x = DataFrame({'x': randn(5)})
        ashtml = x.to_html(bold_rows=False)
        assert('<strong>' not in ashtml)


class TestSeriesFormatting(unittest.TestCase):

    def setUp(self):
        self.ts = tm.makeTimeSeries()

    def test_repr_unicode(self):
        s = Series([u'\u03c3'] * 10)
        repr(s)

    def test_to_string(self):
        from cStringIO import StringIO
        buf = StringIO()

        s = self.ts.to_string()

        retval = self.ts.to_string(buf=buf)
        self.assert_(retval is None)
        self.assertEqual(buf.getvalue().strip(), s)

        # pass float_format
        format = '%.4f'.__mod__
        result = self.ts.to_string(float_format=format)
        result = [x.split()[1] for x in result.split('\n')]
        expected = [format(x) for x in self.ts]
        self.assertEqual(result, expected)

        # empty string
        result = self.ts[:0].to_string()
        self.assertEqual(result, '')

        result = self.ts[:0].to_string(length=0)
        self.assertEqual(result, '')

        # name and length
        cp = self.ts.copy()
        cp.name = 'foo'
        result = cp.to_string(length=True, name=True)
        last_line = result.split('\n')[-1].strip()
        self.assertEqual(last_line, "Name: foo, Length: %d" % len(cp))

    def test_to_string_mixed(self):
        s = Series(['foo', np.nan, -1.23, 4.56])
        result = s.to_string()
        expected = ('0     foo\n'
                    '1     NaN\n'
                    '2   -1.23\n'
                    '3    4.56')
        self.assertEqual(result, expected)

        # but don't count NAs as floats
        s = Series(['foo', np.nan, 'bar', 'baz'])
        result = s.to_string()
        expected = ('0    foo\n'
                    '1    NaN\n'
                    '2    bar\n'
                    '3    baz')
        self.assertEqual(result, expected)

        s = Series(['foo', 5, 'bar', 'baz'])
        result = s.to_string()
        expected = ('0    foo\n'
                    '1      5\n'
                    '2    bar\n'
                    '3    baz')
        self.assertEqual(result, expected)

    def test_to_string_float_na_spacing(self):
        s = Series([0., 1.5678, 2., -3., 4.])
        s[::2] = np.nan

        result = s.to_string()
        expected = ('0       NaN\n'
                    '1    1.5678\n'
                    '2       NaN\n'
                    '3   -3.0000\n'
                    '4       NaN')
        self.assertEqual(result, expected)

class TestEngFormatter(unittest.TestCase):

    def test_eng_float_formatter(self):
        df = DataFrame({'A' : [1.41, 141., 14100, 1410000.]})

        fmt.set_eng_float_format()
        result = df.to_string()
        expected = ('             A\n'
                    '0    1.410E+00\n'
                    '1  141.000E+00\n'
                    '2   14.100E+03\n'
                    '3    1.410E+06')
        self.assertEqual(result, expected)

        fmt.set_eng_float_format(use_eng_prefix=True)
        result = df.to_string()
        expected = ('         A\n'
                    '0    1.410\n'
                    '1  141.000\n'
                    '2  14.100k\n'
                    '3   1.410M')
        self.assertEqual(result, expected)

        fmt.set_eng_float_format(accuracy=0)
        result = df.to_string()
        expected = ('         A\n'
                    '0    1E+00\n'
                    '1  141E+00\n'
                    '2   14E+03\n'
                    '3    1E+06')
        self.assertEqual(result, expected)

        fmt.reset_printoptions()

    def compare(self, formatter, input, output):
        formatted_input = formatter(input)
        msg = ("formatting of %s results in '%s', expected '%s'"
               % (str(input), formatted_input, output))
        self.assertEqual(formatted_input, output, msg)

    def compare_all(self, formatter, in_out):
        """
        Parameters:
        -----------
        formatter: EngFormatter under test
        in_out: list of tuples. Each tuple = (number, expected_formatting)

        It is tested if 'formatter(number) == expected_formatting'.
        *number* should be >= 0 because formatter(-number) == fmt is also
        tested. *fmt* is derived from *expected_formatting*
        """
        for input, output in in_out:
            self.compare(formatter, input, output)
            self.compare(formatter, -input, "-" + output[1:])

    def test_exponents_with_eng_prefix(self):
        formatter = fmt.EngFormatter(accuracy=3, use_eng_prefix=True)
        f = np.sqrt(2)
        in_out = [(f * 10 ** -24, " 1.414y"),
                  (f * 10 ** -23, " 14.142y"),
                  (f * 10 ** -22, " 141.421y"),
                  (f * 10 ** -21, " 1.414z"),
                  (f * 10 ** -20, " 14.142z"),
                  (f * 10 ** -19, " 141.421z"),
                  (f * 10 ** -18, " 1.414a"),
                  (f * 10 ** -17, " 14.142a"),
                  (f * 10 ** -16, " 141.421a"),
                  (f * 10 ** -15, " 1.414f"),
                  (f * 10 ** -14, " 14.142f"),
                  (f * 10 ** -13, " 141.421f"),
                  (f * 10 ** -12, " 1.414p"),
                  (f * 10 ** -11, " 14.142p"),
                  (f * 10 ** -10, " 141.421p"),
                  (f * 10 ** -9, " 1.414n"),
                  (f * 10 ** -8, " 14.142n"),
                  (f * 10 ** -7, " 141.421n"),
                  (f * 10 ** -6, " 1.414u"),
                  (f * 10 ** -5, " 14.142u"),
                  (f * 10 ** -4, " 141.421u"),
                  (f * 10 ** -3, " 1.414m"),
                  (f * 10 ** -2, " 14.142m"),
                  (f * 10 ** -1, " 141.421m"),
                  (f * 10 ** 0, " 1.414"),
                  (f * 10 ** 1, " 14.142"),
                  (f * 10 ** 2, " 141.421"),
                  (f * 10 ** 3, " 1.414k"),
                  (f * 10 ** 4, " 14.142k"),
                  (f * 10 ** 5, " 141.421k"),
                  (f * 10 ** 6, " 1.414M"),
                  (f * 10 ** 7, " 14.142M"),
                  (f * 10 ** 8, " 141.421M"),
                  (f * 10 ** 9, " 1.414G"),
                  (f * 10 ** 10, " 14.142G"),
                  (f * 10 ** 11, " 141.421G"),
                  (f * 10 ** 12, " 1.414T"),
                  (f * 10 ** 13, " 14.142T"),
                  (f * 10 ** 14, " 141.421T"),
                  (f * 10 ** 15, " 1.414P"),
                  (f * 10 ** 16, " 14.142P"),
                  (f * 10 ** 17, " 141.421P"),
                  (f * 10 ** 18, " 1.414E"),
                  (f * 10 ** 19, " 14.142E"),
                  (f * 10 ** 20, " 141.421E"),
                  (f * 10 ** 21, " 1.414Z"),
                  (f * 10 ** 22, " 14.142Z"),
                  (f * 10 ** 23, " 141.421Z"),
                  (f * 10 ** 24, " 1.414Y"),
                  (f * 10 ** 25, " 14.142Y"),
                  (f * 10 ** 26, " 141.421Y")]
        self.compare_all(formatter, in_out)

    def test_exponents_without_eng_prefix(self):
        formatter = fmt.EngFormatter(accuracy=4, use_eng_prefix=False)
        f = np.pi
        in_out = [(f * 10 ** -24, " 3.1416E-24"),
                  (f * 10 ** -23, " 31.4159E-24"),
                  (f * 10 ** -22, " 314.1593E-24"),
                  (f * 10 ** -21, " 3.1416E-21"),
                  (f * 10 ** -20, " 31.4159E-21"),
                  (f * 10 ** -19, " 314.1593E-21"),
                  (f * 10 ** -18, " 3.1416E-18"),
                  (f * 10 ** -17, " 31.4159E-18"),
                  (f * 10 ** -16, " 314.1593E-18"),
                  (f * 10 ** -15, " 3.1416E-15"),
                  (f * 10 ** -14, " 31.4159E-15"),
                  (f * 10 ** -13, " 314.1593E-15"),
                  (f * 10 ** -12, " 3.1416E-12"),
                  (f * 10 ** -11, " 31.4159E-12"),
                  (f * 10 ** -10, " 314.1593E-12"),
                  (f * 10 ** -9, " 3.1416E-09"),
                  (f * 10 ** -8, " 31.4159E-09"),
                  (f * 10 ** -7, " 314.1593E-09"),
                  (f * 10 ** -6, " 3.1416E-06"),
                  (f * 10 ** -5, " 31.4159E-06"),
                  (f * 10 ** -4, " 314.1593E-06"),
                  (f * 10 ** -3, " 3.1416E-03"),
                  (f * 10 ** -2, " 31.4159E-03"),
                  (f * 10 ** -1, " 314.1593E-03"),
                  (f * 10 ** 0, " 3.1416E+00"),
                  (f * 10 ** 1, " 31.4159E+00"),
                  (f * 10 ** 2, " 314.1593E+00"),
                  (f * 10 ** 3, " 3.1416E+03"),
                  (f * 10 ** 4, " 31.4159E+03"),
                  (f * 10 ** 5, " 314.1593E+03"),
                  (f * 10 ** 6, " 3.1416E+06"),
                  (f * 10 ** 7, " 31.4159E+06"),
                  (f * 10 ** 8, " 314.1593E+06"),
                  (f * 10 ** 9, " 3.1416E+09"),
                  (f * 10 ** 10, " 31.4159E+09"),
                  (f * 10 ** 11, " 314.1593E+09"),
                  (f * 10 ** 12, " 3.1416E+12"),
                  (f * 10 ** 13, " 31.4159E+12"),
                  (f * 10 ** 14, " 314.1593E+12"),
                  (f * 10 ** 15, " 3.1416E+15"),
                  (f * 10 ** 16, " 31.4159E+15"),
                  (f * 10 ** 17, " 314.1593E+15"),
                  (f * 10 ** 18, " 3.1416E+18"),
                  (f * 10 ** 19, " 31.4159E+18"),
                  (f * 10 ** 20, " 314.1593E+18"),
                  (f * 10 ** 21, " 3.1416E+21"),
                  (f * 10 ** 22, " 31.4159E+21"),
                  (f * 10 ** 23, " 314.1593E+21"),
                  (f * 10 ** 24, " 3.1416E+24"),
                  (f * 10 ** 25, " 31.4159E+24"),
                  (f * 10 ** 26, " 314.1593E+24")]
        self.compare_all(formatter, in_out)

    def test_rounding(self):
        formatter = fmt.EngFormatter(accuracy=3, use_eng_prefix=True)
        in_out = [(5.55555, ' 5.556'),
                  (55.5555, ' 55.556'),
                  (555.555, ' 555.555'),
                  (5555.55, ' 5.556k'),
                  (55555.5, ' 55.556k'),
                  (555555, ' 555.555k')]
        self.compare_all(formatter, in_out)

        formatter = fmt.EngFormatter(accuracy=1, use_eng_prefix=True)
        in_out = [(5.55555, ' 5.6'),
                  (55.5555, ' 55.6'),
                  (555.555, ' 555.6'),
                  (5555.55, ' 5.6k'),
                  (55555.5, ' 55.6k'),
                  (555555, ' 555.6k')]
        self.compare_all(formatter, in_out)

        formatter = fmt.EngFormatter(accuracy=0, use_eng_prefix=True)
        in_out = [(5.55555, ' 6'),
                  (55.5555, ' 56'),
                  (555.555, ' 556'),
                  (5555.55, ' 6k'),
                  (55555.5, ' 56k'),
                  (555555, ' 556k')]
        self.compare_all(formatter, in_out)

        formatter = fmt.EngFormatter(accuracy=3, use_eng_prefix=True)
        result = formatter(0)
        self.assertEqual(result, u' 0.000')


class TestFloatArrayFormatter(unittest.TestCase):

    def test_misc(self):
        obj = fmt.FloatArrayFormatter(np.array([], dtype=np.float64))
        result = obj.get_result()

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

