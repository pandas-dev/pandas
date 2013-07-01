# -*- coding: utf-8 -*-

try:
    from StringIO import StringIO
except:
    from io import StringIO

import os
import sys
import unittest
from textwrap import dedent
import warnings

from numpy import nan
from numpy.random import randn
import numpy as np

from pandas import DataFrame, Series, Index
from pandas.util.py3compat import lzip, PY3

import pandas.core.format as fmt
import pandas.util.testing as tm
from pandas.util.terminal import get_terminal_size
import pandas
import pandas as pd
from pandas.core.config import (set_option, get_option,
                                option_context, reset_option)

_frame = DataFrame(tm.getSeriesData())


def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth

def has_info_repr(df):
    r = repr(df)
    return r.split('\n')[0].startswith("<class")

def has_expanded_repr(df):
    r = repr(df)
    for line in r.split('\n'):
        if line.endswith('\\'):
            return True
    return False


class TestDataFrameFormatting(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.warn_filters = warnings.filters
        warnings.filterwarnings('ignore',
                                category=FutureWarning,
                                module=".*format")

        self.frame = _frame.copy()

    def tearDown(self):
        warnings.filters = self.warn_filters

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

        fmt.reset_option('^display.')

    def test_repr_tuples(self):
        buf = StringIO()

        df = DataFrame({'tups': zip(range(10), range(10))})
        repr(df)
        df.to_string(col_space=10, buf=buf)

    def test_repr_truncation(self):
        max_len = 20
        with option_context("display.max_colwidth", max_len):
            df = DataFrame({'A': np.random.randn(10),
                            'B': [tm.rands(np.random.randint(max_len - 1,
                                                             max_len + 1)) for i in range(10)]})
            r = repr(df)
            r = r[r.find('\n') + 1:]

            _strlen = fmt._strlen_func()

            for line, value in zip(r.split('\n'), df['B']):
                if _strlen(value) + 1 > max_len:
                    self.assert_('...' in line)
                else:
                    self.assert_('...' not in line)

        with option_context("display.max_colwidth", 999999):
            self.assert_('...' not in repr(df))

        with option_context("display.max_colwidth", max_len + 2):
            self.assert_('...' not in repr(df))

    def test_repr_chop_threshold(self):
        df = DataFrame([[0.1, 0.5],[0.5, -0.1]])
        pd.reset_option("display.chop_threshold") # default None
        self.assertEqual(repr(df), '     0    1\n0  0.1  0.5\n1  0.5 -0.1')

        with option_context("display.chop_threshold", 0.2 ):
            self.assertEqual(repr(df), '     0    1\n0  0.0  0.5\n1  0.5  0.0')

        with option_context("display.chop_threshold", 0.6 ):
            self.assertEqual(repr(df), '   0  1\n0  0  0\n1  0  0')

        with option_context("display.chop_threshold", None ):
            self.assertEqual(repr(df),  '     0    1\n0  0.1  0.5\n1  0.5 -0.1')

    def test_repr_obeys_max_seq_limit(self):
        import pandas.core.common as com

        #unlimited
        reset_option("display.max_seq_items")
        self.assertTrue(len(com.pprint_thing(range(1000)))> 2000)

        with option_context("display.max_seq_items",5):
            self.assertTrue(len(com.pprint_thing(range(1000)))< 100)

    def test_repr_is_valid_construction_code(self):
        import pandas as pd

        # for the case of Index, where the repr is traditional rather then stylized
        idx = pd.Index(['a','b'])
        res = eval("pd."+repr(idx))
        tm.assert_series_equal(Series(res),Series(idx))

    def test_repr_should_return_str(self):
        # http://docs.python.org/py3k/reference/datamodel.html#object.__repr__
        # http://docs.python.org/reference/datamodel.html#object.__repr__
        # "...The return value must be a string object."

        # (str on py2.x, str (unicode) on py3)


        data = [8, 5, 3, 5]
        index1 = [u"\u03c3", u"\u03c4", u"\u03c5", u"\u03c6"]
        cols = [u"\u03c8"]
        df = DataFrame(data, columns=cols, index=index1)
        self.assertTrue(type(df.__repr__() == str))  # both py2 / 3

    def test_repr_no_backslash(self):
        with option_context('mode.sim_interactive', True):
            df = DataFrame(np.random.randn(10, 4))
            self.assertTrue('\\' not in repr(df))

    def test_expand_frame_repr(self):
        df_small = DataFrame('hello', [0], [0])
        df_wide = DataFrame('hello', [0], range(10))
        df_tall = DataFrame('hello', range(30), range(5))

        with option_context('mode.sim_interactive', True):
            with option_context('display.max_columns', 10,
                                'display.width',20,
                                'display.max_rows', 20):
                with option_context('display.expand_frame_repr', True):
                    self.assertFalse(has_info_repr(df_small))
                    self.assertFalse(has_expanded_repr(df_small))
                    self.assertFalse(has_info_repr(df_wide))
                    self.assertTrue(has_expanded_repr(df_wide))
                    self.assertTrue(has_info_repr(df_tall))
                    self.assertFalse(has_expanded_repr(df_tall))

                with option_context('display.expand_frame_repr', False):
                    self.assertFalse(has_info_repr(df_small))
                    self.assertFalse(has_expanded_repr(df_small))
                    self.assertTrue(has_info_repr(df_wide))
                    self.assertFalse(has_expanded_repr(df_wide))
                    self.assertTrue(has_info_repr(df_tall))
                    self.assertFalse(has_expanded_repr(df_tall))

    def test_repr_non_interactive(self):
        # in non interactive mode, there can be no dependency on the
        # result of terminal auto size detection
        df = DataFrame('hello', range(1000), range(5))

        with option_context('mode.sim_interactive', False,
                            'display.width', 0,
                            'display.height', 0,
                            'display.max_rows',5000):
            self.assertFalse(has_info_repr(df))
            self.assertFalse(has_expanded_repr(df))

    def test_repr_max_columns_max_rows(self):
        term_width, term_height = get_terminal_size()
        if term_width < 10 or term_height < 10:
            raise nose.SkipTest

        def mkframe(n):
            index = ['%05d' % i for i in range(n)]
            return DataFrame(0, index, index)

        df6 = mkframe(6)
        df10 = mkframe(10)
        with option_context('mode.sim_interactive', True):
            with option_context('display.width', term_width * 2):
                with option_context('display.max_rows', 5,
                                    'display.max_columns', 5):
                    self.assertFalse(has_expanded_repr(mkframe(4)))
                    self.assertFalse(has_expanded_repr(mkframe(5)))
                    self.assertFalse(has_expanded_repr(df6))
                    self.assertTrue(has_info_repr(df6))

                with option_context('display.max_rows', 20,
                                    'display.max_columns', 10):
                    # Out off max_columns boundary, but no extending
                    # since not exceeding width
                    self.assertFalse(has_expanded_repr(df6))
                    self.assertFalse(has_info_repr(df6))

                with option_context('display.max_rows', 9,
                                    'display.max_columns', 10):
                    # out vertical bounds can not result in exanded repr
                    self.assertFalse(has_expanded_repr(df10))
                    self.assertTrue(has_info_repr(df10))

            # width=None in terminal, auto detection
            with option_context('display.max_columns', 100,
                                'display.max_rows', term_width * 20,
                                'display.width', None):
                df = mkframe((term_width // 7) - 2)
                self.assertFalse(has_expanded_repr(df))
                df = mkframe((term_width // 7) + 2)
                print( df._repr_fits_horizontal_())
                self.assertTrue(has_expanded_repr(df))

    def test_to_string_repr_unicode(self):
        buf = StringIO()

        unicode_values = [u'\u03c3'] * 10
        unicode_values = np.array(unicode_values, dtype=object)
        df = DataFrame({'unicode': unicode_values})
        df.to_string(col_space=10, buf=buf)

        # it works!
        repr(df)

        idx = Index(['abc', u'\u03c3a', 'aegdvg'])
        ser = Series(np.random.randn(len(idx)), idx)
        rs = repr(ser).split('\n')
        line_len = len(rs[0])
        for line in rs[1:]:
            try:
                line = line.decode(get_option("display.encoding"))
            except:
                pass
            if not line.startswith('dtype:'):
                self.assert_(len(line) == line_len)

        # it works even if sys.stdin in None
        _stdin= sys.stdin
        try:
            sys.stdin = None
            repr(df)
        finally:
            sys.stdin = _stdin

    def test_to_string_unicode_columns(self):
        df = DataFrame({u'\u03c3': np.arange(10.)})

        buf = StringIO()
        df.to_string(buf=buf)
        buf.getvalue()

        buf = StringIO()
        df.info(buf=buf)
        buf.getvalue()

        result = self.frame.to_string()
        self.assert_(isinstance(result, unicode))

    def test_to_string_utf8_columns(self):
        n = u"\u05d0".encode('utf-8')

        with option_context('display.max_rows', 1):
            df = pd.DataFrame([1, 2], columns=[n])
            repr(df)

    def test_to_string_unicode_two(self):
        dm = DataFrame({u'c/\u03c3': []})
        buf = StringIO()
        dm.to_string(buf)

    def test_to_string_unicode_three(self):
        dm = DataFrame(['\xc2'])
        buf = StringIO()
        dm.to_string(buf)

    def test_to_string_with_formatters(self):
        df = DataFrame({'int': [1, 2, 3],
                        'float': [1.0, 2.0, 3.0],
                        'object': [(1, 2), True, False]},
                       columns=['int', 'float', 'object'])

        formatters = [('int', lambda x: '0x%x' % x),
                      ('float', lambda x: '[% 4.1f]' % x),
                      ('object', lambda x: '-%s-' % str(x))]
        result = df.to_string(formatters=dict(formatters))
        result2 = df.to_string(formatters=lzip(*formatters)[1])
        self.assertEqual(result, ('  int  float    object\n'
                                  '0 0x1 [ 1.0]  -(1, 2)-\n'
                                  '1 0x2 [ 2.0]    -True-\n'
                                  '2 0x3 [ 3.0]   -False-'))
        self.assertEqual(result, result2)

    def test_to_string_with_formatters_unicode(self):
        df = DataFrame({u'c/\u03c3': [1, 2, 3]})
        result = df.to_string(formatters={u'c/\u03c3': lambda x: '%s' % x})
        self.assertEqual(result, (u'  c/\u03c3\n'
                                  '0   1\n'
                                  '1   2\n'
                                  '2   3'))

    def test_to_string_buffer_all_unicode(self):
        buf = StringIO()

        empty = DataFrame({u'c/\u03c3': Series()})
        nonempty = DataFrame({u'c/\u03c3': Series([1, 2, 3])})

        print >>buf, empty
        print >>buf, nonempty

        # this should work
        buf.getvalue()

    def test_to_string_with_col_space(self):
        df = DataFrame(np.random.random(size=(1, 3)))
        c10 = len(df.to_string(col_space=10).split("\n")[1])
        c20 = len(df.to_string(col_space=20).split("\n")[1])
        c30 = len(df.to_string(col_space=30).split("\n")[1])
        self.assertTrue(c10 < c20 < c30)

    def test_to_html_with_col_space(self):
        def check_with_width(df, col_space):
            import re
            # check that col_space affects HTML generation
            # and be very brittle about it.
            html = df.to_html(col_space=col_space)
            hdrs = [x for x in html.split("\n") if re.search("<th[>\s]", x)]
            self.assertTrue(len(hdrs) > 0)
            for h in hdrs:
                self.assertTrue("min-width" in h)
                self.assertTrue(str(col_space) in h)

        df = DataFrame(np.random.random(size=(1, 3)))

        check_with_width(df, 30)
        check_with_width(df, 50)

    def test_to_html_with_empty_string_label(self):
        # GH3547, to_html regards empty string labels as repeated labels
        data = {'c1': ['a', 'b'], 'c2': ['a', ''], 'data': [1, 2]}
        df = DataFrame(data).set_index(['c1', 'c2'])
        res = df.to_html()
        self.assertTrue("rowspan" not in res)

    def test_to_html_unicode(self):
        # it works!
        df = DataFrame({u'\u03c3': np.arange(10.)})
        df.to_html()
        df = DataFrame({'A': [u'\u03c3']})
        df.to_html()

    def test_to_html_escaped(self):
        a = 'str<ing1 &amp;'
        b = 'stri>ng2 &amp;'

        test_dict = {'co<l1': {a: "<type 'str'>",
                               b: "<type 'str'>"},
                     'co>l2':{a: "<type 'str'>",
                              b: "<type 'str'>"}}
        rs = pd.DataFrame(test_dict).to_html()
        xp = """<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>co&lt;l1</th>
      <th>co&gt;l2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>str&lt;ing1 &amp;amp;</th>
      <td> &lt;type 'str'&gt;</td>
      <td> &lt;type 'str'&gt;</td>
    </tr>
    <tr>
      <th>stri&gt;ng2 &amp;amp;</th>
      <td> &lt;type 'str'&gt;</td>
      <td> &lt;type 'str'&gt;</td>
    </tr>
  </tbody>
</table>"""
        self.assertEqual(xp, rs)

    def test_to_html_escape_disabled(self):
        a = 'str<ing1 &amp;'
        b = 'stri>ng2 &amp;'

        test_dict = {'co<l1': {a: "<b>bold</b>",
                               b: "<b>bold</b>"},
                     'co>l2': {a: "<b>bold</b>",
                               b: "<b>bold</b>"}}
        rs = pd.DataFrame(test_dict).to_html(escape=False)
        xp = """<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>co<l1</th>
      <th>co>l2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>str<ing1 &amp;</th>
      <td> <b>bold</b></td>
      <td> <b>bold</b></td>
    </tr>
    <tr>
      <th>stri>ng2 &amp;</th>
      <td> <b>bold</b></td>
      <td> <b>bold</b></td>
    </tr>
  </tbody>
</table>"""
        self.assertEqual(xp, rs)

    def test_to_html_multiindex_sparsify_false_multi_sparse(self):
        with option_context('display.multi_sparse', False):
            index = pd.MultiIndex.from_arrays([[0, 0, 1, 1], [0, 1, 0, 1]],
                                              names=['foo', None])

            df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]], index=index)

            result = df.to_html()
            expected = """\
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>0</th>
      <th>1</th>
    </tr>
    <tr>
      <th>foo</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <th>0</th>
      <td> 0</td>
      <td> 1</td>
    </tr>
    <tr>
      <th>0</th>
      <th>1</th>
      <td> 2</td>
      <td> 3</td>
    </tr>
    <tr>
      <th>1</th>
      <th>0</th>
      <td> 4</td>
      <td> 5</td>
    </tr>
    <tr>
      <th>1</th>
      <th>1</th>
      <td> 6</td>
      <td> 7</td>
    </tr>
  </tbody>
</table>"""
            self.assertEquals(result, expected)

            df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]],
                           columns=index[::2], index=index)

            result = df.to_html()
            expected = """\
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th></th>
      <th>foo</th>
      <th>0</th>
      <th>1</th>
    </tr>
    <tr>
      <th></th>
      <th></th>
      <th>0</th>
      <th>0</th>
    </tr>
    <tr>
      <th>foo</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <th>0</th>
      <td> 0</td>
      <td> 1</td>
    </tr>
    <tr>
      <th>0</th>
      <th>1</th>
      <td> 2</td>
      <td> 3</td>
    </tr>
    <tr>
      <th>1</th>
      <th>0</th>
      <td> 4</td>
      <td> 5</td>
    </tr>
    <tr>
      <th>1</th>
      <th>1</th>
      <td> 6</td>
      <td> 7</td>
    </tr>
  </tbody>
</table>"""
            self.assertEquals(result, expected)

    def test_to_html_multiindex_sparsify(self):
        index = pd.MultiIndex.from_arrays([[0, 0, 1, 1], [0, 1, 0, 1]],
                                          names=['foo', None])

        df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]], index=index)

        result = df.to_html()
        expected = """<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>0</th>
      <th>1</th>
    </tr>
    <tr>
      <th>foo</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">0</th>
      <th>0</th>
      <td> 0</td>
      <td> 1</td>
    </tr>
    <tr>
      <th>1</th>
      <td> 2</td>
      <td> 3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">1</th>
      <th>0</th>
      <td> 4</td>
      <td> 5</td>
    </tr>
    <tr>
      <th>1</th>
      <td> 6</td>
      <td> 7</td>
    </tr>
  </tbody>
</table>"""
        self.assertEquals(result, expected)

        df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]],
                       columns=index[::2], index=index)

        result = df.to_html()
        expected = """\
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th></th>
      <th>foo</th>
      <th>0</th>
      <th>1</th>
    </tr>
    <tr>
      <th></th>
      <th></th>
      <th>0</th>
      <th>0</th>
    </tr>
    <tr>
      <th>foo</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">0</th>
      <th>0</th>
      <td> 0</td>
      <td> 1</td>
    </tr>
    <tr>
      <th>1</th>
      <td> 2</td>
      <td> 3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">1</th>
      <th>0</th>
      <td> 4</td>
      <td> 5</td>
    </tr>
    <tr>
      <th>1</th>
      <td> 6</td>
      <td> 7</td>
    </tr>
  </tbody>
</table>"""
        self.assertEquals(result, expected)

    def test_to_html_index_formatter(self):
        df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]],
                       columns=['foo', None], index=range(4))

        f = lambda x: 'abcd'[x]
        result = df.to_html(formatters={'__index__': f})
        expected = """\
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>foo</th>
      <th>None</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>a</th>
      <td> 0</td>
      <td> 1</td>
    </tr>
    <tr>
      <th>b</th>
      <td> 2</td>
      <td> 3</td>
    </tr>
    <tr>
      <th>c</th>
      <td> 4</td>
      <td> 5</td>
    </tr>
    <tr>
      <th>d</th>
      <td> 6</td>
      <td> 7</td>
    </tr>
  </tbody>
</table>"""
        self.assertEquals(result, expected)

    def test_nonunicode_nonascii_alignment(self):
        df = DataFrame([["aa\xc3\xa4\xc3\xa4", 1], ["bbbb", 2]])
        rep_str = df.to_string()
        lines = rep_str.split('\n')
        self.assert_(len(lines[1]) == len(lines[2]))

    def test_unicode_problem_decoding_as_ascii(self):
        dm = DataFrame({u'c/\u03c3': Series({'test': np.NaN})})
        unicode(dm.to_string())

    def test_string_repr_encoding(self):
        filepath = tm.get_data_path('unicode_series.csv')
        df = pandas.read_csv(filepath, header=None, encoding='latin1')
        repr(df)
        repr(df[1])

    def test_repr_corner(self):
        # representing infs poses no problems
        df = DataFrame({'foo': np.inf * np.empty(10)})
        foo = repr(df)

    def test_frame_info_encoding(self):
        index = ['\'Til There Was You (1997)',
                 'ldum klaka (Cold Fever) (1994)']
        fmt.set_option('display.max_rows', 1)
        df = DataFrame(columns=['a', 'b', 'c'], index=index)
        repr(df)
        repr(df.T)
        fmt.set_option('display.max_rows', 200)

    def test_large_frame_repr(self):
        def wrap_rows_options(f):
            def _f(*args, **kwargs):
                old_max_rows = pd.get_option('display.max_rows')
                old_max_info_rows = pd.get_option('display.max_info_rows')
                o = f(*args, **kwargs)
                pd.set_option('display.max_rows', old_max_rows)
                pd.set_option('display.max_info_rows', old_max_info_rows)
                return o
            return _f

        @wrap_rows_options
        def test_setting(value, nrows=3, ncols=2):
            if value is None:
                expected_difference = 0
            elif isinstance(value, int):
                expected_difference = ncols
            else:
                raise ValueError("'value' must be int or None")

            with option_context('mode.sim_interactive', True):
                pd.set_option('display.max_rows', nrows - 1)
                pd.set_option('display.max_info_rows', value)

                smallx = DataFrame(np.random.rand(nrows, ncols))
                repr_small = repr(smallx)

                bigx = DataFrame(np.random.rand(nrows + 1, ncols))
                repr_big = repr(bigx)

                diff = len(repr_small.splitlines()) - len(repr_big.splitlines())

                # the difference in line count is the number of columns
                self.assertEqual(diff, expected_difference)

        test_setting(None)
        test_setting(3)
        self.assertRaises(ValueError, test_setting, 'string')

    def test_pprint_thing(self):
        import nose
        from pandas.core.common import pprint_thing as pp_t

        if PY3:
            raise nose.SkipTest()

        self.assertEquals(pp_t('a') , u'a')
        self.assertEquals(pp_t(u'a') , u'a')
        self.assertEquals(pp_t(None) , 'None')
        self.assertEquals(pp_t(u'\u05d0',quote_strings=True) , u"u'\u05d0'")
        self.assertEquals(pp_t(u'\u05d0',quote_strings=False) , u'\u05d0')
        self.assertEquals(pp_t((u'\u05d0', u'\u05d1'),quote_strings=True) ,
                          u"(u'\u05d0', u'\u05d1')")
        self.assertEquals(pp_t((u'\u05d0', (u'\u05d1', u'\u05d2')),quote_strings=True) ,
               u"(u'\u05d0', (u'\u05d1', u'\u05d2'))")
        self.assertEquals(pp_t(('foo', u'\u05d0', (u'\u05d0', u'\u05d0')),quote_strings=True)
                          , u"(u'foo', u'\u05d0', (u'\u05d0', u'\u05d0'))")

        # escape embedded tabs in string
        # GH #2038
        self.assertTrue(not "\t" in pp_t("a\tb", escape_chars=("\t",)))

    def test_wide_repr(self):
        with option_context('mode.sim_interactive', True):
            col = lambda l, k: [tm.rands(k) for _ in xrange(l)]
            max_cols = get_option('display.max_columns')
            df = DataFrame([col(max_cols-1, 25) for _ in range(10)])
            set_option('display.expand_frame_repr', False)
            rep_str = repr(df)
            set_option('display.expand_frame_repr', True)
            wide_repr = repr(df)
            self.assert_(rep_str != wide_repr)

            with option_context('display.width', 120):
                wider_repr = repr(df)
                self.assert_(len(wider_repr) < len(wide_repr))

        reset_option('display.expand_frame_repr')

    def test_wide_repr_wide_columns(self):
        with option_context('mode.sim_interactive', True):
            df = DataFrame(randn(5, 3), columns=['a' * 90, 'b' * 90, 'c' * 90])
            rep_str = repr(df)

            self.assert_(len(rep_str.splitlines()) == 20)

    def test_wide_repr_named(self):
        with option_context('mode.sim_interactive', True):
            col = lambda l, k: [tm.rands(k) for _ in xrange(l)]
            max_cols = get_option('display.max_columns')
            df = DataFrame([col(max_cols-1, 25) for _ in range(10)])
            df.index.name = 'DataFrame Index'
            set_option('display.expand_frame_repr', False)

            rep_str = repr(df)
            set_option('display.expand_frame_repr', True)
            wide_repr = repr(df)
            self.assert_(rep_str != wide_repr)

            with option_context('display.width', 150):
                wider_repr = repr(df)
                self.assert_(len(wider_repr) < len(wide_repr))

            for line in wide_repr.splitlines()[1::13]:
                self.assert_('DataFrame Index' in line)

        reset_option('display.expand_frame_repr')

    def test_wide_repr_multiindex(self):
        with option_context('mode.sim_interactive', True):
            col = lambda l, k: [tm.rands(k) for _ in xrange(l)]
            midx = pandas.MultiIndex.from_arrays([np.array(col(10, 5)),
                                                  np.array(col(10, 5))])
            max_cols = get_option('display.max_columns')
            df = DataFrame([col(max_cols-1, 25) for _ in range(10)],
                           index=midx)
            df.index.names = ['Level 0', 'Level 1']
            set_option('display.expand_frame_repr', False)
            rep_str = repr(df)
            set_option('display.expand_frame_repr', True)
            wide_repr = repr(df)
            self.assert_(rep_str != wide_repr)

            with option_context('display.width', 150):
                wider_repr = repr(df)
                self.assert_(len(wider_repr) < len(wide_repr))

            for line in wide_repr.splitlines()[1::13]:
                self.assert_('Level 0 Level 1' in line)

        reset_option('display.expand_frame_repr')

    def test_wide_repr_multiindex_cols(self):
        with option_context('mode.sim_interactive', True):
            max_cols = get_option('display.max_columns')
            col = lambda l, k: [tm.rands(k) for _ in xrange(l)]
            midx = pandas.MultiIndex.from_arrays([np.array(col(10, 5)),
                                                  np.array(col(10, 5))])
            mcols = pandas.MultiIndex.from_arrays([np.array(col(max_cols-1, 3)),
                                                   np.array(col(max_cols-1, 3))])
            df = DataFrame([col(max_cols-1, 25) for _ in range(10)],
                           index=midx, columns=mcols)
            df.index.names = ['Level 0', 'Level 1']
            set_option('display.expand_frame_repr', False)
            rep_str = repr(df)
            set_option('display.expand_frame_repr', True)
            wide_repr = repr(df)
            self.assert_(rep_str != wide_repr)

        with option_context('display.width', 150):
            wider_repr = repr(df)
            self.assert_(len(wider_repr) < len(wide_repr))

        reset_option('display.expand_frame_repr')

    def test_wide_repr_unicode(self):
        with option_context('mode.sim_interactive', True):
            col = lambda l, k: [tm.randu(k) for _ in xrange(l)]
            max_cols = get_option('display.max_columns')
            df = DataFrame([col(max_cols-1, 25) for _ in range(10)])
            set_option('display.expand_frame_repr', False)
            rep_str = repr(df)
            set_option('display.expand_frame_repr', True)
            wide_repr = repr(df)
            self.assert_(rep_str != wide_repr)

            with option_context('display.width', 150):
                wider_repr = repr(df)
                self.assert_(len(wider_repr) < len(wide_repr))

        reset_option('display.expand_frame_repr')

    def test_wide_repr_wide_long_columns(self):
        with option_context('mode.sim_interactive', True):
            df = DataFrame(
                {'a': ['a' * 30, 'b' * 30], 'b': ['c' * 70, 'd' * 80]})

            result = repr(df)
            self.assertTrue('ccccc' in result)
            self.assertTrue('ddddd' in result)

    def test_long_series(self):
        n = 1000
        s = Series(np.random.randint(-50,50,n),index=['s%04d' % x for x in xrange(n)], dtype='int64')

        import re
        str_rep = str(s)
        nmatches = len(re.findall('dtype',str_rep))
        self.assert_(nmatches == 1)

    def test_index_with_nan(self):
        #  GH 2850
        df = DataFrame({'id1': {0: '1a3', 1: '9h4'}, 'id2': {0: np.nan, 1: 'd67'},
                        'id3': {0: '78d', 1: '79d'}, 'value': {0: 123, 1: 64}})

        # multi-index
        y = df.set_index(['id1', 'id2', 'id3'])
        result = y.to_string()
        expected = u'             value\nid1 id2 id3       \n1a3 NaN 78d    123\n9h4 d67 79d     64'
        self.assert_(result == expected)

        # index
        y = df.set_index('id2')
        result = y.to_string()
        expected = u'     id1  id3  value\nid2                 \nNaN  1a3  78d    123\nd67  9h4  79d     64'
        self.assert_(result == expected)

        # all-nan in mi
        df2 = df.copy()
        df2.ix[:,'id2'] = np.nan
        y = df2.set_index('id2')
        result = y.to_string()
        expected = u'     id1  id3  value\nid2                 \nNaN  1a3  78d    123\nNaN  9h4  79d     64'
        self.assert_(result == expected)

        # partial nan in mi
        df2 = df.copy()
        df2.ix[:,'id2'] = np.nan
        y = df2.set_index(['id2','id3'])
        result = y.to_string()
        expected = u'         id1  value\nid2 id3            \nNaN 78d  1a3    123\n    79d  9h4     64'
        self.assert_(result == expected)

        df = DataFrame({'id1': {0: np.nan, 1: '9h4'}, 'id2': {0: np.nan, 1: 'd67'},
                        'id3': {0: np.nan, 1: '79d'}, 'value': {0: 123, 1: 64}})

        y = df.set_index(['id1','id2','id3'])
        result = y.to_string()
        expected = u'             value\nid1 id2 id3       \nNaN NaN NaN    123\n9h4 d67 79d     64'
        self.assert_(result == expected)

    def test_to_string(self):
        from pandas import read_table
        import re

        # big mixed
        biggie = DataFrame({'A': randn(200),
                            'B': tm.makeStringIndex(200)},
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
        recons = read_table(StringIO(joined), names=header,
                            header=None, sep=' ')
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
                         formatters={'A': lambda x: '%.1f' % x})

        biggie.to_string(columns=['B', 'A'], float_format=str)
        biggie.to_string(columns=['B', 'A'], col_space=12,
                         float_format=str)

        frame = DataFrame(index=np.arange(200))
        frame.to_string()

    def test_to_string_no_header(self):
        df = DataFrame({'x': [1, 2, 3],
                        'y': [4, 5, 6]})

        df_s = df.to_string(header=False)
        expected = "0  1  4\n1  2  5\n2  3  6"

        assert(df_s == expected)

    def test_to_string_no_index(self):
        df = DataFrame({'x': [1, 2, 3],
                        'y': [4, 5, 6]})

        df_s = df.to_string(index=False)
        expected = " x  y\n 1  4\n 2  5\n 3  6"

        assert(df_s == expected)

    def test_to_string_float_formatting(self):
        fmt.reset_option('^display.')
        fmt.set_option('display.precision', 6, 'display.column_space',
                       12, 'display.notebook_repr_html', False)

        df = DataFrame({'x': [0, 0.25, 3456.000, 12e+45, 1.64e+6,
                              1.7e+8, 1.253456, np.pi, -1e6]})

        df_s = df.to_string()

        # Python 2.5 just wants me to be sad. And debian 32-bit
        # sys.version_info[0] == 2 and sys.version_info[1] < 6:
        if _three_digit_exp():
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

        df = DataFrame({'x': [3234, 0.253]})
        df_s = df.to_string()

        expected = ('          x\n'
                    '0  3234.000\n'
                    '1     0.253')
        assert(df_s == expected)

        fmt.reset_option('^display.')
        self.assertEqual(get_option("display.precision"), 7)

        df = DataFrame({'x': [1e9, 0.2512]})
        df_s = df.to_string()
        # Python 2.5 just wants me to be sad. And debian 32-bit
        # sys.version_info[0] == 2 and sys.version_info[1] < 6:
        if _three_digit_exp():
            expected = ('               x\n'
                        '0  1.000000e+009\n'
                        '1  2.512000e-001')
        else:
            expected = ('              x\n'
                        '0  1.000000e+09\n'
                        '1  2.512000e-01')
        assert(df_s == expected)

    def test_to_string_small_float_values(self):
        df = DataFrame({'a': [1.5, 1e-17, -5.5e-7]})

        result = df.to_string()
        # sadness per above
        if '%.4g' % 1.7e8 == '1.7e+008':
            expected = ('               a\n'
                        '0  1.500000e+000\n'
                        '1  1.000000e-017\n'
                        '2 -5.500000e-007')
        else:
            expected = ('              a\n'
                        '0  1.500000e+00\n'
                        '1  1.000000e-17\n'
                        '2 -5.500000e-07')
        self.assertEqual(result, expected)

        # but not all exactly zero
        df = df * 0
        result = df.to_string()
        expected = ('   0\n'
                    '0  0\n'
                    '1  0\n'
                    '2 -0')

    def test_to_string_float_index(self):
        index = Index([1.5, 2, 3, 4, 5])
        df = DataFrame(range(5), index=index)

        result = df.to_string()
        expected = ('     0\n'
                    '1.5  0\n'
                    '2    1\n'
                    '3    2\n'
                    '4    3\n'
                    '5    4')
        self.assertEqual(result, expected)

    def test_to_string_ascii_error(self):
        data = [('0  ',
                 u'                        .gitignore ',
                 u'     5 ',
                 ' \xe2\x80\xa2\xe2\x80\xa2\xe2\x80'
                 '\xa2\xe2\x80\xa2\xe2\x80\xa2')]
        df = DataFrame(data)

        # it works!
        repr(df)

    def test_to_string_int_formatting(self):
        df = DataFrame({'x': [-15, 20, 25, -35]})
        self.assert_(issubclass(df['x'].dtype.type, np.integer))

        output = df.to_string()
        expected = ('    x\n'
                    '0 -15\n'
                    '1  20\n'
                    '2  25\n'
                    '3 -35')
        self.assertEqual(output, expected)

    def test_to_string_index_formatter(self):
        df = DataFrame([range(5), range(5, 10), range(10, 15)])

        rs = df.to_string(formatters={'__index__': lambda x: 'abc'[x]})

        xp = """\
    0   1   2   3   4
a   0   1   2   3   4
b   5   6   7   8   9
c  10  11  12  13  14\
"""
        self.assertEqual(rs, xp)

    def test_to_string_left_justify_cols(self):
        fmt.reset_option('^display.')
        df = DataFrame({'x': [3234, 0.253]})
        df_s = df.to_string(justify='left')
        expected = ('   x       \n'
                    '0  3234.000\n'
                    '1     0.253')
        assert(df_s == expected)

    def test_to_string_format_na(self):
        fmt.reset_option('^display.')
        df = DataFrame({'A': [np.nan, -1, -2.1234, 3, 4],
                        'B': [np.nan, 'foo', 'foooo', 'fooooo', 'bar']})
        result = df.to_string()

        expected = ('        A       B\n'
                    '0     NaN     NaN\n'
                    '1 -1.0000     foo\n'
                    '2 -2.1234   foooo\n'
                    '3  3.0000  fooooo\n'
                    '4  4.0000     bar')
        self.assertEqual(result, expected)

        df = DataFrame({'A': [np.nan, -1., -2., 3., 4.],
                        'B': [np.nan, 'foo', 'foooo', 'fooooo', 'bar']})
        result = df.to_string()

        expected = ('    A       B\n'
                    '0 NaN     NaN\n'
                    '1  -1     foo\n'
                    '2  -2   foooo\n'
                    '3   3  fooooo\n'
                    '4   4     bar')
        self.assertEqual(result, expected)

    def test_to_string_line_width(self):
        df = pd.DataFrame(123, range(10, 15), range(30))
        s = df.to_string(line_width=80)
        self.assertEqual(max(len(l) for l in s.split('\n')), 80)

    def test_to_html(self):
        # big mixed
        biggie = DataFrame({'A': randn(200),
                            'B': tm.makeStringIndex(200)},
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
                       formatters={'A': lambda x: '%.1f' % x})

        biggie.to_html(columns=['B', 'A'], float_format=str)
        biggie.to_html(columns=['B', 'A'], col_space=12,
                       float_format=str)

        frame = DataFrame(index=np.arange(200))
        frame.to_html()

    def test_to_html_filename(self):
        biggie = DataFrame({'A': randn(200),
                            'B': tm.makeStringIndex(200)},
                           index=range(200))

        biggie['A'][:20] = nan
        biggie['B'][:20] = nan
        with tm.ensure_clean('test.html') as path:
            biggie.to_html(path)
            with open(path, 'r') as f:
                s = biggie.to_html()
                s2 = f.read()
                self.assertEqual(s, s2)

        frame = DataFrame(index=np.arange(200))
        with tm.ensure_clean('test.html') as path:
            frame.to_html(path)
            with open(path, 'r') as f:
                self.assertEqual(frame.to_html(), f.read())

    def test_to_html_with_no_bold(self):
        x = DataFrame({'x': randn(5)})
        ashtml = x.to_html(bold_rows=False)
        assert('<strong>' not in ashtml[ashtml.find('</thead>')])

    def test_to_html_columns_arg(self):
        result = self.frame.to_html(columns=['A'])
        self.assert_('<th>B</th>' not in result)

    def test_to_html_multiindex(self):
        columns = pandas.MultiIndex.from_tuples(zip(np.arange(2).repeat(2),
                                                    np.mod(range(4), 2)),
                                                names=['CL0', 'CL1'])
        df = pandas.DataFrame([list('abcd'), list('efgh')], columns=columns)
        result = df.to_html(justify='left')
        expected = ('<table border="1" class="dataframe">\n'
                    '  <thead>\n'
                    '    <tr>\n'
                    '      <th>CL0</th>\n'
                    '      <th colspan="2" halign="left">0</th>\n'
                    '      <th colspan="2" halign="left">1</th>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>CL1</th>\n'
                    '      <th>0</th>\n'
                    '      <th>1</th>\n'
                    '      <th>0</th>\n'
                    '      <th>1</th>\n'
                    '    </tr>\n'
                    '  </thead>\n'
                    '  <tbody>\n'
                    '    <tr>\n'
                    '      <th>0</th>\n'
                    '      <td> a</td>\n'
                    '      <td> b</td>\n'
                    '      <td> c</td>\n'
                    '      <td> d</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>1</th>\n'
                    '      <td> e</td>\n'
                    '      <td> f</td>\n'
                    '      <td> g</td>\n'
                    '      <td> h</td>\n'
                    '    </tr>\n'
                    '  </tbody>\n'
                    '</table>')

        self.assertEqual(result, expected)

        columns = pandas.MultiIndex.from_tuples(zip(range(4),
                                                    np.mod(range(4), 2)))
        df = pandas.DataFrame([list('abcd'), list('efgh')], columns=columns)

        result = df.to_html(justify='right')
        expected = ('<table border="1" class="dataframe">\n'
                    '  <thead>\n'
                    '    <tr>\n'
                    '      <th></th>\n'
                    '      <th>0</th>\n'
                    '      <th>1</th>\n'
                    '      <th>2</th>\n'
                    '      <th>3</th>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th></th>\n'
                    '      <th>0</th>\n'
                    '      <th>1</th>\n'
                    '      <th>0</th>\n'
                    '      <th>1</th>\n'
                    '    </tr>\n'
                    '  </thead>\n'
                    '  <tbody>\n'
                    '    <tr>\n'
                    '      <th>0</th>\n'
                    '      <td> a</td>\n'
                    '      <td> b</td>\n'
                    '      <td> c</td>\n'
                    '      <td> d</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>1</th>\n'
                    '      <td> e</td>\n'
                    '      <td> f</td>\n'
                    '      <td> g</td>\n'
                    '      <td> h</td>\n'
                    '    </tr>\n'
                    '  </tbody>\n'
                    '</table>')

        self.assertEqual(result, expected)

    def test_to_html_justify(self):
        df = pandas.DataFrame({'A': [6, 30000, 2],
                               'B': [1, 2, 70000],
                               'C': [223442, 0, 1]},
                              columns=['A', 'B', 'C'])
        result = df.to_html(justify='left')
        expected = ('<table border="1" class="dataframe">\n'
                    '  <thead>\n'
                    '    <tr style="text-align: left;">\n'
                    '      <th></th>\n'
                    '      <th>A</th>\n'
                    '      <th>B</th>\n'
                    '      <th>C</th>\n'
                    '    </tr>\n'
                    '  </thead>\n'
                    '  <tbody>\n'
                    '    <tr>\n'
                    '      <th>0</th>\n'
                    '      <td>     6</td>\n'
                    '      <td>     1</td>\n'
                    '      <td> 223442</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>1</th>\n'
                    '      <td> 30000</td>\n'
                    '      <td>     2</td>\n'
                    '      <td>      0</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>2</th>\n'
                    '      <td>     2</td>\n'
                    '      <td> 70000</td>\n'
                    '      <td>      1</td>\n'
                    '    </tr>\n'
                    '  </tbody>\n'
                    '</table>')

        self.assertEqual(result, expected)

        result = df.to_html(justify='right')
        expected = ('<table border="1" class="dataframe">\n'
                    '  <thead>\n'
                    '    <tr style="text-align: right;">\n'
                    '      <th></th>\n'
                    '      <th>A</th>\n'
                    '      <th>B</th>\n'
                    '      <th>C</th>\n'
                    '    </tr>\n'
                    '  </thead>\n'
                    '  <tbody>\n'
                    '    <tr>\n'
                    '      <th>0</th>\n'
                    '      <td>     6</td>\n'
                    '      <td>     1</td>\n'
                    '      <td> 223442</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>1</th>\n'
                    '      <td> 30000</td>\n'
                    '      <td>     2</td>\n'
                    '      <td>      0</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>2</th>\n'
                    '      <td>     2</td>\n'
                    '      <td> 70000</td>\n'
                    '      <td>      1</td>\n'
                    '    </tr>\n'
                    '  </tbody>\n'
                    '</table>')
        self.assertEqual(result, expected)

    def test_to_html_index(self):
        index = ['foo', 'bar', 'baz']
        df = pandas.DataFrame({'A': [1, 2, 3],
                               'B': [1.2, 3.4, 5.6],
                               'C': ['one', 'two', np.NaN]},
                              columns=['A', 'B', 'C'],
                              index=index)
        result = df.to_html(index=False)
        for i in index:
            self.assert_(i not in result)

        tuples = [('foo', 'car'), ('foo', 'bike'), ('bar', 'car')]
        df.index = pandas.MultiIndex.from_tuples(tuples)
        result = df.to_html(index=False)
        for i in ['foo', 'bar', 'car', 'bike']:
            self.assert_(i not in result)

    def test_repr_html(self):
        self.frame._repr_html_()

        fmt.set_option('display.max_rows', 1, 'display.max_columns', 1)
        self.frame._repr_html_()

        fmt.set_option('display.notebook_repr_html', False)
        self.frame._repr_html_()

        fmt.reset_option('^display.')

    def test_fake_qtconsole_repr_html(self):
        def get_ipython():
            return {'config':
                   {'KernelApp':
                   {'parent_appname': 'ipython-qtconsole'}}}

        repstr = self.frame._repr_html_()
        self.assert_(repstr is not None)

        fmt.set_option('display.max_rows', 5, 'display.max_columns', 2)
        repstr = self.frame._repr_html_()
        self.assert_('class' in repstr)  # info fallback

        fmt.reset_option('^display.')

    def test_to_html_with_classes(self):
        df = pandas.DataFrame()
        result = df.to_html(classes="sortable draggable")
        expected = dedent("""

            <table border="1" class="dataframe sortable draggable">
              <tbody>
                <tr>
                  <td>Index([], dtype=object)</td>
                  <td>Empty DataFrame</td>
                </tr>
              </tbody>
            </table>

        """).strip()
        self.assertEqual(result, expected)

        result = df.to_html(classes=["sortable", "draggable"])
        self.assertEqual(result, expected)

    def test_pprint_pathological_object(self):
        """
        if the test fails, the stack will overflow and nose crash,
        but it won't hang.
        """
        class A:
            def __getitem__(self, key):
                return 3 # obviously simplified
        df = pandas.DataFrame([A()])
        repr(df) # just don't dine

    def test_float_trim_zeros(self):
        vals = [2.08430917305e+10, 3.52205017305e+10, 2.30674817305e+10,
                2.03954217305e+10, 5.59897817305e+10]
        skip = True
        for line in repr(DataFrame({'A': vals})).split('\n'):
            if line.startswith('dtype:'):
                continue
            if _three_digit_exp():
                self.assert_(('+010' in line) or skip)
            else:
                self.assert_(('+10' in line) or skip)
            skip = False

    def test_dict_entries(self):
        df = DataFrame({'A': [{'a': 1, 'b': 2}]})

        val = df.to_string()
        self.assertTrue("'a': 1" in val)
        self.assertTrue("'b': 2" in val)

    def test_to_latex_filename(self):
        with tm.ensure_clean('test.tex') as path:
            self.frame.to_latex(path)

            with open(path, 'r') as f:
                self.assertEqual(self.frame.to_latex(), f.read())

    def test_to_latex(self):
        # it works!
        self.frame.to_latex()

        df = DataFrame({'a': [1, 2],
                        'b': ['b1', 'b2']})
        withindex_result = df.to_latex()
        withindex_expected = r"""\begin{tabular}{lrl}
\toprule
{} &  a &   b \\
\midrule
0 &  1 &  b1 \\
1 &  2 &  b2 \\
\bottomrule
\end{tabular}
"""
        self.assertEqual(withindex_result, withindex_expected)

        withoutindex_result = df.to_latex(index=False)
        withoutindex_expected = r"""\begin{tabular}{rl}
\toprule
 a &   b \\
\midrule
 1 &  b1 \\
 2 &  b2 \\
\bottomrule
\end{tabular}
"""
        self.assertEqual(withoutindex_result, withoutindex_expected)

class TestSeriesFormatting(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.ts = tm.makeTimeSeries()

    def test_repr_unicode(self):
        s = Series([u'\u03c3'] * 10)
        repr(s)

        a = Series([u"\u05d0"] * 1000)
        a.name = 'title1'
        repr(a)

    def test_to_string(self):
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
        result = cp.to_string(length=True, name=True, dtype=True)
        last_line = result.split('\n')[-1].strip()
        self.assertEqual(last_line, "Freq: B, Name: foo, Length: %d, dtype: float64" % len(cp))

    def test_freq_name_separation(self):
        s = Series(np.random.randn(10),
                   index=pd.date_range('1/1/2000', periods=10), name=0)

        result = repr(s)
        self.assertTrue('Freq: D, Name: 0' in result)

    def test_to_string_mixed(self):
        s = Series(['foo', np.nan, -1.23, 4.56])
        result = s.to_string()
        expected = (u'0     foo\n'
                    u'1     NaN\n'
                    u'2   -1.23\n'
                    u'3    4.56')
        self.assertEqual(result, expected)

        # but don't count NAs as floats
        s = Series(['foo', np.nan, 'bar', 'baz'])
        result = s.to_string()
        expected = (u'0    foo\n'
                    '1    NaN\n'
                    '2    bar\n'
                    '3    baz')
        self.assertEqual(result, expected)

        s = Series(['foo', 5, 'bar', 'baz'])
        result = s.to_string()
        expected = (u'0    foo\n'
                    '1      5\n'
                    '2    bar\n'
                    '3    baz')
        self.assertEqual(result, expected)

    def test_to_string_float_na_spacing(self):
        s = Series([0., 1.5678, 2., -3., 4.])
        s[::2] = np.nan

        result = s.to_string()
        expected = (u'0       NaN\n'
                    '1    1.5678\n'
                    '2       NaN\n'
                    '3   -3.0000\n'
                    '4       NaN')
        self.assertEqual(result, expected)

    def test_unicode_name_in_footer(self):
        s = Series([1, 2], name=u'\u05e2\u05d1\u05e8\u05d9\u05ea')
        sf = fmt.SeriesFormatter(s, name=u'\u05e2\u05d1\u05e8\u05d9\u05ea')
        sf._get_footer()  # should not raise exception

    def test_float_trim_zeros(self):
        vals = [2.08430917305e+10, 3.52205017305e+10, 2.30674817305e+10,
                2.03954217305e+10, 5.59897817305e+10]
        for line in repr(Series(vals)).split('\n'):
            if line.startswith('dtype:'):
                continue
            if _three_digit_exp():
                self.assert_('+010' in line)
            else:
                self.assert_('+10' in line)

    def test_datetimeindex(self):

        from pandas import date_range, NaT, Timestamp
        index = date_range('20130102',periods=6)
        s = Series(1,index=index)
        result = s.to_string()
        self.assertTrue('2013-01-02' in result)

        # nat in index
        s2 = Series(2, index=[ Timestamp('20130111'), NaT ])
        s = s2.append(s)
        result = s.to_string()
        self.assertTrue('NaT' in result)

        # nat in summary
        result = str(s2.index)
        self.assertTrue('NaT' in result)

    def test_timedelta64(self):

        from pandas import date_range
        from datetime import datetime, timedelta

        Series(np.array([1100, 20], dtype='timedelta64[s]')).to_string()

        s = Series(date_range('2012-1-1', periods=3, freq='D'))

        # GH2146

        # adding NaTs
        y = s-s.shift(1)
        result = y.to_string()
        self.assertTrue('1 days, 00:00:00' in result)
        self.assertTrue('NaT' in result)

        # with frac seconds
        o = Series([datetime(2012,1,1,microsecond=150)]*3)
        y = s-o
        result = y.to_string()
        self.assertTrue('-00:00:00.000150' in result)

        # rounding?
        o = Series([datetime(2012,1,1,1)]*3)
        y = s-o
        result = y.to_string()
        self.assertTrue('-01:00:00' in result)
        self.assertTrue('1 days, 23:00:00' in result)

        o = Series([datetime(2012,1,1,1,1)]*3)
        y = s-o
        result = y.to_string()
        self.assertTrue('-01:01:00' in result)
        self.assertTrue('1 days, 22:59:00' in result)

        o = Series([datetime(2012,1,1,1,1,microsecond=150)]*3)
        y = s-o
        result = y.to_string()
        self.assertTrue('-01:01:00.000150' in result)
        self.assertTrue('1 days, 22:58:59.999850' in result)

        # neg time
        td = timedelta(minutes=5,seconds=3)
        s2 = Series(date_range('2012-1-1', periods=3, freq='D')) + td
        y = s - s2
        result = y.to_string()
        self.assertTrue('-00:05:03' in result)

        td = timedelta(microseconds=550)
        s2 = Series(date_range('2012-1-1', periods=3, freq='D')) + td
        y = s - td
        result = y.to_string()
        self.assertTrue('2012-01-01 23:59:59.999450' in result)

    def test_mixed_datetime64(self):
        df = DataFrame({'A': [1, 2],
                        'B': ['2012-01-01', '2012-01-02']})
        df['B'] = pd.to_datetime(df.B)

        result = repr(df.ix[0])
        self.assertTrue('2012-01-01' in result)


class TestEngFormatter(unittest.TestCase):
    _multiprocess_can_split_ = True

    def test_eng_float_formatter(self):
        df = DataFrame({'A': [1.41, 141., 14100, 1410000.]})

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

        fmt.reset_option('^display.')

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


def _three_digit_exp():
    return '%.4g' % 1.7e8 == '1.7e+008'


class TestFloatArrayFormatter(unittest.TestCase):

    def test_misc(self):
        obj = fmt.FloatArrayFormatter(np.array([], dtype=np.float64))
        result = obj.get_result()

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
