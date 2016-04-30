# -*- coding: utf-8 -*-

# TODO(wesm): lots of issues making flake8 hard
# flake8: noqa

from __future__ import print_function
from distutils.version import LooseVersion
import re

from pandas.compat import (range, zip, lrange, StringIO, PY3,
                           u, lzip, is_platform_windows,
                           is_platform_32bit)
import pandas.compat as compat
import itertools
from operator import methodcaller
import os
import sys
from textwrap import dedent
import warnings

from numpy import nan
from numpy.random import randn
import numpy as np

import codecs

div_style = ''
try:
    import IPython
    if IPython.__version__ < LooseVersion('3.0.0'):
        div_style = ' style="max-width:1500px;overflow:auto;"'
except (ImportError, AttributeError):
    pass

from pandas import DataFrame, Series, Index, Timestamp, MultiIndex, date_range, NaT

import pandas.formats.format as fmt
import pandas.util.testing as tm
import pandas.core.common as com
import pandas.formats.printing as printing
from pandas.util.terminal import get_terminal_size
import pandas as pd
from pandas.core.config import (set_option, get_option, option_context,
                                reset_option)
from datetime import datetime

import nose

use_32bit_repr = is_platform_windows() or is_platform_32bit()

_frame = DataFrame(tm.getSeriesData())


def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth


def has_info_repr(df):
    r = repr(df)
    c1 = r.split('\n')[0].startswith("<class")
    c2 = r.split('\n')[0].startswith(r"&lt;class")  # _repr_html_
    return c1 or c2


def has_non_verbose_info_repr(df):
    has_info = has_info_repr(df)
    r = repr(df)
    nv = len(r.split(
        '\n')) == 6  # 1. <class>, 2. Index, 3. Columns, 4. dtype, 5. memory usage, 6. trailing newline
    return has_info and nv


def has_horizontally_truncated_repr(df):
    try:  # Check header row
        fst_line = np.array(repr(df).splitlines()[0].split())
        cand_col = np.where(fst_line == '...')[0][0]
    except:
        return False
    # Make sure each row has this ... in the same place
    r = repr(df)
    for ix, l in enumerate(r.splitlines()):
        if not r.split()[cand_col] == '...':
            return False
    return True


def has_vertically_truncated_repr(df):
    r = repr(df)
    only_dot_row = False
    for row in r.splitlines():
        if re.match('^[\.\ ]+$', row):
            only_dot_row = True
    return only_dot_row


def has_truncated_repr(df):
    return has_horizontally_truncated_repr(
        df) or has_vertically_truncated_repr(df)


def has_doubly_truncated_repr(df):
    return has_horizontally_truncated_repr(
        df) and has_vertically_truncated_repr(df)


def has_expanded_repr(df):
    r = repr(df)
    for line in r.split('\n'):
        if line.endswith('\\'):
            return True
    return False


class TestDataFrameFormatting(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.warn_filters = warnings.filters
        warnings.filterwarnings('ignore', category=FutureWarning,
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
        repr(self.frame)

        fmt.set_eng_float_format(use_eng_prefix=True)
        repr(self.frame)

        fmt.set_eng_float_format(accuracy=0)
        repr(self.frame)
        self.reset_display_options()

    def test_show_null_counts(self):

        df = DataFrame(1, columns=range(10), index=range(10))
        df.iloc[1, 1] = np.nan

        def check(null_counts, result):
            buf = StringIO()
            df.info(buf=buf, null_counts=null_counts)
            self.assertTrue(('non-null' in buf.getvalue()) is result)

        with option_context('display.max_info_rows', 20,
                            'display.max_info_columns', 20):
            check(None, True)
            check(True, True)
            check(False, False)

        with option_context('display.max_info_rows', 5,
                            'display.max_info_columns', 5):
            check(None, False)
            check(True, False)
            check(False, False)

    def test_repr_tuples(self):
        buf = StringIO()

        df = DataFrame({'tups': lzip(range(10), range(10))})
        repr(df)
        df.to_string(col_space=10, buf=buf)

    def test_repr_truncation(self):
        max_len = 20
        with option_context("display.max_colwidth", max_len):
            df = DataFrame({'A': np.random.randn(10),
                            'B': [tm.rands(np.random.randint(
                                max_len - 1, max_len + 1)) for i in range(10)
            ]})
            r = repr(df)
            r = r[r.find('\n') + 1:]

            adj = fmt._get_adjustment()

            for line, value in lzip(r.split('\n'), df['B']):
                if adj.len(value) + 1 > max_len:
                    self.assertIn('...', line)
                else:
                    self.assertNotIn('...', line)

        with option_context("display.max_colwidth", 999999):
            self.assertNotIn('...', repr(df))

        with option_context("display.max_colwidth", max_len + 2):
            self.assertNotIn('...', repr(df))

    def test_repr_chop_threshold(self):
        df = DataFrame([[0.1, 0.5], [0.5, -0.1]])
        pd.reset_option("display.chop_threshold")  # default None
        self.assertEqual(repr(df), '     0    1\n0  0.1  0.5\n1  0.5 -0.1')

        with option_context("display.chop_threshold", 0.2):
            self.assertEqual(repr(df), '     0    1\n0  0.0  0.5\n1  0.5  0.0')

        with option_context("display.chop_threshold", 0.6):
            self.assertEqual(repr(df), '     0    1\n0  0.0  0.0\n1  0.0  0.0')

        with option_context("display.chop_threshold", None):
            self.assertEqual(repr(df), '     0    1\n0  0.1  0.5\n1  0.5 -0.1')

    def test_repr_obeys_max_seq_limit(self):
        with option_context("display.max_seq_items", 2000):
            self.assertTrue(len(printing.pprint_thing(lrange(1000))) > 1000)

        with option_context("display.max_seq_items", 5):
            self.assertTrue(len(printing.pprint_thing(lrange(1000))) < 100)

    def test_repr_set(self):
        self.assertEqual(printing.pprint_thing(set([1])), '{1}')

    def test_repr_is_valid_construction_code(self):
        # for the case of Index, where the repr is traditional rather then
        # stylized
        idx = Index(['a', 'b'])
        res = eval("pd." + repr(idx))
        tm.assert_series_equal(Series(res), Series(idx))

    def test_repr_should_return_str(self):
        # http://docs.python.org/py3k/reference/datamodel.html#object.__repr__
        # http://docs.python.org/reference/datamodel.html#object.__repr__
        # "...The return value must be a string object."

        # (str on py2.x, str (unicode) on py3)

        data = [8, 5, 3, 5]
        index1 = [u("\u03c3"), u("\u03c4"), u("\u03c5"), u("\u03c6")]
        cols = [u("\u03c8")]
        df = DataFrame(data, columns=cols, index=index1)
        self.assertTrue(type(df.__repr__()) == str)  # both py2 / 3

    def test_repr_no_backslash(self):
        with option_context('mode.sim_interactive', True):
            df = DataFrame(np.random.randn(10, 4))
            self.assertTrue('\\' not in repr(df))

    def test_expand_frame_repr(self):
        df_small = DataFrame('hello', [0], [0])
        df_wide = DataFrame('hello', [0], lrange(10))
        df_tall = DataFrame('hello', lrange(30), lrange(5))

        with option_context('mode.sim_interactive', True):
            with option_context('display.max_columns', 10, 'display.width', 20,
                                'display.max_rows', 20,
                                'display.show_dimensions', True):
                with option_context('display.expand_frame_repr', True):
                    self.assertFalse(has_truncated_repr(df_small))
                    self.assertFalse(has_expanded_repr(df_small))
                    self.assertFalse(has_truncated_repr(df_wide))
                    self.assertTrue(has_expanded_repr(df_wide))
                    self.assertTrue(has_vertically_truncated_repr(df_tall))
                    self.assertTrue(has_expanded_repr(df_tall))

                with option_context('display.expand_frame_repr', False):
                    self.assertFalse(has_truncated_repr(df_small))
                    self.assertFalse(has_expanded_repr(df_small))
                    self.assertFalse(has_horizontally_truncated_repr(df_wide))
                    self.assertFalse(has_expanded_repr(df_wide))
                    self.assertTrue(has_vertically_truncated_repr(df_tall))
                    self.assertFalse(has_expanded_repr(df_tall))

    def test_repr_non_interactive(self):
        # in non interactive mode, there can be no dependency on the
        # result of terminal auto size detection
        df = DataFrame('hello', lrange(1000), lrange(5))

        with option_context('mode.sim_interactive', False, 'display.width', 0,
                            'display.height', 0, 'display.max_rows', 5000):
            self.assertFalse(has_truncated_repr(df))
            self.assertFalse(has_expanded_repr(df))

    def test_repr_max_columns_max_rows(self):
        term_width, term_height = get_terminal_size()
        if term_width < 10 or term_height < 10:
            raise nose.SkipTest("terminal size too small, "
                                "{0} x {1}".format(term_width, term_height))

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
                    self.assertTrue(has_doubly_truncated_repr(df6))

                with option_context('display.max_rows', 20,
                                    'display.max_columns', 10):
                    # Out off max_columns boundary, but no extending
                    # since not exceeding width
                    self.assertFalse(has_expanded_repr(df6))
                    self.assertFalse(has_truncated_repr(df6))

                with option_context('display.max_rows', 9,
                                    'display.max_columns', 10):
                    # out vertical bounds can not result in exanded repr
                    self.assertFalse(has_expanded_repr(df10))
                    self.assertTrue(has_vertically_truncated_repr(df10))

            # width=None in terminal, auto detection
            with option_context('display.max_columns', 100, 'display.max_rows',
                                term_width * 20, 'display.width', None):
                df = mkframe((term_width // 7) - 2)
                self.assertFalse(has_expanded_repr(df))
                df = mkframe((term_width // 7) + 2)
                printing.pprint_thing(df._repr_fits_horizontal_())
                self.assertTrue(has_expanded_repr(df))

    def test_str_max_colwidth(self):
        # GH 7856
        df = pd.DataFrame([{'a': 'foo',
                            'b': 'bar',
                            'c': 'uncomfortably long line with lots of stuff',
                            'd': 1}, {'a': 'foo',
                                      'b': 'bar',
                                      'c': 'stuff',
                                      'd': 1}])
        df.set_index(['a', 'b', 'c'])
        self.assertTrue(
            str(df) ==
            '     a    b                                           c  d\n'
            '0  foo  bar  uncomfortably long line with lots of stuff  1\n'
            '1  foo  bar                                       stuff  1')
        with option_context('max_colwidth', 20):
            self.assertTrue(str(df) == '     a    b                    c  d\n'
                            '0  foo  bar  uncomfortably lo...  1\n'
                            '1  foo  bar                stuff  1')

    def test_auto_detect(self):
        term_width, term_height = get_terminal_size()
        fac = 1.05  # Arbitrary large factor to exceed term widht
        cols = range(int(term_width * fac))
        index = range(10)
        df = DataFrame(index=index, columns=cols)
        with option_context('mode.sim_interactive', True):
            with option_context('max_rows', None):
                with option_context('max_columns', None):
                    # Wrap around with None
                    self.assertTrue(has_expanded_repr(df))
            with option_context('max_rows', 0):
                with option_context('max_columns', 0):
                    # Truncate with auto detection.
                    self.assertTrue(has_horizontally_truncated_repr(df))

            index = range(int(term_height * fac))
            df = DataFrame(index=index, columns=cols)
            with option_context('max_rows', 0):
                with option_context('max_columns', None):
                    # Wrap around with None
                    self.assertTrue(has_expanded_repr(df))
                    # Truncate vertically
                    self.assertTrue(has_vertically_truncated_repr(df))

            with option_context('max_rows', None):
                with option_context('max_columns', 0):
                    self.assertTrue(has_horizontally_truncated_repr(df))

    def test_to_string_repr_unicode(self):
        buf = StringIO()

        unicode_values = [u('\u03c3')] * 10
        unicode_values = np.array(unicode_values, dtype=object)
        df = DataFrame({'unicode': unicode_values})
        df.to_string(col_space=10, buf=buf)

        # it works!
        repr(df)

        idx = Index(['abc', u('\u03c3a'), 'aegdvg'])
        ser = Series(np.random.randn(len(idx)), idx)
        rs = repr(ser).split('\n')
        line_len = len(rs[0])
        for line in rs[1:]:
            try:
                line = line.decode(get_option("display.encoding"))
            except:
                pass
            if not line.startswith('dtype:'):
                self.assertEqual(len(line), line_len)

        # it works even if sys.stdin in None
        _stdin = sys.stdin
        try:
            sys.stdin = None
            repr(df)
        finally:
            sys.stdin = _stdin

    def test_to_string_unicode_columns(self):
        df = DataFrame({u('\u03c3'): np.arange(10.)})

        buf = StringIO()
        df.to_string(buf=buf)
        buf.getvalue()

        buf = StringIO()
        df.info(buf=buf)
        buf.getvalue()

        result = self.frame.to_string()
        tm.assertIsInstance(result, compat.text_type)

    def test_to_string_utf8_columns(self):
        n = u("\u05d0").encode('utf-8')

        with option_context('display.max_rows', 1):
            df = DataFrame([1, 2], columns=[n])
            repr(df)

    def test_to_string_unicode_two(self):
        dm = DataFrame({u('c/\u03c3'): []})
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
        df = DataFrame({u('c/\u03c3'): [1, 2, 3]})
        result = df.to_string(formatters={u('c/\u03c3'): lambda x: '%s' % x})
        self.assertEqual(result, u('  c/\u03c3\n') + '0   1\n1   2\n2   3')

    def test_east_asian_unicode_frame(self):
        if PY3:
            _rep = repr
        else:
            _rep = unicode

        # not alighned properly because of east asian width

        # mid col
        df = DataFrame({'a': [u'あ', u'いいい', u'う', u'ええええええ'],
                        'b': [1, 222, 33333, 4]},
                       index=['a', 'bb', 'c', 'ddd'])
        expected = (u"          a      b\na         あ      1\n"
                    u"bb      いいい    222\nc         う  33333\n"
                    u"ddd  ええええええ      4")
        self.assertEqual(_rep(df), expected)

        # last col
        df = DataFrame({'a': [1, 222, 33333, 4],
                        'b': [u'あ', u'いいい', u'う', u'ええええええ']},
                       index=['a', 'bb', 'c', 'ddd'])
        expected = (u"         a       b\na        1       あ\n"
                    u"bb     222     いいい\nc    33333       う\n"
                    u"ddd      4  ええええええ")
        self.assertEqual(_rep(df), expected)

        # all col
        df = DataFrame({'a': [u'あああああ', u'い', u'う', u'えええ'],
                        'b': [u'あ', u'いいい', u'う', u'ええええええ']},
                       index=['a', 'bb', 'c', 'ddd'])
        expected = (u"         a       b\na    あああああ       あ\n"
                    u"bb       い     いいい\nc        う       う\n"
                    u"ddd    えええ  ええええええ")
        self.assertEqual(_rep(df), expected)

        # column name
        df = DataFrame({u'あああああ': [1, 222, 33333, 4],
                        'b': [u'あ', u'いいい', u'う', u'ええええええ']},
                       index=['a', 'bb', 'c', 'ddd'])
        expected = (u"          b  あああああ\na         あ      1\n"
                    u"bb      いいい    222\nc         う  33333\n"
                    u"ddd  ええええええ      4")
        self.assertEqual(_rep(df), expected)

        # index
        df = DataFrame({'a': [u'あああああ', u'い', u'う', u'えええ'],
                        'b': [u'あ', u'いいい', u'う', u'ええええええ']},
                       index=[u'あああ', u'いいいいいい', u'うう', u'え'])
        expected = (u"            a       b\nあああ     あああああ       あ\n"
                    u"いいいいいい      い     いいい\nうう          う       う\n"
                    u"え         えええ  ええええええ")
        self.assertEqual(_rep(df), expected)

        # index name
        df = DataFrame({'a': [u'あああああ', u'い', u'う', u'えええ'],
                        'b': [u'あ', u'いいい', u'う', u'ええええええ']},
                       index=pd.Index([u'あ', u'い', u'うう', u'え'], name=u'おおおお'))
        expected = (u"          a       b\nおおおお               \nあ     あああああ       あ\n"
                    u"い         い     いいい\nうう        う       う\nえ       えええ  ええええええ"
                    )
        self.assertEqual(_rep(df), expected)

        # all
        df = DataFrame({u'あああ': [u'あああ', u'い', u'う', u'えええええ'],
                        u'いいいいい': [u'あ', u'いいい', u'う', u'ええ']},
                       index=pd.Index([u'あ', u'いいい', u'うう', u'え'], name=u'お'))
        expected = (u"       あああ いいいいい\nお               \nあ      あああ     あ\n"
                    u"いいい      い   いいい\nうう       う     う\nえ    えええええ    ええ")
        self.assertEqual(_rep(df), expected)

        # MultiIndex
        idx = pd.MultiIndex.from_tuples([(u'あ', u'いい'), (u'う', u'え'), (
            u'おおお', u'かかかか'), (u'き', u'くく')])
        df = DataFrame({'a': [u'あああああ', u'い', u'う', u'えええ'],
                        'b': [u'あ', u'いいい', u'う', u'ええええええ']}, index=idx)
        expected = (u"              a       b\nあ   いい    あああああ       あ\n"
                    u"う   え         い     いいい\nおおお かかかか      う       う\n"
                    u"き   くく      えええ  ええええええ")
        self.assertEqual(_rep(df), expected)

        # truncate
        with option_context('display.max_rows', 3, 'display.max_columns', 3):
            df = pd.DataFrame({'a': [u'あああああ', u'い', u'う', u'えええ'],
                               'b': [u'あ', u'いいい', u'う', u'ええええええ'],
                               'c': [u'お', u'か', u'ききき', u'くくくくくく'],
                               u'ああああ': [u'さ', u'し', u'す', u'せ']},
                              columns=['a', 'b', 'c', u'ああああ'])

            expected = (u"        a ...  ああああ\n0   あああああ ...     さ\n"
                        u"..    ... ...   ...\n3     えええ ...     せ\n"
                        u"\n[4 rows x 4 columns]")
            self.assertEqual(_rep(df), expected)

            df.index = [u'あああ', u'いいいい', u'う', 'aaa']
            expected = (u"         a ...  ああああ\nあああ  あああああ ...     さ\n"
                        u"..     ... ...   ...\naaa    えええ ...     せ\n"
                        u"\n[4 rows x 4 columns]")
            self.assertEqual(_rep(df), expected)

        # Emable Unicode option -----------------------------------------
        with option_context('display.unicode.east_asian_width', True):

            # mid col
            df = DataFrame({'a': [u'あ', u'いいい', u'う', u'ええええええ'],
                            'b': [1, 222, 33333, 4]},
                           index=['a', 'bb', 'c', 'ddd'])
            expected = (u"                a      b\na              あ      1\n"
                        u"bb         いいい    222\nc              う  33333\n"
                        u"ddd  ええええええ      4")
            self.assertEqual(_rep(df), expected)

            # last col
            df = DataFrame({'a': [1, 222, 33333, 4],
                            'b': [u'あ', u'いいい', u'う', u'ええええええ']},
                           index=['a', 'bb', 'c', 'ddd'])
            expected = (u"         a             b\na        1            あ\n"
                        u"bb     222        いいい\nc    33333            う\n"
                        u"ddd      4  ええええええ")
            self.assertEqual(_rep(df), expected)

            # all col
            df = DataFrame({'a': [u'あああああ', u'い', u'う', u'えええ'],
                            'b': [u'あ', u'いいい', u'う', u'ええええええ']},
                           index=['a', 'bb', 'c', 'ddd'])
            expected = (u"              a             b\na    あああああ            あ\n"
                        u"bb           い        いいい\nc            う            う\n"
                        u"ddd      えええ  ええええええ"
                        "")
            self.assertEqual(_rep(df), expected)

            # column name
            df = DataFrame({u'あああああ': [1, 222, 33333, 4],
                            'b': [u'あ', u'いいい', u'う', u'ええええええ']},
                           index=['a', 'bb', 'c', 'ddd'])
            expected = (u"                b  あああああ\na              あ           1\n"
                        u"bb         いいい         222\nc              う       33333\n"
                        u"ddd  ええええええ           4")
            self.assertEqual(_rep(df), expected)

            # index
            df = DataFrame({'a': [u'あああああ', u'い', u'う', u'えええ'],
                            'b': [u'あ', u'いいい', u'う', u'ええええええ']},
                           index=[u'あああ', u'いいいいいい', u'うう', u'え'])
            expected = (u"                       a             b\nあああ        あああああ            あ\n"
                        u"いいいいいい          い        いいい\nうう                  う            う\n"
                        u"え                えええ  ええええええ")
            self.assertEqual(_rep(df), expected)

            # index name
            df = DataFrame({'a': [u'あああああ', u'い', u'う', u'えええ'],
                            'b': [u'あ', u'いいい', u'う', u'ええええええ']},
                           index=pd.Index([u'あ', u'い', u'うう', u'え'], name=u'おおおお'))
            expected = (u"                   a             b\nおおおお                          \n"
                        u"あ        あああああ            あ\nい                い        いいい\n"
                        u"うう              う            う\nえ            えええ  ええええええ"
                        )
            self.assertEqual(_rep(df), expected)

            # all
            df = DataFrame({u'あああ': [u'あああ', u'い', u'う', u'えええええ'],
                            u'いいいいい': [u'あ', u'いいい', u'う', u'ええ']},
                           index=pd.Index([u'あ', u'いいい', u'うう', u'え'], name=u'お'))
            expected = (u"            あああ いいいいい\nお                           \n"
                        u"あ          あああ         あ\nいいい          い     いいい\n"
                        u"うう            う         う\nえ      えええええ       ええ")
            self.assertEqual(_rep(df), expected)

            # MultiIndex
            idx = pd.MultiIndex.from_tuples([(u'あ', u'いい'), (u'う', u'え'), (
                u'おおお', u'かかかか'), (u'き', u'くく')])
            df = DataFrame({'a': [u'あああああ', u'い', u'う', u'えええ'],
                            'b': [u'あ', u'いいい', u'う', u'ええええええ']}, index=idx)
            expected = (u"                          a             b\nあ     いい      あああああ            あ\n"
                        u"う     え                い        いいい\nおおお かかかか          う            う\n"
                        u"き     くく          えええ  ええええええ")
            self.assertEqual(_rep(df), expected)

            # truncate
            with option_context('display.max_rows', 3, 'display.max_columns',
                                3):

                df = pd.DataFrame({'a': [u'あああああ', u'い', u'う', u'えええ'],
                                   'b': [u'あ', u'いいい', u'う', u'ええええええ'],
                                   'c': [u'お', u'か', u'ききき', u'くくくくくく'],
                                   u'ああああ': [u'さ', u'し', u'す', u'せ']},
                                  columns=['a', 'b', 'c', u'ああああ'])

                expected = (u"             a   ...    ああああ\n0   あああああ   ...          さ\n"
                            u"..         ...   ...         ...\n3       えええ   ...          せ\n"
                            u"\n[4 rows x 4 columns]")
                self.assertEqual(_rep(df), expected)

                df.index = [u'あああ', u'いいいい', u'う', 'aaa']
                expected = (u"                 a   ...    ああああ\nあああ  あああああ   ...          さ\n"
                            u"...            ...   ...         ...\naaa         えええ   ...          せ\n"
                            u"\n[4 rows x 4 columns]")
                self.assertEqual(_rep(df), expected)

            # ambiguous unicode
            df = DataFrame({u'あああああ': [1, 222, 33333, 4],
                            'b': [u'あ', u'いいい', u'¡¡', u'ええええええ']},
                           index=['a', 'bb', 'c', '¡¡¡'])
            expected = (u"                b  あああああ\na              あ           1\n"
                        u"bb         いいい         222\nc              ¡¡       33333\n"
                        u"¡¡¡  ええええええ           4")
            self.assertEqual(_rep(df), expected)

    def test_to_string_buffer_all_unicode(self):
        buf = StringIO()

        empty = DataFrame({u('c/\u03c3'): Series()})
        nonempty = DataFrame({u('c/\u03c3'): Series([1, 2, 3])})

        print(empty, file=buf)
        print(nonempty, file=buf)

        # this should work
        buf.getvalue()

    def test_to_string_with_col_space(self):
        df = DataFrame(np.random.random(size=(1, 3)))
        c10 = len(df.to_string(col_space=10).split("\n")[1])
        c20 = len(df.to_string(col_space=20).split("\n")[1])
        c30 = len(df.to_string(col_space=30).split("\n")[1])
        self.assertTrue(c10 < c20 < c30)

        # GH 8230
        # col_space wasn't being applied with header=False
        with_header = df.to_string(col_space=20)
        with_header_row1 = with_header.splitlines()[1]
        no_header = df.to_string(col_space=20, header=False)
        self.assertEqual(len(with_header_row1), len(no_header))

    def test_to_string_truncate_indices(self):
        for index in [tm.makeStringIndex, tm.makeUnicodeIndex, tm.makeIntIndex,
                      tm.makeDateIndex, tm.makePeriodIndex]:
            for column in [tm.makeStringIndex]:
                for h in [10, 20]:
                    for w in [10, 20]:
                        with option_context("display.expand_frame_repr",
                                            False):
                            df = DataFrame(index=index(h), columns=column(w))
                            with option_context("display.max_rows", 15):
                                if h == 20:
                                    self.assertTrue(
                                        has_vertically_truncated_repr(df))
                                else:
                                    self.assertFalse(
                                        has_vertically_truncated_repr(df))
                            with option_context("display.max_columns", 15):
                                if w == 20:
                                    self.assertTrue(
                                        has_horizontally_truncated_repr(df))
                                else:
                                    self.assertFalse(
                                        has_horizontally_truncated_repr(df))
                            with option_context("display.max_rows", 15,
                                                "display.max_columns", 15):
                                if h == 20 and w == 20:
                                    self.assertTrue(has_doubly_truncated_repr(
                                        df))
                                else:
                                    self.assertFalse(has_doubly_truncated_repr(
                                        df))

    def test_to_string_truncate_multilevel(self):
        arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                  ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        df = DataFrame(index=arrays, columns=arrays)
        with option_context("display.max_rows", 7, "display.max_columns", 7):
            self.assertTrue(has_doubly_truncated_repr(df))

    def test_truncate_with_different_dtypes(self):

        # 11594, 12045
        # when truncated the dtypes of the splits can differ

        # 11594
        import datetime
        s = Series([datetime.datetime(2012, 1, 1)]*10 + [datetime.datetime(1012,1,2)] + [datetime.datetime(2012, 1, 3)]*10)

        with pd.option_context('display.max_rows', 8):
            result = str(s)
            self.assertTrue('object' in result)

        # 12045
        df = DataFrame({'text': ['some words'] + [None]*9})

        with pd.option_context('display.max_rows', 8, 'display.max_columns', 3):
            result = str(df)
            self.assertTrue('None' in result)
            self.assertFalse('NaN' in result)

    def test_datetimelike_frame(self):

        # GH 12211
        df = DataFrame({'date' : [pd.Timestamp('20130101').tz_localize('UTC')] + [pd.NaT]*5})

        with option_context("display.max_rows", 5):
            result = str(df)
            self.assertTrue('2013-01-01 00:00:00+00:00' in result)
            self.assertTrue('NaT' in result)
            self.assertTrue('...' in result)
            self.assertTrue('[6 rows x 1 columns]' in result)

        dts = [pd.Timestamp('2011-01-01', tz='US/Eastern')] * 5 + [pd.NaT] * 5
        df = pd.DataFrame({"dt": dts,
                           "x": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]})
        with option_context('display.max_rows', 5):
            expected = ('                          dt   x\n'
                        '0  2011-01-01 00:00:00-05:00   1\n'
                        '1  2011-01-01 00:00:00-05:00   2\n'
                        '..                       ...  ..\n'
                        '8                        NaT   9\n'
                        '9                        NaT  10\n\n'
                        '[10 rows x 2 columns]')
            self.assertEqual(repr(df), expected)

        dts = [pd.NaT] * 5 + [pd.Timestamp('2011-01-01', tz='US/Eastern')] * 5
        df = pd.DataFrame({"dt": dts,
                           "x": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]})
        with option_context('display.max_rows', 5):
            expected = ('                          dt   x\n'
                        '0                        NaT   1\n'
                        '1                        NaT   2\n'
                        '..                       ...  ..\n'
                        '8  2011-01-01 00:00:00-05:00   9\n'
                        '9  2011-01-01 00:00:00-05:00  10\n\n'
                        '[10 rows x 2 columns]')
            self.assertEqual(repr(df), expected)

        dts = ([pd.Timestamp('2011-01-01', tz='Asia/Tokyo')] * 5 +
               [pd.Timestamp('2011-01-01', tz='US/Eastern')] * 5)
        df = pd.DataFrame({"dt": dts,
                           "x": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]})
        with option_context('display.max_rows', 5):
            expected = ('                           dt   x\n'
                        '0   2011-01-01 00:00:00+09:00   1\n'
                        '1   2011-01-01 00:00:00+09:00   2\n'
                        '..                        ...  ..\n'
                        '8   2011-01-01 00:00:00-05:00   9\n'
                        '9   2011-01-01 00:00:00-05:00  10\n\n'
                        '[10 rows x 2 columns]')
            self.assertEqual(repr(df), expected)

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
        df = DataFrame({u('\u03c3'): np.arange(10.)})
        expected = u'<table border="1" class="dataframe">\n  <thead>\n    <tr style="text-align: right;">\n      <th></th>\n      <th>\u03c3</th>\n    </tr>\n  </thead>\n  <tbody>\n    <tr>\n      <th>0</th>\n      <td>0.0</td>\n    </tr>\n    <tr>\n      <th>1</th>\n      <td>1.0</td>\n    </tr>\n    <tr>\n      <th>2</th>\n      <td>2.0</td>\n    </tr>\n    <tr>\n      <th>3</th>\n      <td>3.0</td>\n    </tr>\n    <tr>\n      <th>4</th>\n      <td>4.0</td>\n    </tr>\n    <tr>\n      <th>5</th>\n      <td>5.0</td>\n    </tr>\n    <tr>\n      <th>6</th>\n      <td>6.0</td>\n    </tr>\n    <tr>\n      <th>7</th>\n      <td>7.0</td>\n    </tr>\n    <tr>\n      <th>8</th>\n      <td>8.0</td>\n    </tr>\n    <tr>\n      <th>9</th>\n      <td>9.0</td>\n    </tr>\n  </tbody>\n</table>'
        self.assertEqual(df.to_html(), expected)
        df = DataFrame({'A': [u('\u03c3')]})
        expected = u'<table border="1" class="dataframe">\n  <thead>\n    <tr style="text-align: right;">\n      <th></th>\n      <th>A</th>\n    </tr>\n  </thead>\n  <tbody>\n    <tr>\n      <th>0</th>\n      <td>\u03c3</td>\n    </tr>\n  </tbody>\n</table>'
        self.assertEqual(df.to_html(), expected)

    def test_to_html_decimal(self):
        # GH 12031
        df = DataFrame({'A': [6.0, 3.1, 2.2]})
        result = df.to_html(decimal=',')
        expected = ('<table border="1" class="dataframe">\n'
                    '  <thead>\n'
                    '    <tr style="text-align: right;">\n'
                    '      <th></th>\n'
                    '      <th>A</th>\n'
                    '    </tr>\n'
                    '  </thead>\n'
                    '  <tbody>\n'
                    '    <tr>\n'
                    '      <th>0</th>\n'
                    '      <td>6,0</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>1</th>\n'
                    '      <td>3,1</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>2</th>\n'
                    '      <td>2,2</td>\n'
                    '    </tr>\n'
                    '  </tbody>\n'
                    '</table>')
        self.assertEqual(result, expected)

    def test_to_html_escaped(self):
        a = 'str<ing1 &amp;'
        b = 'stri>ng2 &amp;'

        test_dict = {'co<l1': {a: "<type 'str'>",
                               b: "<type 'str'>"},
                     'co>l2': {a: "<type 'str'>",
                               b: "<type 'str'>"}}
        rs = DataFrame(test_dict).to_html()
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
      <td>&lt;type 'str'&gt;</td>
      <td>&lt;type 'str'&gt;</td>
    </tr>
    <tr>
      <th>stri&gt;ng2 &amp;amp;</th>
      <td>&lt;type 'str'&gt;</td>
      <td>&lt;type 'str'&gt;</td>
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
        rs = DataFrame(test_dict).to_html(escape=False)
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
      <td><b>bold</b></td>
      <td><b>bold</b></td>
    </tr>
    <tr>
      <th>stri>ng2 &amp;</th>
      <td><b>bold</b></td>
      <td><b>bold</b></td>
    </tr>
  </tbody>
</table>"""

        self.assertEqual(xp, rs)

    def test_to_html_multiindex_index_false(self):
        # issue 8452
        df = DataFrame({
            'a': range(2),
            'b': range(3, 5),
            'c': range(5, 7),
            'd': range(3, 5)
        })
        df.columns = MultiIndex.from_product([['a', 'b'], ['c', 'd']])
        result = df.to_html(index=False)
        expected = """\
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th colspan="2" halign="left">a</th>
      <th colspan="2" halign="left">b</th>
    </tr>
    <tr>
      <th>c</th>
      <th>d</th>
      <th>c</th>
      <th>d</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>0</td>
      <td>3</td>
      <td>5</td>
      <td>3</td>
    </tr>
    <tr>
      <td>1</td>
      <td>4</td>
      <td>6</td>
      <td>4</td>
    </tr>
  </tbody>
</table>"""

        self.assertEqual(result, expected)

        df.index = Index(df.index.values, name='idx')
        result = df.to_html(index=False)
        self.assertEqual(result, expected)

    def test_to_html_multiindex_sparsify_false_multi_sparse(self):
        with option_context('display.multi_sparse', False):
            index = MultiIndex.from_arrays([[0, 0, 1, 1], [0, 1, 0, 1]],
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
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>0</th>
      <th>1</th>
      <td>2</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <th>0</th>
      <td>4</td>
      <td>5</td>
    </tr>
    <tr>
      <th>1</th>
      <th>1</th>
      <td>6</td>
      <td>7</td>
    </tr>
  </tbody>
</table>"""

            self.assertEqual(result, expected)

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
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>0</th>
      <th>1</th>
      <td>2</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <th>0</th>
      <td>4</td>
      <td>5</td>
    </tr>
    <tr>
      <th>1</th>
      <th>1</th>
      <td>6</td>
      <td>7</td>
    </tr>
  </tbody>
</table>"""

            self.assertEqual(result, expected)

    def test_to_html_multiindex_sparsify(self):
        index = MultiIndex.from_arrays([[0, 0, 1, 1], [0, 1, 0, 1]],
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
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">1</th>
      <th>0</th>
      <td>4</td>
      <td>5</td>
    </tr>
    <tr>
      <th>1</th>
      <td>6</td>
      <td>7</td>
    </tr>
  </tbody>
</table>"""

        self.assertEqual(result, expected)

        df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]], columns=index[::2],
                       index=index)

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
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">1</th>
      <th>0</th>
      <td>4</td>
      <td>5</td>
    </tr>
    <tr>
      <th>1</th>
      <td>6</td>
      <td>7</td>
    </tr>
  </tbody>
</table>"""

        self.assertEqual(result, expected)

    def test_to_html_index_formatter(self):
        df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]], columns=['foo', None],
                       index=lrange(4))

        f = lambda x: 'abcd' [x]
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
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>b</th>
      <td>2</td>
      <td>3</td>
    </tr>
    <tr>
      <th>c</th>
      <td>4</td>
      <td>5</td>
    </tr>
    <tr>
      <th>d</th>
      <td>6</td>
      <td>7</td>
    </tr>
  </tbody>
</table>"""

        self.assertEqual(result, expected)

    def test_to_html_regression_GH6098(self):
        df = DataFrame({u('clé1'): [u('a'), u('a'), u('b'), u('b'), u('a')],
                        u('clé2'): [u('1er'), u('2ème'), u('1er'), u('2ème'),
                                    u('1er')],
                        'données1': np.random.randn(5),
                        'données2': np.random.randn(5)})
        # it works
        df.pivot_table(index=[u('clé1')], columns=[u('clé2')])._repr_html_()

    def test_to_html_truncate(self):
        raise nose.SkipTest("unreliable on travis")
        index = pd.DatetimeIndex(start='20010101', freq='D', periods=20)
        df = DataFrame(index=index, columns=range(20))
        fmt.set_option('display.max_rows', 8)
        fmt.set_option('display.max_columns', 4)
        result = df._repr_html_()
        expected = '''\
<div{0}>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>...</th>
      <th>18</th>
      <th>19</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2001-01-01</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2001-01-02</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2001-01-03</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2001-01-04</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>2001-01-17</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2001-01-18</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2001-01-19</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2001-01-20</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>20 rows × 20 columns</p>
</div>'''.format(div_style)
        if compat.PY2:
            expected = expected.decode('utf-8')
        self.assertEqual(result, expected)

    def test_to_html_truncate_multi_index(self):
        raise nose.SkipTest("unreliable on travis")
        arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                  ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        df = DataFrame(index=arrays, columns=arrays)
        fmt.set_option('display.max_rows', 7)
        fmt.set_option('display.max_columns', 7)
        result = df._repr_html_()
        expected = '''\
<div{0}>
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th></th>
      <th></th>
      <th colspan="2" halign="left">bar</th>
      <th>baz</th>
      <th>...</th>
      <th>foo</th>
      <th colspan="2" halign="left">qux</th>
    </tr>
    <tr>
      <th></th>
      <th></th>
      <th>one</th>
      <th>two</th>
      <th>one</th>
      <th>...</th>
      <th>two</th>
      <th>one</th>
      <th>two</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">bar</th>
      <th>one</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>two</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>baz</th>
      <th>one</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>...</th>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>foo</th>
      <th>two</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">qux</th>
      <th>one</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>two</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>8 rows × 8 columns</p>
</div>'''.format(div_style)
        if compat.PY2:
            expected = expected.decode('utf-8')
        self.assertEqual(result, expected)

    def test_to_html_truncate_multi_index_sparse_off(self):
        raise nose.SkipTest("unreliable on travis")
        arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                  ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        df = DataFrame(index=arrays, columns=arrays)
        fmt.set_option('display.max_rows', 7)
        fmt.set_option('display.max_columns', 7)
        fmt.set_option('display.multi_sparse', False)
        result = df._repr_html_()
        expected = '''\
<div{0}>
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th></th>
      <th></th>
      <th>bar</th>
      <th>bar</th>
      <th>baz</th>
      <th>...</th>
      <th>foo</th>
      <th>qux</th>
      <th>qux</th>
    </tr>
    <tr>
      <th></th>
      <th></th>
      <th>one</th>
      <th>two</th>
      <th>one</th>
      <th>...</th>
      <th>two</th>
      <th>one</th>
      <th>two</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>bar</th>
      <th>one</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>bar</th>
      <th>two</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>baz</th>
      <th>one</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>foo</th>
      <th>two</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>qux</th>
      <th>one</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>qux</th>
      <th>two</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>8 rows × 8 columns</p>
</div>'''.format(div_style)
        if compat.PY2:
            expected = expected.decode('utf-8')
        self.assertEqual(result, expected)

    def test_nonunicode_nonascii_alignment(self):
        df = DataFrame([["aa\xc3\xa4\xc3\xa4", 1], ["bbbb", 2]])
        rep_str = df.to_string()
        lines = rep_str.split('\n')
        self.assertEqual(len(lines[1]), len(lines[2]))

    def test_unicode_problem_decoding_as_ascii(self):
        dm = DataFrame({u('c/\u03c3'): Series({'test': np.NaN})})
        compat.text_type(dm.to_string())

    def test_string_repr_encoding(self):
        filepath = tm.get_data_path('unicode_series.csv')
        df = pd.read_csv(filepath, header=None, encoding='latin1')
        repr(df)
        repr(df[1])

    def test_repr_corner(self):
        # representing infs poses no problems
        df = DataFrame({'foo': np.inf * np.empty(10)})
        repr(df)

    def test_frame_info_encoding(self):
        index = ['\'Til There Was You (1997)',
                 'ldum klaka (Cold Fever) (1994)']
        fmt.set_option('display.max_rows', 1)
        df = DataFrame(columns=['a', 'b', 'c'], index=index)
        repr(df)
        repr(df.T)
        fmt.set_option('display.max_rows', 200)

    def test_pprint_thing(self):
        from pandas.formats.printing import pprint_thing as pp_t

        if PY3:
            raise nose.SkipTest("doesn't work on Python 3")

        self.assertEqual(pp_t('a'), u('a'))
        self.assertEqual(pp_t(u('a')), u('a'))
        self.assertEqual(pp_t(None), 'None')
        self.assertEqual(pp_t(u('\u05d0'), quote_strings=True), u("u'\u05d0'"))
        self.assertEqual(pp_t(u('\u05d0'), quote_strings=False), u('\u05d0'))
        self.assertEqual(pp_t((u('\u05d0'),
                               u('\u05d1')), quote_strings=True),
                         u("(u'\u05d0', u'\u05d1')"))
        self.assertEqual(pp_t((u('\u05d0'), (u('\u05d1'),
                                             u('\u05d2'))),
                              quote_strings=True),
                         u("(u'\u05d0', (u'\u05d1', u'\u05d2'))"))
        self.assertEqual(pp_t(('foo', u('\u05d0'), (u('\u05d0'),
                                                    u('\u05d0'))),
                              quote_strings=True),
                         u("(u'foo', u'\u05d0', (u'\u05d0', u'\u05d0'))"))

        # escape embedded tabs in string
        # GH #2038
        self.assertTrue(not "\t" in pp_t("a\tb", escape_chars=("\t", )))

    def test_wide_repr(self):
        with option_context('mode.sim_interactive', True,
                            'display.show_dimensions', True):
            max_cols = get_option('display.max_columns')
            df = DataFrame(tm.rands_array(25, size=(10, max_cols - 1)))
            set_option('display.expand_frame_repr', False)
            rep_str = repr(df)

            assert "10 rows x %d columns" % (max_cols - 1) in rep_str
            set_option('display.expand_frame_repr', True)
            wide_repr = repr(df)
            self.assertNotEqual(rep_str, wide_repr)

            with option_context('display.width', 120):
                wider_repr = repr(df)
                self.assertTrue(len(wider_repr) < len(wide_repr))

        reset_option('display.expand_frame_repr')

    def test_wide_repr_wide_columns(self):
        with option_context('mode.sim_interactive', True):
            df = DataFrame(randn(5, 3), columns=['a' * 90, 'b' * 90, 'c' * 90])
            rep_str = repr(df)

            self.assertEqual(len(rep_str.splitlines()), 20)

    def test_wide_repr_named(self):
        with option_context('mode.sim_interactive', True):
            max_cols = get_option('display.max_columns')
            df = DataFrame(tm.rands_array(25, size=(10, max_cols - 1)))
            df.index.name = 'DataFrame Index'
            set_option('display.expand_frame_repr', False)

            rep_str = repr(df)
            set_option('display.expand_frame_repr', True)
            wide_repr = repr(df)
            self.assertNotEqual(rep_str, wide_repr)

            with option_context('display.width', 150):
                wider_repr = repr(df)
                self.assertTrue(len(wider_repr) < len(wide_repr))

            for line in wide_repr.splitlines()[1::13]:
                self.assertIn('DataFrame Index', line)

        reset_option('display.expand_frame_repr')

    def test_wide_repr_multiindex(self):
        with option_context('mode.sim_interactive', True):
            midx = MultiIndex.from_arrays(tm.rands_array(5, size=(2, 10)))
            max_cols = get_option('display.max_columns')
            df = DataFrame(tm.rands_array(25, size=(10, max_cols - 1)),
                           index=midx)
            df.index.names = ['Level 0', 'Level 1']
            set_option('display.expand_frame_repr', False)
            rep_str = repr(df)
            set_option('display.expand_frame_repr', True)
            wide_repr = repr(df)
            self.assertNotEqual(rep_str, wide_repr)

            with option_context('display.width', 150):
                wider_repr = repr(df)
                self.assertTrue(len(wider_repr) < len(wide_repr))

            for line in wide_repr.splitlines()[1::13]:
                self.assertIn('Level 0 Level 1', line)

        reset_option('display.expand_frame_repr')

    def test_wide_repr_multiindex_cols(self):
        with option_context('mode.sim_interactive', True):
            max_cols = get_option('display.max_columns')
            midx = MultiIndex.from_arrays(tm.rands_array(5, size=(2, 10)))
            mcols = MultiIndex.from_arrays(tm.rands_array(3, size=(2, max_cols
                                                                   - 1)))
            df = DataFrame(tm.rands_array(25, (10, max_cols - 1)),
                           index=midx, columns=mcols)
            df.index.names = ['Level 0', 'Level 1']
            set_option('display.expand_frame_repr', False)
            rep_str = repr(df)
            set_option('display.expand_frame_repr', True)
            wide_repr = repr(df)
            self.assertNotEqual(rep_str, wide_repr)

        with option_context('display.width', 150):
            wider_repr = repr(df)
            self.assertTrue(len(wider_repr) < len(wide_repr))

        reset_option('display.expand_frame_repr')

    def test_wide_repr_unicode(self):
        with option_context('mode.sim_interactive', True):
            max_cols = get_option('display.max_columns')
            df = DataFrame(tm.rands_array(25, size=(10, max_cols - 1)))
            set_option('display.expand_frame_repr', False)
            rep_str = repr(df)
            set_option('display.expand_frame_repr', True)
            wide_repr = repr(df)
            self.assertNotEqual(rep_str, wide_repr)

            with option_context('display.width', 150):
                wider_repr = repr(df)
                self.assertTrue(len(wider_repr) < len(wide_repr))

        reset_option('display.expand_frame_repr')

    def test_wide_repr_wide_long_columns(self):
        with option_context('mode.sim_interactive', True):
            df = DataFrame({'a': ['a' * 30, 'b' * 30],
                            'b': ['c' * 70, 'd' * 80]})

            result = repr(df)
            self.assertTrue('ccccc' in result)
            self.assertTrue('ddddd' in result)

    def test_long_series(self):
        n = 1000
        s = Series(
            np.random.randint(-50, 50, n),
            index=['s%04d' % x for x in range(n)], dtype='int64')

        import re
        str_rep = str(s)
        nmatches = len(re.findall('dtype', str_rep))
        self.assertEqual(nmatches, 1)

    def test_index_with_nan(self):
        #  GH 2850
        df = DataFrame({'id1': {0: '1a3',
                                1: '9h4'},
                        'id2': {0: np.nan,
                                1: 'd67'},
                        'id3': {0: '78d',
                                1: '79d'},
                        'value': {0: 123,
                                  1: 64}})

        # multi-index
        y = df.set_index(['id1', 'id2', 'id3'])
        result = y.to_string()
        expected = u(
            '             value\nid1 id2 id3       \n1a3 NaN 78d    123\n9h4 d67 79d     64')
        self.assertEqual(result, expected)

        # index
        y = df.set_index('id2')
        result = y.to_string()
        expected = u(
            '     id1  id3  value\nid2                 \nNaN  1a3  78d    123\nd67  9h4  79d     64')
        self.assertEqual(result, expected)

        # with append (this failed in 0.12)
        y = df.set_index(['id1', 'id2']).set_index('id3', append=True)
        result = y.to_string()
        expected = u(
            '             value\nid1 id2 id3       \n1a3 NaN 78d    123\n9h4 d67 79d     64')
        self.assertEqual(result, expected)

        # all-nan in mi
        df2 = df.copy()
        df2.ix[:, 'id2'] = np.nan
        y = df2.set_index('id2')
        result = y.to_string()
        expected = u(
            '     id1  id3  value\nid2                 \nNaN  1a3  78d    123\nNaN  9h4  79d     64')
        self.assertEqual(result, expected)

        # partial nan in mi
        df2 = df.copy()
        df2.ix[:, 'id2'] = np.nan
        y = df2.set_index(['id2', 'id3'])
        result = y.to_string()
        expected = u(
            '         id1  value\nid2 id3            \nNaN 78d  1a3    123\n    79d  9h4     64')
        self.assertEqual(result, expected)

        df = DataFrame({'id1': {0: np.nan,
                                1: '9h4'},
                        'id2': {0: np.nan,
                                1: 'd67'},
                        'id3': {0: np.nan,
                                1: '79d'},
                        'value': {0: 123,
                                  1: 64}})

        y = df.set_index(['id1', 'id2', 'id3'])
        result = y.to_string()
        expected = u(
            '             value\nid1 id2 id3       \nNaN NaN NaN    123\n9h4 d67 79d     64')
        self.assertEqual(result, expected)

    def test_to_string(self):
        from pandas import read_table
        import re

        # big mixed
        biggie = DataFrame({'A': randn(200),
                            'B': tm.makeStringIndex(200)},
                           index=lrange(200))

        biggie.loc[:20, 'A'] = nan
        biggie.loc[:20, 'B'] = nan
        s = biggie.to_string()

        buf = StringIO()
        retval = biggie.to_string(buf=buf)
        self.assertIsNone(retval)
        self.assertEqual(buf.getvalue(), s)

        tm.assertIsInstance(s, compat.string_types)

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
        self.assertTrue((np.abs(recons['A'].dropna() - biggie['A'].dropna()) <
                         0.1).all())

        # expected = ['B', 'A']
        # self.assertEqual(header, expected)

        result = biggie.to_string(columns=['A'], col_space=17)
        header = result.split('\n')[0].strip().split()
        expected = ['A']
        self.assertEqual(header, expected)

        biggie.to_string(columns=['B', 'A'],
                         formatters={'A': lambda x: '%.1f' % x})

        biggie.to_string(columns=['B', 'A'], float_format=str)
        biggie.to_string(columns=['B', 'A'], col_space=12, float_format=str)

        frame = DataFrame(index=np.arange(200))
        frame.to_string()

    def test_to_string_no_header(self):
        df = DataFrame({'x': [1, 2, 3], 'y': [4, 5, 6]})

        df_s = df.to_string(header=False)
        expected = "0  1  4\n1  2  5\n2  3  6"

        self.assertEqual(df_s, expected)

    def test_to_string_no_index(self):
        df = DataFrame({'x': [1, 2, 3], 'y': [4, 5, 6]})

        df_s = df.to_string(index=False)
        expected = "x  y\n1  4\n2  5\n3  6"

        self.assertEqual(df_s, expected)

    def test_to_string_float_formatting(self):
        self.reset_display_options()
        fmt.set_option('display.precision', 5, 'display.column_space', 12,
                       'display.notebook_repr_html', False)

        df = DataFrame({'x': [0, 0.25, 3456.000, 12e+45, 1.64e+6, 1.7e+8,
                              1.253456, np.pi, -1e6]})

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
        self.assertEqual(df_s, expected)

        df = DataFrame({'x': [3234, 0.253]})
        df_s = df.to_string()

        expected = ('          x\n' '0  3234.000\n' '1     0.253')
        self.assertEqual(df_s, expected)

        self.reset_display_options()
        self.assertEqual(get_option("display.precision"), 6)

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
        self.assertEqual(df_s, expected)

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
        expected = ('   0\n' '0  0\n' '1  0\n' '2 -0')

    def test_to_string_float_index(self):
        index = Index([1.5, 2, 3, 4, 5])
        df = DataFrame(lrange(5), index=index)

        result = df.to_string()
        expected = ('     0\n'
                    '1.5  0\n'
                    '2.0  1\n'
                    '3.0  2\n'
                    '4.0  3\n'
                    '5.0  4')
        self.assertEqual(result, expected)

    def test_to_string_ascii_error(self):
        data = [('0  ', u('                        .gitignore '), u('     5 '),
                 ' \xe2\x80\xa2\xe2\x80\xa2\xe2\x80'
                 '\xa2\xe2\x80\xa2\xe2\x80\xa2')]
        df = DataFrame(data)

        # it works!
        repr(df)

    def test_to_string_int_formatting(self):
        df = DataFrame({'x': [-15, 20, 25, -35]})
        self.assertTrue(issubclass(df['x'].dtype.type, np.integer))

        output = df.to_string()
        expected = ('    x\n' '0 -15\n' '1  20\n' '2  25\n' '3 -35')
        self.assertEqual(output, expected)

    def test_to_string_index_formatter(self):
        df = DataFrame([lrange(5), lrange(5, 10), lrange(10, 15)])

        rs = df.to_string(formatters={'__index__': lambda x: 'abc' [x]})

        xp = """\
    0   1   2   3   4
a   0   1   2   3   4
b   5   6   7   8   9
c  10  11  12  13  14\
"""

        self.assertEqual(rs, xp)

    def test_to_string_left_justify_cols(self):
        self.reset_display_options()
        df = DataFrame({'x': [3234, 0.253]})
        df_s = df.to_string(justify='left')
        expected = ('   x       \n' '0  3234.000\n' '1     0.253')
        self.assertEqual(df_s, expected)

    def test_to_string_format_na(self):
        self.reset_display_options()
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

        expected = ('     A       B\n'
                    '0  NaN     NaN\n'
                    '1 -1.0     foo\n'
                    '2 -2.0   foooo\n'
                    '3  3.0  fooooo\n'
                    '4  4.0     bar')
        self.assertEqual(result, expected)

    def test_to_string_line_width(self):
        df = DataFrame(123, lrange(10, 15), lrange(30))
        s = df.to_string(line_width=80)
        self.assertEqual(max(len(l) for l in s.split('\n')), 80)

    def test_show_dimensions(self):
        df = DataFrame(123, lrange(10, 15), lrange(30))

        with option_context('display.max_rows', 10, 'display.max_columns', 40,
                            'display.width', 500, 'display.expand_frame_repr',
                            'info', 'display.show_dimensions', True):
            self.assertTrue('5 rows' in str(df))
            self.assertTrue('5 rows' in df._repr_html_())
        with option_context('display.max_rows', 10, 'display.max_columns', 40,
                            'display.width', 500, 'display.expand_frame_repr',
                            'info', 'display.show_dimensions', False):
            self.assertFalse('5 rows' in str(df))
            self.assertFalse('5 rows' in df._repr_html_())
        with option_context('display.max_rows', 2, 'display.max_columns', 2,
                            'display.width', 500, 'display.expand_frame_repr',
                            'info', 'display.show_dimensions', 'truncate'):
            self.assertTrue('5 rows' in str(df))
            self.assertTrue('5 rows' in df._repr_html_())
        with option_context('display.max_rows', 10, 'display.max_columns', 40,
                            'display.width', 500, 'display.expand_frame_repr',
                            'info', 'display.show_dimensions', 'truncate'):
            self.assertFalse('5 rows' in str(df))
            self.assertFalse('5 rows' in df._repr_html_())

    def test_to_html(self):
        # big mixed
        biggie = DataFrame({'A': randn(200),
                            'B': tm.makeStringIndex(200)},
                           index=lrange(200))

        biggie.loc[:20, 'A'] = nan
        biggie.loc[:20, 'B'] = nan
        s = biggie.to_html()

        buf = StringIO()
        retval = biggie.to_html(buf=buf)
        self.assertIsNone(retval)
        self.assertEqual(buf.getvalue(), s)

        tm.assertIsInstance(s, compat.string_types)

        biggie.to_html(columns=['B', 'A'], col_space=17)
        biggie.to_html(columns=['B', 'A'],
                       formatters={'A': lambda x: '%.1f' % x})

        biggie.to_html(columns=['B', 'A'], float_format=str)
        biggie.to_html(columns=['B', 'A'], col_space=12, float_format=str)

        frame = DataFrame(index=np.arange(200))
        frame.to_html()

    def test_to_html_filename(self):
        biggie = DataFrame({'A': randn(200),
                            'B': tm.makeStringIndex(200)},
                           index=lrange(200))

        biggie.loc[:20, 'A'] = nan
        biggie.loc[:20, 'B'] = nan
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
        self.assertFalse('<strong' in ashtml[ashtml.find("</thead>")])

    def test_to_html_columns_arg(self):
        result = self.frame.to_html(columns=['A'])
        self.assertNotIn('<th>B</th>', result)

    def test_to_html_multiindex(self):
        columns = MultiIndex.from_tuples(list(zip(np.arange(2).repeat(2),
                                                  np.mod(lrange(4), 2))),
                                         names=['CL0', 'CL1'])
        df = DataFrame([list('abcd'), list('efgh')], columns=columns)
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
                    '      <td>a</td>\n'
                    '      <td>b</td>\n'
                    '      <td>c</td>\n'
                    '      <td>d</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>1</th>\n'
                    '      <td>e</td>\n'
                    '      <td>f</td>\n'
                    '      <td>g</td>\n'
                    '      <td>h</td>\n'
                    '    </tr>\n'
                    '  </tbody>\n'
                    '</table>')

        self.assertEqual(result, expected)

        columns = MultiIndex.from_tuples(list(zip(
            range(4), np.mod(
                lrange(4), 2))))
        df = DataFrame([list('abcd'), list('efgh')], columns=columns)

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
                    '      <td>a</td>\n'
                    '      <td>b</td>\n'
                    '      <td>c</td>\n'
                    '      <td>d</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>1</th>\n'
                    '      <td>e</td>\n'
                    '      <td>f</td>\n'
                    '      <td>g</td>\n'
                    '      <td>h</td>\n'
                    '    </tr>\n'
                    '  </tbody>\n'
                    '</table>')

        self.assertEqual(result, expected)

    def test_to_html_justify(self):
        df = DataFrame({'A': [6, 30000, 2],
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
                    '      <td>6</td>\n'
                    '      <td>1</td>\n'
                    '      <td>223442</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>1</th>\n'
                    '      <td>30000</td>\n'
                    '      <td>2</td>\n'
                    '      <td>0</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>2</th>\n'
                    '      <td>2</td>\n'
                    '      <td>70000</td>\n'
                    '      <td>1</td>\n'
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
                    '      <td>6</td>\n'
                    '      <td>1</td>\n'
                    '      <td>223442</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>1</th>\n'
                    '      <td>30000</td>\n'
                    '      <td>2</td>\n'
                    '      <td>0</td>\n'
                    '    </tr>\n'
                    '    <tr>\n'
                    '      <th>2</th>\n'
                    '      <td>2</td>\n'
                    '      <td>70000</td>\n'
                    '      <td>1</td>\n'
                    '    </tr>\n'
                    '  </tbody>\n'
                    '</table>')
        self.assertEqual(result, expected)

    def test_to_html_index(self):
        index = ['foo', 'bar', 'baz']
        df = DataFrame({'A': [1, 2, 3],
                        'B': [1.2, 3.4, 5.6],
                        'C': ['one', 'two', np.NaN]},
                       columns=['A', 'B', 'C'],
                       index=index)
        expected_with_index = ('<table border="1" class="dataframe">\n'
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
                               '      <th>foo</th>\n'
                               '      <td>1</td>\n'
                               '      <td>1.2</td>\n'
                               '      <td>one</td>\n'
                               '    </tr>\n'
                               '    <tr>\n'
                               '      <th>bar</th>\n'
                               '      <td>2</td>\n'
                               '      <td>3.4</td>\n'
                               '      <td>two</td>\n'
                               '    </tr>\n'
                               '    <tr>\n'
                               '      <th>baz</th>\n'
                               '      <td>3</td>\n'
                               '      <td>5.6</td>\n'
                               '      <td>NaN</td>\n'
                               '    </tr>\n'
                               '  </tbody>\n'
                               '</table>')
        self.assertEqual(df.to_html(), expected_with_index)

        expected_without_index = ('<table border="1" class="dataframe">\n'
                                  '  <thead>\n'
                                  '    <tr style="text-align: right;">\n'
                                  '      <th>A</th>\n'
                                  '      <th>B</th>\n'
                                  '      <th>C</th>\n'
                                  '    </tr>\n'
                                  '  </thead>\n'
                                  '  <tbody>\n'
                                  '    <tr>\n'
                                  '      <td>1</td>\n'
                                  '      <td>1.2</td>\n'
                                  '      <td>one</td>\n'
                                  '    </tr>\n'
                                  '    <tr>\n'
                                  '      <td>2</td>\n'
                                  '      <td>3.4</td>\n'
                                  '      <td>two</td>\n'
                                  '    </tr>\n'
                                  '    <tr>\n'
                                  '      <td>3</td>\n'
                                  '      <td>5.6</td>\n'
                                  '      <td>NaN</td>\n'
                                  '    </tr>\n'
                                  '  </tbody>\n'
                                  '</table>')
        result = df.to_html(index=False)
        for i in index:
            self.assertNotIn(i, result)
        self.assertEqual(result, expected_without_index)
        df.index = Index(['foo', 'bar', 'baz'], name='idx')
        expected_with_index = ('<table border="1" class="dataframe">\n'
                               '  <thead>\n'
                               '    <tr style="text-align: right;">\n'
                               '      <th></th>\n'
                               '      <th>A</th>\n'
                               '      <th>B</th>\n'
                               '      <th>C</th>\n'
                               '    </tr>\n'
                               '    <tr>\n'
                               '      <th>idx</th>\n'
                               '      <th></th>\n'
                               '      <th></th>\n'
                               '      <th></th>\n'
                               '    </tr>\n'
                               '  </thead>\n'
                               '  <tbody>\n'
                               '    <tr>\n'
                               '      <th>foo</th>\n'
                               '      <td>1</td>\n'
                               '      <td>1.2</td>\n'
                               '      <td>one</td>\n'
                               '    </tr>\n'
                               '    <tr>\n'
                               '      <th>bar</th>\n'
                               '      <td>2</td>\n'
                               '      <td>3.4</td>\n'
                               '      <td>two</td>\n'
                               '    </tr>\n'
                               '    <tr>\n'
                               '      <th>baz</th>\n'
                               '      <td>3</td>\n'
                               '      <td>5.6</td>\n'
                               '      <td>NaN</td>\n'
                               '    </tr>\n'
                               '  </tbody>\n'
                               '</table>')
        self.assertEqual(df.to_html(), expected_with_index)
        self.assertEqual(df.to_html(index=False), expected_without_index)

        tuples = [('foo', 'car'), ('foo', 'bike'), ('bar', 'car')]
        df.index = MultiIndex.from_tuples(tuples)

        expected_with_index = ('<table border="1" class="dataframe">\n'
                               '  <thead>\n'
                               '    <tr style="text-align: right;">\n'
                               '      <th></th>\n'
                               '      <th></th>\n'
                               '      <th>A</th>\n'
                               '      <th>B</th>\n'
                               '      <th>C</th>\n'
                               '    </tr>\n'
                               '  </thead>\n'
                               '  <tbody>\n'
                               '    <tr>\n'
                               '      <th rowspan="2" valign="top">foo</th>\n'
                               '      <th>car</th>\n'
                               '      <td>1</td>\n'
                               '      <td>1.2</td>\n'
                               '      <td>one</td>\n'
                               '    </tr>\n'
                               '    <tr>\n'
                               '      <th>bike</th>\n'
                               '      <td>2</td>\n'
                               '      <td>3.4</td>\n'
                               '      <td>two</td>\n'
                               '    </tr>\n'
                               '    <tr>\n'
                               '      <th>bar</th>\n'
                               '      <th>car</th>\n'
                               '      <td>3</td>\n'
                               '      <td>5.6</td>\n'
                               '      <td>NaN</td>\n'
                               '    </tr>\n'
                               '  </tbody>\n'
                               '</table>')
        self.assertEqual(df.to_html(), expected_with_index)

        result = df.to_html(index=False)
        for i in ['foo', 'bar', 'car', 'bike']:
            self.assertNotIn(i, result)
        # must be the same result as normal index
        self.assertEqual(result, expected_without_index)

        df.index = MultiIndex.from_tuples(tuples, names=['idx1', 'idx2'])
        expected_with_index = ('<table border="1" class="dataframe">\n'
                               '  <thead>\n'
                               '    <tr style="text-align: right;">\n'
                               '      <th></th>\n'
                               '      <th></th>\n'
                               '      <th>A</th>\n'
                               '      <th>B</th>\n'
                               '      <th>C</th>\n'
                               '    </tr>\n'
                               '    <tr>\n'
                               '      <th>idx1</th>\n'
                               '      <th>idx2</th>\n'
                               '      <th></th>\n'
                               '      <th></th>\n'
                               '      <th></th>\n'
                               '    </tr>\n'
                               '  </thead>\n'
                               '  <tbody>\n'
                               '    <tr>\n'
                               '      <th rowspan="2" valign="top">foo</th>\n'
                               '      <th>car</th>\n'
                               '      <td>1</td>\n'
                               '      <td>1.2</td>\n'
                               '      <td>one</td>\n'
                               '    </tr>\n'
                               '    <tr>\n'
                               '      <th>bike</th>\n'
                               '      <td>2</td>\n'
                               '      <td>3.4</td>\n'
                               '      <td>two</td>\n'
                               '    </tr>\n'
                               '    <tr>\n'
                               '      <th>bar</th>\n'
                               '      <th>car</th>\n'
                               '      <td>3</td>\n'
                               '      <td>5.6</td>\n'
                               '      <td>NaN</td>\n'
                               '    </tr>\n'
                               '  </tbody>\n'
                               '</table>')
        self.assertEqual(df.to_html(), expected_with_index)
        self.assertEqual(df.to_html(index=False), expected_without_index)

    def test_repr_html(self):
        self.frame._repr_html_()

        fmt.set_option('display.max_rows', 1, 'display.max_columns', 1)
        self.frame._repr_html_()

        fmt.set_option('display.notebook_repr_html', False)
        self.frame._repr_html_()

        self.reset_display_options()

        df = DataFrame([[1, 2], [3, 4]])
        fmt.set_option('display.show_dimensions', True)
        self.assertTrue('2 rows' in df._repr_html_())
        fmt.set_option('display.show_dimensions', False)
        self.assertFalse('2 rows' in df._repr_html_())

        self.reset_display_options()

    def test_repr_html_wide(self):
        max_cols = get_option('display.max_columns')
        df = DataFrame(tm.rands_array(25, size=(10, max_cols - 1)))
        reg_repr = df._repr_html_()
        assert "..." not in reg_repr

        wide_df = DataFrame(tm.rands_array(25, size=(10, max_cols + 1)))
        wide_repr = wide_df._repr_html_()
        assert "..." in wide_repr

    def test_repr_html_wide_multiindex_cols(self):
        max_cols = get_option('display.max_columns')

        mcols = MultiIndex.from_product([np.arange(max_cols // 2),
                                         ['foo', 'bar']],
                                        names=['first', 'second'])
        df = DataFrame(tm.rands_array(25, size=(10, len(mcols))),
                       columns=mcols)
        reg_repr = df._repr_html_()
        assert '...' not in reg_repr

        mcols = MultiIndex.from_product((np.arange(1 + (max_cols // 2)),
                                         ['foo', 'bar']),
                                        names=['first', 'second'])
        df = DataFrame(tm.rands_array(25, size=(10, len(mcols))),
                       columns=mcols)
        wide_repr = df._repr_html_()
        assert '...' in wide_repr

    def test_repr_html_long(self):
        max_rows = get_option('display.max_rows')
        h = max_rows - 1
        df = DataFrame({'A': np.arange(1, 1 + h), 'B': np.arange(41, 41 + h)})
        reg_repr = df._repr_html_()
        assert '..' not in reg_repr
        assert str(41 + max_rows // 2) in reg_repr

        h = max_rows + 1
        df = DataFrame({'A': np.arange(1, 1 + h), 'B': np.arange(41, 41 + h)})
        long_repr = df._repr_html_()
        assert '..' in long_repr
        assert str(41 + max_rows // 2) not in long_repr
        assert u('%d rows ') % h in long_repr
        assert u('2 columns') in long_repr

    def test_repr_html_float(self):
        max_rows = get_option('display.max_rows')
        h = max_rows - 1
        df = DataFrame({'idx': np.linspace(-10, 10, h),
                        'A': np.arange(1, 1 + h),
                        'B': np.arange(41, 41 + h)}).set_index('idx')
        reg_repr = df._repr_html_()
        assert '..' not in reg_repr
        assert str(40 + h) in reg_repr

        h = max_rows + 1
        df = DataFrame({'idx': np.linspace(-10, 10, h),
                        'A': np.arange(1, 1 + h),
                        'B': np.arange(41, 41 + h)}).set_index('idx')
        long_repr = df._repr_html_()
        assert '..' in long_repr
        assert '31' not in long_repr
        assert u('%d rows ') % h in long_repr
        assert u('2 columns') in long_repr

    def test_repr_html_long_multiindex(self):
        max_rows = get_option('display.max_rows')
        max_L1 = max_rows // 2

        tuples = list(itertools.product(np.arange(max_L1), ['foo', 'bar']))
        idx = MultiIndex.from_tuples(tuples, names=['first', 'second'])
        df = DataFrame(np.random.randn(max_L1 * 2, 2), index=idx,
                       columns=['A', 'B'])
        reg_repr = df._repr_html_()
        assert '...' not in reg_repr

        tuples = list(itertools.product(np.arange(max_L1 + 1), ['foo', 'bar']))
        idx = MultiIndex.from_tuples(tuples, names=['first', 'second'])
        df = DataFrame(np.random.randn((max_L1 + 1) * 2, 2), index=idx,
                       columns=['A', 'B'])
        long_repr = df._repr_html_()
        assert '...' in long_repr

    def test_repr_html_long_and_wide(self):
        max_cols = get_option('display.max_columns')
        max_rows = get_option('display.max_rows')

        h, w = max_rows - 1, max_cols - 1
        df = DataFrame(dict((k, np.arange(1, 1 + h)) for k in np.arange(w)))
        assert '...' not in df._repr_html_()

        h, w = max_rows + 1, max_cols + 1
        df = DataFrame(dict((k, np.arange(1, 1 + h)) for k in np.arange(w)))
        assert '...' in df._repr_html_()

    def test_info_repr(self):
        max_rows = get_option('display.max_rows')
        max_cols = get_option('display.max_columns')
        # Long
        h, w = max_rows + 1, max_cols - 1
        df = DataFrame(dict((k, np.arange(1, 1 + h)) for k in np.arange(w)))
        assert has_vertically_truncated_repr(df)
        with option_context('display.large_repr', 'info'):
            assert has_info_repr(df)

        # Wide
        h, w = max_rows - 1, max_cols + 1
        df = DataFrame(dict((k, np.arange(1, 1 + h)) for k in np.arange(w)))
        assert has_horizontally_truncated_repr(df)
        with option_context('display.large_repr', 'info'):
            assert has_info_repr(df)

    def test_info_repr_max_cols(self):
        # GH #6939
        df = DataFrame(randn(10, 5))
        with option_context('display.large_repr', 'info',
                            'display.max_columns', 1,
                            'display.max_info_columns', 4):
            self.assertTrue(has_non_verbose_info_repr(df))

        with option_context('display.large_repr', 'info',
                            'display.max_columns', 1,
                            'display.max_info_columns', 5):
            self.assertFalse(has_non_verbose_info_repr(df))

        # test verbose overrides
        # fmt.set_option('display.max_info_columns', 4)  # exceeded

    def test_info_repr_html(self):
        max_rows = get_option('display.max_rows')
        max_cols = get_option('display.max_columns')
        # Long
        h, w = max_rows + 1, max_cols - 1
        df = DataFrame(dict((k, np.arange(1, 1 + h)) for k in np.arange(w)))
        assert r'&lt;class' not in df._repr_html_()
        with option_context('display.large_repr', 'info'):
            assert r'&lt;class' in df._repr_html_()

        # Wide
        h, w = max_rows - 1, max_cols + 1
        df = DataFrame(dict((k, np.arange(1, 1 + h)) for k in np.arange(w)))
        assert '<class' not in df._repr_html_()
        with option_context('display.large_repr', 'info'):
            assert '&lt;class' in df._repr_html_()

    def test_fake_qtconsole_repr_html(self):
        def get_ipython():
            return {'config': {'KernelApp':
                               {'parent_appname': 'ipython-qtconsole'}}}

        repstr = self.frame._repr_html_()
        self.assertIsNotNone(repstr)

        fmt.set_option('display.max_rows', 5, 'display.max_columns', 2)
        repstr = self.frame._repr_html_()
        self.assertIn('class', repstr)  # info fallback

        self.reset_display_options()

    def test_to_html_with_classes(self):
        df = DataFrame()
        result = df.to_html(classes="sortable draggable")
        expected = dedent("""

            <table border="1" class="dataframe sortable draggable">
              <thead>
                <tr style="text-align: right;">
                  <th></th>
                </tr>
              </thead>
              <tbody>
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
                return 3  # obviously simplified

        df = DataFrame([A()])
        repr(df)  # just don't dine

    def test_float_trim_zeros(self):
        vals = [2.08430917305e+10, 3.52205017305e+10, 2.30674817305e+10,
                2.03954217305e+10, 5.59897817305e+10]
        skip = True
        for line in repr(DataFrame({'A': vals})).split('\n')[:-2]:
            if line.startswith('dtype:'):
                continue
            if _three_digit_exp():
                self.assertTrue(('+010' in line) or skip)
            else:
                self.assertTrue(('+10' in line) or skip)
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

        # test with utf-8 and encoding option (GH 7061)
        df = DataFrame([[u'au\xdfgangen']])
        with tm.ensure_clean('test.tex') as path:
            df.to_latex(path, encoding='utf-8')
            with codecs.open(path, 'r', encoding='utf-8') as f:
                self.assertEqual(df.to_latex(), f.read())

        # test with utf-8 without encoding option
        if compat.PY3:  # python3 default encoding is utf-8
            with tm.ensure_clean('test.tex') as path:
                df.to_latex(path)
                with codecs.open(path, 'r') as f:
                    self.assertEqual(df.to_latex(), f.read())
        else:
            # python2 default encoding is ascii, so an error should be raised
            with tm.ensure_clean('test.tex') as path:
                self.assertRaises(UnicodeEncodeError, df.to_latex, path)

    def test_to_latex(self):
        # it works!
        self.frame.to_latex()

        df = DataFrame({'a': [1, 2], 'b': ['b1', 'b2']})
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

    def test_to_latex_format(self):
        # GH Bug #9402
        self.frame.to_latex(column_format='ccc')

        df = DataFrame({'a': [1, 2], 'b': ['b1', 'b2']})
        withindex_result = df.to_latex(column_format='ccc')
        withindex_expected = r"""\begin{tabular}{ccc}
\toprule
{} &  a &   b \\
\midrule
0 &  1 &  b1 \\
1 &  2 &  b2 \\
\bottomrule
\end{tabular}
"""

        self.assertEqual(withindex_result, withindex_expected)

    def test_to_latex_multiindex(self):
        df = DataFrame({('x', 'y'): ['a']})
        result = df.to_latex()
        expected = r"""\begin{tabular}{ll}
\toprule
{} &  x \\
{} &  y \\
\midrule
0 &  a \\
\bottomrule
\end{tabular}
"""

        self.assertEqual(result, expected)

        result = df.T.to_latex()
        expected = r"""\begin{tabular}{lll}
\toprule
  &   &  0 \\
\midrule
x & y &  a \\
\bottomrule
\end{tabular}
"""

        self.assertEqual(result, expected)

        df = DataFrame.from_dict({
            ('c1', 0): pd.Series(dict((x, x) for x in range(4))),
            ('c1', 1): pd.Series(dict((x, x + 4) for x in range(4))),
            ('c2', 0): pd.Series(dict((x, x) for x in range(4))),
            ('c2', 1): pd.Series(dict((x, x + 4) for x in range(4))),
            ('c3', 0): pd.Series(dict((x, x) for x in range(4))),
        }).T
        result = df.to_latex()
        expected = r"""\begin{tabular}{llrrrr}
\toprule
   &   &  0 &  1 &  2 &  3 \\
\midrule
c1 & 0 &  0 &  1 &  2 &  3 \\
   & 1 &  4 &  5 &  6 &  7 \\
c2 & 0 &  0 &  1 &  2 &  3 \\
   & 1 &  4 &  5 &  6 &  7 \\
c3 & 0 &  0 &  1 &  2 &  3 \\
\bottomrule
\end{tabular}
"""

        self.assertEqual(result, expected)

        # GH 10660
        df = pd.DataFrame({'a': [0, 0, 1, 1],
                           'b': list('abab'),
                           'c': [1, 2, 3, 4]})
        result = df.set_index(['a', 'b']).to_latex()
        expected = r"""\begin{tabular}{llr}
\toprule
  &   &  c \\
a & b &    \\
\midrule
0 & a &  1 \\
  & b &  2 \\
1 & a &  3 \\
  & b &  4 \\
\bottomrule
\end{tabular}
"""

        self.assertEqual(result, expected)

        result = df.groupby('a').describe().to_latex()
        expected = r"""\begin{tabular}{llr}
\toprule
  &       &         c \\
a & {} &           \\
\midrule
0 & count &  2.000000 \\
  & mean &  1.500000 \\
  & std &  0.707107 \\
  & min &  1.000000 \\
  & 25\% &  1.250000 \\
  & 50\% &  1.500000 \\
  & 75\% &  1.750000 \\
  & max &  2.000000 \\
1 & count &  2.000000 \\
  & mean &  3.500000 \\
  & std &  0.707107 \\
  & min &  3.000000 \\
  & 25\% &  3.250000 \\
  & 50\% &  3.500000 \\
  & 75\% &  3.750000 \\
  & max &  4.000000 \\
\bottomrule
\end{tabular}
"""

        self.assertEqual(result, expected)

    def test_to_latex_escape(self):
        a = 'a'
        b = 'b'

        test_dict = {u('co^l1'): {a: "a",
                                  b: "b"},
                     u('co$e^x$'): {a: "a",
                                    b: "b"}}

        unescaped_result = DataFrame(test_dict).to_latex(escape=False)
        escaped_result = DataFrame(test_dict).to_latex(
        )  # default: escape=True

        unescaped_expected = r'''\begin{tabular}{lll}
\toprule
{} & co$e^x$ & co^l1 \\
\midrule
a &       a &     a \\
b &       b &     b \\
\bottomrule
\end{tabular}
'''

        escaped_expected = r'''\begin{tabular}{lll}
\toprule
{} & co\$e\textasciicircumx\$ & co\textasciicircuml1 \\
\midrule
a &       a &     a \\
b &       b &     b \\
\bottomrule
\end{tabular}
'''

        self.assertEqual(unescaped_result, unescaped_expected)
        self.assertEqual(escaped_result, escaped_expected)

    def test_to_latex_longtable(self):
        self.frame.to_latex(longtable=True)

        df = DataFrame({'a': [1, 2], 'b': ['b1', 'b2']})
        withindex_result = df.to_latex(longtable=True)
        withindex_expected = r"""\begin{longtable}{lrl}
\toprule
{} &  a &   b \\
\midrule
\endhead
\midrule
\multicolumn{3}{r}{{Continued on next page}} \\
\midrule
\endfoot

\bottomrule
\endlastfoot
0 &  1 &  b1 \\
1 &  2 &  b2 \\
\end{longtable}
"""

        self.assertEqual(withindex_result, withindex_expected)

        withoutindex_result = df.to_latex(index=False, longtable=True)
        withoutindex_expected = r"""\begin{longtable}{rl}
\toprule
 a &   b \\
\midrule
\endhead
\midrule
\multicolumn{3}{r}{{Continued on next page}} \\
\midrule
\endfoot

\bottomrule
\endlastfoot
 1 &  b1 \\
 2 &  b2 \\
\end{longtable}
"""

        self.assertEqual(withoutindex_result, withoutindex_expected)

    def test_to_latex_escape_special_chars(self):
        special_characters = ['&', '%', '$', '#', '_', '{', '}', '~', '^',
                              '\\']
        df = DataFrame(data=special_characters)
        observed = df.to_latex()
        expected = r"""\begin{tabular}{ll}
\toprule
{} &  0 \\
\midrule
0 &  \& \\
1 &  \% \\
2 &  \$ \\
3 &  \# \\
4 &  \_ \\
5 &  \{ \\
6 &  \} \\
7 &  \textasciitilde \\
8 &  \textasciicircum \\
9 &  \textbackslash \\
\bottomrule
\end{tabular}
"""

        self.assertEqual(observed, expected)

    def test_to_latex_no_header(self):
        # GH 7124
        df = DataFrame({'a': [1, 2], 'b': ['b1', 'b2']})
        withindex_result = df.to_latex(header=False)
        withindex_expected = r"""\begin{tabular}{lrl}
\toprule
0 &  1 &  b1 \\
1 &  2 &  b2 \\
\bottomrule
\end{tabular}
"""

        self.assertEqual(withindex_result, withindex_expected)

        withoutindex_result = df.to_latex(index=False, header=False)
        withoutindex_expected = r"""\begin{tabular}{rl}
\toprule
 1 &  b1 \\
 2 &  b2 \\
\bottomrule
\end{tabular}
"""

        self.assertEqual(withoutindex_result, withoutindex_expected)

    def test_to_latex_decimal(self):
        # GH 12031
        self.frame.to_latex()
        df = DataFrame({'a': [1.0, 2.1], 'b': ['b1', 'b2']})
        withindex_result = df.to_latex(decimal=',')
        print("WHAT THE")
        withindex_expected = r"""\begin{tabular}{lrl}
\toprule
{} &    a &   b \\
\midrule
0 &  1,0 &  b1 \\
1 &  2,1 &  b2 \\
\bottomrule
\end{tabular}
"""

        self.assertEqual(withindex_result, withindex_expected)

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
                self.assertEqual(f.read(), expected)
        with tm.ensure_clean('test.csv') as path:
            df.to_csv(path, quoting=1, engine='python')
            with open(path, 'r') as f:
                self.assertEqual(f.read(), expected)

        expected = """\
$$,$col$
$0$,$1$
$1$,$2$
"""

        with tm.ensure_clean('test.csv') as path:
            df.to_csv(path, quoting=1, quotechar="$")
            with open(path, 'r') as f:
                self.assertEqual(f.read(), expected)
        with tm.ensure_clean('test.csv') as path:
            df.to_csv(path, quoting=1, quotechar="$", engine='python')
            with open(path, 'r') as f:
                self.assertEqual(f.read(), expected)

        with tm.ensure_clean('test.csv') as path:
            with tm.assertRaisesRegexp(TypeError, 'quotechar'):
                df.to_csv(path, quoting=1, quotechar=None)
        with tm.ensure_clean('test.csv') as path:
            with tm.assertRaisesRegexp(TypeError, 'quotechar'):
                df.to_csv(path, quoting=1, quotechar=None, engine='python')

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
                self.assertEqual(f.read(), expected)
        with tm.ensure_clean('test.csv') as path:
            df.to_csv(path, quoting=1, doublequote=True, engine='python')
            with open(path, 'r') as f:
                self.assertEqual(f.read(), expected)

        from _csv import Error
        with tm.ensure_clean('test.csv') as path:
            with tm.assertRaisesRegexp(Error, 'escapechar'):
                df.to_csv(path, doublequote=False)  # no escapechar set
        with tm.ensure_clean('test.csv') as path:
            with tm.assertRaisesRegexp(Error, 'escapechar'):
                df.to_csv(path, doublequote=False, engine='python')

    def test_to_csv_escapechar(self):
        df = DataFrame({'col': ['a"a', '"bb"']})
        expected = """\
"","col"
"0","a\\"a"
"1","\\"bb\\""
"""

        with tm.ensure_clean('test.csv') as path:  # QUOTE_ALL
            df.to_csv(path, quoting=1, doublequote=False, escapechar='\\')
            with open(path, 'r') as f:
                self.assertEqual(f.read(), expected)
        with tm.ensure_clean('test.csv') as path:
            df.to_csv(path, quoting=1, doublequote=False, escapechar='\\',
                      engine='python')
            with open(path, 'r') as f:
                self.assertEqual(f.read(), expected)

        df = DataFrame({'col': ['a,a', ',bb,']})
        expected = """\
,col
0,a\\,a
1,\\,bb\\,
"""

        with tm.ensure_clean('test.csv') as path:
            df.to_csv(path, quoting=3, escapechar='\\')  # QUOTE_NONE
            with open(path, 'r') as f:
                self.assertEqual(f.read(), expected)
        with tm.ensure_clean('test.csv') as path:
            df.to_csv(path, quoting=3, escapechar='\\', engine='python')
            with open(path, 'r') as f:
                self.assertEqual(f.read(), expected)

    def test_csv_to_string(self):
        df = DataFrame({'col': [1, 2]})
        expected = ',col\n0,1\n1,2\n'
        self.assertEqual(df.to_csv(), expected)

    def test_to_csv_decimal(self):
        # GH 781
        df = DataFrame({'col1': [1], 'col2': ['a'], 'col3': [10.1]})

        expected_default = ',col1,col2,col3\n0,1,a,10.1\n'
        self.assertEqual(df.to_csv(), expected_default)

        expected_european_excel = ';col1;col2;col3\n0;1;a;10,1\n'
        self.assertEqual(
            df.to_csv(decimal=',', sep=';'), expected_european_excel)

        expected_float_format_default = ',col1,col2,col3\n0,1,a,10.10\n'
        self.assertEqual(
            df.to_csv(float_format='%.2f'), expected_float_format_default)

        expected_float_format = ';col1;col2;col3\n0;1;a;10,10\n'
        self.assertEqual(
            df.to_csv(decimal=',', sep=';',
                      float_format='%.2f'), expected_float_format)

        # GH 11553: testing if decimal is taken into account for '0.0'
        df = pd.DataFrame({'a': [0, 1.1], 'b': [2.2, 3.3], 'c': 1})
        expected = 'a,b,c\n0^0,2^2,1\n1^1,3^3,1\n'
        self.assertEqual(df.to_csv(index=False, decimal='^'), expected)

        # same but for an index
        self.assertEqual(df.set_index('a').to_csv(decimal='^'), expected)

        # same for a multi-index
        self.assertEqual(
            df.set_index(['a', 'b']).to_csv(decimal="^"), expected)

    def test_to_csv_float_format(self):
        # testing if float_format is taken into account for the index
        # GH 11553
        df = pd.DataFrame({'a': [0, 1], 'b': [2.2, 3.3], 'c': 1})
        expected = 'a,b,c\n0,2.20,1\n1,3.30,1\n'
        self.assertEqual(
            df.set_index('a').to_csv(float_format='%.2f'), expected)

        # same for a multi-index
        self.assertEqual(
            df.set_index(['a', 'b']).to_csv(float_format='%.2f'), expected)

    def test_to_csv_na_rep(self):
        # testing if NaN values are correctly represented in the index
        # GH 11553
        df = DataFrame({'a': [0, np.NaN], 'b': [0, 1], 'c': [2, 3]})
        expected = "a,b,c\n0.0,0,2\n_,1,3\n"
        self.assertEqual(df.set_index('a').to_csv(na_rep='_'), expected)
        self.assertEqual(df.set_index(['a', 'b']).to_csv(na_rep='_'), expected)

        # now with an index containing only NaNs
        df = DataFrame({'a': np.NaN, 'b': [0, 1], 'c': [2, 3]})
        expected = "a,b,c\n_,0,2\n_,1,3\n"
        self.assertEqual(df.set_index('a').to_csv(na_rep='_'), expected)
        self.assertEqual(df.set_index(['a', 'b']).to_csv(na_rep='_'), expected)

        # check if na_rep parameter does not break anything when no NaN
        df = DataFrame({'a': 0, 'b': [0, 1], 'c': [2, 3]})
        expected = "a,b,c\n0,0,2\n0,1,3\n"
        self.assertEqual(df.set_index('a').to_csv(na_rep='_'), expected)
        self.assertEqual(df.set_index(['a', 'b']).to_csv(na_rep='_'), expected)

    def test_to_csv_date_format(self):
        # GH 10209
        df_sec = DataFrame({'A': pd.date_range('20130101', periods=5, freq='s')
                            })
        df_day = DataFrame({'A': pd.date_range('20130101', periods=5, freq='d')
                            })

        expected_default_sec = ',A\n0,2013-01-01 00:00:00\n1,2013-01-01 00:00:01\n2,2013-01-01 00:00:02' + \
                               '\n3,2013-01-01 00:00:03\n4,2013-01-01 00:00:04\n'
        self.assertEqual(df_sec.to_csv(), expected_default_sec)

        expected_ymdhms_day = ',A\n0,2013-01-01 00:00:00\n1,2013-01-02 00:00:00\n2,2013-01-03 00:00:00' + \
                              '\n3,2013-01-04 00:00:00\n4,2013-01-05 00:00:00\n'
        self.assertEqual(
            df_day.to_csv(
                date_format='%Y-%m-%d %H:%M:%S'), expected_ymdhms_day)

        expected_ymd_sec = ',A\n0,2013-01-01\n1,2013-01-01\n2,2013-01-01\n3,2013-01-01\n4,2013-01-01\n'
        self.assertEqual(
            df_sec.to_csv(date_format='%Y-%m-%d'), expected_ymd_sec)

        expected_default_day = ',A\n0,2013-01-01\n1,2013-01-02\n2,2013-01-03\n3,2013-01-04\n4,2013-01-05\n'
        self.assertEqual(df_day.to_csv(), expected_default_day)
        self.assertEqual(
            df_day.to_csv(date_format='%Y-%m-%d'), expected_default_day)

        # testing if date_format parameter is taken into account for
        # multi-indexed dataframes (GH 7791)
        df_sec['B'] = 0
        df_sec['C'] = 1
        expected_ymd_sec = 'A,B,C\n2013-01-01,0,1\n'
        df_sec_grouped = df_sec.groupby([pd.Grouper(key='A', freq='1h'), 'B'])
        self.assertEqual(df_sec_grouped.mean().to_csv(date_format='%Y-%m-%d'),
                         expected_ymd_sec)

    # deprecation GH11274
    def test_to_csv_engine_kw_deprecation(self):
        with tm.assert_produces_warning(FutureWarning):
            df = DataFrame({'col1': [1], 'col2': ['a'], 'col3': [10.1]})
            df.to_csv(engine='python')

    def test_period(self):
        # GH 12615
        df = pd.DataFrame({'A': pd.period_range('2013-01',
                                                periods=4, freq='M'),
                           'B': [pd.Period('2011-01', freq='M'),
                                 pd.Period('2011-02-01', freq='D'),
                                 pd.Period('2011-03-01 09:00', freq='H'),
                                 pd.Period('2011-04', freq='M')],
                           'C': list('abcd')})
        exp = ("        A                B  C\n0 2013-01          2011-01  a\n"
               "1 2013-02       2011-02-01  b\n2 2013-03 2011-03-01 09:00  c\n"
               "3 2013-04          2011-04  d")
        self.assertEqual(str(df), exp)


class TestSeriesFormatting(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.ts = tm.makeTimeSeries()

    def test_repr_unicode(self):
        s = Series([u('\u03c3')] * 10)
        repr(s)

        a = Series([u("\u05d0")] * 1000)
        a.name = 'title1'
        repr(a)

    def test_to_string(self):
        buf = StringIO()

        s = self.ts.to_string()

        retval = self.ts.to_string(buf=buf)
        self.assertIsNone(retval)
        self.assertEqual(buf.getvalue().strip(), s)

        # pass float_format
        format = '%.4f'.__mod__
        result = self.ts.to_string(float_format=format)
        result = [x.split()[1] for x in result.split('\n')[:-1]]
        expected = [format(x) for x in self.ts]
        self.assertEqual(result, expected)

        # empty string
        result = self.ts[:0].to_string()
        self.assertEqual(result, 'Series([], Freq: B)')

        result = self.ts[:0].to_string(length=0)
        self.assertEqual(result, 'Series([], Freq: B)')

        # name and length
        cp = self.ts.copy()
        cp.name = 'foo'
        result = cp.to_string(length=True, name=True, dtype=True)
        last_line = result.split('\n')[-1].strip()
        self.assertEqual(last_line,
                         "Freq: B, Name: foo, Length: %d, dtype: float64" %
                         len(cp))

    def test_freq_name_separation(self):
        s = Series(np.random.randn(10),
                   index=date_range('1/1/2000', periods=10), name=0)

        result = repr(s)
        self.assertTrue('Freq: D, Name: 0' in result)

    def test_to_string_mixed(self):
        s = Series(['foo', np.nan, -1.23, 4.56])
        result = s.to_string()
        expected = (u('0     foo\n') + u('1     NaN\n') + u('2   -1.23\n') +
                    u('3    4.56'))
        self.assertEqual(result, expected)

        # but don't count NAs as floats
        s = Series(['foo', np.nan, 'bar', 'baz'])
        result = s.to_string()
        expected = (u('0    foo\n') + '1    NaN\n' + '2    bar\n' + '3    baz')
        self.assertEqual(result, expected)

        s = Series(['foo', 5, 'bar', 'baz'])
        result = s.to_string()
        expected = (u('0    foo\n') + '1      5\n' + '2    bar\n' + '3    baz')
        self.assertEqual(result, expected)

    def test_to_string_float_na_spacing(self):
        s = Series([0., 1.5678, 2., -3., 4.])
        s[::2] = np.nan

        result = s.to_string()
        expected = (u('0       NaN\n') + '1    1.5678\n' + '2       NaN\n' +
                    '3   -3.0000\n' + '4       NaN')
        self.assertEqual(result, expected)

    def test_to_string_without_index(self):
        # GH 11729 Test index=False option
        s = Series([1, 2, 3, 4])
        result = s.to_string(index=False)
        expected = (u('1\n') + '2\n' + '3\n' + '4')
        self.assertEqual(result, expected)

    def test_unicode_name_in_footer(self):
        s = Series([1, 2], name=u('\u05e2\u05d1\u05e8\u05d9\u05ea'))
        sf = fmt.SeriesFormatter(s, name=u('\u05e2\u05d1\u05e8\u05d9\u05ea'))
        sf._get_footer()  # should not raise exception

    def test_east_asian_unicode_series(self):
        if PY3:
            _rep = repr
        else:
            _rep = unicode
        # not alighned properly because of east asian width

        # unicode index
        s = Series(['a', 'bb', 'CCC', 'D'],
                   index=[u'あ', u'いい', u'ううう', u'ええええ'])
        expected = (u"あ         a\nいい       bb\nううう     CCC\n"
                    u"ええええ      D\ndtype: object")
        self.assertEqual(_rep(s), expected)

        # unicode values
        s = Series([u'あ', u'いい', u'ううう', u'ええええ'],
                   index=['a', 'bb', 'c', 'ddd'])
        expected = (u"a         あ\nbb       いい\nc       ううう\n"
                    u"ddd    ええええ\ndtype: object")
        self.assertEqual(_rep(s), expected)

        # both
        s = Series([u'あ', u'いい', u'ううう', u'ええええ'],
                   index=[u'ああ', u'いいいい', u'う', u'えええ'])
        expected = (u"ああ         あ\nいいいい      いい\nう        ううう\n"
                    u"えええ     ええええ\ndtype: object")
        self.assertEqual(_rep(s), expected)

        # unicode footer
        s = Series([u'あ', u'いい', u'ううう', u'ええええ'],
                   index=[u'ああ', u'いいいい', u'う', u'えええ'], name=u'おおおおおおお')
        expected = (u"ああ         あ\nいいいい      いい\nう        ううう\n"
                    u"えええ     ええええ\nName: おおおおおおお, dtype: object")
        self.assertEqual(_rep(s), expected)

        # MultiIndex
        idx = pd.MultiIndex.from_tuples([(u'あ', u'いい'), (u'う', u'え'), (
            u'おおお', u'かかかか'), (u'き', u'くく')])
        s = Series([1, 22, 3333, 44444], index=idx)
        expected = (u"あ    いい          1\nう    え          22\nおおお  かかかか     3333\n"
                    u"き    くく      44444\ndtype: int64")
        self.assertEqual(_rep(s), expected)

        # object dtype, shorter than unicode repr
        s = Series([1, 22, 3333, 44444], index=[1, 'AB', np.nan, u'あああ'])
        expected = (u"1          1\nAB        22\nNaN     3333\n"
                    u"あああ    44444\ndtype: int64")
        self.assertEqual(_rep(s), expected)

        # object dtype, longer than unicode repr
        s = Series([1, 22, 3333, 44444],
                   index=[1, 'AB', pd.Timestamp('2011-01-01'), u'あああ'])
        expected = (u"1                          1\nAB                        22\n"
                    u"2011-01-01 00:00:00     3333\nあああ                    44444\ndtype: int64"
                    )
        self.assertEqual(_rep(s), expected)

        # truncate
        with option_context('display.max_rows', 3):
            s = Series([u'あ', u'いい', u'ううう', u'ええええ'], name=u'おおおおおおお')

            expected = (u"0       あ\n     ... \n"
                        u"3    ええええ\nName: おおおおおおお, dtype: object")
            self.assertEqual(_rep(s), expected)

            s.index = [u'ああ', u'いいいい', u'う', u'えええ']
            expected = (u"ああ        あ\n       ... \n"
                        u"えええ    ええええ\nName: おおおおおおお, dtype: object")
            self.assertEqual(_rep(s), expected)

        # Emable Unicode option -----------------------------------------
        with option_context('display.unicode.east_asian_width', True):

            # unicode index
            s = Series(['a', 'bb', 'CCC', 'D'],
                       index=[u'あ', u'いい', u'ううう', u'ええええ'])
            expected = (u"あ            a\nいい         bb\nううう      CCC\n"
                        u"ええええ      D\ndtype: object")
            self.assertEqual(_rep(s), expected)

            # unicode values
            s = Series([u'あ', u'いい', u'ううう', u'ええええ'],
                       index=['a', 'bb', 'c', 'ddd'])
            expected = (u"a            あ\nbb         いい\nc        ううう\n"
                        u"ddd    ええええ\ndtype: object")
            self.assertEqual(_rep(s), expected)

            # both
            s = Series([u'あ', u'いい', u'ううう', u'ええええ'],
                       index=[u'ああ', u'いいいい', u'う', u'えええ'])
            expected = (u"ああ              あ\nいいいい        いい\nう            ううう\n"
                        u"えええ      ええええ\ndtype: object")
            self.assertEqual(_rep(s), expected)

            # unicode footer
            s = Series([u'あ', u'いい', u'ううう', u'ええええ'],
                       index=[u'ああ', u'いいいい', u'う', u'えええ'], name=u'おおおおおおお')
            expected = (u"ああ              あ\nいいいい        いい\nう            ううう\n"
                        u"えええ      ええええ\nName: おおおおおおお, dtype: object")
            self.assertEqual(_rep(s), expected)

            # MultiIndex
            idx = pd.MultiIndex.from_tuples([(u'あ', u'いい'), (u'う', u'え'), (
                u'おおお', u'かかかか'), (u'き', u'くく')])
            s = Series([1, 22, 3333, 44444], index=idx)
            expected = (u"あ      いい            1\nう      え             22\nおおお  かかかか     3333\n"
                        u"き      くく        44444\ndtype: int64")
            self.assertEqual(_rep(s), expected)

            # object dtype, shorter than unicode repr
            s = Series([1, 22, 3333, 44444], index=[1, 'AB', np.nan, u'あああ'])
            expected = (u"1             1\nAB           22\nNaN        3333\n"
                        u"あああ    44444\ndtype: int64")
            self.assertEqual(_rep(s), expected)

            # object dtype, longer than unicode repr
            s = Series([1, 22, 3333, 44444],
                       index=[1, 'AB', pd.Timestamp('2011-01-01'), u'あああ'])
            expected = (u"1                          1\nAB                        22\n"
                        u"2011-01-01 00:00:00     3333\nあああ                 44444\ndtype: int64"
                        )
            self.assertEqual(_rep(s), expected)

            # truncate
            with option_context('display.max_rows', 3):
                s = Series([u'あ', u'いい', u'ううう', u'ええええ'], name=u'おおおおおおお')
                expected = (u"0          あ\n       ...   \n"
                            u"3    ええええ\nName: おおおおおおお, dtype: object")
                self.assertEqual(_rep(s), expected)

                s.index = [u'ああ', u'いいいい', u'う', u'えええ']
                expected = (u"ああ            あ\n            ...   \n"
                            u"えええ    ええええ\nName: おおおおおおお, dtype: object")
                self.assertEqual(_rep(s), expected)

            # ambiguous unicode
            s = Series([u'¡¡', u'い¡¡', u'ううう', u'ええええ'],
                       index=[u'ああ', u'¡¡¡¡いい', u'¡¡', u'えええ'])
            expected = (u"ああ              ¡¡\n¡¡¡¡いい        い¡¡\n¡¡            ううう\n"
                        u"えええ      ええええ\ndtype: object")
            self.assertEqual(_rep(s), expected)

    def test_float_trim_zeros(self):
        vals = [2.08430917305e+10, 3.52205017305e+10, 2.30674817305e+10,
                2.03954217305e+10, 5.59897817305e+10]
        for line in repr(Series(vals)).split('\n'):
            if line.startswith('dtype:'):
                continue
            if _three_digit_exp():
                self.assertIn('+010', line)
            else:
                self.assertIn('+10', line)

    def test_datetimeindex(self):

        index = date_range('20130102', periods=6)
        s = Series(1, index=index)
        result = s.to_string()
        self.assertTrue('2013-01-02' in result)

        # nat in index
        s2 = Series(2, index=[Timestamp('20130111'), NaT])
        s = s2.append(s)
        result = s.to_string()
        self.assertTrue('NaT' in result)

        # nat in summary
        result = str(s2.index)
        self.assertTrue('NaT' in result)

    def test_timedelta64(self):

        from datetime import datetime, timedelta

        Series(np.array([1100, 20], dtype='timedelta64[ns]')).to_string()

        s = Series(date_range('2012-1-1', periods=3, freq='D'))

        # GH2146

        # adding NaTs
        y = s - s.shift(1)
        result = y.to_string()
        self.assertTrue('1 days' in result)
        self.assertTrue('00:00:00' not in result)
        self.assertTrue('NaT' in result)

        # with frac seconds
        o = Series([datetime(2012, 1, 1, microsecond=150)] * 3)
        y = s - o
        result = y.to_string()
        self.assertTrue('-1 days +23:59:59.999850' in result)

        # rounding?
        o = Series([datetime(2012, 1, 1, 1)] * 3)
        y = s - o
        result = y.to_string()
        self.assertTrue('-1 days +23:00:00' in result)
        self.assertTrue('1 days 23:00:00' in result)

        o = Series([datetime(2012, 1, 1, 1, 1)] * 3)
        y = s - o
        result = y.to_string()
        self.assertTrue('-1 days +22:59:00' in result)
        self.assertTrue('1 days 22:59:00' in result)

        o = Series([datetime(2012, 1, 1, 1, 1, microsecond=150)] * 3)
        y = s - o
        result = y.to_string()
        self.assertTrue('-1 days +22:58:59.999850' in result)
        self.assertTrue('0 days 22:58:59.999850' in result)

        # neg time
        td = timedelta(minutes=5, seconds=3)
        s2 = Series(date_range('2012-1-1', periods=3, freq='D')) + td
        y = s - s2
        result = y.to_string()
        self.assertTrue('-1 days +23:54:57' in result)

        td = timedelta(microseconds=550)
        s2 = Series(date_range('2012-1-1', periods=3, freq='D')) + td
        y = s - td
        result = y.to_string()
        self.assertTrue('2012-01-01 23:59:59.999450' in result)

        # no boxing of the actual elements
        td = Series(pd.timedelta_range('1 days', periods=3))
        result = td.to_string()
        self.assertEqual(result, u("0   1 days\n1   2 days\n2   3 days"))

    def test_mixed_datetime64(self):
        df = DataFrame({'A': [1, 2], 'B': ['2012-01-01', '2012-01-02']})
        df['B'] = pd.to_datetime(df.B)

        result = repr(df.ix[0])
        self.assertTrue('2012-01-01' in result)

    def test_period(self):
        # GH 12615
        index = pd.period_range('2013-01', periods=6, freq='M')
        s = Series(np.arange(6, dtype='int64'), index=index)
        exp = ("2013-01    0\n2013-02    1\n2013-03    2\n2013-04    3\n"
               "2013-05    4\n2013-06    5\nFreq: M, dtype: int64")
        self.assertEqual(str(s), exp)

        s = Series(index)
        exp = ("0   2013-01\n1   2013-02\n2   2013-03\n3   2013-04\n"
               "4   2013-05\n5   2013-06\ndtype: object")
        self.assertEqual(str(s), exp)

        # periods with mixed freq
        s = Series([pd.Period('2011-01', freq='M'),
                    pd.Period('2011-02-01', freq='D'),
                    pd.Period('2011-03-01 09:00', freq='H')])
        exp = ("0            2011-01\n1         2011-02-01\n"
               "2   2011-03-01 09:00\ndtype: object")
        self.assertEqual(str(s), exp)

    def test_max_multi_index_display(self):
        # GH 7101

        # doc example (indexing.rst)

        # multi-index
        arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                  ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        tuples = list(zip(*arrays))
        index = MultiIndex.from_tuples(tuples, names=['first', 'second'])
        s = Series(randn(8), index=index)

        with option_context("display.max_rows", 10):
            self.assertEqual(len(str(s).split('\n')), 10)
        with option_context("display.max_rows", 3):
            self.assertEqual(len(str(s).split('\n')), 5)
        with option_context("display.max_rows", 2):
            self.assertEqual(len(str(s).split('\n')), 5)
        with option_context("display.max_rows", 1):
            self.assertEqual(len(str(s).split('\n')), 4)
        with option_context("display.max_rows", 0):
            self.assertEqual(len(str(s).split('\n')), 10)

        # index
        s = Series(randn(8), None)

        with option_context("display.max_rows", 10):
            self.assertEqual(len(str(s).split('\n')), 9)
        with option_context("display.max_rows", 3):
            self.assertEqual(len(str(s).split('\n')), 4)
        with option_context("display.max_rows", 2):
            self.assertEqual(len(str(s).split('\n')), 4)
        with option_context("display.max_rows", 1):
            self.assertEqual(len(str(s).split('\n')), 3)
        with option_context("display.max_rows", 0):
            self.assertEqual(len(str(s).split('\n')), 9)

    # Make sure #8532 is fixed
    def test_consistent_format(self):
        s = pd.Series([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.9999, 1, 1] * 10)
        with option_context("display.max_rows", 10):
            res = repr(s)
        exp = ('0      1.0000\n1      1.0000\n2      1.0000\n3      '
               '1.0000\n4      1.0000\n        ...  \n125    '
               '1.0000\n126    1.0000\n127    0.9999\n128    '
               '1.0000\n129    1.0000\ndtype: float64')
        self.assertEqual(res, exp)

    @staticmethod
    def gen_test_series():
        s1 = pd.Series(['a'] * 100)
        s2 = pd.Series(['ab'] * 100)
        s3 = pd.Series(['a', 'ab', 'abc', 'abcd', 'abcde', 'abcdef'])
        s4 = s3[::-1]
        test_sers = {'onel': s1, 'twol': s2, 'asc': s3, 'desc': s4}
        return test_sers

    def chck_ncols(self, s):
        with option_context("display.max_rows", 10):
            res = repr(s)
        lines = res.split('\n')
        lines = [line for line in repr(s).split('\n')
                 if not re.match('[^\.]*\.+', line)][:-1]
        ncolsizes = len(set(len(line.strip()) for line in lines))
        self.assertEqual(ncolsizes, 1)

    def test_format_explicit(self):
        test_sers = self.gen_test_series()
        with option_context("display.max_rows", 4):
            res = repr(test_sers['onel'])
            exp = '0     a\n1     a\n     ..\n98    a\n99    a\ndtype: object'
            self.assertEqual(exp, res)
            res = repr(test_sers['twol'])
            exp = ('0     ab\n1     ab\n      ..\n98    ab\n99    ab\ndtype:'
                   ' object')
            self.assertEqual(exp, res)
            res = repr(test_sers['asc'])
            exp = ('0         a\n1        ab\n      ...  \n4     abcde\n5'
                   '    abcdef\ndtype: object')
            self.assertEqual(exp, res)
            res = repr(test_sers['desc'])
            exp = ('5    abcdef\n4     abcde\n      ...  \n1        ab\n0'
                   '         a\ndtype: object')
            self.assertEqual(exp, res)

    def test_ncols(self):
        test_sers = self.gen_test_series()
        for s in test_sers.values():
            self.chck_ncols(s)

    def test_max_rows_eq_one(self):
        s = Series(range(10), dtype='int64')
        with option_context("display.max_rows", 1):
            strrepr = repr(s).split('\n')
        exp1 = ['0', '0']
        res1 = strrepr[0].split()
        self.assertEqual(exp1, res1)
        exp2 = ['..']
        res2 = strrepr[1].split()
        self.assertEqual(exp2, res2)

    def test_truncate_ndots(self):
        def getndots(s):
            return len(re.match('[^\.]*(\.*)', s).groups()[0])

        s = Series([0, 2, 3, 6])
        with option_context("display.max_rows", 2):
            strrepr = repr(s).replace('\n', '')
        self.assertEqual(getndots(strrepr), 2)

        s = Series([0, 100, 200, 400])
        with option_context("display.max_rows", 2):
            strrepr = repr(s).replace('\n', '')
        self.assertEqual(getndots(strrepr), 3)

    def test_to_string_name(self):
        s = Series(range(100), dtype='int64')
        s.name = 'myser'
        res = s.to_string(max_rows=2, name=True)
        exp = '0      0\n      ..\n99    99\nName: myser'
        self.assertEqual(res, exp)
        res = s.to_string(max_rows=2, name=False)
        exp = '0      0\n      ..\n99    99'
        self.assertEqual(res, exp)

    def test_to_string_dtype(self):
        s = Series(range(100), dtype='int64')
        res = s.to_string(max_rows=2, dtype=True)
        exp = '0      0\n      ..\n99    99\ndtype: int64'
        self.assertEqual(res, exp)
        res = s.to_string(max_rows=2, dtype=False)
        exp = '0      0\n      ..\n99    99'
        self.assertEqual(res, exp)

    def test_to_string_length(self):
        s = Series(range(100), dtype='int64')
        res = s.to_string(max_rows=2, length=True)
        exp = '0      0\n      ..\n99    99\nLength: 100'
        self.assertEqual(res, exp)

    def test_to_string_na_rep(self):
        s = pd.Series(index=range(100))
        res = s.to_string(na_rep='foo', max_rows=2)
        exp = '0    foo\n      ..\n99   foo'
        self.assertEqual(res, exp)

    def test_to_string_float_format(self):
        s = pd.Series(range(10), dtype='float64')
        res = s.to_string(float_format=lambda x: '{0:2.1f}'.format(x),
                          max_rows=2)
        exp = '0   0.0\n     ..\n9   9.0'
        self.assertEqual(res, exp)

    def test_to_string_header(self):
        s = pd.Series(range(10), dtype='int64')
        s.index.name = 'foo'
        res = s.to_string(header=True, max_rows=2)
        exp = 'foo\n0    0\n    ..\n9    9'
        self.assertEqual(res, exp)
        res = s.to_string(header=False, max_rows=2)
        exp = '0    0\n    ..\n9    9'
        self.assertEqual(res, exp)

    def test_sparse_max_row(self):
        s = pd.Series([1, np.nan, np.nan, 3, np.nan]).to_sparse()
        result = repr(s)
        dtype = '' if use_32bit_repr else ', dtype=int32'
        exp = ("0    1.0\n1    NaN\n2    NaN\n3    3.0\n"
               "4    NaN\ndtype: float64\nBlockIndex\n"
               "Block locations: array([0, 3]{0})\n"
               "Block lengths: array([1, 1]{0})".format(dtype))
        self.assertEqual(result, exp)

        with option_context("display.max_rows", 3):
            # GH 10560
            result = repr(s)
            exp = ("0    1.0\n    ... \n4    NaN\n"
                   "dtype: float64\nBlockIndex\n"
                   "Block locations: array([0, 3]{0})\n"
                   "Block lengths: array([1, 1]{0})".format(dtype))
            self.assertEqual(result, exp)


class TestEngFormatter(tm.TestCase):
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

        self.reset_display_options()

    def compare(self, formatter, input, output):
        formatted_input = formatter(input)
        msg = ("formatting of %s results in '%s', expected '%s'" %
               (str(input), formatted_input, output))
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
        in_out = [(f * 10 ** -24, " 1.414y"), (f * 10 ** -23, " 14.142y"),
                  (f * 10 ** -22, " 141.421y"), (f * 10 ** -21, " 1.414z"),
                  (f * 10 ** -20, " 14.142z"), (f * 10 ** -19, " 141.421z"),
                  (f * 10 ** -18, " 1.414a"), (f * 10 ** -17, " 14.142a"),
                  (f * 10 ** -16, " 141.421a"), (f * 10 ** -15, " 1.414f"),
                  (f * 10 ** -14, " 14.142f"), (f * 10 ** -13, " 141.421f"),
                  (f * 10 ** -12, " 1.414p"), (f * 10 ** -11, " 14.142p"),
                  (f * 10 ** -10, " 141.421p"), (f * 10 ** -9, " 1.414n"),
                  (f * 10 ** -8, " 14.142n"), (f * 10 ** -7, " 141.421n"),
                  (f * 10 ** -6, " 1.414u"), (f * 10 ** -5, " 14.142u"),
                  (f * 10 ** -4, " 141.421u"), (f * 10 ** -3, " 1.414m"),
                  (f * 10 ** -2, " 14.142m"), (f * 10 ** -1, " 141.421m"),
                  (f * 10 ** 0, " 1.414"), (f * 10 ** 1, " 14.142"),
                  (f * 10 ** 2, " 141.421"), (f * 10 ** 3, " 1.414k"),
                  (f * 10 ** 4, " 14.142k"), (f * 10 ** 5, " 141.421k"),
                  (f * 10 ** 6, " 1.414M"), (f * 10 ** 7, " 14.142M"),
                  (f * 10 ** 8, " 141.421M"), (f * 10 ** 9, " 1.414G"), (
                      f * 10 ** 10, " 14.142G"), (f * 10 ** 11, " 141.421G"),
                  (f * 10 ** 12, " 1.414T"), (f * 10 ** 13, " 14.142T"), (
                      f * 10 ** 14, " 141.421T"), (f * 10 ** 15, " 1.414P"), (
                          f * 10 ** 16, " 14.142P"), (f * 10 ** 17, " 141.421P"), (
                              f * 10 ** 18, " 1.414E"), (f * 10 ** 19, " 14.142E"),
                  (f * 10 ** 20, " 141.421E"), (f * 10 ** 21, " 1.414Z"), (
                      f * 10 ** 22, " 14.142Z"), (f * 10 ** 23, " 141.421Z"), (
                          f * 10 ** 24, " 1.414Y"), (f * 10 ** 25, " 14.142Y"), (
                              f * 10 ** 26, " 141.421Y")]
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
                  (f * 10 ** -9, " 3.1416E-09"), (f * 10 ** -8, " 31.4159E-09"),
                  (f * 10 ** -7, " 314.1593E-09"), (f * 10 ** -6, " 3.1416E-06"),
                  (f * 10 ** -5, " 31.4159E-06"), (f * 10 ** -4,
                                                   " 314.1593E-06"),
                  (f * 10 ** -3, " 3.1416E-03"), (f * 10 ** -2, " 31.4159E-03"),
                  (f * 10 ** -1, " 314.1593E-03"), (f * 10 ** 0, " 3.1416E+00"), (
                      f * 10 ** 1, " 31.4159E+00"), (f * 10 ** 2, " 314.1593E+00"),
                  (f * 10 ** 3, " 3.1416E+03"), (f * 10 ** 4, " 31.4159E+03"), (
                      f * 10 ** 5, " 314.1593E+03"), (f * 10 ** 6, " 3.1416E+06"),
                  (f * 10 ** 7, " 31.4159E+06"), (f * 10 ** 8, " 314.1593E+06"), (
                      f * 10 ** 9, " 3.1416E+09"), (f * 10 ** 10, " 31.4159E+09"),
                  (f * 10 ** 11, " 314.1593E+09"), (f * 10 ** 12, " 3.1416E+12"),
                  (f * 10 ** 13, " 31.4159E+12"), (f * 10 ** 14, " 314.1593E+12"),
                  (f * 10 ** 15, " 3.1416E+15"), (f * 10 ** 16, " 31.4159E+15"),
                  (f * 10 ** 17, " 314.1593E+15"), (f * 10 ** 18, " 3.1416E+18"),
                  (f * 10 ** 19, " 31.4159E+18"), (f * 10 ** 20, " 314.1593E+18"),
                  (f * 10 ** 21, " 3.1416E+21"), (f * 10 ** 22, " 31.4159E+21"),
                  (f * 10 ** 23, " 314.1593E+21"), (f * 10 ** 24, " 3.1416E+24"),
                  (f * 10 ** 25, " 31.4159E+24"), (f * 10 ** 26, " 314.1593E+24")]
        self.compare_all(formatter, in_out)

    def test_rounding(self):
        formatter = fmt.EngFormatter(accuracy=3, use_eng_prefix=True)
        in_out = [(5.55555, ' 5.556'), (55.5555, ' 55.556'),
                  (555.555, ' 555.555'), (5555.55, ' 5.556k'),
                  (55555.5, ' 55.556k'), (555555, ' 555.555k')]
        self.compare_all(formatter, in_out)

        formatter = fmt.EngFormatter(accuracy=1, use_eng_prefix=True)
        in_out = [(5.55555, ' 5.6'), (55.5555, ' 55.6'), (555.555, ' 555.6'),
                  (5555.55, ' 5.6k'), (55555.5, ' 55.6k'), (555555, ' 555.6k')]
        self.compare_all(formatter, in_out)

        formatter = fmt.EngFormatter(accuracy=0, use_eng_prefix=True)
        in_out = [(5.55555, ' 6'), (55.5555, ' 56'), (555.555, ' 556'),
                  (5555.55, ' 6k'), (55555.5, ' 56k'), (555555, ' 556k')]
        self.compare_all(formatter, in_out)

        formatter = fmt.EngFormatter(accuracy=3, use_eng_prefix=True)
        result = formatter(0)
        self.assertEqual(result, u(' 0.000'))


def _three_digit_exp():
    return '%.4g' % 1.7e8 == '1.7e+008'


class TestFloatArrayFormatter(tm.TestCase):

    def test_misc(self):
        obj = fmt.FloatArrayFormatter(np.array([], dtype=np.float64))
        result = obj.get_result()
        self.assertTrue(len(result) == 0)

    def test_format(self):
        obj = fmt.FloatArrayFormatter(np.array([12, 0], dtype=np.float64))
        result = obj.get_result()
        self.assertEqual(result[0], " 12.0")
        self.assertEqual(result[1], "  0.0")

    def test_output_significant_digits(self):
        # Issue #9764

        # In case default display precision changes:
        with pd.option_context('display.precision', 6):
            # DataFrame example from issue #9764
            d = pd.DataFrame(
                {'col1': [9.999e-8, 1e-7, 1.0001e-7, 2e-7, 4.999e-7, 5e-7,
                          5.0001e-7, 6e-7, 9.999e-7, 1e-6, 1.0001e-6, 2e-6,
                          4.999e-6, 5e-6, 5.0001e-6, 6e-6]})

            expected_output = {
                (0, 6):
                '           col1\n0  9.999000e-08\n1  1.000000e-07\n2  1.000100e-07\n3  2.000000e-07\n4  4.999000e-07\n5  5.000000e-07',
                (1, 6):
                '           col1\n1  1.000000e-07\n2  1.000100e-07\n3  2.000000e-07\n4  4.999000e-07\n5  5.000000e-07',
                (1, 8):
                '           col1\n1  1.000000e-07\n2  1.000100e-07\n3  2.000000e-07\n4  4.999000e-07\n5  5.000000e-07\n6  5.000100e-07\n7  6.000000e-07',
                (8, 16):
                '            col1\n8   9.999000e-07\n9   1.000000e-06\n10  1.000100e-06\n11  2.000000e-06\n12  4.999000e-06\n13  5.000000e-06\n14  5.000100e-06\n15  6.000000e-06',
                (9, 16):
                '        col1\n9   0.000001\n10  0.000001\n11  0.000002\n12  0.000005\n13  0.000005\n14  0.000005\n15  0.000006'
            }

            for (start, stop), v in expected_output.items():
                self.assertEqual(str(d[start:stop]), v)

    def test_too_long(self):
        # GH 10451
        with pd.option_context('display.precision', 4):
            # need both a number > 1e6 and something that normally formats to
            # having length > display.precision + 6
            df = pd.DataFrame(dict(x=[12345.6789]))
            self.assertEqual(str(df), '            x\n0  12345.6789')
            df = pd.DataFrame(dict(x=[2e6]))
            self.assertEqual(str(df), '           x\n0  2000000.0')
            df = pd.DataFrame(dict(x=[12345.6789, 2e6]))
            self.assertEqual(
                str(df), '            x\n0  1.2346e+04\n1  2.0000e+06')


class TestRepr_timedelta64(tm.TestCase):

    def test_none(self):
        delta_1d = pd.to_timedelta(1, unit='D')
        delta_0d = pd.to_timedelta(0, unit='D')
        delta_1s = pd.to_timedelta(1, unit='s')
        delta_500ms = pd.to_timedelta(500, unit='ms')

        drepr = lambda x: x._repr_base()
        self.assertEqual(drepr(delta_1d), "1 days")
        self.assertEqual(drepr(-delta_1d), "-1 days")
        self.assertEqual(drepr(delta_0d), "0 days")
        self.assertEqual(drepr(delta_1s), "0 days 00:00:01")
        self.assertEqual(drepr(delta_500ms), "0 days 00:00:00.500000")
        self.assertEqual(drepr(delta_1d + delta_1s), "1 days 00:00:01")
        self.assertEqual(
            drepr(delta_1d + delta_500ms), "1 days 00:00:00.500000")

    def test_even_day(self):
        delta_1d = pd.to_timedelta(1, unit='D')
        delta_0d = pd.to_timedelta(0, unit='D')
        delta_1s = pd.to_timedelta(1, unit='s')
        delta_500ms = pd.to_timedelta(500, unit='ms')

        drepr = lambda x: x._repr_base(format='even_day')
        self.assertEqual(drepr(delta_1d), "1 days")
        self.assertEqual(drepr(-delta_1d), "-1 days")
        self.assertEqual(drepr(delta_0d), "0 days")
        self.assertEqual(drepr(delta_1s), "0 days 00:00:01")
        self.assertEqual(drepr(delta_500ms), "0 days 00:00:00.500000")
        self.assertEqual(drepr(delta_1d + delta_1s), "1 days 00:00:01")
        self.assertEqual(
            drepr(delta_1d + delta_500ms), "1 days 00:00:00.500000")

    def test_sub_day(self):
        delta_1d = pd.to_timedelta(1, unit='D')
        delta_0d = pd.to_timedelta(0, unit='D')
        delta_1s = pd.to_timedelta(1, unit='s')
        delta_500ms = pd.to_timedelta(500, unit='ms')

        drepr = lambda x: x._repr_base(format='sub_day')
        self.assertEqual(drepr(delta_1d), "1 days")
        self.assertEqual(drepr(-delta_1d), "-1 days")
        self.assertEqual(drepr(delta_0d), "00:00:00")
        self.assertEqual(drepr(delta_1s), "00:00:01")
        self.assertEqual(drepr(delta_500ms), "00:00:00.500000")
        self.assertEqual(drepr(delta_1d + delta_1s), "1 days 00:00:01")
        self.assertEqual(
            drepr(delta_1d + delta_500ms), "1 days 00:00:00.500000")

    def test_long(self):
        delta_1d = pd.to_timedelta(1, unit='D')
        delta_0d = pd.to_timedelta(0, unit='D')
        delta_1s = pd.to_timedelta(1, unit='s')
        delta_500ms = pd.to_timedelta(500, unit='ms')

        drepr = lambda x: x._repr_base(format='long')
        self.assertEqual(drepr(delta_1d), "1 days 00:00:00")
        self.assertEqual(drepr(-delta_1d), "-1 days +00:00:00")
        self.assertEqual(drepr(delta_0d), "0 days 00:00:00")
        self.assertEqual(drepr(delta_1s), "0 days 00:00:01")
        self.assertEqual(drepr(delta_500ms), "0 days 00:00:00.500000")
        self.assertEqual(drepr(delta_1d + delta_1s), "1 days 00:00:01")
        self.assertEqual(
            drepr(delta_1d + delta_500ms), "1 days 00:00:00.500000")

    def test_all(self):
        delta_1d = pd.to_timedelta(1, unit='D')
        delta_0d = pd.to_timedelta(0, unit='D')
        delta_1ns = pd.to_timedelta(1, unit='ns')

        drepr = lambda x: x._repr_base(format='all')
        self.assertEqual(drepr(delta_1d), "1 days 00:00:00.000000000")
        self.assertEqual(drepr(delta_0d), "0 days 00:00:00.000000000")
        self.assertEqual(drepr(delta_1ns), "0 days 00:00:00.000000001")


class TestTimedelta64Formatter(tm.TestCase):

    def test_days(self):
        x = pd.to_timedelta(list(range(5)) + [pd.NaT], unit='D')
        result = fmt.Timedelta64Formatter(x, box=True).get_result()
        self.assertEqual(result[0].strip(), "'0 days'")
        self.assertEqual(result[1].strip(), "'1 days'")

        result = fmt.Timedelta64Formatter(x[1:2], box=True).get_result()
        self.assertEqual(result[0].strip(), "'1 days'")

        result = fmt.Timedelta64Formatter(x, box=False).get_result()
        self.assertEqual(result[0].strip(), "0 days")
        self.assertEqual(result[1].strip(), "1 days")

        result = fmt.Timedelta64Formatter(x[1:2], box=False).get_result()
        self.assertEqual(result[0].strip(), "1 days")

    def test_days_neg(self):
        x = pd.to_timedelta(list(range(5)) + [pd.NaT], unit='D')
        result = fmt.Timedelta64Formatter(-x, box=True).get_result()
        self.assertEqual(result[0].strip(), "'0 days'")
        self.assertEqual(result[1].strip(), "'-1 days'")

    def test_subdays(self):
        y = pd.to_timedelta(list(range(5)) + [pd.NaT], unit='s')
        result = fmt.Timedelta64Formatter(y, box=True).get_result()
        self.assertEqual(result[0].strip(), "'00:00:00'")
        self.assertEqual(result[1].strip(), "'00:00:01'")

    def test_subdays_neg(self):
        y = pd.to_timedelta(list(range(5)) + [pd.NaT], unit='s')
        result = fmt.Timedelta64Formatter(-y, box=True).get_result()
        self.assertEqual(result[0].strip(), "'00:00:00'")
        self.assertEqual(result[1].strip(), "'-1 days +23:59:59'")

    def test_zero(self):
        x = pd.to_timedelta(list(range(1)) + [pd.NaT], unit='D')
        result = fmt.Timedelta64Formatter(x, box=True).get_result()
        self.assertEqual(result[0].strip(), "'0 days'")

        x = pd.to_timedelta(list(range(1)), unit='D')
        result = fmt.Timedelta64Formatter(x, box=True).get_result()
        self.assertEqual(result[0].strip(), "'0 days'")


class TestDatetime64Formatter(tm.TestCase):

    def test_mixed(self):
        x = Series([datetime(2013, 1, 1), datetime(2013, 1, 1, 12), pd.NaT])
        result = fmt.Datetime64Formatter(x).get_result()
        self.assertEqual(result[0].strip(), "2013-01-01 00:00:00")
        self.assertEqual(result[1].strip(), "2013-01-01 12:00:00")

    def test_dates(self):
        x = Series([datetime(2013, 1, 1), datetime(2013, 1, 2), pd.NaT])
        result = fmt.Datetime64Formatter(x).get_result()
        self.assertEqual(result[0].strip(), "2013-01-01")
        self.assertEqual(result[1].strip(), "2013-01-02")

    def test_date_nanos(self):
        x = Series([Timestamp(200)])
        result = fmt.Datetime64Formatter(x).get_result()
        self.assertEqual(result[0].strip(), "1970-01-01 00:00:00.000000200")

    def test_dates_display(self):

        # 10170
        # make sure that we are consistently display date formatting
        x = Series(date_range('20130101 09:00:00', periods=5, freq='D'))
        x.iloc[1] = np.nan
        result = fmt.Datetime64Formatter(x).get_result()
        self.assertEqual(result[0].strip(), "2013-01-01 09:00:00")
        self.assertEqual(result[1].strip(), "NaT")
        self.assertEqual(result[4].strip(), "2013-01-05 09:00:00")

        x = Series(date_range('20130101 09:00:00', periods=5, freq='s'))
        x.iloc[1] = np.nan
        result = fmt.Datetime64Formatter(x).get_result()
        self.assertEqual(result[0].strip(), "2013-01-01 09:00:00")
        self.assertEqual(result[1].strip(), "NaT")
        self.assertEqual(result[4].strip(), "2013-01-01 09:00:04")

        x = Series(date_range('20130101 09:00:00', periods=5, freq='ms'))
        x.iloc[1] = np.nan
        result = fmt.Datetime64Formatter(x).get_result()
        self.assertEqual(result[0].strip(), "2013-01-01 09:00:00.000")
        self.assertEqual(result[1].strip(), "NaT")
        self.assertEqual(result[4].strip(), "2013-01-01 09:00:00.004")

        x = Series(date_range('20130101 09:00:00', periods=5, freq='us'))
        x.iloc[1] = np.nan
        result = fmt.Datetime64Formatter(x).get_result()
        self.assertEqual(result[0].strip(), "2013-01-01 09:00:00.000000")
        self.assertEqual(result[1].strip(), "NaT")
        self.assertEqual(result[4].strip(), "2013-01-01 09:00:00.000004")

        x = Series(date_range('20130101 09:00:00', periods=5, freq='N'))
        x.iloc[1] = np.nan
        result = fmt.Datetime64Formatter(x).get_result()
        self.assertEqual(result[0].strip(), "2013-01-01 09:00:00.000000000")
        self.assertEqual(result[1].strip(), "NaT")
        self.assertEqual(result[4].strip(), "2013-01-01 09:00:00.000000004")


class TestNaTFormatting(tm.TestCase):

    def test_repr(self):
        self.assertEqual(repr(pd.NaT), "NaT")

    def test_str(self):
        self.assertEqual(str(pd.NaT), "NaT")


class TestDatetimeIndexFormat(tm.TestCase):

    def test_datetime(self):
        formatted = pd.to_datetime([datetime(2003, 1, 1, 12), pd.NaT]).format()
        self.assertEqual(formatted[0], "2003-01-01 12:00:00")
        self.assertEqual(formatted[1], "NaT")

    def test_date(self):
        formatted = pd.to_datetime([datetime(2003, 1, 1), pd.NaT]).format()
        self.assertEqual(formatted[0], "2003-01-01")
        self.assertEqual(formatted[1], "NaT")

    def test_date_tz(self):
        formatted = pd.to_datetime([datetime(2013, 1, 1)], utc=True).format()
        self.assertEqual(formatted[0], "2013-01-01 00:00:00+00:00")

        formatted = pd.to_datetime(
            [datetime(2013, 1, 1), pd.NaT], utc=True).format()
        self.assertEqual(formatted[0], "2013-01-01 00:00:00+00:00")

    def test_date_explict_date_format(self):
        formatted = pd.to_datetime([datetime(2003, 2, 1), pd.NaT]).format(
            date_format="%m-%d-%Y", na_rep="UT")
        self.assertEqual(formatted[0], "02-01-2003")
        self.assertEqual(formatted[1], "UT")


class TestDatetimeIndexUnicode(tm.TestCase):

    def test_dates(self):
        text = str(pd.to_datetime([datetime(2013, 1, 1), datetime(2014, 1, 1)
                                   ]))
        self.assertTrue("['2013-01-01'," in text)
        self.assertTrue(", '2014-01-01']" in text)

    def test_mixed(self):
        text = str(pd.to_datetime([datetime(2013, 1, 1), datetime(
            2014, 1, 1, 12), datetime(2014, 1, 1)]))
        self.assertTrue("'2013-01-01 00:00:00'," in text)
        self.assertTrue("'2014-01-01 00:00:00']" in text)


class TestStringRepTimestamp(tm.TestCase):

    def test_no_tz(self):
        dt_date = datetime(2013, 1, 2)
        self.assertEqual(str(dt_date), str(Timestamp(dt_date)))

        dt_datetime = datetime(2013, 1, 2, 12, 1, 3)
        self.assertEqual(str(dt_datetime), str(Timestamp(dt_datetime)))

        dt_datetime_us = datetime(2013, 1, 2, 12, 1, 3, 45)
        self.assertEqual(str(dt_datetime_us), str(Timestamp(dt_datetime_us)))

        ts_nanos_only = Timestamp(200)
        self.assertEqual(str(ts_nanos_only), "1970-01-01 00:00:00.000000200")

        ts_nanos_micros = Timestamp(1200)
        self.assertEqual(str(ts_nanos_micros), "1970-01-01 00:00:00.000001200")

    def test_tz_pytz(self):
        tm._skip_if_no_pytz()

        import pytz

        dt_date = datetime(2013, 1, 2, tzinfo=pytz.utc)
        self.assertEqual(str(dt_date), str(Timestamp(dt_date)))

        dt_datetime = datetime(2013, 1, 2, 12, 1, 3, tzinfo=pytz.utc)
        self.assertEqual(str(dt_datetime), str(Timestamp(dt_datetime)))

        dt_datetime_us = datetime(2013, 1, 2, 12, 1, 3, 45, tzinfo=pytz.utc)
        self.assertEqual(str(dt_datetime_us), str(Timestamp(dt_datetime_us)))

    def test_tz_dateutil(self):
        tm._skip_if_no_dateutil()
        import dateutil
        utc = dateutil.tz.tzutc()

        dt_date = datetime(2013, 1, 2, tzinfo=utc)
        self.assertEqual(str(dt_date), str(Timestamp(dt_date)))

        dt_datetime = datetime(2013, 1, 2, 12, 1, 3, tzinfo=utc)
        self.assertEqual(str(dt_datetime), str(Timestamp(dt_datetime)))

        dt_datetime_us = datetime(2013, 1, 2, 12, 1, 3, 45, tzinfo=utc)
        self.assertEqual(str(dt_datetime_us), str(Timestamp(dt_datetime_us)))

    def test_nat_representations(self):
        for f in (str, repr, methodcaller('isoformat')):
            self.assertEqual(f(pd.NaT), 'NaT')


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
