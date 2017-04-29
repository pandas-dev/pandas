# -*- coding: utf-8 -*-
import pytest

import numpy as np
import pandas as pd

from pandas import compat
from pandas.errors import UnserializableWarning
import pandas.io.formats.printing as printing
import pandas.io.formats.format as fmt
import pandas.util.testing as tm
import pandas.core.config as cf


def test_adjoin():
    data = [['a', 'b', 'c'], ['dd', 'ee', 'ff'], ['ggg', 'hhh', 'iii']]
    expected = 'a  dd  ggg\nb  ee  hhh\nc  ff  iii'

    adjoined = printing.adjoin(2, *data)

    assert (adjoined == expected)


def test_repr_binary_type():
    import string
    letters = string.ascii_letters
    btype = compat.binary_type
    try:
        raw = btype(letters, encoding=cf.get_option('display.encoding'))
    except TypeError:
        raw = btype(letters)
    b = compat.text_type(compat.bytes_to_str(raw))
    res = printing.pprint_thing(b, quote_strings=True)
    assert res == repr(b)
    res = printing.pprint_thing(b, quote_strings=False)
    assert res == b


class TestFormattBase(tm.TestCase):

    def test_adjoin(self):
        data = [['a', 'b', 'c'], ['dd', 'ee', 'ff'], ['ggg', 'hhh', 'iii']]
        expected = 'a  dd  ggg\nb  ee  hhh\nc  ff  iii'

        adjoined = printing.adjoin(2, *data)

        assert adjoined == expected

    def test_adjoin_unicode(self):
        data = [[u'あ', 'b', 'c'], ['dd', u'ええ', 'ff'], ['ggg', 'hhh', u'いいい']]
        expected = u'あ  dd  ggg\nb  ええ  hhh\nc  ff  いいい'
        adjoined = printing.adjoin(2, *data)
        assert adjoined == expected

        adj = fmt.EastAsianTextAdjustment()

        expected = u"""あ  dd    ggg
b   ええ  hhh
c   ff    いいい"""

        adjoined = adj.adjoin(2, *data)
        assert adjoined == expected
        cols = adjoined.split('\n')
        assert adj.len(cols[0]) == 13
        assert adj.len(cols[1]) == 13
        assert adj.len(cols[2]) == 16

        expected = u"""あ       dd         ggg
b        ええ       hhh
c        ff         いいい"""

        adjoined = adj.adjoin(7, *data)
        assert adjoined == expected
        cols = adjoined.split('\n')
        assert adj.len(cols[0]) == 23
        assert adj.len(cols[1]) == 23
        assert adj.len(cols[2]) == 26

    def test_justify(self):
        adj = fmt.EastAsianTextAdjustment()

        def just(x, *args, **kwargs):
            # wrapper to test single str
            return adj.justify([x], *args, **kwargs)[0]

        assert just('abc', 5, mode='left') == 'abc  '
        assert just('abc', 5, mode='center') == ' abc '
        assert just('abc', 5, mode='right') == '  abc'
        assert just(u'abc', 5, mode='left') == 'abc  '
        assert just(u'abc', 5, mode='center') == ' abc '
        assert just(u'abc', 5, mode='right') == '  abc'

        assert just(u'パンダ', 5, mode='left') == u'パンダ'
        assert just(u'パンダ', 5, mode='center') == u'パンダ'
        assert just(u'パンダ', 5, mode='right') == u'パンダ'

        assert just(u'パンダ', 10, mode='left') == u'パンダ    '
        assert just(u'パンダ', 10, mode='center') == u'  パンダ  '
        assert just(u'パンダ', 10, mode='right') == u'    パンダ'

    def test_east_asian_len(self):
        adj = fmt.EastAsianTextAdjustment()

        assert adj.len('abc') == 3
        assert adj.len(u'abc') == 3

        assert adj.len(u'パンダ') == 6
        assert adj.len(u'ﾊﾟﾝﾀﾞ') == 5
        assert adj.len(u'パンダpanda') == 11
        assert adj.len(u'ﾊﾟﾝﾀﾞpanda') == 10

    def test_ambiguous_width(self):
        adj = fmt.EastAsianTextAdjustment()
        assert adj.len(u'¡¡ab') == 4

        with cf.option_context('display.unicode.ambiguous_as_wide', True):
            adj = fmt.EastAsianTextAdjustment()
            assert adj.len(u'¡¡ab') == 6

        data = [[u'あ', 'b', 'c'], ['dd', u'ええ', 'ff'],
                ['ggg', u'¡¡ab', u'いいい']]
        expected = u'あ  dd    ggg \nb   ええ  ¡¡ab\nc   ff    いいい'
        adjoined = adj.adjoin(2, *data)
        assert adjoined == expected


class TestTableSchemaRepr(tm.TestCase):

    @classmethod
    def setUpClass(cls):
        pytest.importorskip('IPython')
        try:
            import mock
        except ImportError:
            try:
                from unittest import mock
            except ImportError:
                pytest.skip("Mock is not installed")
        cls.mock = mock

    def test_publishes(self):
        df = pd.DataFrame({"A": [1, 2]})
        objects = [df['A'], df, df]  # dataframe / series
        expected_keys = [
            {'text/plain', 'application/vnd.dataresource+json'},
            {'text/plain', 'text/html', 'application/vnd.dataresource+json'},
        ]

        make_patch = self.mock.patch('IPython.display.display')
        opt = pd.option_context('display.html.table_schema', True)
        for obj, expected in zip(objects, expected_keys):
            with opt, make_patch as mock_display:
                handle = obj._ipython_display_()
                assert mock_display.call_count == 1
                assert handle is None
                args, kwargs = mock_display.call_args
                arg, = args  # just one argument

            assert kwargs == {"raw": True}
            assert set(arg.keys()) == expected

        with_latex = pd.option_context('display.latex.repr', True)

        with opt, with_latex, make_patch as mock_display:
            handle = obj._ipython_display_()
            args, kwargs = mock_display.call_args
            arg, = args

        expected = {'text/plain', 'text/html', 'text/latex',
                    'application/vnd.dataresource+json'}
        assert set(arg.keys()) == expected

    def test_publishes_not_implemented(self):
        # column MultiIndex
        # GH 15996
        midx = pd.MultiIndex.from_product([['A', 'B'], ['a', 'b', 'c']])
        df = pd.DataFrame(np.random.randn(5, len(midx)), columns=midx)

        make_patch = self.mock.patch('IPython.display.display')
        opt = pd.option_context('display.html.table_schema', True)
        with opt, make_patch as mock_display:
            with pytest.warns(UnserializableWarning) as record:
                df._ipython_display_()
                args, _ = mock_display.call_args
                arg, = args  # just one argument

        expected = {'text/plain', 'text/html'}
        assert set(arg.keys()) == expected
        assert "orient='table' is not supported for MultiIndex" in (
            record[-1].message.args[0])

    def test_config_on(self):
        df = pd.DataFrame({"A": [1, 2]})
        with pd.option_context("display.html.table_schema", True):
            result = df._repr_table_schema_()

        assert result is not None

    def test_config_default_off(self):
        df = pd.DataFrame({"A": [1, 2]})
        with pd.option_context("display.html.table_schema", False):
            result = df._repr_table_schema_()

        assert result is None


# TODO: fix this broken test

# def test_console_encode():
#     """
#     On Python 2, if sys.stdin.encoding is None (IPython with zmq frontend)
#     common.console_encode should encode things as utf-8.
#     """
#     if compat.PY3:
#         pytest.skip

#     with tm.stdin_encoding(encoding=None):
#         result = printing.console_encode(u"\u05d0")
#         expected = u"\u05d0".encode('utf-8')
#         assert (result == expected)
