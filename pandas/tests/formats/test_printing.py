# -*- coding: utf-8 -*-
import nose
from pandas import compat
import pandas.formats.printing as printing
import pandas.formats.format as fmt
import pandas.util.testing as tm
import pandas.core.config as cf

_multiprocess_can_split_ = True


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
    tm.assert_equal(res, repr(b))
    res = printing.pprint_thing(b, quote_strings=False)
    tm.assert_equal(res, b)


class TestFormattBase(tm.TestCase):

    def test_adjoin(self):
        data = [['a', 'b', 'c'], ['dd', 'ee', 'ff'], ['ggg', 'hhh', 'iii']]
        expected = 'a  dd  ggg\nb  ee  hhh\nc  ff  iii'

        adjoined = printing.adjoin(2, *data)

        self.assertEqual(adjoined, expected)

    def test_adjoin_unicode(self):
        data = [[u'あ', 'b', 'c'], ['dd', u'ええ', 'ff'], ['ggg', 'hhh', u'いいい']]
        expected = u'あ  dd  ggg\nb  ええ  hhh\nc  ff  いいい'
        adjoined = printing.adjoin(2, *data)
        self.assertEqual(adjoined, expected)

        adj = fmt.EastAsianTextAdjustment()

        expected = u"""あ  dd    ggg
b   ええ  hhh
c   ff    いいい"""

        adjoined = adj.adjoin(2, *data)
        self.assertEqual(adjoined, expected)
        cols = adjoined.split('\n')
        self.assertEqual(adj.len(cols[0]), 13)
        self.assertEqual(adj.len(cols[1]), 13)
        self.assertEqual(adj.len(cols[2]), 16)

        expected = u"""あ       dd         ggg
b        ええ       hhh
c        ff         いいい"""

        adjoined = adj.adjoin(7, *data)
        self.assertEqual(adjoined, expected)
        cols = adjoined.split('\n')
        self.assertEqual(adj.len(cols[0]), 23)
        self.assertEqual(adj.len(cols[1]), 23)
        self.assertEqual(adj.len(cols[2]), 26)

    def test_justify(self):
        adj = fmt.EastAsianTextAdjustment()

        def just(x, *args, **kwargs):
            # wrapper to test single str
            return adj.justify([x], *args, **kwargs)[0]

        self.assertEqual(just('abc', 5, mode='left'), 'abc  ')
        self.assertEqual(just('abc', 5, mode='center'), ' abc ')
        self.assertEqual(just('abc', 5, mode='right'), '  abc')
        self.assertEqual(just(u'abc', 5, mode='left'), 'abc  ')
        self.assertEqual(just(u'abc', 5, mode='center'), ' abc ')
        self.assertEqual(just(u'abc', 5, mode='right'), '  abc')

        self.assertEqual(just(u'パンダ', 5, mode='left'), u'パンダ')
        self.assertEqual(just(u'パンダ', 5, mode='center'), u'パンダ')
        self.assertEqual(just(u'パンダ', 5, mode='right'), u'パンダ')

        self.assertEqual(just(u'パンダ', 10, mode='left'), u'パンダ    ')
        self.assertEqual(just(u'パンダ', 10, mode='center'), u'  パンダ  ')
        self.assertEqual(just(u'パンダ', 10, mode='right'), u'    パンダ')

    def test_east_asian_len(self):
        adj = fmt.EastAsianTextAdjustment()

        self.assertEqual(adj.len('abc'), 3)
        self.assertEqual(adj.len(u'abc'), 3)

        self.assertEqual(adj.len(u'パンダ'), 6)
        self.assertEqual(adj.len(u'ﾊﾟﾝﾀﾞ'), 5)
        self.assertEqual(adj.len(u'パンダpanda'), 11)
        self.assertEqual(adj.len(u'ﾊﾟﾝﾀﾞpanda'), 10)

    def test_ambiguous_width(self):
        adj = fmt.EastAsianTextAdjustment()
        self.assertEqual(adj.len(u'¡¡ab'), 4)

        with cf.option_context('display.unicode.ambiguous_as_wide', True):
            adj = fmt.EastAsianTextAdjustment()
            self.assertEqual(adj.len(u'¡¡ab'), 6)

        data = [[u'あ', 'b', 'c'], ['dd', u'ええ', 'ff'],
                ['ggg', u'¡¡ab', u'いいい']]
        expected = u'あ  dd    ggg \nb   ええ  ¡¡ab\nc   ff    いいい'
        adjoined = adj.adjoin(2, *data)
        self.assertEqual(adjoined, expected)


# TODO: fix this broken test

# def test_console_encode():
#     """
#     On Python 2, if sys.stdin.encoding is None (IPython with zmq frontend)
#     common.console_encode should encode things as utf-8.
#     """
#     if compat.PY3:
#         raise nose.SkipTest

#     with tm.stdin_encoding(encoding=None):
#         result = printing.console_encode(u"\u05d0")
#         expected = u"\u05d0".encode('utf-8')
#         assert (result == expected)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
