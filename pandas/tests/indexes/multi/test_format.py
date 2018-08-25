import warnings

import pytest

import pandas as pd
from pandas import MultiIndex
import pandas.util.testing as tm


def test_dtype_str(indices):
    dtype = indices.dtype_str
    assert isinstance(dtype, str)
    assert dtype == str(indices.dtype)


def test_format(idx):
    idx.format()
    idx[:0].format()


def test_format_integer_names():
    index = MultiIndex(levels=[[0, 1], [0, 1]],
                       codes=[[0, 0, 1, 1], [0, 1, 0, 1]], names=[0, 1])
    index.format(names=True)


def test_format_sparse_config(idx):
    warn_filters = warnings.filters
    warnings.filterwarnings('ignore', category=FutureWarning,
                            module=".*format")
    # GH1538
    pd.set_option('display.multi_sparse', False)

    result = idx.format()
    assert result[1] == 'foo  two'

    tm.reset_display_options()

    warnings.filters = warn_filters


def test_format_sparse_display():
    index = MultiIndex(levels=[[0, 1], [0, 1], [0, 1], [0]],
                       codes=[[0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 0, 1],
                              [0, 1, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0]])

    result = index.format()
    assert result[3] == '1  0  0  0'


def test_repr_with_unicode_data():
    with pd.option_context("display.encoding", 'UTF-8'):
        d = {"a": ["\u05d0", 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
        index = pd.DataFrame(d).set_index(["a", "b"]).index
        assert "\\" not in repr(index)  # we don't want unicode-escaped


@pytest.mark.skip(reason="#22511 will remove this test")
def test_repr_roundtrip():

    mi = MultiIndex.from_product([list('ab'), range(3)],
                                 names=['first', 'second'])
    str(mi)

    tm.assert_index_equal(eval(repr(mi)), mi, exact=True)

    mi_u = MultiIndex.from_product(
        [list('ab'), range(3)], names=['first', 'second'])
    result = eval(repr(mi_u))
    tm.assert_index_equal(result, mi_u, exact=True)

    # formatting
    str(mi)

    # long format
    mi = MultiIndex.from_product([list('abcdefg'), range(10)],
                                 names=['first', 'second'])

    tm.assert_index_equal(eval(repr(mi)), mi, exact=True)

    result = eval(repr(mi_u))
    tm.assert_index_equal(result, mi_u, exact=True)


def test_unicode_string_with_unicode():
    d = {"a": ["\u05d0", 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
    idx = pd.DataFrame(d).set_index(["a", "b"]).index
    str(idx)


def test_repr_max_seq_item_setting(idx):
    # GH10182
    idx = idx.repeat(50)
    with pd.option_context("display.max_seq_items", None):
        repr(idx)
        assert '...' not in str(idx)


class TestRepr(object):

    def setup_class(self):
        n = 1000
        ci = pd.CategoricalIndex(list('a' * n) + (['abc'] * n))
        dti = pd.date_range('2000-01-01', freq='s', periods=n * 2)
        self.narrow_mi = pd.MultiIndex.from_arrays([ci, ci.codes + 9, dti],
                                                   names=['a', 'b', 'dti'])

        levels = [ci, ci.codes + 9, dti, dti, dti]
        names = ['a', 'b', 'dti_1', 'dti_2', 'dti_3']
        self.wide_mi = pd.MultiIndex.from_arrays(levels, names=names)

    def test_repr(self, idx):
        result = idx[:1].__repr__()
        expected = """MultiIndex([('foo', 'one')],
           dtype='object', names=['first', 'second'])"""
        assert result == expected

        result = idx.__repr__()
        expected = """MultiIndex([('foo', 'one'),
            ('foo', 'two'),
            ('bar', 'one'),
            ('baz', 'two'),
            ('qux', 'one'),
            ('qux', 'two')],
           dtype='object', names=['first', 'second'])"""
        assert result == expected

        with pd.option_context('display.max_seq_items', 5):
            result = idx.__repr__()
            expected = """MultiIndex([('foo', 'one'),
            ('foo', 'two'),
            ...
            ('qux', 'one'),
            ('qux', 'two')],
           dtype='object', names=['first', 'second'], length=6)"""
            assert result == expected

    def test_rjust(self):
        result = self.narrow_mi[:1].__repr__()
        expected = """\
MultiIndex([('a', 9, '2000-01-01 00:00:00')],
           dtype='object', names=['a', 'b', 'dti'])"""
        assert result == expected

        result = self.narrow_mi[::500].__repr__()
        expected = """\
MultiIndex([(  'a',  9, '2000-01-01 00:00:00'),
            (  'a',  9, '2000-01-01 00:08:20'),
            ('abc', 10, '2000-01-01 00:16:40'),
            ('abc', 10, '2000-01-01 00:25:00')],
           dtype='object', names=['a', 'b', 'dti'])"""
        assert result == expected

        result = self.narrow_mi.__repr__()
        expected = """\
MultiIndex([(  'a',  9, '2000-01-01 00:00:00'),
            (  'a',  9, '2000-01-01 00:00:01'),
            (  'a',  9, '2000-01-01 00:00:02'),
            (  'a',  9, '2000-01-01 00:00:03'),
            (  'a',  9, '2000-01-01 00:00:04'),
            (  'a',  9, '2000-01-01 00:00:05'),
            (  'a',  9, '2000-01-01 00:00:06'),
            (  'a',  9, '2000-01-01 00:00:07'),
            (  'a',  9, '2000-01-01 00:00:08'),
            (  'a',  9, '2000-01-01 00:00:09'),
            ...
            ('abc', 10, '2000-01-01 00:33:10'),
            ('abc', 10, '2000-01-01 00:33:11'),
            ('abc', 10, '2000-01-01 00:33:12'),
            ('abc', 10, '2000-01-01 00:33:13'),
            ('abc', 10, '2000-01-01 00:33:14'),
            ('abc', 10, '2000-01-01 00:33:15'),
            ('abc', 10, '2000-01-01 00:33:16'),
            ('abc', 10, '2000-01-01 00:33:17'),
            ('abc', 10, '2000-01-01 00:33:18'),
            ('abc', 10, '2000-01-01 00:33:19')],
           dtype='object', names=['a', 'b', 'dti'], length=2000)"""
        assert result == expected

    def test_tuple_width(self):
        result = self.wide_mi[:1].__repr__()
        expected = """MultiIndex([('a', 9, '2000-01-01 00:00:00', '2000-01-01 00:00:00', ...)],
           dtype='object', names=['a', 'b', 'dti_1', 'dti_2', 'dti_3'])"""
        assert result == expected

        result = self.wide_mi[:10].__repr__()
        expected = """\
MultiIndex([('a', 9, '2000-01-01 00:00:00', '2000-01-01 00:00:00', ...),
            ('a', 9, '2000-01-01 00:00:01', '2000-01-01 00:00:01', ...),
            ('a', 9, '2000-01-01 00:00:02', '2000-01-01 00:00:02', ...),
            ('a', 9, '2000-01-01 00:00:03', '2000-01-01 00:00:03', ...),
            ('a', 9, '2000-01-01 00:00:04', '2000-01-01 00:00:04', ...),
            ('a', 9, '2000-01-01 00:00:05', '2000-01-01 00:00:05', ...),
            ('a', 9, '2000-01-01 00:00:06', '2000-01-01 00:00:06', ...),
            ('a', 9, '2000-01-01 00:00:07', '2000-01-01 00:00:07', ...),
            ('a', 9, '2000-01-01 00:00:08', '2000-01-01 00:00:08', ...),
            ('a', 9, '2000-01-01 00:00:09', '2000-01-01 00:00:09', ...)],
           dtype='object', names=['a', 'b', 'dti_1', 'dti_2', 'dti_3'])"""
        assert result == expected

        result = self.wide_mi.__repr__()
        expected = """\
MultiIndex([(  'a',  9, '2000-01-01 00:00:00', '2000-01-01 00:00:00', ...),
            (  'a',  9, '2000-01-01 00:00:01', '2000-01-01 00:00:01', ...),
            (  'a',  9, '2000-01-01 00:00:02', '2000-01-01 00:00:02', ...),
            (  'a',  9, '2000-01-01 00:00:03', '2000-01-01 00:00:03', ...),
            (  'a',  9, '2000-01-01 00:00:04', '2000-01-01 00:00:04', ...),
            (  'a',  9, '2000-01-01 00:00:05', '2000-01-01 00:00:05', ...),
            (  'a',  9, '2000-01-01 00:00:06', '2000-01-01 00:00:06', ...),
            (  'a',  9, '2000-01-01 00:00:07', '2000-01-01 00:00:07', ...),
            (  'a',  9, '2000-01-01 00:00:08', '2000-01-01 00:00:08', ...),
            (  'a',  9, '2000-01-01 00:00:09', '2000-01-01 00:00:09', ...),
            ...
            ('abc', 10, '2000-01-01 00:33:10', '2000-01-01 00:33:10', ...),
            ('abc', 10, '2000-01-01 00:33:11', '2000-01-01 00:33:11', ...),
            ('abc', 10, '2000-01-01 00:33:12', '2000-01-01 00:33:12', ...),
            ('abc', 10, '2000-01-01 00:33:13', '2000-01-01 00:33:13', ...),
            ('abc', 10, '2000-01-01 00:33:14', '2000-01-01 00:33:14', ...),
            ('abc', 10, '2000-01-01 00:33:15', '2000-01-01 00:33:15', ...),
            ('abc', 10, '2000-01-01 00:33:16', '2000-01-01 00:33:16', ...),
            ('abc', 10, '2000-01-01 00:33:17', '2000-01-01 00:33:17', ...),
            ('abc', 10, '2000-01-01 00:33:18', '2000-01-01 00:33:18', ...),
            ('abc', 10, '2000-01-01 00:33:19', '2000-01-01 00:33:19', ...)],
           dtype='object', names=['a', 'b', 'dti_1', 'dti_2', 'dti_3'], length=2000)"""  # noqa
        assert result == expected
