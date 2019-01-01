# -*- coding: utf-8 -*-


import warnings

import pytest

from pandas.compat import PY3, range, u

import pandas as pd
from pandas import MultiIndex, compat
import pandas.util.testing as tm


def test_dtype_str(indices):
    dtype = indices.dtype_str
    assert isinstance(dtype, compat.string_types)
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
    with pd.core.config.option_context("display.encoding", 'UTF-8'):
        d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
        index = pd.DataFrame(d).set_index(["a", "b"]).index
        assert "\\u" not in repr(index)  # we don't want unicode-escaped


@pytest.mark.skip(reason="#22511 will remove this test")
def test_repr_roundtrip():

    mi = MultiIndex.from_product([list('ab'), range(3)],
                                 names=['first', 'second'])
    str(mi)

    if PY3:
        tm.assert_index_equal(eval(repr(mi)), mi, exact=True)
    else:
        result = eval(repr(mi))
        # string coerces to unicode
        tm.assert_index_equal(result, mi, exact=False)
        assert mi.get_level_values('first').inferred_type == 'string'
        assert result.get_level_values('first').inferred_type == 'unicode'

    mi_u = MultiIndex.from_product(
        [list(u'ab'), range(3)], names=['first', 'second'])
    result = eval(repr(mi_u))
    tm.assert_index_equal(result, mi_u, exact=True)

    # formatting
    if PY3:
        str(mi)
    else:
        compat.text_type(mi)

    # long format
    mi = MultiIndex.from_product([list('abcdefg'), range(10)],
                                 names=['first', 'second'])

    if PY3:
        tm.assert_index_equal(eval(repr(mi)), mi, exact=True)
    else:
        result = eval(repr(mi))
        # string coerces to unicode
        tm.assert_index_equal(result, mi, exact=False)
        assert mi.get_level_values('first').inferred_type == 'string'
        assert result.get_level_values('first').inferred_type == 'unicode'

    result = eval(repr(mi_u))
    tm.assert_index_equal(result, mi_u, exact=True)


def test_unicode_string_with_unicode():
    d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
    idx = pd.DataFrame(d).set_index(["a", "b"]).index

    if PY3:
        str(idx)
    else:
        compat.text_type(idx)


def test_bytestring_with_unicode():
    d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
    idx = pd.DataFrame(d).set_index(["a", "b"]).index

    if PY3:
        bytes(idx)
    else:
        str(idx)


def test_repr_max_seq_item_setting(idx):
    # GH10182
    idx = idx.repeat(50)
    with pd.option_context("display.max_seq_items", None):
        repr(idx)
        assert '...' not in str(idx)
