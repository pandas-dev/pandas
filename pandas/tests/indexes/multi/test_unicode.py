# -*- coding: utf-8 -*-

import pandas as pd

from pandas import (MultiIndex, compat)
from pandas.compat import PY3, range, u

import pandas.util.testing as tm

from pandas.tests.indexes.common import Base

class TestUnicode(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def test_repr_roundtrip(self):

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

    def test_repr_with_unicode_data(self):
        with pd.core.config.option_context("display.encoding", 'UTF-8'):
            d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
            index = pd.DataFrame(d).set_index(["a", "b"]).index
            assert "\\u" not in repr(index)  # we don't want unicode-escaped

    def test_unicode_string_with_unicode(self):
        d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
        idx = pd.DataFrame(d).set_index(["a", "b"]).index

        if PY3:
            str(idx)
        else:
            compat.text_type(idx)

    def test_bytestring_with_unicode(self):
        d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
        idx = pd.DataFrame(d).set_index(["a", "b"]).index

        if PY3:
            bytes(idx)
        else:
            str(idx)
