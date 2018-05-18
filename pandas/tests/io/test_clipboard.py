# -*- coding: utf-8 -*-
import numpy as np
from numpy.random import randint
from textwrap import dedent

import pytest
import pandas as pd

from pandas import DataFrame
from pandas import read_clipboard
from pandas import get_option
from pandas.util import testing as tm
from pandas.util.testing import makeCustomDataframe as mkdf
from pandas.io.clipboard.exceptions import PyperclipException
from pandas.io.clipboard import clipboard_set, clipboard_get


try:
    DataFrame({'A': [1, 2]}).to_clipboard()
    _DEPS_INSTALLED = 1
except (PyperclipException, RuntimeError):
    _DEPS_INSTALLED = 0


@pytest.mark.single
@pytest.mark.skipif(not _DEPS_INSTALLED,
                    reason="clipboard primitives not installed")
class TestClipboard(object):

    @classmethod
    def setup_class(cls):
        cls.data = {}
        cls.data['string'] = mkdf(5, 3, c_idx_type='s', r_idx_type='i',
                                  c_idx_names=[None], r_idx_names=[None])
        cls.data['int'] = mkdf(5, 3, data_gen_f=lambda *args: randint(2),
                               c_idx_type='s', r_idx_type='i',
                               c_idx_names=[None], r_idx_names=[None])
        cls.data['float'] = mkdf(5, 3,
                                 data_gen_f=lambda r, c: float(r) + 0.01,
                                 c_idx_type='s', r_idx_type='i',
                                 c_idx_names=[None], r_idx_names=[None])
        cls.data['mixed'] = DataFrame({'a': np.arange(1.0, 6.0) + 0.01,
                                       'b': np.arange(1, 6),
                                       'c': list('abcde')})

        # Test columns exceeding "max_colwidth" (GH8305)
        _cw = get_option('display.max_colwidth') + 1
        cls.data['colwidth'] = mkdf(5, 3, data_gen_f=lambda *args: 'x' * _cw,
                                    c_idx_type='s', r_idx_type='i',
                                    c_idx_names=[None], r_idx_names=[None])
        # Test GH-5346
        max_rows = get_option('display.max_rows')
        cls.data['longdf'] = mkdf(max_rows + 1, 3,
                                  data_gen_f=lambda *args: randint(2),
                                  c_idx_type='s', r_idx_type='i',
                                  c_idx_names=[None], r_idx_names=[None])
        # Test for non-ascii text: GH9263
        cls.data['nonascii'] = pd.DataFrame({'en': 'in English'.split(),
                                             'es': 'en español'.split()})
        # unicode round trip test for GH 13747, GH 12529
        cls.data['utf8'] = pd.DataFrame({'a': ['µasd', 'Ωœ∑´'],
                                         'b': ['øπ∆˚¬', 'œ∑´®']})
        # Test for quotes and common delimiters in text
        cls.data['delim_symbols'] = pd.DataFrame({'a': ['"a,\t"b|c', 'd\tef´'],
                                                  'b': ['hi\'j', 'k\'\'lm']})
        cls.data_types = list(cls.data.keys())

    @classmethod
    def teardown_class(cls):
        del cls.data_types, cls.data

    def check_round_trip_frame(self, data_type, excel=None, sep=None,
                               encoding=None):
        data = self.data[data_type]
        if excel in [None, True] and sep is not None and len(sep) > 1:
            with tm.assert_produces_warning():
                data.to_clipboard(excel=excel, sep=sep, encoding=encoding)
        else:
            data.to_clipboard(excel=excel, sep=sep, encoding=encoding)

            if excel in [None, True] and sep is not None:
                # Expect Excel
                result = read_clipboard(sep=sep, index_col=0, encoding=encoding)
            elif excel in [None, True] and sep is None:
                # Expect Excel with tabs
                result = read_clipboard(sep='\t', index_col=0, encoding=encoding)
            else:
                # Expect df.__repr__ format
                result = read_clipboard(encoding=encoding)

            if excel in [None, True]:
                tm.assert_frame_equal(data, result, check_dtype=False)
            else:
                assert data.to_string() == result.to_string()

    def test_round_trip_frame_sep(self):
        for dt in self.data_types:
            self.check_round_trip_frame(dt, sep=',')
            self.check_round_trip_frame(dt, sep=r'\s+')
            self.check_round_trip_frame(dt, sep='|')
            self.check_round_trip_frame(dt, sep='\t')

    def test_round_trip_frame_string(self):
        for dt in self.data_types:
            self.check_round_trip_frame(dt, excel=False)

    def test_round_trip_frame(self):
        for dt in self.data_types:
            self.check_round_trip_frame(dt)

    def test_read_clipboard_infer_excel(self):
        # gh-19010: avoid warnings
        clip_kwargs = dict(engine="python")

        text = dedent("""
            John James	Charlie Mingus
            1	2
            4	Harry Carney
            """.strip())
        clipboard_set(text)
        df = pd.read_clipboard(**clip_kwargs)

        # excel data is parsed correctly
        assert df.iloc[1][1] == 'Harry Carney'

        # having diff tab counts doesn't trigger it
        text = dedent("""
            a\t b
            1  2
            3  4
            """.strip())
        clipboard_set(text)
        res = pd.read_clipboard(**clip_kwargs)

        text = dedent("""
            a  b
            1  2
            3  4
            """.strip())
        clipboard_set(text)
        exp = pd.read_clipboard(**clip_kwargs)

        tm.assert_frame_equal(res, exp)

    def test_excel_clipboard_tabs(self):
        for dt in self.data_types:
            for sep in ['\t', None, 'default', ',', '|']:
                for excel in [True, None, 'default', False]:
                    # Function default should be to to produce tab delimited
                    kwargs = {}
                    if excel != 'default':
                        kwargs['excel'] = excel
                    if sep != 'default':
                        kwargs['sep'] = sep
                    data = self.data[dt]
                    data.to_clipboard(**kwargs)
                    if sep in ['\t', None, 'default'] and excel is not False:
                        # Expect tab delimited
                        result = read_clipboard(sep='\t', index_col=0)
                        tm.assert_frame_equal(data, result, check_dtype=False)
                        assert clipboard_get().count('\t') > 0
                    elif excel is False:
                        # Expect spaces (ie. df.__repr__() default)
                        result = read_clipboard(sep=r'\s+')
                        assert result.to_string() == data.to_string()
                    else:
                        # Expect other delimited ',' and '|'
                        result = read_clipboard(sep=sep, index_col=0)
                        tm.assert_frame_equal(data, result, check_dtype=False)
                        assert clipboard_get().count(sep) > 0

    def test_invalid_encoding(self):
        # test case for testing invalid encoding
        data = self.data['string']
        with pytest.raises(ValueError):
            data.to_clipboard(encoding='ascii')
        with pytest.raises(NotImplementedError):
            pd.read_clipboard(encoding='ascii')

    def test_round_trip_valid_encodings(self):
        for enc in ['UTF-8', 'utf-8', 'utf8']:
            for dt in self.data_types:
                self.check_round_trip_frame(dt, encoding=enc)
