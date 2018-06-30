# -*- coding: utf-8 -*-
import numpy as np
from numpy.random import randint
from textwrap import dedent

import pytest
import pandas as pd

from pandas import DataFrame
from pandas import read_clipboard
from pandas import get_option
from pandas.compat import PY2
from pandas.util import testing as tm
from pandas.util.testing import makeCustomDataframe as mkdf
from pandas.io.clipboard.exceptions import PyperclipException
from pandas.io.clipboard import clipboard_set, clipboard_get


try:
    DataFrame({'A': [1, 2]}).to_clipboard()
    _DEPS_INSTALLED = 1
except (PyperclipException, RuntimeError):
    _DEPS_INSTALLED = 0


def build_kwargs(sep, excel):
    kwargs = {}
    if excel != 'default':
        kwargs['excel'] = excel
    if sep != 'default':
        kwargs['sep'] = sep
    return kwargs


@pytest.fixture(params=['delims', 'utf8', 'string', 'long', 'nonascii',
                        'colwidth', 'mixed', 'float', 'int'])
def df(request):
    data_type = request.param

    if data_type == 'delims':
        return pd.DataFrame({'a': ['"a,\t"b|c', 'd\tef´'],
                             'b': ['hi\'j', 'k\'\'lm']})
    elif data_type == 'utf8':
        return pd.DataFrame({'a': ['µasd', 'Ωœ∑´'],
                             'b': ['øπ∆˚¬', 'œ∑´®']})
    elif data_type == 'string':
        return mkdf(5, 3, c_idx_type='s', r_idx_type='i',
                    c_idx_names=[None], r_idx_names=[None])
    elif data_type == 'long':
        max_rows = get_option('display.max_rows')
        return mkdf(max_rows + 1, 3,
                    data_gen_f=lambda *args: randint(2),
                    c_idx_type='s', r_idx_type='i',
                    c_idx_names=[None], r_idx_names=[None])
    elif data_type == 'nonascii':
        return pd.DataFrame({'en': 'in English'.split(),
                             'es': 'en español'.split()})
    elif data_type == 'colwidth':
        _cw = get_option('display.max_colwidth') + 1
        return mkdf(5, 3, data_gen_f=lambda *args: 'x' * _cw,
                    c_idx_type='s', r_idx_type='i',
                    c_idx_names=[None], r_idx_names=[None])
    elif data_type == 'mixed':
        return DataFrame({'a': np.arange(1.0, 6.0) + 0.01,
                          'b': np.arange(1, 6),
                          'c': list('abcde')})
    elif data_type == 'float':
        return mkdf(5, 3, data_gen_f=lambda r, c: float(r) + 0.01,
                    c_idx_type='s', r_idx_type='i',
                    c_idx_names=[None], r_idx_names=[None])
    elif data_type == 'int':
        return mkdf(5, 3, data_gen_f=lambda *args: randint(2),
                    c_idx_type='s', r_idx_type='i',
                    c_idx_names=[None], r_idx_names=[None])
    else:
        raise ValueError


@pytest.mark.single
@pytest.mark.skipif(not _DEPS_INSTALLED,
                    reason="clipboard primitives not installed")
class TestClipboard(object):
    def check_round_trip_frame(self, data, excel=None, sep=None,
                               encoding=None):
        data.to_clipboard(excel=excel, sep=sep, encoding=encoding)
        result = read_clipboard(sep=sep or '\t', index_col=0,
                                encoding=encoding)
        tm.assert_frame_equal(data, result, check_dtype=False)

    # Test that default arguments copy as tab delimited
    def test_round_trip_frame(self, df):
        self.check_round_trip_frame(df)

    # Test that explicit delimiters are respected
    @pytest.mark.parametrize('sep', ['\t', ',', '|'])
    def test_round_trip_frame_sep(self, df, sep):
        self.check_round_trip_frame(df, sep=sep)

    # Test white space separator
    def test_round_trip_frame_string(self, df):
        df.to_clipboard(excel=False, sep=None)
        result = read_clipboard()
        assert df.to_string() == result.to_string()
        assert df.shape == result.shape

    # Two character separator is not supported in to_clipboard
    # Test that multi-character separators are not silently passed
    def test_excel_sep_warning(self, df):
        with tm.assert_produces_warning():
            df.to_clipboard(excel=True, sep=r'\t')

    # Separator is ignored when excel=False and should produce a warning
    def test_copy_delim_warning(self, df):
        with tm.assert_produces_warning():
            df.to_clipboard(excel=False, sep='\t')

    # Tests that the default behavior of to_clipboard is tab
    # delimited and excel="True"
    @pytest.mark.parametrize('sep', ['\t', None, 'default'])
    @pytest.mark.parametrize('excel', [True, None, 'default'])
    def test_clipboard_copy_tabs_default(self, sep, excel, df):
        kwargs = build_kwargs(sep, excel)
        df.to_clipboard(**kwargs)
        if PY2:
            # to_clipboard copies unicode, to_csv produces bytes. This is
            # expected behavior
            assert clipboard_get().encode('utf-8') == df.to_csv(sep='\t')
        else:
            assert clipboard_get() == df.to_csv(sep='\t')

    # Tests reading of white space separated tables
    @pytest.mark.parametrize('sep', [None, 'default'])
    @pytest.mark.parametrize('excel', [False])
    def test_clipboard_copy_strings(self, sep, excel, df):
        kwargs = build_kwargs(sep, excel)
        df.to_clipboard(**kwargs)
        result = read_clipboard(sep=r'\s+')
        assert result.to_string() == df.to_string()
        assert df.shape == result.shape

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

    def test_invalid_encoding(self, df):
        # test case for testing invalid encoding
        with pytest.raises(ValueError):
            df.to_clipboard(encoding='ascii')
        with pytest.raises(NotImplementedError):
            pd.read_clipboard(encoding='ascii')

    @pytest.mark.parametrize('enc', ['UTF-8', 'utf-8', 'utf8'])
    def test_round_trip_valid_encodings(self, enc, df):
        self.check_round_trip_frame(df, encoding=enc)
