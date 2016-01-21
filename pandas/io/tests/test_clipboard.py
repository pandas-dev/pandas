# -*- coding: utf-8 -*-
import numpy as np
from numpy.random import randint

import nose
import pandas as pd

from pandas import DataFrame
from pandas import read_clipboard
from pandas import get_option
from pandas.util import testing as tm
from pandas.util.testing import makeCustomDataframe as mkdf, disabled


try:
    import pandas.util.clipboard  # noqa
except OSError:
    raise nose.SkipTest("no clipboard found")


@disabled
class TestClipboard(tm.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestClipboard, cls).setUpClass()
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
        cls.data_types = list(cls.data.keys())

    @classmethod
    def tearDownClass(cls):
        super(TestClipboard, cls).tearDownClass()
        del cls.data_types, cls.data

    def check_round_trip_frame(self, data_type, excel=None, sep=None):
        data = self.data[data_type]
        data.to_clipboard(excel=excel, sep=sep)
        if sep is not None:
            result = read_clipboard(sep=sep, index_col=0)
        else:
            result = read_clipboard()
        tm.assert_frame_equal(data, result, check_dtype=False)

    def test_round_trip_frame_sep(self):
        for dt in self.data_types:
            self.check_round_trip_frame(dt, sep=',')

    def test_round_trip_frame_string(self):
        for dt in self.data_types:
            self.check_round_trip_frame(dt, excel=False)

    def test_round_trip_frame(self):
        for dt in self.data_types:
            self.check_round_trip_frame(dt)

    def test_read_clipboard_infer_excel(self):
        from textwrap import dedent
        from pandas.util.clipboard import clipboard_set

        text = dedent("""
            John James	Charlie Mingus
            1	2
            4	Harry Carney
            """.strip())
        clipboard_set(text)
        df = pd.read_clipboard()

        # excel data is parsed correctly
        self.assertEqual(df.iloc[1][1], 'Harry Carney')

        # having diff tab counts doesn't trigger it
        text = dedent("""
            a\t b
            1  2
            3  4
            """.strip())
        clipboard_set(text)
        res = pd.read_clipboard()

        text = dedent("""
            a  b
            1  2
            3  4
            """.strip())
        clipboard_set(text)
        exp = pd.read_clipboard()

        tm.assert_frame_equal(res, exp)
