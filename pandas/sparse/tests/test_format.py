# -*- coding: utf-8 -*-
from __future__ import print_function

import numpy as np
import pandas as pd

import pandas.util.testing as tm
from pandas.compat import (is_platform_windows,
                           is_platform_32bit)
from pandas.core.config import option_context


use_32bit_repr = is_platform_windows() or is_platform_32bit()


class TestSeriesFormatting(tm.TestCase):

    _multiprocess_can_split_ = True

    @property
    def dtype_format_for_platform(self):
        return '' if use_32bit_repr else ', dtype=int32'

    def test_sparse_max_row(self):
        s = pd.Series([1, np.nan, np.nan, 3, np.nan]).to_sparse()
        result = repr(s)
        dfm = self.dtype_format_for_platform
        exp = ("0    1.0\n1    NaN\n2    NaN\n3    3.0\n"
               "4    NaN\ndtype: float64\nBlockIndex\n"
               "Block locations: array([0, 3]{0})\n"
               "Block lengths: array([1, 1]{0})".format(dfm))
        self.assertEqual(result, exp)

        with option_context("display.max_rows", 3):
            # GH 10560
            result = repr(s)
            exp = ("0    1.0\n    ... \n4    NaN\n"
                   "dtype: float64\nBlockIndex\n"
                   "Block locations: array([0, 3]{0})\n"
                   "Block lengths: array([1, 1]{0})".format(dfm))
            self.assertEqual(result, exp)

    def test_sparse_mi_max_row(self):
        idx = pd.MultiIndex.from_tuples([('A', 0), ('A', 1), ('B', 0),
                                         ('C', 0), ('C', 1), ('C', 2)])
        s = pd.Series([1, np.nan, np.nan, 3, np.nan, np.nan],
                      index=idx).to_sparse()
        result = repr(s)
        dfm = self.dtype_format_for_platform
        exp = ("A  0    1.0\n   1    NaN\nB  0    NaN\n"
               "C  0    3.0\n   1    NaN\n   2    NaN\n"
               "dtype: float64\nBlockIndex\n"
               "Block locations: array([0, 3]{0})\n"
               "Block lengths: array([1, 1]{0})".format(dfm))
        self.assertEqual(result, exp)

        with option_context("display.max_rows", 3):
            # GH 13144
            result = repr(s)
            exp = ("A  0    1.0\n       ... \nC  2    NaN\n"
                   "dtype: float64\nBlockIndex\n"
                   "Block locations: array([0, 3]{0})\n"
                   "Block lengths: array([1, 1]{0})".format(dfm))
            self.assertEqual(result, exp)
