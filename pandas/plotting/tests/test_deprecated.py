# coding: utf-8

import nose
import string

import pandas as pd
import pandas.util.testing as tm
from pandas.util.testing import slow

import numpy as np
from numpy.random import randn

import pandas.tools.plotting as plotting

from pandas.plotting.tests.common import TestPlotBase


"""
Test cases for plot functions imported from deprecated
pandas.tools.plotting
"""


@tm.mplskip
class TestDeprecatedNameSpace(TestPlotBase):

    @slow
    def test_scatter_plot_legacy(self):
        tm._skip_if_no_scipy()

        df = pd.DataFrame(randn(100, 2))

        with tm.assert_produces_warning(FutureWarning):
            plotting.scatter_matrix(df)

        with tm.assert_produces_warning(FutureWarning):
            pd.scatter_matrix(df)

    @slow
    def test_boxplot_deprecated(self):
        df = pd.DataFrame(randn(6, 4),
                          index=list(string.ascii_letters[:6]),
                          columns=['one', 'two', 'three', 'four'])
        df['indic'] = ['foo', 'bar'] * 3

        with tm.assert_produces_warning(FutureWarning):
            plotting.boxplot(df, column=['one', 'two'],
                             by='indic')

    @slow
    def test_grouped_hist_legacy(self):
        df = pd.DataFrame(randn(500, 2), columns=['A', 'B'])
        df['C'] = np.random.randint(0, 4, 500)
        df['D'] = ['X'] * 500

        with tm.assert_produces_warning(FutureWarning):
            plotting.grouped_hist(df.A, by=df.C)

    @slow
    def test_radviz_deprecated(self):
        df = self.iris
        with tm.assert_produces_warning(FutureWarning):
            plotting.radviz(frame=df, class_column='Name')


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
