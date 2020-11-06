""" Test cases for DataFrame.plot """

import numpy as np

import pandas.util._test_decorators as td

from pandas import DataFrame
import pandas._testing as tm
from pandas.tests.plotting.common import TestPlotBase


@td.skip_if_no_mpl
class TestDataFramePlotsGroupby(TestPlotBase):
    def setup_method(self, method):
        TestPlotBase.setup_method(self, method)
        import matplotlib as mpl

        mpl.rcdefaults()

        self.tdf = tm.makeTimeDataFrame()
        self.hexbin_df = DataFrame(
            {
                "A": np.random.uniform(size=20),
                "B": np.random.uniform(size=20),
                "C": np.arange(20) + np.random.uniform(size=20),
            }
        )

    def _assert_ytickslabels_visibility(self, axes, expected):
        for ax, exp in zip(axes, expected):
            self._check_visible(ax.get_yticklabels(), visible=exp)

    def _assert_xtickslabels_visibility(self, axes, expected):
        for ax, exp in zip(axes, expected):
            self._check_visible(ax.get_xticklabels(), visible=exp)

    def test_groupby_boxplot_sharey(self):
        # https://github.com/pandas-dev/pandas/issues/20968
        # sharey can now be switched check whether the right
        # pair of axes is turned on or off

        df = DataFrame(
            {
                "a": [-1.43, -0.15, -3.70, -1.43, -0.14],
                "b": [0.56, 0.84, 0.29, 0.56, 0.85],
                "c": [0, 1, 2, 3, 1],
            },
            index=[0, 1, 2, 3, 4],
        )

        # behavior without keyword
        axes = df.groupby("c").boxplot()
        expected = [True, False, True, False]
        self._assert_ytickslabels_visibility(axes, expected)

        # set sharey=True should be identical
        axes = df.groupby("c").boxplot(sharey=True)
        expected = [True, False, True, False]
        self._assert_ytickslabels_visibility(axes, expected)

        # sharey=False, all yticklabels should be visible
        axes = df.groupby("c").boxplot(sharey=False)
        expected = [True, True, True, True]
        self._assert_ytickslabels_visibility(axes, expected)

    def test_groupby_boxplot_sharex(self):
        # https://github.com/pandas-dev/pandas/issues/20968
        # sharex can now be switched check whether the right
        # pair of axes is turned on or off

        df = DataFrame(
            {
                "a": [-1.43, -0.15, -3.70, -1.43, -0.14],
                "b": [0.56, 0.84, 0.29, 0.56, 0.85],
                "c": [0, 1, 2, 3, 1],
            },
            index=[0, 1, 2, 3, 4],
        )

        # behavior without keyword
        axes = df.groupby("c").boxplot()
        expected = [True, True, True, True]
        self._assert_xtickslabels_visibility(axes, expected)

        # set sharex=False should be identical
        axes = df.groupby("c").boxplot(sharex=False)
        expected = [True, True, True, True]
        self._assert_xtickslabels_visibility(axes, expected)

        # sharex=True, yticklabels should be visible
        # only for bottom plots
        axes = df.groupby("c").boxplot(sharex=True)
        expected = [False, False, True, True]
        self._assert_xtickslabels_visibility(axes, expected)
