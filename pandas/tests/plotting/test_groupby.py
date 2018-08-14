# coding: utf-8

""" Test cases for GroupBy.plot """


from pandas import Series, DataFrame
import pandas.util.testing as tm
import pandas.util._test_decorators as td

import numpy as np

from pandas.tests.plotting.common import TestPlotBase
import pytest


@td.skip_if_no_mpl
class TestDataFrameGroupByPlots(TestPlotBase):

    @pytest.mark.parametrize('equal_bins, bins, expected1, expected2', (
        (True, 10,
         [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
         [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4]),
        (True, None,
         [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
         [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4]),
        (True, [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4, 3.],
         [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
         [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4]),
        (False, [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4, 3.],
         [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
         [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4]),
        (False, 10,
         [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
         [-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8]),
        (False, None,
         [-3., -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
         [-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8])
    ))
    def test_hist_bins_match(self, equal_bins, bins, expected1, expected2):
        # GH-22222
        df = DataFrame(np.append(np.linspace(-3, 3, 100),
                                 np.linspace(-1, 1, 100)),
                       columns=['value'])
        df['group'] = [0] * 100 + [1] * 100
        g = df.groupby('group')['value']
        ax = g.hist(bins=bins, alpha=0.7, equal_bins=equal_bins)[0]

        # x-axis (leftmost) hist bar coordinates:
        # group0: points[:num_bins]
        # group1: points[num_bins:]
        num_bins = len(expected1)
        points = np.array([patch.get_bbox().get_points()
                           for patch in ax.patches])[:, 0, 0]

        tm.assert_almost_equal(points[:num_bins], np.array(expected1))
        tm.assert_almost_equal(points[num_bins:], np.array(expected2))
        tm.close()

    def test_series_groupby_plotting_nominally_works(self):
        n = 10
        weight = Series(np.random.normal(166, 20, size=n))
        height = Series(np.random.normal(60, 10, size=n))
        with tm.RNGContext(42):
            gender = np.random.choice(['male', 'female'], size=n)

        weight.groupby(gender).plot()
        tm.close()
        height.groupby(gender).hist()
        tm.close()
        # Regression test for GH8733
        height.groupby(gender).plot(alpha=0.5)
        tm.close()

    def test_plotting_with_float_index_works(self):
        # GH 7025
        df = DataFrame({'def': [1, 1, 1, 2, 2, 2, 3, 3, 3],
                        'val': np.random.randn(9)},
                       index=[1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0])

        df.groupby('def')['val'].plot()
        tm.close()
        df.groupby('def')['val'].apply(lambda x: x.plot())
        tm.close()

    def test_hist_single_row(self):
        # GH10214
        bins = np.arange(80, 100 + 2, 1)
        df = DataFrame({"Name": ["AAA", "BBB"],
                        "ByCol": [1, 2],
                        "Mark": [85, 89]})
        df["Mark"].hist(by=df["ByCol"], bins=bins)
        df = DataFrame({"Name": ["AAA"], "ByCol": [1], "Mark": [85]})
        df["Mark"].hist(by=df["ByCol"], bins=bins)

    def test_plot_submethod_works(self):
        df = DataFrame({'x': [1, 2, 3, 4, 5],
                        'y': [1, 2, 3, 2, 1],
                        'z': list('ababa')})
        df.groupby('z').plot.scatter('x', 'y')
        tm.close()
        df.groupby('z')['x'].plot.line()
        tm.close()

    def test_plot_kwargs(self):

        df = DataFrame({'x': [1, 2, 3, 4, 5],
                        'y': [1, 2, 3, 2, 1],
                        'z': list('ababa')})

        res = df.groupby('z').plot(kind='scatter', x='x', y='y')
        # check that a scatter plot is effectively plotted: the axes should
        # contain a PathCollection from the scatter plot (GH11805)
        assert len(res['a'].collections) == 1

        res = df.groupby('z').plot.scatter(x='x', y='y')
        assert len(res['a'].collections) == 1
