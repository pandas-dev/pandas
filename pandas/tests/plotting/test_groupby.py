# coding: utf-8

""" Test cases for GroupBy.plot """


from pandas import Series, DataFrame
import pandas.util.testing as tm
import pandas.util._test_decorators as td

import numpy as np

from pandas.tests.plotting.common import TestPlotBase

import matplotlib.pyplot as plt
import pandas as pd
import pytest
from pandas._libs import reduction


@pytest.mark.mpl_image_compare
def test_first_plot_only_once():
    class multiplot:
        def __init__(self):
            self.fig, self.axes = plt.subplots(3, 1)
            plt.tight_layout()
            self.axidx = 0

        def plot(self, group, **kwargs):
            self.axes[self.axidx].scatter(group['x'], group['y'], **kwargs)
            self.axidx += 1
            return self.fig

    df = pd.DataFrame([[1, 2], [-3, -4], [5, 6], [-7, -8], [9, 10]],
                      index=['A', 'B', 'C', 'D', 'E'], columns=['x', 'y'])
    df['cat'] = [1, 2, 1, 2, 1]
    sdata = df.sort_values('cat')
    mlplt = multiplot()
    f = mlplt.plot
    names = pd.Int64Index([1, 2], dtype='int64', name='cat')
    starts = np.array([0, 3])
    ends = np.array([3, 5])

    results, mutated = reduction.apply_frame_axis0(
        sdata, f, names, starts, ends)

    return results[0]


@td.skip_if_no_mpl
class TestDataFrameGroupByPlots(TestPlotBase):

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
