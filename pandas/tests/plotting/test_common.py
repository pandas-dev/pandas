import pytest
import numpy as np
from pandas import DataFrame
from pandas import unique
from pandas.tests.plotting.common import (
    _check_plot_works,
    _check_ticks_props,
    _gen_two_subplots,
)

plt = pytest.importorskip("matplotlib.pyplot")


class TestCommon:
    def test__check_ticks_props(self):
        # GH 34768
        df = DataFrame({"b": [0, 1, 0], "a": [1, 2, 3]})
        ax = _check_plot_works(df.plot, rot=30)
        ax.yaxis.set_tick_params(rotation=30)
        msg = "expected 0.00000 but got "
        with pytest.raises(AssertionError, match=msg):
            _check_ticks_props(ax, xrot=0)
        with pytest.raises(AssertionError, match=msg):
            _check_ticks_props(ax, xlabelsize=0)
        with pytest.raises(AssertionError, match=msg):
            _check_ticks_props(ax, yrot=0)
        with pytest.raises(AssertionError, match=msg):
            _check_ticks_props(ax, ylabelsize=0)

    def test__gen_two_subplots_with_ax(self):
        fig = plt.gcf()
        gen = _gen_two_subplots(f=lambda **kwargs: None, fig=fig, ax="test")
        # On the first yield, no subplot should be added since ax was passed
        next(gen)
        assert fig.get_axes() == []
        # On the second, the one axis should match fig.subplot(2, 1, 2)
        next(gen)
        axes = fig.get_axes()
        assert len(axes) == 1
        subplot_geometry = list(axes[0].get_subplotspec().get_geometry()[:-1])
        subplot_geometry[-1] += 1
        assert subplot_geometry == [2, 1, 2]

    def test_colorbar_layout(self):
        fig = plt.figure()

        axes = fig.subplot_mosaic(
            """
            AB
            CC
            """
        )

        x = [1, 2, 3]
        y = [1, 2, 3]

        cs0 = axes["A"].scatter(x, y)
        axes["B"].scatter(x, y)

        fig.colorbar(cs0, ax=[axes["A"], axes["B"]], location="right")
        DataFrame(x).plot(ax=axes["C"])
    
    def test_bar_subplot_stacking(self):
        #GH Issue 61018
        #Extracts height and location data 
        test_data = np.random.default_rng(3).integers(0,100,5)
        df = DataFrame({"a": test_data, "b": test_data[::-1]})
        ax = _check_plot_works(df.plot, subplots= [('a','b')], kind="bar", stacked=True)

        #get xy and height of squares that represent the data graphed from the df
        #we would expect the height value of A to be reflected in the Y coord of B
        data_from_plot_mat = [(x.get_x(), x.get_y(), x.get_height()) for x in ax[0].findobj(plt.Rectangle) if x.get_height() in test_data]
        data_from_plot_df = DataFrame(data = data_from_plot_mat, columns = ["x_coord", "y_coord", "height"])
        unique_x_loc = unique(data_from_plot_df["x_coord"])

        plot_a_df = data_from_plot_df.iloc[:len(test_data)]
        plot_b_df = data_from_plot_df.iloc[len(test_data):].reset_index()
        total_bar_height = plot_a_df["height"].add(plot_b_df["height"])

        #check number of bars matches the number of data plotted
        assert len(unique_x_loc) == len(test_data)

        #checks that the first set of bars are the correct height and that the second one starts at the top of the first, additional checks the combined height of the bars are correct
        assert (plot_a_df["height"] == test_data).all()
        assert (plot_b_df["y_coord"] == test_data).all()
        assert (total_bar_height == test_data + test_data[::-1]).all()



