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
        test_data = np.random.default_rng(3).integers(0,100,5)
        df = DataFrame({"A": test_data, "B": test_data[::-1], "C": test_data[0]})
        ax = df.plot(subplots= [('A','B')], kind="bar", stacked=True)

        #finds all the rectangles that represent the values from both subplots
        data_from_subplots = [[(x.get_x(), x.get_y(), x.get_height()) for x in ax[i].findobj(plt.Rectangle) if x.get_height() in test_data] for i in range(0,2)]

        #get xy and height of squares that represent the data graphed from the df
        #we would expect the height value of A to be reflected in the Y coord of B in subplot 1
        subplot_data_df_list = []
        unique_x_loc_list = []
        for i in range(0,len(data_from_subplots)):
            subplot_data_df= DataFrame(data = data_from_subplots[i], columns = ["x_coord", "y_coord", "height"])
            unique_x_loc = unique(subplot_data_df["x_coord"])

            subplot_data_df_list.append(subplot_data_df)
            unique_x_loc_list.append(unique_x_loc)

        #Checks subplot 1
        plot_A_df = subplot_data_df_list[0].iloc[:len(test_data)]
        plot_B_df = subplot_data_df_list[0].iloc[len(test_data):].reset_index()
        total_bar_height = plot_A_df["height"].add(plot_B_df["height"])
        #check number of bars matches the number of data plotted
        assert len(unique_x_loc_list[0]) == len(test_data)
        #checks that the first set of bars are the correct height and that the second one starts at the top of the first, additional checks the combined height of the bars are correct
        assert (plot_A_df["height"] == test_data).all()
        assert (plot_B_df["y_coord"] == test_data).all()
        assert (total_bar_height == test_data + test_data[::-1]).all()

        #Checks subplot 2
        plot_C_df = subplot_data_df_list[1].iloc[:len(test_data)]
        #check number of bars matches the number of data plotted
        assert len(unique_x_loc_list[1]) == len(test_data)
        #checks that all the bars start at zero and are the correct height
        assert (plot_C_df["height"] == test_data[0]).all()
        assert (plot_C_df["y_coord"] == 0).all()




