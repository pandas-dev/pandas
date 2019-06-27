from pandas._config import get_option

from pandas.plotting._matplotlib.boxplot import (
    BoxPlot, boxplot, boxplot_frame, boxplot_frame_groupby)
from pandas.plotting._matplotlib.converter import deregister, register
from pandas.plotting._matplotlib.core import (
    AreaPlot, BarhPlot, BarPlot, HexBinPlot, LinePlot, PiePlot, ScatterPlot)
from pandas.plotting._matplotlib.hist import (
    HistPlot, KdePlot, hist_frame, hist_series)
from pandas.plotting._matplotlib.misc import (
    andrews_curves, autocorrelation_plot, bootstrap_plot, lag_plot,
    parallel_coordinates, radviz, scatter_matrix)
from pandas.plotting._matplotlib.timeseries import tsplot
from pandas.plotting._matplotlib.tools import table

if get_option("plotting.matplotlib.register_converters"):
    register(explicit=False)


__all__ = ['LinePlot', 'BarPlot', 'BarhPlot', 'HistPlot', 'BoxPlot', 'KdePlot',
           'AreaPlot', 'PiePlot', 'ScatterPlot', 'HexBinPlot', 'hist_series',
           'hist_frame', 'boxplot', 'boxplot_frame', 'boxplot_frame_groupby',
           'tsplot', 'table', 'andrews_curves', 'autocorrelation_plot',
           'bootstrap_plot', 'lag_plot', 'parallel_coordinates', 'radviz',
           'scatter_matrix', 'register', 'deregister']
