"""
Plotting api
"""

# flake8: noqa

try:  # mpl optional
    from pandas.plotting import converter
    converter.register()  # needs to override so set_xlim works with str/number
except ImportError:
    pass

from pandas.plotting.misc import (scatter_matrix, radviz,
                                  andrews_curves, bootstrap_plot,
                                  parallel_coordinates, lag_plot,
                                  autocorrelation_plot)
from pandas.plotting.core import (boxplot, scatter_plot, grouped_hist,
                                  hist_frame, hist_series)
from pandas.plotting.style import plot_params
from pandas.plotting.tools import table
