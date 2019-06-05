"""
Plotting api
"""

# flake8: noqa

from pandas.plotting._core import boxplot
from pandas.plotting._misc import (
    andrews_curves, autocorrelation_plot, bootstrap_plot, lag_plot,
    parallel_coordinates, radviz, scatter_matrix)
from pandas.plotting._style import plot_params
from pandas.plotting._tools import table

try:
    from pandas.plotting._converter import (
        register as register_matplotlib_converters)
    from pandas.plotting._converter import (
        deregister as deregister_matplotlib_converters)
except ImportError:
    pass
