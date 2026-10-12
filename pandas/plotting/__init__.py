"""
Plotting public API.

Authors of third-party plotting backends should implement a module with a
public ``plot(data, kind, **kwargs)``. The parameter `data` will contain
the data structure and can be a `Series` or a `DataFrame`. For example,
for ``df.plot()`` the parameter `data` will contain the DataFrame `df`.
In some cases, the data structure is transformed before being sent to
the backend (see PlotAccessor.__call__ in pandas/plotting/_core.py for
the exact transformations).

The parameter `kind` will be one of:

- line
- bar
- barh
- box
- hist
- kde
- area
- pie
- scatter
- hexbin

See the pandas API reference for documentation on each kind of plot.

Any other keyword argument is currently assumed to be backend specific,
but some parameters may be unified and added to the signature in the
future (e.g. `title` which should be useful for any backend).

Besides ``plot``, the selected backend is used for these functions:

- hist_series (for `Series.hist`)
- hist_frame (for `DataFrame.hist`)
- boxplot_frame (for `DataFrame.boxplot`)
- boxplot_frame_groupby (for `DataFrameGroupBy.boxplot`)

The other functions in `pandas.plotting` (e.g. `pandas.plotting.boxplot`,
`scatter_matrix`, `register_matplotlib_converters`) always use Matplotlib,
regardless of the selected backend.

Backends are found through the ``pandas_plotting_backends`` entry point,
falling back to importing the backend name as a module.

Use the code in pandas/plotting/_matplotlib/ and
https://github.com/holoviz/hvplot as a reference on how to write a backend.
"""

from pandas.plotting._core import (
    PlotAccessor,
    boxplot,
    boxplot_frame,
    boxplot_frame_groupby,
    hist_frame,
    hist_series,
)
from pandas.plotting._misc import (
    andrews_curves,
    autocorrelation_plot,
    bootstrap_plot,
    deregister as deregister_matplotlib_converters,
    lag_plot,
    parallel_coordinates,
    plot_params,
    radviz,
    register as register_matplotlib_converters,
    scatter_matrix,
    table,
)

__all__ = [
    "PlotAccessor",
    "andrews_curves",
    "autocorrelation_plot",
    "bootstrap_plot",
    "boxplot",
    "boxplot_frame",
    "boxplot_frame_groupby",
    "deregister_matplotlib_converters",
    "hist_frame",
    "hist_series",
    "lag_plot",
    "parallel_coordinates",
    "plot_params",
    "radviz",
    "register_matplotlib_converters",
    "scatter_matrix",
    "table",
]
