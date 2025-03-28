"""
Use this module directly:
    import xarray.plot as xplt

Or use the methods on a DataArray or Dataset:
    DataArray.plot._____
    Dataset.plot._____
"""

from xarray.plot.dataarray_plot import (
    contour,
    contourf,
    hist,
    imshow,
    line,
    pcolormesh,
    plot,
    step,
    surface,
)
from xarray.plot.dataset_plot import scatter
from xarray.plot.facetgrid import FacetGrid

__all__ = [
    "plot",
    "line",
    "step",
    "contour",
    "contourf",
    "hist",
    "imshow",
    "pcolormesh",
    "FacetGrid",
    "scatter",
    "surface",
]
